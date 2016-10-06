#####################
### open packages ###
#####################
library(dlnm)
library(splines)
library(tidyverse)
library(nlme)
library(lme4)

#################
### data prep ###
#################
## open files ##
data_path = "F:/Users/Jason/Desktop/Health_data/DLNMs/data_files"

dane_admits = read_csv(paste(data_path,"admissions_daily_countywide.csv",sep="/"))
zip_temps = read_csv(paste(data_path,"temperatures_ZCTA_daily.csv",sep="/"))
zip_admits = read_csv(paste(data_path,"admissions_daily_ZCTA.csv",sep="/"))
zip_census = read_csv(paste(data_path,"census_ZCTA.csv",sep="/"))

## prep zip_df ##
zip_df = merge(zip_admits, zip_temps, by=c("DATE","ZIP"), all.x=T)
zip_df = merge(zip_df, zip_census, by="ZIP", all.x=T)
zip_df$DATE = as.Date(zip_df$DATE, format="%m-%d-%y")
zip_df = zip_df[order(zip_df$DATE),]
zip_df$DOW = weekdays(zip_df$DATE)
zip_df$DOY = as.integer(strftime(zip_df$DATE, format = "%j"))
zip_df$WEEK = (zip_df$DOY-1)%/%7 + 1
zip_df$MONTH = as.integer(strftime(zip_df$DATE, format = "%m")) 
zip_df$YEAR = as.integer(strftime(zip_df$DATE, format = "%y"))
zip_df$WEEKEND = ifelse(zip_df$DOW == "Saturday" | zip_df$DOW == "Sunday", 1, 0)

## population weighted ATs for countywide analysis ##
w1 = zip_df
w1$wATMAX = w1$ATMAX*w1$POP/sum(zip_census$POP)
w1$wATMIN = w1$ATMIN*w1$POP/sum(zip_census$POP)
dane_temps <- aggregate(c(w1["wATMAX"],w1["wATMIN"]), by=w1[c("DATE")], FUN=sum, na.rm=T)
colnames(dane_temps) = c("DATE","ATMAX","ATMIN")

## prep dane_df
dane_admits$DATE = as.Date(dane_admits$DATE, format="%m-%d-%y")
dane_df = merge(dane_admits, dane_temps, by="DATE", all.x=T)
dane_df = dane_df[order(dane_df$DATE),]
dane_df$DOW = weekdays(dane_df$DATE)
dane_df$DOY = as.integer(strftime(dane_df$DATE, format = "%j"))
dane_df$WEEK = (dane_df$DOY-1)%/%7 + 1
dane_df$MONTH = as.integer(strftime(dane_df$DATE, format = "%m")) 
dane_df$YEAR = as.integer(strftime(dane_df$DATE, format = "%y"))
dane_df$WEEKEND = ifelse(dane_df$DOW == "Saturday" | dane_df$DOW == "Sunday", 1, 0)
dane_df$obsno = seq(1,nrow(dane_df),1)

obsno = data.frame(dane_df[,c("DATE","obsno")])
zip_df = merge(zip_df,obsno,by="DATE")

########################
### define functions ###
########################

### overdispersion testing function for poisson models ###
overdisp_fun <- function(model) {
	## number of variance parameters in 
	##   an n-by-n variance-covariance matrix
	vpars <- function(m) {
		nrow(m)*(nrow(m)+1)/2
		}
	model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
	rdf <- nrow(model.frame(model))-model.df
	rp <- residuals(model,type="pearson")
	Pearson.chisq <- sum(rp^2)
	prat <- Pearson.chisq/rdf
	pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
	c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

#############################
### countywide admissions ###
#############################
varknots=with(dane_df, equalknots(ATMAX,fun="ns",df=4))
lagknots=logknots(5,1)
cb1=with(dane_df, crossbasis(ATMAX,lag=7,argvar=list(fun="ns",knots=varknots,cen=20),arglag=list(knots=lagknots)))
mdl1=with(dane_df, glmer(ADMITS~cb1+ns(obsno,3)+DOW+(1|YEAR/WEEK/DOW),family=poisson())); AIC(mdl1)

overdisp_fun(mdl1)
pred1=crosspred(cb1,mdl1,by=1)

plot(pred1,xlab="ATmax",zlab="RR",theta=250,phi=30,lphi=30,main="3D graph of temperature effect")
plot(pred1,"contour",xlab="ATmax",key.title=title("RR"),plot.title=title("20C baseline",xlab="ATmax",ylab="Lag"))

plot(pred1,"slices",var=20,ci="n",col=1,ylim=c(0.95,1.1),lwd=1.5,main="Lag response curves for different temperatures; reference temp 20C")
for(i in 1:3) lines(pred1,"slices",var=c(19,30,35)[i],col=i,lwd=1.5)
legend("topright",paste("ATmax=",c(19,30,35)),col=1:3,lwd=1.5)

plot(pred1,"slices",var=c(15,25,30,35),col=2,ci.arg=list(density=40,col=grey(0.7)))

##extract cumulative effect of 32C (relative to 20C) over 5 days of lag, plus 95% CI:
pred1$allRRfit["32"]	##cumulative
cbind(pred1$allRRlow,pred1$allRRhigh)["32",]	##95% CI
cbind(sum(zip_census$POP),mean(dane_df$ADMITS))


############################################
######## individual zip admissions #########
############################################
zips = levels(factor(zip_census$ZIP))

for (i in 1:length(zips)){
   tryCatch({
	z1 = zips[i]
	df1 = subset(zip_df,ZIP==z1)
	varknots=with(df1, equalknots(ATMAX,fun="ns",df=4))
	lagknots=logknots(5,1)
	cb1=with(df1, crossbasis(ATMAX,lag=7,argvar=list(fun="ns",knots=varknots,cen=20),arglag=list(knots=lagknots)))
	mdl1=with(df1, glmer(ADMITS~cb1+ns(obsno,3)+DOW+(1|YEAR/WEEK/DOW),family=poisson())); AIC(mdl1)

	pred1=crosspred(cb1,mdl1,by=1)

	#overdisp_fun(mdl1)
	#plot(pred1,"contour",xlab="ATmax",key.title=title("RR"),plot.title=title("20C baseline",xlab="ATmax",ylab="Lag"))
	#plot(pred1,"slices",var=20,ci="n",col=1,ylim=range(pred1$matRRfit),lwd=1.5,main="Lag response curves for different temperatures; reference temp 20C")
	#for(i in 1:3) lines(pred1,"slices",var=c(19,30,35)[i],col=i,lwd=1.5)
	#legend("topright",paste("ATmax=",c(19,30,35)),col=1:3,lwd=1.5)

	outname = paste(z1,"curves.tiff",sep="_")
	tiff(outname,height=2400,width=700,res=300)
	par(mar=c(0.1,0.1,0.1,0.1), mai = c(0.1, 0.1, 0.1, 0.1))
	plot(pred1,"slices",var=c(15,25,30,35),col=2,ci.arg=list(density=40,col=grey(0.7)),xlab='n')
	dev.off()

	##extract cumulative effect of 32C (relative to 20C) over 5 days of lag, plus 95% CI:
	pred1$allRRfit["32"]	##cumulative
	cbind(pred1$allRRlow,pred1$allRRhigh)["32",]	##95% CI
	cbind(mean(df1$POP),mean(df1$ADMITS))
   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#######################################
######## inter-zip admissions #########
#######################################
varknots = with(zip_df, equalknots(ATMAX,fun="ns",df=5))
lagknots = logknots(5,1)
zip_df = scale(zip_df)
cb1 = with(zip_df, crossbasis(ATMAX,lag=7,argvar=list(fun="ns",knots=varknots,cen=20),arglag=list(knots=lagknots)))
md11 = with(zip_df, glmer(ADMITS ~ offset(log(POP)) + scale(MEDIAN_INC) + scale(PNONWHITE) + scale(PDISABLED) + DOW + cb1 + ns(obsno,3) + (1|ZIP/YEAR/DOY), family=poisson))
mdl2 = with(zip_df, glmer(ADMITS ~ offset(log(POP)) + scale(MEDIAN_INC) + scale(PNONWHITE) + scale(PDISABLED) + DOW + cb1 + (1|ZIP/YEAR/DOY), family=poisson))

pred1=crosspred(cb1,mdl1,by=1)

plot(pred1,xlab="Temperature",zlab="RR",theta=250,phi=30,lphi=30,main="3D graph of temperature effect")
plot(pred1,"slices",var=20,ci="n",col=1,ylim=c(0.95,1.1),lwd=1.5,main="Lag response curves for different temperatures; reference temp 20C")
for(i in 1:3) lines(pred1,"slices",var=c(19,30,35)[i],col=i,lwd=1.5)
legend("topright",paste("Temperature=",c(19,30,35)),col=1:3,lwd=1.5)

plot(pred1,"slices",var=c(15,25,30,35),col=2,ci.arg=list(density=40,col=grey(0.7)))

##extract cumulative effect of 32C (relative to 20C) over 5 days of lag, plus 95% CI:
pred3$allRRfit["32"]	##cumulative
cbind(pred3$allRRlow,pred3$allRRhigh)["32",]	##95% CI
#cbind(mean(data3$pop),mean(data3$admissions))


dd=read.csv("Rdata.csv",header=T)
dd=subset(dd,month>4)
hist(dd$Atmax,breaks=15,xlab=expression(AT[max]~"("*degree*"C)"),ylab="Days")

#################################
###### Model comparisons ########
#################################
varknots=equalknots(Tmax,fun="ns",df=3)
cb3=crossbasis(Tmax,lag=10,argvar=list(fun="ns",knots=varknots,cen=20),arglag=list(fun="ns",df=3))
model1=glmer(admissions~ offset(log(pop))+Per.capita.income + p.nonwhite + p.disabled +dow + cb3 +ns(obsno,3) + (1|doy/zip),family=poisson)

pred3=crosspred(cb3,model1,by=1)
















##########################
#### reducing a DLNM #####
##########################
cb4=crossbasis(Atmax,lag=3,argvar=list(df=5,cen=18.33),arglag=list(fun="strata",breaks=1))
model4=glm(admissions~cb4+ns(doy,3)+dow,family=quasipoisson(),data)
pred4=crosspred(cb4,model4,by=1)

redall=crossreduce(cb4,model4)	#overall cumulative
redlag=crossreduce(cb4,model4,type="lag",value=5)	#lag specific exposure/response relationship
redvar=crossreduce(cb4,model4,type="var",value=33)	#temperature specific lag/response relationship

plot(pred4,"overall",xlab="Temperature",ylab="RR",ylim=c(0.8,1.6),main="Overall cumulative association")
lines(redall,ci="lines",col=4,lty=2)
legend("top",c("Original","Reduced"),col=c(2,4),lty=1:2,ins=0.1)

b4=onebasis(0:40,knots=attributes(cb4)$arglag$knots,int=TRUE,cen=FALSE)
pred4b=crosspred(b4,coef=coef(redvar),vcov=vcov(redvar),model.link="log",by=1)

plot(pred4,"slices",var=10,ylab="RR",ylim=c(0.95,1.05),main="Predictor-specific association at 33C")
lines(redvar,ci="lines",col=4,lty=2)
points(pred4b,col=1,pch=19,cex=0.6)
legend("top",c("Original","Reduced","Reconstructed"),col=c(2,4,1),lty=c(1:2,NA),pch=c(NA,NA,19),pt.cex=0.6,ins=0.1)

cb1=crossbasis(Atmin,lag=15,argvar=list(fun="lin",cen=18.66),arglag=list(fun="poly",degree=4))
summary(cb1)
model1=glm(admissions~cb1+ns(doy,7)+dow,family=quasipoisson(),data)
pred1=crosspred(cb1,model1,at=0:40,bylag=0.2,cumul=TRUE)

plot(pred1,"slices",var=30,col=3,ylab="RR",ci.arg=list(density=15,lwd=2),main="Temp v admissions")
plot(pred1,"slices",var=29,cumul=TRUE,ylab="Cumulative RR",main="Temp v admissions (cumulative)")
	
	#var=10 sets plot to display the effect of a 10 unit increase in X on Y given the specified set of lags

plot(pred1,"overall",xlab="temp",ci="lines",ylim=c(0.9,1.3),lwd=2,ci.arg=list(col=1,lty=3),main="Overall cumulative association for 5 lags")
	#overall cumulative exposure-response relationship


##extract cumulative effect of 10 unit increase over 15 days of lag, plus 95% CI:
pred1$allRRfit["10"]	##cumulative
cbind(pred1$allRRlow,pred1$allRRhigh)["10",]	##95% CI


