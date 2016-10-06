
library(spdep)
library(sp)
library(rgdal)
library(INLA)
library(maptools)
library(rgeos)
library(reshape2)

df.hosp = read.csv("F:/Users/Jason/Desktop/Health_data/INLA_models/Data/admits_daily_ZCTA.csv",header=T)
df.temps = read.csv("F:/Users/Jason/Desktop/Health_data/INLA_models/Data/zip_temperatures_daily.csv",header=T)
df.daily = merge(df.hosp,df.temps,by=c("DATE","ZIP"),all.x=T)

zips = readShapePoly("F:/Users/Jason/Desktop/Health_data/INLA_models/Data/Shapefiles/dane_zip_data.shp")
temp = poly2nb(zips)
nb2INLA("LDN.graph",temp)
adj1 = paste(getwd(),"/LDN.graph",sep="")
H = inla.read.graph(filename="LDN.graph")

data1 = zips@data; data1$ID = as.numeric(data1$ID); data1$MEDIAN_INC = as.numeric(as.character(data1$MEDIAN_INC))
data.final = merge(df.daily,data1,by="ZIP",all.x=T)
Nareas <- length(data1[,1])

summary(data1)

data1$rate = round(data1$ADMISSIONS/data1$POP * 100000)
formula1 = rate ~ 1 + PDISABLED + PNONWHITE + TMIN + f(ID, model="bym", graph=adj1)
mod1 = inla(formula1,family="poisson",data=data1,control.comput=list(dic=TRUE))
round(mod1$summary.fixed,3)
head(round(mod1$summary.random$ID,3))

inla.emarginal(exp,mod1$marginals.fixed[[1]])		#dane county mean rates (9.65%)
inla.qmarginal(c(0.025,0.975),inla.tmarginal(exp,mod1$marginals.fixed[[1]]))	#CI for dane county rates

csi = mod1$marginals.random$ID[1:Nareas]
zeta = lapply(csi,function(x) inla.emarginal(exp,x))
zeta.cutoff = quantile(as.numeric(zeta))
cat.zeta = cut(unlist(zeta),breaks=zeta.cutoff,include.lowest=T)
maps.cat.zeta = data.frame(SP_ID=data1$ID,cat.zeta=cat.zeta)
data.zips = attr(zips,"data")
attr(zips,"data") = merge(data.zips,maps.cat.zeta,by.x="ID",by.y="SP_ID")
spplot(obj=zips, zcol="cat.zeta",col.regions=gray(seq(0.9,0.1,length=4)),asp=1)

a = 0
prob.csi = lapply(csi, function(x) {1-inla.pmarginal(a,x)})
prob.csi.cutoff = quantile(as.numeric(prob.csi))
cat.prob.csi = cut(unlist(prob.csi),breaks=prob.csi.cutoff,include.lowest=T)
maps.cat.prob.csi = data.frame(SP_ID=data1$ID,cat.prob.csi=cat.prob.csi)
data.zips = attr(zips,"data")
attr(zips,"data") = merge(data.zips,maps.cat.prob.csi,by.x="ID",by.y="SP_ID")
spplot(obj=zips, zcol="cat.prob.csi.y",col.regions=gray(seq(0.9,0.1,length=4)),asp=1)

mat.marg = matrix(NA,nrow=Nareas, ncol=100000)
m = mod1$marginals.random$ID
for (i in 1:Nareas){
	u = m[[Nareas+i]]
	mat.marg[i,] = inla.rmarginal(100000,u)
}
var.u = apply(mat.marg,2,var)
var.v = inla.rmarginal(100000, inla.tmarginal(function(x) 1/x,
	mod1$marginals.hyper$"Precision for SP_ID (iid component)"))
perc.var.u = mean(var.u/(var.u+var.v)); perc.var.u

marg.hyper = inla.hyperpar.sample(100000,mod1)
perc.var.u1 = mean(marg.hyper[,1]/(marg.hyper[,1]+marg.hyper[,2])); perc.var.u1

df.hosp = read.csv("F:/Users/Jason/Desktop/Health_data/INLA_models/Data/admits_daily_ZCTA.csv",header=T)
df.temps = read.csv("F:/Users/Jason/Desktop/Health_data/INLA_models/Data/zip_temperatures_daily.csv",header=T)
df.daily = merge(df.hosp,df.temps,by=c("DATE","ZIP"),all.x=T)

zips = readShapePoly("F:/Users/Jason/Desktop/Health_data/INLA_models/Data/Shapefiles/dane_zip_data.shp")
temp = poly2nb(zips)
nb2INLA("LDN.graph",temp)
adj1 = paste(getwd(),"/LDN.graph",sep="")
H = inla.read.graph(filename="LDN.graph")

data1 = zips@data; data1$ID = as.numeric(data1$ID); data1$MEDIAN_INC = as.numeric(as.character(data1$MEDIAN_INC))
data1$TMIN = NULL; data1$TMAX = NULL; data1$RATE = NULL; data1$JUL_MAX = NULL; data1$JUL_MIN = NULL; data1$GOEGRAPHY = NULL; data1$ID2 = NULL
data.final = merge(df.daily,data1,by="ZIP",all.x=T)
Nareas <- length(data1[,1])

dates = data.frame(levels(factor(data.final$DATE))); colnames(dates)="DATE"
dates$DATE = as.Date(dates$DATE,"%m-%d-%y")
dates$DATE = dates[order(dates$DATE),]
dates$DAY = seq(1,nrow(dates),1)
data.final$DATE = as.Date(data.final$DATE,"%m-%d-%y")
data.final = merge(data.final,dates,by="DATE",all.x=T)
data.final$ID1 = data.final$ID

##linear temporal trend
formula1 = ADMITS ~ TMIN + PDISABLED + PNONWHITE + f(ID, model="bym", graph=H, constr=TRUE) + 
              f(ID1, DAY, model="iid", constr=TRUE) + 
              DAY

mdl1 = inla(formula1,family="poisson",data=data.final,E=EXPECTED,
            control.predictor=list(compute=TRUE),
            control.compute=list(dic=TRUE,cpo=TRUE))

##fixed effects
round(mdl1$summary.fixed[,1:5],3)

##posterior mean of spatial effect
m = mdl1$marginals.random[[1]][1:64]
prob.csi = unlist(lapply(m,function(x) inla.emarginal(exp,x)))
    
##map relative risks 
prob.csi.cutoff = quantile(as.numeric(prob.csi))
cat.prob.csi = cut(unlist(prob.csi),breaks=prob.csi.cutoff,include.lowest=T)
maps.cat.prob.csi = data.frame(SP_ID=data1$ID,cat.prob.csi=cat.prob.csi)
data.zips = attr(zips,"data")
attr(zips,"data") = merge(data.zips, maps.cat.prob.csi, by.x = "ID", by.y="SP_ID")
spplot(obj=zips, zcol="cat.prob.csi.x",col.regions=gray(seq(0.9,0.1,length=4)),asp=1)
