
library(lattice)
library(sp)
library(gstat)
library(geoR)
library(raster)
library(reshape)
library(rgdal)
library(maptools)
library(plyr)
library(data.table)
library(foreach)
library(doParallel)
library(doSNOW)

registerDoParallel(3)  #sets number of cores to be used for operation (parallel processing)

##Specify "AT-max" or "AT-min"
obj1 = "AT-min"

##Set working directory
setwd("F:/Users/Jason/Desktop/Health_data/Daily_rasters")

coords.weights = read.csv("bg_cell_weights.csv", header=T)     #these are the weights each cell in the raster receives in calculating block group mean temperatures
#coords.weights = read.csv("zip_cell_weights.csv", header=T)   #same as above, but for zips

##Temp_data
temp=read.csv(paste(obj1,"csv",sep="."),header=T)
temp=melt(temp)
names(temp)=c("Date","SID","Temp")

##Station_locations
locs=read.csv("Station_covariates_2011.csv",header=T)
locations=locs
SID=as.data.frame(locations[,1])
names(SID)="SID"
locations=locations[,2:3]
locations=SpatialPointsDataFrame(locations,SID)
proj4string(locations)<-CRS("+init=epsg:3746")

##merged data
data2=merge(temp,locs,by=c("SID","SID"))

##Grid
grid=read.csv("Dane_grid_2011.csv",header=T)
coordinates(grid)=~X+Y
gridded(grid)=T
grid = as(grid, "SpatialPixelsDataFrame") # to full grid

dates=levels(factor(data2$Date))
list1 = list()

cl<-makeCluster(3) #change to your number of CPU cores (usually minus 1, so computer doesn't freeze up)
registerDoSNOW(cl)

ptm <- proc.time()
list1 = foreach(i=1:length(dates), .packages=c('geoR','sp','raster')) %dopar% { 
	tdate = dates[i]
	data4 = data2[data2$Date %in% tdate,]
	data4 = subset(data4, !is.na(Temp))

	geoRdata=as.data.frame(data4)
	geodata=as.geodata(geoRdata,4:5,3)
	IMP=geoRdata$Imp500
	LAKE=geoRdata$Lake.noscale.quarterexp
	geoRdata$Topo<-2/(1+exp((-geoRdata$ELEV)/6))-1	##logistic relationship
	ELEV=geoRdata$Topo

	#grid covariates
	grids=as.data.frame(grid)
	grids$Imp500[is.na(grids$Imp500)] <- 0
	grid2=data.frame(grids$X,grids$Y)
	names(grid2)=c("X","Y")
	IMP2=grids$Imp500
	LAKE2=grids$Lake.noscale.quarterexp
	grids$Topo<-2/(1+exp((-grids$ELEV)/6))-1	##logistic relationship
	ELEV2=grids$Topo

	## temperature fit
	vario=variog(geodata,estimator.type="modulus",max=20000)
	sill1=median(vario$v[5:11])
	#readline(prompt="Press [enter] to continue")

	#maximum likelihood fit for model
	sill=sill1
	range=7000

	ml.1=likfit(geodata,trend=~IMP+LAKE+ELEV,ini.cov.pars=c(sill,range),
		cov.model="exponential")
	ml.10=likfit(geodata,trend=~IMP+LAKE+ELEV,ini.cov.pars=c(sill,range),
		cov.model="spherical")
	covmod = ifelse(abs(ml.1$BIC)<abs(ml.10$BIC),"ml.1","ml.10")

	#Universal kriging (w/covariates)
	kc=krige.control(trend.d=~IMP+LAKE+ELEV,
		trend.l=~IMP2+LAKE2+ELEV2,obj.model=get(covmod))
	uk=krige.conv(geodata,loc=grid2,krige=kc)

	#Temps
	a=uk$predict
	uk.df=data.frame(a,grid2$X,grid2$Y)
	names(uk.df)=c("UK","X","Y")
	coordinates(uk.df)=~X+Y
	gridded(uk.df)=T
	uk.spdf = as(uk.df, "SpatialPixelsDataFrame")
	uk.r=raster(uk.spdf)
	
	matrix1 = t(as.matrix(uk.r))
	coords.weights$val = matrix1[coords.weights[,2]]
	coords.weights$wms = coords.weights$val * coords.weights$weight
	result1 = aggregate(cbind(wms,weight)~ID, data=coords.weights, FUN=sum)
	list1[[i]] = data.frame(t(result1$wms/result1$weight))
}
proc.time() - ptm
stopCluster(cl)

x=rbindlist(list1)
x$DATE = dates
write.csv(x, paste(obj1,"results.csv",sep="_"))


