#This R code is not for running (because all original climate inputs were prepared in author's desktop only).
#Therefore we submitted .R rather than .Rmd file here.
#Overall, the work flow is:
#(1)

library(spgwr)
library(sp)

#load("E:/R code/Elevation_Global_CN_project.Rdata")

######

sites <- read.csv(file="E:/C-N cycling/Carbon allocation/GWR/final2_LAI.csv") #lat, lon, z, site...
sites <- sites[,1:4]
dim(sites) # now, do third analysis 448-

tmn <- list.files(path = "E:/climate/tmn/",full.names = T) #2000-2016 (1-17),1981-1999 (18-36)
tmx <- list.files(path = "E:/climate/tmx/",full.names = T) #2000-2016 (1-17),1981-1999 (18-36)
vap <- list.files(path = "E:/climate/vap/",full.names = T) #2000-2016 (1-17),1981-1999 (18-36)
rad <- list.files(path = "E:/climate/radiation/",full.names = T) #1981-1999 (1-19),2000-2016 (20-36)
pre <- list.files(path = "E:/climate/pre/",full.names = T) #1981-1999 (1-19),2000-2016 (20-36)

summary(pre)

#calculate monthly year data
x <- 1 #which year of data
i <- 1 #which points
a <- 1.5 # which degree of grid

final1 <- data.frame(matrix(NA)) # for final sites climate data
climatelist <- vector(mode = "list", length = length(tmx)) # for climate list
mylist <- vector(mode = "list", length = nrow(sites)) # for specifying gridded data at certain degrees
mylist2 <- vector(mode = "list", length = nrow(sites)) # for sites

#first round for CRU: 1:295,303:322,324:364,366:388; first round for WFDEI: 1:295,303:322,324:364,366:388
#for tmn, tmx, vap, pew
for (x in 1:51){ 
  print(x)
  climatelist[[x]] <- read.csv(file=tmn[x])
  names(climatelist[[x]]) <- c("z","lon","lat","a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12")
  for (i in c(448:672)){
    mylist[[i]]<- subset(climatelist[[x]],lon>(sites[i,2]-a)&lon<(sites[i,2]+a)&lat>(sites[i,1]-a)&lat<(sites[i,1]+a))
    coordinates(mylist[[i]]) <- c("lon","lat")
    gridded(mylist[[i]]) <- TRUE
    
    mylist2[[i]] <- sites[i,1:3]
    coordinates(mylist2[[i]]) <- c("lon","lat")
    final1[i,1+12*(x-1)] <- (gwr(a1 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,2+12*(x-1)] <- (gwr(a2 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,3+12*(x-1)] <- (gwr(a3 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,4+12*(x-1)] <- (gwr(a4 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,5+12*(x-1)] <- (gwr(a5 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,6+12*(x-1)] <- (gwr(a6 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,7+12*(x-1)] <- (gwr(a7 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,8+12*(x-1)] <- (gwr(a8 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,9+12*(x-1)] <- (gwr(a9 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,10+12*(x-1)] <- (gwr(a10 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,11+12*(x-1)] <- (gwr(a11 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,12+12*(x-1)] <- (gwr(a12 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  } # let's do site step by step
}

tmn_final <- final1
csvfile <- paste("E:/C-N cycling/Carbon allocation/","final1_tmn3.csv",sep = "")
write.csv(tmn_final, csvfile, row.names = TRUE)

for (x in 21:21){ 
  print(x)
  climatelist[[x]] <- read.csv(file=vap[x])
  names(climatelist[[x]]) <- c("z","lon","lat","a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12")
  for (i in c(448:672)){
  mylist[[i]]<- subset(climatelist[[x]],lon>(sites[i,2]-a)&lon<(sites[i,2]+a)&lat>(sites[i,1]-a)&lat<(sites[i,1]+a))
  coordinates(mylist[[i]]) <- c("lon","lat")
  gridded(mylist[[i]]) <- TRUE
  
  mylist2[[i]] <- sites[i,1:3]
  coordinates(mylist2[[i]]) <- c("lon","lat")
  final1[i,1+12*(x-1)] <- (gwr(a1 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  final1[i,2+12*(x-1)] <- (gwr(a2 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  final1[i,3+12*(x-1)] <- (gwr(a3 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  final1[i,4+12*(x-1)] <- (gwr(a4 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  final1[i,5+12*(x-1)] <- (gwr(a5 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  final1[i,6+12*(x-1)] <- (gwr(a6 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  final1[i,7+12*(x-1)] <- (gwr(a7 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  final1[i,8+12*(x-1)] <- (gwr(a8 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  final1[i,9+12*(x-1)] <- (gwr(a9 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  final1[i,10+12*(x-1)] <- (gwr(a10 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  final1[i,11+12*(x-1)] <- (gwr(a11 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  final1[i,12+12*(x-1)] <- (gwr(a12 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  } # let's do site step by step
}
vap_final <- final1
csvfile <- paste("E:/C-N cycling/Carbon allocation/","final1_vap3.csv",sep = "")
write.csv(vap_final, csvfile, row.names = FALSE)

for (x in 1:51){ 
  print(x)
  climatelist[[x]] <- read.csv(file=tmx[x])
  names(climatelist[[x]]) <- c("z","lon","lat","a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12")
  for (i in c(448:672)){
    mylist[[i]]<- subset(climatelist[[x]],lon>(sites[i,2]-a)&lon<(sites[i,2]+a)&lat>(sites[i,1]-a)&lat<(sites[i,1]+a))
    coordinates(mylist[[i]]) <- c("lon","lat")
    gridded(mylist[[i]]) <- TRUE
    
    mylist2[[i]] <- sites[i,1:3]
    coordinates(mylist2[[i]]) <- c("lon","lat")
    final1[i,1+12*(x-1)] <- (gwr(a1 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,2+12*(x-1)] <- (gwr(a2 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,3+12*(x-1)] <- (gwr(a3 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,4+12*(x-1)] <- (gwr(a4 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,5+12*(x-1)] <- (gwr(a5 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,6+12*(x-1)] <- (gwr(a6 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,7+12*(x-1)] <- (gwr(a7 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,8+12*(x-1)] <- (gwr(a8 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,9+12*(x-1)] <- (gwr(a9 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,10+12*(x-1)] <- (gwr(a10 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,11+12*(x-1)] <- (gwr(a11 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,12+12*(x-1)] <- (gwr(a12 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  } # let's do site step by step
}
tmx_final <- final1
csvfile <- paste("E:/C-N cycling/Carbon allocation/","final1_tmx3.csv",sep = "")
write.csv(tmx_final, csvfile, row.names = FALSE)

for (x in 1:51){ 
  print(x)
  climatelist[[x]] <- read.csv(file=pre[x])
  names(climatelist[[x]]) <- c("z","lon","lat","a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12")
  for (i in c(448:672)){
    mylist[[i]]<- subset(climatelist[[x]],lon>(sites[i,2]-a)&lon<(sites[i,2]+a)&lat>(sites[i,1]-a)&lat<(sites[i,1]+a))
    coordinates(mylist[[i]]) <- c("lon","lat")
    gridded(mylist[[i]]) <- TRUE
    
    mylist2[[i]] <- sites[i,1:3]
    coordinates(mylist2[[i]]) <- c("lon","lat")
    final1[i,1+12*(x-1)] <- (gwr(a1 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,2+12*(x-1)] <- (gwr(a2 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,3+12*(x-1)] <- (gwr(a3 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,4+12*(x-1)] <- (gwr(a4 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,5+12*(x-1)] <- (gwr(a5 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,6+12*(x-1)] <- (gwr(a6 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,7+12*(x-1)] <- (gwr(a7 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,8+12*(x-1)] <- (gwr(a8 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,9+12*(x-1)] <- (gwr(a9 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,10+12*(x-1)] <- (gwr(a10 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,11+12*(x-1)] <- (gwr(a11 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,12+12*(x-1)] <- (gwr(a12 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  } # let's do site step by step
}
pre_final <- final1
csvfile <- paste("E:/C-N cycling/Carbon allocation/","final1_pre3.csv",sep = "")
write.csv(pre_final, csvfile, row.names = FALSE)

#for rad
for (x in 1:31){ 
  print(x)
  climatelist[[x]] <- read.csv(file=rad[x])
  names(climatelist[[x]]) <- c("lon","lat","z","a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12")
  for (i in c(448:672)){
    mylist[[i]]<- subset(climatelist[[x]],lon>(sites[i,2]-a)&lon<(sites[i,2]+a)&lat>(sites[i,1]-a)&lat<(sites[i,1]+a))
    coordinates(mylist[[i]]) <- c("lon","lat")
    gridded(mylist[[i]]) <- TRUE
    
    mylist2[[i]] <- sites[i,1:3]
    coordinates(mylist2[[i]]) <- c("lon","lat")
    final1[i,1+12*(x-1)] <- (gwr(a1 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,2+12*(x-1)] <- (gwr(a2 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,3+12*(x-1)] <- (gwr(a3 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,4+12*(x-1)] <- (gwr(a4 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,5+12*(x-1)] <- (gwr(a5 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,6+12*(x-1)] <- (gwr(a6 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,7+12*(x-1)] <- (gwr(a7 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,8+12*(x-1)] <- (gwr(a8 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,9+12*(x-1)] <- (gwr(a9 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,10+12*(x-1)] <- (gwr(a10 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,11+12*(x-1)] <- (gwr(a11 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
    final1[i,12+12*(x-1)] <- (gwr(a12 ~ z, mylist[[i]], bandwidth = 1.06, fit.points = mylist2[[i]],predictions=TRUE))$SDF$pred
  } # let's do site step by step
}
rad_final <- final1
csvfile <- paste("E:/C-N cycling/Carbon allocation/","final1_rad3.csv",sep = "")
write.csv(rad_final, csvfile, row.names = FALSE)

tmx1 <- tmx_final
tmn1 <- tmn_final
vap1 <- vap_final
rad1 <- rad_final 
pre1 <- pre_final

#now, calculate average at specfic year, and also input LAI
sites <- read.csv(file="E:/C-N cycling/Carbon allocation/GWR/final2_LAI.csv") #lat, lon, z, site...FAPAR,LAI

#select data within years
tmxlist <- vector(mode = "list", length = nrow(sites)) # for specifying gridded data at certain degrees
tmnlist <- vector(mode = "list", length = nrow(sites)) # for specifying gridded data at certain degrees
vaplist <- vector(mode = "list", length = nrow(sites)) # for specifying gridded data at certain degrees
radlist <- vector(mode = "list", length = nrow(sites)) # for specifying gridded data at certain degrees
prelist <- vector(mode = "list", length = nrow(sites)) # for specifying gridded data at certain degrees

#calculate output datafrmae after averages
tmxlist0 <- tmxlist1
tmnlist0 <- tmnlist1
radlist0 <- radlist1
vaplist0 <- vaplist1

tmxlist1 <- data.frame(matrix(NA))
tmnlist1 <- data.frame(matrix(NA))
vaplist1 <- data.frame(matrix(NA))
radlist1 <- data.frame(matrix(NA))
prelist1 <- data.frame(matrix(NA))

i <-1
m <-1

#select specific year and calculate average
for (i in 448:672){
  tmxlist[[i]] <- tmx1[i,((sites$Begin_year[i]-1961)*12+1):((sites$End_year[i]-1961)*12+12)]
  tmnlist[[i]] <- tmn1[i,((sites$Begin_year[i]-1961)*12+1):((sites$End_year[i]-1961)*12+12)]
  vaplist[[i]] <- vap1[i,((sites$Begin_year[i]-1961)*12+1):((sites$End_year[i]-1961)*12+12)]
  radlist[[i]] <- rad1[i,((sites$Begin_year_radi[i]-1981)*12+1):((sites$End_year_radi[i]-1981)*12+12)]
  if (length(tmxlist[[i]]) == 12) {
    tmxlist1[i,] <- tmxlist[[i]]
    tmnlist1[i,] <- tmnlist[[i]]
    vaplist1[i,] <- vaplist[[i]]
    radlist1[i,] <- radlist[[i]]} else {for (m in 1:12){
      tmxlist1[i,m] <- rowMeans(tmxlist[[i]][,c(seq(m, length(tmxlist[[i]])-12+m, by=12))],na.rm=TRUE)
      tmnlist1[i,m] <- rowMeans(tmnlist[[i]][,c(seq(m, length(tmnlist[[i]])-12+m, by=12))],na.rm=TRUE)
      vaplist1[i,m] <- rowMeans(vaplist[[i]][,c(seq(m, length(vaplist[[i]])-12+m, by=12))],na.rm=TRUE)
      radlist1[i,m] <- rowMeans(radlist[[i]][,c(seq(m, length(radlist[[i]])-12+m, by=12))],na.rm=TRUE)
    }}
}

for (i in 448:672){
  prelist[[i]] <- pre1[i,((sites$Begin_year[i]-1961)*12+1):((sites$End_year[i]-1961)*12+12)]
  if (length(prelist[[i]]) == 12) {
    prelist1[i,] <- prelist[[i]]} else {for (m in 1:12){
      prelist1[i,m] <- rowMeans(prelist[[i]][,c(seq(m, length(prelist[[i]])-12+m, by=12))],na.rm=TRUE)
    }}
}


#output:
summary(tmxlist1)
summary(tmnlist1)
summary(vaplist1)
summary(radlist1) 
summary(prelist1)

#final
tmxlist3 <- rbind(tmxlist2,tmxlist1[448:672,])
tmnlist3 <- rbind(tmnlist2,tmnlist1[448:672,])
vaplist3 <- rbind(vaplist2,vaplist1[448:672,])
radlist3 <- rbind(radlist2,radlist1[448:672,])
prelist3 <- rbind(prelist2,prelist1[448:672,])


#now, calculate splash
tc <- Tg#see below
sw_in <- radlist3
pn <- prelist3
lat_all <-sites$lat
lon_all <-sites$lon
elev_all <-sites$z

m <- 1961 #year
j <- 1 #site

mydata <- data.frame(matrix(NA,672,672))

splash<-function(sw_in, tc, pn, lat,elev,outfolder=getwd()){
  source("E:/Master/splashv1/new/splashv1/const.R")
  source("E:/Master/splashv1/new/splashv1/evap.R")
  source("E:/Master/splashv1/new/splashv1/solar.R")
  source("E:/Master/splashv1/new/splashv1/splash.R")
  
  require(compiler)
  enableJIT(3)
  require(xts)
  # Extract time info from data
  y<-as.numeric(unique(format(time(pn),'%Y')))
  ny <- julian_day(y + 1, 1, 1) - julian_day(y, 1, 1)
  ztime<-time(pn)
  time.freq<-abs(as.numeric(ztime[1]-ztime[2], units = "days"))
  
  
  if (time.freq<2){		
    
    if (length(y)==1){
      initial<-spin_up(lat,elev, sw_in, tc, pn, y[1])
      
      result<-run_one_year(lat,elev, sw_in, tc, pn,initial,y[1])
      result<-xts(result,ztime)
    }
    else if(length(y)>1){
      
      end<-cumsum(ny)
      start<-end+1
      result<-list()
      initial<-spin_up(lat,elev, sw_in[1:ny[1]], tc[1:ny[1]], pn[1:ny[1]], y[1])
      
      result[[1]]<-run_one_year(lat,elev, sw_in[1:ny[1]], tc[1:ny[1]], pn[1:ny[1]],initial,y[1])
      
      for (i in 2:length(y)){
        stidx<-i-1
        # correct for leap years	
        if(ny[stidx]<ny[i]){
          result[[i]]<-run_one_year(lat,elev, sw_in[start[stidx]:end[i]], tc[start[stidx]:end[i]], pn[start[stidx]:end[i]],c(result[[stidx]]$wn,result[[stidx]]$wn[1]),y[i])
          
        }
        else if(ny[stidx]>ny[i]){
          result[[i]]<-run_one_year(lat,elev, sw_in[start[stidx]:end[i]], tc[start[stidx]:end[i]], pn[start[stidx]:end[i]],result[[stidx]]$wn[1:365],y[i])
        }
        else if(ny[stidx]==ny[i]){
          result[[i]]<-run_one_year(lat,elev, sw_in[start[stidx]:end[i]], tc[start[stidx]:end[i]], pn[start[stidx]:end[i]],result[[stidx]]$wn,y[i])
        }	
        
      }
      # order results as time series
      result<-Reduce(rbind,result)
      result<-xts(result,ztime)	
      
      
    }
    
    
    
    
    
  }
  
  else if (time.freq>20){		
    ztime.days<-seq(as.Date(paste(y[1],1,sep="-"),format="%Y-%j"),as.Date(paste(y[length(y)],ny[length(y)],sep="-"),format="%Y-%j"), by="day")
    if (length(y)==1){
      initial<-spin_up(lat,elev, sw_in, tc, pn, y[1])
      
      result<-run_one_year(lat,elev, sw_in, tc, pn,initial,y[1])
      result<-xts(result,ztime)
    }
    else if(length(y)>1){
      nm <- rep(12,length(y))
      end<-cumsum(nm)
      start<-end-11
      
      result<-list()
      initial<-spin_up(lat,elev, sw_in[1:end[1]], tc[1:end[1]], pn[1:end[1]], y[1])
      
      result[[1]]<-run_one_year(lat,elev, sw_in[1:end[1]], tc[1:end[1]], pn[1:end[1]],initial,y[1])
      
      for (i in 2:length(y)){
        stidx<-i-1
        # correct for leap years	
        if(ny[stidx]<ny[i]){
          result[[i]]<-run_one_year(lat,elev, sw_in[start[i]:end[i]], tc[start[i]:end[i]], pn[start[i]:end[i]],c(result[[stidx]]$wn,result[[stidx]]$wn[1]),y[i])
          
        }
        else if(ny[stidx]>ny[i]){
          result[[i]]<-run_one_year(lat,elev, sw_in[start[i]:end[i]], tc[start[i]:end[i]], pn[start[i]:end[i]],result[[stidx]]$wn[1:365],y[i])
        }
        else if(ny[stidx]==ny[i]){
          result[[i]]<-run_one_year(lat,elev, sw_in[start[i]:end[i]], tc[start[i]:end[i]], pn[start[i]:end[i]],result[[stidx]]$wn,y[i])
        }	
        
      }
      # order results as time series
      result<-Reduce(rbind,result)
      result<-xts(result,ztime.days)	
      
      
    }
    
    
    
    
    
  }
  return(result)	
}
source("E:/Master/splashv1/new/splashv1/splash_Main_point.R")
buildts<-function(data,y){
  
  require(xts)
  
  var.ts<-xts(data, seq(from=as.Date(paste0(y,"-","01"),format="%Y-%j"),to=as.Date(paste0(y,"-","365"),format="%Y-%j"),by="month"))
  
}

#for alpha
for (j in c(1:672)){
  for (m in sites$Begin_year[j]:sites$End_year[j]){
    pn1<-buildts(as.numeric(pn[j,]),m) #specify timeline 
    
    yunk1<- splash(sw_in=as.numeric(swin[j,]),tc=as.numeric(tc[j,]),pn=pn1,lat_all[j],elev_all[j])
    
    yunk2 <- apply.monthly(yunk1,mean)
    
    yunk3 <- as.numeric(yunk2[,9]/yunk2[,8]) # aet, pet
    
    mydata[j,c(((m-1961)*12+1):((m-1961)*12+12))]<- yunk3
  }
  print(j)
}
alpha <- rowMeans(mydata,na.rm = TRUE) 
sites$alpha <- alpha

#for cwd
for (j in c(366:672)){
  for (m in sites$Begin_year[j]:sites$End_year[j]){
    pn1<-buildts(as.numeric(pn[j,]),m) #specify timeline 
    
    yunk1<- splash(sw_in=as.numeric(swin[j,]),tc=as.numeric(tc[j,]),pn=pn1,lat_all[j],elev_all[j])
    
    yunk2 <- apply.monthly(yunk1,mean)
    
    yunk3 <- as.numeric(yunk2[,8]-yunk2[,9]) # aet, pet
    
    mydata[j,c(((m-1961)*12+1):((m-1961)*12+12))]<- yunk3
  }
  print(j)
}

cwd <- rowMeans(mydata,na.rm = TRUE) #mm/hr
sites$cwd_mm_h <- cwd

summary(lm(cwd~sites$cwd))


#Last but not least, calculate Ca at each site
Ca <- read.csv(file="E:/C-N cycling/Carbon allocation/GWR/Ca.csv") 
Caa <- Ca$Ca[1:56] # select Ca: 1961-2016

Calist1 <- data.frame(matrix(NA))

for (i in 1:nrow(sites)){
  Calist1[i,1] <- mean(Caa[((sites$Begin_year[i]-1961)+1):((sites$End_year[i]-1961)+1)])
}

#output
Caaaa <- Calist1$matrix.NA.

