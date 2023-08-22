## Local regression downscale using locally fitted models from 2bRegressionForests.R

library(rgdal)
library(fields)
library(raster)
library(ranger) 

ncols <- 14752
nrows <- 9712

#extent(matrix index boundaries) of nNA data in Landsat
bottom.left.col <- 2794
bottom.left.row <- 9479
top.left.col <- 4363
top.left.row <- 3504
bottom.right.col <- 8866
bottom.right.row <- 10732

block_size <- 750

#partitioning sequence of Landsat 
##into 750*750 block and get their start row and column index
rowBeg <- seq(top.left.row, bottom.right.row, block_size)
rowBeg <- rowBeg[-length(rowBeg)] # remove last item of this vector
blockBeg <- seq(bottom.left.col, nrows, block_size)
blockBeg <- blockBeg[-length(blockBeg)] # remove last item of this vector

bl.size <- block_size-1

#regression forests path
Rdata_path <- "/pl/active/rittger_esp/SierraBighorn/Rdata/localForests/regression_forest_750/forest/"

#path to save regression downscales
Downscaled_path <- "/pl/active/rittger_esp/SierraBighorn/downscaledv5_localForests/downscaled_regression_750/"

#Landsat data path
landsat.path <- "/pl/active/rittger_esp/SierraBighorn/Rdata/localForests/SierraBighornData/Landsat/SC/"
landsat.sc.file.names <- list.files(landsat.path)

#MODIS data and getting the dates
modis.path <- "/pl/active/rittger_esp/SierraBighorn/Rdata/localForests/SierraBighornData/MODIS/v03/"
modis.file.names <- list.files(modis.path)
mod.date <- NULL
for(i in 1:length(modis.file.names)){
  mod.date[i] <- regmatches(modis.file.names[i], regexpr("\\d{4}\\d{2}\\d{2}",modis.file.names[i]))
}
mod.date <- as.integer(mod.date)
mod.date <- strptime(x=as.character(mod.date),format="%Y%m%d")
mod.da <- mod.date$yday
mod.yr <- as.integer(format(mod.date,"%Y"))
these <- ((mod.yr %% 4) == 0) & (mod.da >= 59)
mod.da[these] <- mod.da[these] - 1
rm(these,mod.yr)
mod.da <- mod.da + 1
nday <- length(mod.da)

#predictors 
predictors.path <- "/pl/active/rittger_esp/SierraBighorn/predictors/"
predictors.file.names <- list.files(predictors.path)

out <- raster(paste0(predictors.path,"SouthernSierraNevada_Elevation.tif"))
elev <- t(matrix(values(out),nc=ncols,nr=nrows))
elev[elev < (-100)] <- NA
rm(out)

out <- raster(paste0(predictors.path,"SouthernSierraNevada_Slope.tif"))
slope <- t(matrix(values(out),nc=ncols,nr=nrows))
slope[slope < 0] <- NA
rm(out)

out <- raster(paste0(predictors.path,"SouthernSierraNevada_Aspect.tif"))
asp <- t(matrix(values(out),nc=ncols,nr=nrows))
asp[asp < (-2)] <- NA
rm(out)

out <- raster(paste0(predictors.path,"downscaled_s_sierra_winds_dec_april_climatology_nldas2.tif"))
windspeed <- t(matrix(values(out),nc=ncols,nr=nrows))
rm(out)

out <- raster(paste0(predictors.path, "SouthernSierraNevada_NorthWestBarrierDistance.tif"))
nw.barrierdist <- t(matrix(values(out),nc=ncols,nr=nrows))
rm(out)

out <- raster(paste0(predictors.path, "SouthernSierraNevada_SouthWestBarrierDistance.tif"))
sw.barrierdist <- t(matrix(values(out),nc=ncols,nr=nrows))
rm(out)

out <- raster(paste0(predictors.path, "SouthernSierraNevada_WestBarrierDistance.tif"))
w.barrierdist <- t(matrix(values(out),nc=ncols,nr=nrows))
rm(out)

out <- raster(paste0(predictors.path, "SouthernSierraNevada_SouthWestDistanceToWater.tif"))
sw.waterdist <- t(matrix(values(out),nc=ncols,nr=nrows))
rm(out)

out <- raster(paste0(predictors.path, "SouthernSierraNevada_WestDistanceToWater.tif"))
w.waterdist <- t(matrix(values(out),nc=ncols,nr=nrows))
rm(out)

## Standardize predictors
min.elev <- min(elev[-which(is.na(elev))])
max.elev <- max(elev[-which(is.na(elev))])
elev <- (elev-min.elev)/(max.elev - min.elev)
rm(min.elev, max.elev)

min.slope <- min(slope[-which(is.na(slope))])
max.slope <- max(slope[-which(is.na(slope))])
slope <- (slope-min.slope)/(max.slope - min.slope)
rm(min.slope, max.slope)

windspeed <- (windspeed-min(windspeed))/(max(windspeed) - min(windspeed))
nw.barrierdist <- nw.barrierdist/1000
sw.barrierdist <- sw.barrierdist/1000
w.barrierdist <- w.barrierdist/1000
sw.waterdist <- (sw.waterdist-min(sw.waterdist))/(max(sw.waterdist) - min(sw.waterdist))
w.waterdist <- (w.waterdist-min(w.waterdist))/(max(w.waterdist) - min(w.waterdist))

## Outer loop to downscale each day
for (pred.day in 1:length(mod.da)) {
  
  print(paste0("Starting prediction day ",pred.day))
  Pred.Start.time <- Sys.time()
  pred <- array(255,dim = c(ncols, nrows))
  
  ## Inner loop to downscale each block on a given day
  for (i in rowBeg) {
    
    for (j in blockBeg) {
      
      if(!file.exists(paste0(Rdata_path,"ranger.regression.750.",i,".",j,".RData"))){
        next
      }
      
      elev.block <- elev[i:(i+bl.size),j:(j+bl.size)]
      slope.block <- slope[i:(i+bl.size),j:(j+bl.size)]
      asp.block <- asp[i:(i+bl.size),j:(j+bl.size)]
      windspeed.block <- windspeed[i:(i+bl.size),j:(j+bl.size)]
      nw.barrierdist.block <- nw.barrierdist[i:(i+bl.size),j:(j+bl.size)]
      sw.barrierdist.block <- sw.barrierdist[i:(i+bl.size),j:(j+bl.size)]
      w.barrierdist.block <- w.barrierdist[i:(i+bl.size),j:(j+bl.size)]
      sw.waterdist.block <- sw.waterdist[i:(i+bl.size),j:(j+bl.size)]
      w.waterdist.block <- w.waterdist[i:(i+bl.size),j:(j+bl.size)]
      
      #asp and slope has missing data in the same places as elev.
      #nw.barrierdiast is missing in the same places as w.barrierdist
      theseNA <- union(which(is.na(elev.block)),which(is.na(w.barrierdist.block)))
      
      MOD.ras <- raster(paste0(modis.path, modis.file.names[pred.day]))
      MOD <- t(matrix(values(MOD.ras),nc=ncols ,nr=nrows))
      MOD.block <- MOD[i:(i+bl.size),j:(j+bl.size)]
      
      if(length(theseNA) == 0){
        daydat <- data.frame(elev=elev.block[1:(block_size*block_size)],
                             slope=slope.block[1:(block_size*block_size)],
                             asp=asp.block[1:(block_size*block_size)],
                             windspeed=windspeed.block[1:(block_size*block_size)],
                             nw.barrierdist=nw.barrierdist.block[1:(block_size*block_size)],
                             sw.barrierdist=sw.barrierdist.block[1:(block_size*block_size)],
                             w.barrierdist=w.barrierdist.block[1:(block_size*block_size)],
                             sw.waterdist=sw.waterdist.block[1:(block_size*block_size)],
                             w.waterdist=w.waterdist.block[1:(block_size*block_size)])
        daydat$mod <- MOD.block[1:(block_size*block_size)]
      }else{
        daydat <- data.frame(elev=elev.block[-theseNA],
                             slope=slope.block[-theseNA],
                             asp=asp.block[-theseNA],
                             windspeed=windspeed.block[-theseNA],
                             nw.barrierdist=nw.barrierdist.block[-theseNA],
                             sw.barrierdist=sw.barrierdist.block[-theseNA],
                             w.barrierdist=w.barrierdist.block[-theseNA],
                             sw.waterdist=sw.waterdist.block[-theseNA],
                             w.waterdist=w.waterdist.block[-theseNA])
        daydat$mod <- MOD.block[-theseNA]
      }
      
      daydat$da <- mod.da[pred.day]
      
      rm(MOD.ras, MOD, MOD.block)   
      
      ## Starting predictions
      print(paste0("Loading block ","i= ",i," j= ",j))
      load(file = paste0(Rdata_path,"ranger.regression.750.",i,".",j,".RData"))
      
      ranger.temp <- predict(ranger.regressor,data=daydat,num.threads=128)
      
      
      downscaled <- rep(255,times=block_size*block_size)
      if(length(theseNA) == 0){
        downscaled <- ranger.temp$predictions
      }else{
        downscaled[-theseNA] <- ranger.temp$predictions
      }
      
      pred[i:(i+bl.size),j:(j+bl.size)] <- downscaled
      
      rm(ranger.regressor, daydat, ranger.temp, downscaled, theseNA)
    }
  }
  
  pred.rast <- raster(paste0(landsat.path,landsat.sc.file.names[pred.day]))
  pred <- as.numeric(pred)
  values(pred.rast) <- c(t(matrix(pred,nr=ncols,nc=nrows)))
  
  print(paste("Saving Raster downscaled for day",pred.day))
  writeRaster(pred.rast, filename=paste0(Downscaled_path,"downscaled/",
    "downscaled.regression.750.", format(mod.date[pred.day],"%Y%m%d"),".tif"), 
    format="GTiff", option="COMPRESS=LZW", datatype="INT1U", overwrite=TRUE)
  
  Pred.End.time <- Sys.time()
  
  print(paste0(difftime(Pred.End.time, Pred.Start.time, units = "mins")," mins"))
  rm(pred.rast, pred)
}

