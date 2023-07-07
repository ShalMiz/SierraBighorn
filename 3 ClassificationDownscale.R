###############################################################################################
# Local Classification Downscale Script
###############################################################################################

#___________________________________
#  Loading Libraries
#___________________________________

library(rgdal)
library(fields)
library(raster)
library(ranger) 

ncols <- 14752
nrows <- 9712

#___________________________________
#  Importing saved models 
#___________________________________

#extent(matrix index boundaries) of nNA data in Landsat
bottom.left.col <- 2794
bottom.left.row <- 9479
top.left.col <- 4363
top.left.row <- 3504
bottom.right.col <- 8866
bottom.right.row <- 10732

block_size <- 450

#___________________________________
#  Partitioning sequence of Landsat into 450*450 blocks
#  and get their start row and column index
#___________________________________

rowBeg <- seq(top.left.row, bottom.right.row, block_size)
rowBeg <- rowBeg[-length(rowBeg)] # remove last item of this vector
blockBeg <- seq(bottom.left.col, nrows, block_size)
blockBeg <- blockBeg[-length(blockBeg)] # remove last item of this vector

bl.size <- block_size-1

#classification forests path
Rdata_path <- "~/SierraBighorn/Rdata/localForests/classification_forest/forest"

#path to save classification downscales
Downscaled_path <- "~/SierraBighorn/downscaledv5_localForests/downscaled_classification/"

#Landsat data path
landsat.path <- "~/SierraBighorn/Rdata/localForests/SierraBighornData/Landsat/SC/"
landsat.sc.file.names <- list.files(landsat.path)

#MODIS data and getting the dates
modis.path <- "~/SierraBighorn/Rdata/localForests/SierraBighornData/MODIS/v03/"
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
predictors.path <- "~/SierraBighorn/predictors/"
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


#___________________________________
#  Data Pre-processing
#___________________________________

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

#___________________________________
#  Downscaling
#___________________________________

#outer loop to downscale each day
for (pred.day in 1:length(mod.da)) {
  
  print(Sys.time())
  print(paste0("Starting prediction day ",pred.day))
  pred <- array(255,dim = c(ncols, nrows))
  pred.prob <- array(255,dim = c(ncols, nrows,3))
  
  #inner loop to downscale each block on a given day
  for (i in rowBeg) {
    
    for (j in blockBeg) {
      
      #move to downscale next index if such forest is not available 
      if(!file.exists(paste0(Rdata_path,"ranger.classifier.",i,".",j,".Rda"))){
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
      
      #starting predictions
      load(file = paste0(Rdata_path,"ranger.classifier.",i,".",j,".Rda"))
      
      ranger.temp <- predict(ranger.classifier,data=daydat,num.threads=128)
      
      idx <- array(dim=c(3,2))
      cat.names <- c("btwn","hundred","zero")
      
      #some models only have two classes out of "btwn","hundred","zero"
      if(length(dimnames(ranger.temp$predictions)[[2]])==2){
        pred.cats <- intersect(cat.names,dimnames(ranger.temp$predictions)[[2]])
        N.pred.cat <- setdiff(cat.names,dimnames(ranger.temp$predictions)[[2]])
        idx[1,1] <- which(dimnames(ranger.temp$predictions)[[2]]==pred.cats[1])
        idx[2,1] <- which(dimnames(ranger.temp$predictions)[[2]]==pred.cats[2])
        idx[3,1] <- 0
        idx[1,2] <- which(cat.names==pred.cats[1])
        idx[2,2] <- which(cat.names==pred.cats[2])
        idx[3,2] <- which(cat.names==N.pred.cat)
      }
      
      #some blocks have NAs on predictors
      downscaled <- rep(255,times=block_size*block_size)
      if(length(theseNA) == 0){
        downscaled <- colnames(ranger.temp$predictions)[max.col(ranger.temp$predictions,"first")]
        pred[i:(i+bl.size),j:(j+bl.size)] <- downscaled
        #adjusting for when there are only two classes
        if(length(dimnames(ranger.temp$predictions)[[2]])==2){
          #saving probabilities
          downscaled <- rep(255,times=block_size*block_size)
          downscaled <- as.integer(round(100*ranger.temp$predictions[,idx[1,1]]))
          pred.prob[i:(i+bl.size),j:(j+bl.size),idx[1,2]] <- downscaled
          downscaled <- as.integer(round(100*ranger.temp$predictions[,idx[2,1]]))
          pred.prob[i:(i+bl.size),j:(j+bl.size),idx[2,2]] <- downscaled
          pred.prob[i:(i+bl.size),j:(j+bl.size),idx[3,2]] <- 0
        }else{
          #saving probabilities
          downscaled <- rep(255,times=block_size*block_size)
          downscaled <- as.integer(round(100*ranger.temp$predictions[,1]))
          pred.prob[i:(i+bl.size),j:(j+bl.size),1] <- downscaled
          downscaled <- as.integer(round(100*ranger.temp$predictions[,2]))
          pred.prob[i:(i+bl.size),j:(j+bl.size),2] <- downscaled
          downscaled <- as.integer(round(100*ranger.temp$predictions[,3]))
          pred.prob[i:(i+bl.size),j:(j+bl.size),3] <- downscaled
        }
      }else{
        downscaled <- rep(255,times=block_size*block_size)
        downscaled[-theseNA] <- colnames(ranger.temp$predictions)[max.col(ranger.temp$predictions,"first")]
        pred[i:(i+bl.size),j:(j+bl.size)] <- downscaled
        if(length(dimnames(ranger.temp$predictions)[[2]])==2){
          downscaled <- rep(255,times=block_size*block_size)
          #saving probabilities
          downscaled[-theseNA] <- as.integer(round(100*ranger.temp$predictions[,idx[1,1]]))
          pred.prob[i:(i+bl.size),j:(j+bl.size),idx[1,2]] <- downscaled
          downscaled[-theseNA] <- as.integer(round(100*ranger.temp$predictions[,idx[2,1]]))
          pred.prob[i:(i+bl.size),j:(j+bl.size),idx[2,2]] <- downscaled
          downscaled[-theseNA] <- 0
          pred.prob[i:(i+bl.size),j:(j+bl.size),idx[3,2]] <- downscaled
        }else{
          downscaled <- rep(255,times=block_size*block_size)
          downscaled[-theseNA] <- colnames(ranger.temp$predictions)[max.col(ranger.temp$predictions,"first")]
          pred[i:(i+bl.size),j:(j+bl.size)] <- downscaled
          #saving probabilities
          downscaled <- rep(255,times=block_size*block_size)
          downscaled[-theseNA] <- as.integer(round(100*ranger.temp$predictions[,1]))
          pred.prob[i:(i+bl.size),j:(j+bl.size),1] <- downscaled
          downscaled[-theseNA] <- as.integer(round(100*ranger.temp$predictions[,2]))
          pred.prob[i:(i+bl.size),j:(j+bl.size),2] <- downscaled
          downscaled[-theseNA] <- as.integer(round(100*ranger.temp$predictions[,3]))
          pred.prob[i:(i+bl.size),j:(j+bl.size),3] <- downscaled
        }
      }
      
      rm(ranger.classifier, daydat, ranger.temp, downscaled, theseNA)
    }
  }
  
  pred.rast <- raster(paste0(landsat.path,landsat.sc.file.names[pred.day]))
  pred[which(pred=="hundred")] <- 2
  pred[which(pred=="btwn")] <- 1
  pred[which(pred=="zero")] <- 0
  pred <- as.numeric(pred)
  
  values(pred.rast) <- c(t(matrix(pred,nr=ncols,nc=nrows)))
  print("Saving Raster downscaled")
  writeRaster(pred.rast, filename=paste0(Downscaled_path,"downscaled/","downscaled.", format(mod.date[pred.day],"%Y%m%d"),".tif"),     format="GTiff",option="COMPRESS=LZW",datatype="INT1U",overwrite=TRUE)
  
  values(pred.rast) <- c(t(matrix(pred.prob[,,1],nr=ncols,nc=nrows)))
  print("Saving Raster prob.btwn")
  writeRaster(pred.rast, filename=paste0(Downscaled_path,"prob.btwn/","prob.btwn.", format(mod.date[pred.day],"%Y%m%d"),".tif"),     format="GTiff",option="COMPRESS=LZW",datatype="INT1U",overwrite=TRUE)
  
  values(pred.rast) <- c(t(matrix(pred.prob[,,2],nr=ncols,nc=nrows)))
  print("Saving Raster prob.btwn")
  writeRaster(pred.rast, filename=paste0(Downscaled_path,"prob.hundred/","prob.hundred.", format(mod.date[pred.day],"%Y%m%d"),".tif"),     format="GTiff",option="COMPRESS=LZW",datatype="INT1U",overwrite=TRUE)
  
  values(pred.rast) <- c(t(matrix(pred.prob[,,3],nr=ncols,nc=nrows)))
  print("Saving Raster prob.btwn")
  writeRaster(pred.rast, filename=paste0(Downscaled_path,"prob.zero/","prob.zero.", format(mod.date[pred.day],"%Y%m%d"),".tif"),     format="GTiff",option="COMPRESS=LZW",datatype="INT1U",overwrite=TRUE)
  
  rm(pred.rast, pred, pred.prob)
}


