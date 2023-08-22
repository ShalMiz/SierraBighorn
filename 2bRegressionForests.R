
# Local regression forest training

library(rgdal)
library(fields)
library(raster)
library(ranger) 

set.seed(999)

#data locations
landsat.path <- "/pl/active/rittger_esp/SierraBighorn/Rdata/localForests/SierraBighornData/Landsat/SC/"
landsat.sc.file.names <- list.files(landsat.path)
sat.mask.path <- "/pl/active/rittger_esp/SierraBighorn/Rdata/localForests/SierraBighornData/Landsat/sat.mask/"
sat.mask.file.names <- list.files(sat.mask.path)
modis.path <- "/pl/active/rittger_esp/SierraBighorn/Rdata/localForests/SierraBighornData/MODIS/v03/"
modis.file.names <- list.files(modis.path)
predictors.path <- "/pl/active/rittger_esp/SierraBighorn/predictors/"
predictors.file.names <- list.files(predictors.path)
Rdata_path <- "/pl/active/rittger_esp/SierraBighorn/Rdata/localForests/regression_forest/"

ncols <- 14752
nrows <- 9712

#predictors 
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

## Rescaling elevation, slope, wind, barrier dists
# note that some blocks have elevation and slope NAs
elev_min <- min(elev)
elev_max <- max(elev)
if(length(which(is.na(elev)))>0){
  elev_min <- min(elev[-which(is.na(elev))])
  elev_max <- max(elev[-which(is.na(elev))])
}
elev <- (elev-elev_min)/(elev_max - elev_min)
#asp <- (asp-mean(asp))/sd(asp)
rm(elev_min,elev_max)

slope_min <- min(slope)
slope_max <- max(slope)
if(length(which(is.na(slope)))>0){
  slope_min <- min(slope[-which(is.na(slope))])
  slope_max <- max(slope[-which(is.na(slope))])
}
slope <- (slope-slope_min)/(slope_max - slope_min)
rm(slope_min,slope_max)

windspeed <- (windspeed-min(windspeed))/(max(windspeed) - min(windspeed))
nw.barrierdist <- nw.barrierdist/1000
sw.barrierdist <- sw.barrierdist/1000
w.barrierdist <- w.barrierdist/1000
sw.waterdist <- (sw.waterdist-min(sw.waterdist))/(max(sw.waterdist) - min(sw.waterdist))
w.waterdist <- (w.waterdist-min(w.waterdist))/(max(w.waterdist) - min(w.waterdist))

# mod.da creation
mod.date <- NULL
for(i in 1:length(modis.file.names)){
  mod.date[i] <- regmatches(modis.file.names[i], regexpr("\\d{4}\\d{2}\\d{2}",modis.file.names[i]))
}

mod.date <- as.integer(mod.date)
mod.date <- strptime(x=as.character(mod.date),format="%Y%m%d")

# Fix count for leap years for every day past Feb 29 = day 59 of year
mod.da <- mod.date$yday
mod.yr <- as.integer(format(mod.date,"%Y"))

these <- ((mod.yr %% 4) == 0) & (mod.da >= 59)
mod.da[these] <- mod.da[these] - 1
rm(these,mod.yr,mod.date)

mod.da <- mod.da + 1
nday <- length(mod.da)

## Store variable importances for each local RF
counter <- 1
var.importance <- list()
var.importance.id <- list()
#var.importance <- array(dim = c(260,13)) # 260 = 20 * 13 = ceiling(14752*9712/750^2)

## extent(matrix index boundaries) of nNA data in Landsat
bottom.left.col <- 2794
bottom.left.row <- 9479
top.left.col <- 4363
top.left.row <- 3504
bottom.right.col <- 8866
bottom.right.row <- 10732

## Block size for each local RF
block_size <- 750 # was 450 in SM's original code
train.pct <- 0.5

## Partitioning sequence of Landsat 
## Into 750*750 block and get their start row and column index
rowBeg <- seq(top.left.row, bottom.right.row, block_size)
rowBeg <- rowBeg[-length(rowBeg)] # remove last item of this vector
blockBeg <- seq(bottom.left.col, nrows, block_size)
blockBeg <- blockBeg[-length(blockBeg)] # remove last item of this vector

bl.size <- block_size-1

counter <- 1

set.seed(56210)
for(i in rowBeg){
  for(j in blockBeg){
    print(paste0("Starting model training block: ",i," ",j," at: ",Sys.time()))
    
    ## Array to save training indexes for training chosen from each day from the block
    train.idx <- array(dim=c(nday,block_size*block_size*train.pct*1.1)) # 10% buffer
    test.idx <- array(dim=c(nday,block_size*block_size*(1-train.pct)*1.1)) # 10% buffer
    
    ## Select predictor blocks
    elev.block <- elev[i:(i+bl.size),j:(j+bl.size)]
    slope.block <- slope[i:(i+bl.size),j:(j+bl.size)]
    asp.block <- asp[i:(i+bl.size),j:(j+bl.size)]
    windspeed.block <- windspeed[i:(i+bl.size),j:(j+bl.size)]
    nw.barrierdist.block <- nw.barrierdist[i:(i+bl.size),j:(j+bl.size)]
    sw.barrierdist.block <- sw.barrierdist[i:(i+bl.size),j:(j+bl.size)]
    w.barrierdist.block <- w.barrierdist[i:(i+bl.size),j:(j+bl.size)]
    sw.waterdist.block <- sw.waterdist[i:(i+bl.size),j:(j+bl.size)]
    w.waterdist.block <- w.waterdist[i:(i+bl.size),j:(j+bl.size)]
    
    ## indexes of NAs in elevation in the selected block for later to be removed from training 
    elev.NA <- which(is.na(elev.block))
    
    ## list to save training data from each day
    train.LST <- train.MOD <- train.slope <- train.asp  <- train.elev <- train.windspeed <-
      train.nw.barrierdist <- train.sw.barrierdist <- train.w.barrierdist <-
      train.sw.waterdist <- train.w.waterdist <- list()
    test.LST <- test.MOD <- test.slope <- test.asp  <- test.elev <- test.windspeed <-
      test.nw.barrierdist <- test.sw.barrierdist <- test.w.barrierdist <-
      test.sw.waterdist <- test.w.waterdist <- list()

    train.day <- test.day <- list()
    
    ## Selecting training data from each day
    for(day in 1:nday){
      LST.ras <- raster(paste0(landsat.path,landsat.sc.file.names[day]))
      LST <- t(matrix(values(LST.ras),nc=ncols,nr=nrows))
      LST.block <- LST[i:(i+bl.size),j:(j+bl.size)]
      LST.NA <- which(LST.block==255)
      rm(LST.ras, LST)
      
      sat.mask.ras <- raster(paste0(sat.mask.path,sat.mask.file.names[day]))
      sat.mask <- t(matrix(values(sat.mask.ras),nc=ncols,nr=nrows))
      sat.block <- sat.mask[i:(i+bl.size),j:(j+bl.size)]
      LST.sat <- which(sat.block==1)
      rm(sat.mask.ras, sat.mask, sat.block)

      LST.zero.hundred <- which( (LST.block == 0) | (LST.block >= 100) )
      
      # pool of indexes to sample training data from
      these.good.day <- (1:(block_size*block_size))[-c(LST.zero.hundred, LST.sat, LST.NA,
        elev.NA)]

      ## Reqiure more than 100 of (0,100) values, then use half of the samples for training
      if(length(these.good.day) <= 100){
        next
      }

      ## Sample each day
      train.size <- floor(length(these.good.day)*train.pct)
      good.ind.rand <- sample(x=these.good.day,size=length(these.good.day),replace=FALSE)
      train.indices.day <- good.ind.rand[1:train.size]
      test.indices.day <- good.ind.rand[(train.size+1):length(good.ind.rand)]
      rm(good.ind.rand)
      
      ## Store training/testing indices
      if(length(train.indices.day)!=0){
        train.idx[day,1:length(train.indices.day)] <- train.indices.day
        test.idx[day,1:length(test.indices.day)] <- test.indices.day
      }
      
      ## Grab training/testing Landsat
      train.LST[[day]] <- LST.block[train.indices.day]
      test.LST[[day]] <- LST.block[test.indices.day]
      rm(LST.block, LST.NA)
      
      ## Grab training/testing  MODIS
      MOD.ras <- raster(paste0(modis.path, modis.file.names[day]))
      MOD <- t(matrix(values(MOD.ras),nc=ncols ,nr=nrows))
      MOD.block <- MOD[i:(i+bl.size),j:(j+bl.size)]

      train.MOD[[day]] <- MOD.block[train.indices.day]
      test.MOD[[day]] <- MOD.block[test.indices.day]
      rm(MOD.ras, MOD, MOD.block)
      
      train.slope[[day]] <- slope.block[train.indices.day]
      train.asp[[day]] <- asp.block[train.indices.day]
      train.elev[[day]] <- elev.block[train.indices.day]
      train.windspeed[[day]] <- windspeed.block[train.indices.day]
      train.nw.barrierdist[[day]] <- nw.barrierdist.block[train.indices.day]
      train.sw.barrierdist[[day]] <- sw.barrierdist.block[train.indices.day]
      train.w.barrierdist[[day]] <- w.barrierdist.block[train.indices.day]
      train.sw.waterdist[[day]] <- sw.waterdist.block[train.indices.day]
      train.w.waterdist[[day]] <- w.waterdist.block[train.indices.day]
      train.day[[day]] <- length(train.indices.day)

      test.slope[[day]] <- slope.block[test.indices.day]
      test.asp[[day]] <- asp.block[test.indices.day]
      test.elev[[day]] <- elev.block[test.indices.day]
      test.windspeed[[day]] <- windspeed.block[test.indices.day]
      test.nw.barrierdist[[day]] <- nw.barrierdist.block[test.indices.day]
      test.sw.barrierdist[[day]] <- sw.barrierdist.block[test.indices.day]
      test.w.barrierdist[[day]] <- w.barrierdist.block[test.indices.day]
      test.sw.waterdist[[day]] <- sw.waterdist.block[test.indices.day]
      test.w.waterdist[[day]] <- w.waterdist.block[test.indices.day]
      test.day[[day]] <- length(test.indices.day)
    }

    ## If fewer than 100 data points available, don't train
    if(length(unlist(train.day)) < 100){
      next
    }
    
    train.da.nNA <- unlist(lapply(train.LST,length))
    test.da.nNA <- unlist(lapply(test.LST,length))

    train.dat <- data.frame(lst=unlist(train.LST),
                            elev=unlist(train.elev),
                            slope=unlist(train.slope),
                            asp=unlist(train.asp),
                            windspeed= unlist(train.windspeed),
                            nw.barrierdist= unlist(train.nw.barrierdist),
                            sw.barrierdist= unlist(train.sw.barrierdist),
                            w.barrierdist= unlist(train.w.barrierdist),
                            sw.waterdist= unlist(train.sw.waterdist),
                            w.waterdist= unlist(train.w.waterdist),
                            mod= unlist(train.MOD),
                            da=rep(mod.da[which(train.da.nNA>0)],
                              times=train.da.nNA[which(train.da.nNA>0)])
                            )
  
    test.dat <- data.frame(lst=unlist(test.LST),
                           elev=unlist(test.elev),
                           slope=unlist(test.slope),
                           asp=unlist(test.asp),
                           windspeed= unlist(test.windspeed),
                           nw.barrierdist= unlist(test.nw.barrierdist),
                           sw.barrierdist= unlist(test.sw.barrierdist),
                           w.barrierdist= unlist(test.w.barrierdist),
                           sw.waterdist= unlist(test.sw.waterdist),
                           w.waterdist= unlist(test.w.waterdist),
                           mod= unlist(test.MOD),
                           da=rep(mod.da[which(test.da.nNA>0)],
                             times=test.da.nNA[which(test.da.nNA>0)])
                           )

    rm(train.LST, train.elev, train.slope, train.asp, train.MOD, train.windspeed,
      train.nw.barrierdist, train.sw.barrierdist, train.w.barrierdist, train.sw.waterdist,
      train.w.waterdist)
    rm(test.LST, test.elev, test.slope, test.asp, test.MOD, test.windspeed,
        test.nw.barrierdist, test.sw.barrierdist, test.w.barrierdist, test.sw.waterdist,
        test.w.waterdist)
    rm(asp.block,elev.block, slope.block,windspeed.block, nw.barrierdist.block,
       sw.barrierdist.block, w.barrierdist.block, w.waterdist.block)

    if(length(which(is.na(test.dat$elev)))!=0){
      test.dat <- test.dat[-which(is.na(test.dat$elev)),]
    }

    ## Regression random forest training and testing
    print(paste("Training started:",Sys.time()))
    set.seed(55337)
    ranger.regressor <- ranger(lst ~ elev+ slope + mod + asp + windspeed +
      nw.barrierdist + sw.barrierdist + w.barrierdist + sw.waterdist + w.waterdist + da,
      data=train.dat, num.trees=100, importance="impurity")
    print(paste("Training ended:",Sys.time()))
  #  ranger.temp <- predict(ranger.regressor, data=test.dat, num.threads=128)
  #  print(paste("Prediction ended:",Sys.time()))

    ## Don't save out predicted values though!
    rm(train.dat, train.da.nNA, test.dat, test.da.nNA)

    var.importance[[counter]] <- ranger.regressor$variable.importance
    var.importance.id[[counter]] <- c(i,j)
    counter <- counter + 1

    save(ranger.regressor,file=paste0("/pl/active/rittger_esp/SierraBighorn/Rdata/localForests/regression_forest_750/forest/ranger.regression.",block_size,".",i,".",j,".RData"))
    rm(ranger.regressor)

    save(var.importance,file=paste0("/pl/active/rittger_esp/SierraBighorn/Rdata/localForests/regression_forest_750/variable.importance.regression.RData"))
    save(var.importance.id,file=paste0("/pl/active/rittger_esp/SierraBighorn/Rdata/localForests/regression_forest_750/variable.importance.regression.id.RData"))

    save(train.idx,file=paste0("/pl/active/rittger_esp/SierraBighorn/Rdata/localForests/regression_forest_750/train_indices_locations/train.indices.regression.",block_size,".",i,".",j,".RData"))
    save(test.idx,file=paste0("/pl/active/rittger_esp/SierraBighorn/Rdata/localForests/regression_forest_750/test_indices_locations/test.indices.regression.",block_size,".",i,".",j,".RData"))

    rm(train.idx,test.idx)

    print(paste0("Finished model training block: ",i," ",j," at: ",Sys.time()))
  }
}

