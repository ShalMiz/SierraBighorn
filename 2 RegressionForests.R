###############################################################################################
# Local Regression Random Forest Training Script
###############################################################################################

#___________________________________
#  Loading Libraries
#___________________________________

library(rgdal)
library(fields)
library(raster)
library(ranger) 

set.seed(999)

#___________________________________
#  Importing Dataset
#___________________________________

landsat.path <- "~/SierraBighornData/Landsat/SC/"
landsat.sc.file.names <- list.files(landsat.path)
sat.mask.path <- "~/SierraBighornData/Landsat/sat.mask/"
sat.mask.file.names <- list.files(sat.mask.path)
modis.path <- "~/SierraBighornData/MODIS/v03/"
modis.file.names <- list.files(modis.path)
predictors.path <- "~/SierraBighorn/predictors/"
predictors.file.names <- list.files(predictors.path)
Rdata_path <- "~/regression_forest/"

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

#___________________________________
#  Model training 
#___________________________________

# store variable importance for each local RF
var.importance <- array(dim = c(250,13))
model.time <- array(dim = c(250,6))
colnames(model.time) <- c("i", "j", "Block time", "Model time", "Train size", "Model size")

#extent(matrix index boundaries) of nNA data in Landsat
bottom.left.col <- 2794
bottom.left.row <- 9479
top.left.col <- 4363
top.left.row <- 3504
bottom.right.col <- 8866
bottom.right.row <- 10732

#Need k to keep track of the local RF number to update var.importance array
k <- 1
#block size for each local RF
block_size <- 450
train.pct <- 0.35

#partitioning sequence of Landsat 
##into 450*450 block and get their start row and column index
rowBeg <- seq(top.left.row, bottom.right.row, block_size)
rowBeg <- rowBeg[-length(rowBeg)] # remove last item of this vector
blockBeg <- seq(bottom.left.col, nrows, block_size)
blockBeg <- blockBeg[-length(blockBeg)] # remove last item of this vector

bl.size <- block_size-1


for (i in rowBeg) {
  
  for (j in blockBeg) {
    
    print(paste0("Starting model training block: ",i," ",j," at: ",Sys.time()))
    
    # array to save training indexes for training chosen from each day from the block
    train.idx <- array(dim = c(nday,block_size*block_size*train.pct))
    
    Block.Start.time <- Sys.time()
    
    # select predictor blocks
    elev.block <- elev[i:(i+bl.size),j:(j+bl.size)]
    slope.block <- slope[i:(i+bl.size),j:(j+bl.size)]
    asp.block <- asp[i:(i+bl.size),j:(j+bl.size)]
    windspeed.block <- windspeed[i:(i+bl.size),j:(j+bl.size)]
    nw.barrierdist.block <- nw.barrierdist[i:(i+bl.size),j:(j+bl.size)]
    sw.barrierdist.block <- sw.barrierdist[i:(i+bl.size),j:(j+bl.size)]
    w.barrierdist.block <- w.barrierdist[i:(i+bl.size),j:(j+bl.size)]
    sw.waterdist.block <- sw.waterdist[i:(i+bl.size),j:(j+bl.size)]
    w.waterdist.block <- w.waterdist[i:(i+bl.size),j:(j+bl.size)]
    
    # indexes of NAs in elevation in the selected block for later to be removed from training 
    elev.NA <- which(is.na(elev.block))
    
    # list to save training data from each day
    train.LST <- train.MOD <- train.slope <- train.asp  <- train.elev <- train.windspeed <- train.nw.barrierdist <-train.sw.barrierdist <-train.w.barrierdist <-train.sw.waterdist <-train.w.waterdist <- list()
    
    # selecting training data from each day
    for(day in 1:nday){ 
      
      LST.ras <- raster(paste0(landsat.path,landsat.sc.file.names[day]))
      LST <- t(matrix(values(LST.ras),nc= ncols, nr= nrows))
      LST.block <- LST[i:(i+bl.size),j:(j+bl.size)]
      LST.NA <- which(LST.block==255)
      rm(LST.ras, LST)
      
      sat.mask.ras <- raster(paste0(sat.mask.path,sat.mask.file.names[day]))
      sat.mask <- t(matrix(values(sat.mask.ras),nc= ncols, nr= nrows))
      sat.block <- sat.mask[i:(i+bl.size),j:(j+bl.size)]
      LST.sat <- which(sat.block==1)
      rm(sat.mask.ras, sat.mask, sat.block)
      
      # pool of indexes to sample training data from
      these.good.day <- (1:(block_size*block_size))[-c(LST.sat,LST.NA,elev.NA)]
      
      # number of samples to be taken from each day
      train.size <- length(these.good.day)*train.pct
      
      # sampling training data
      train.indices.day <- sample(these.good.day,train.size)
      
      # saving sampled indexs for future reference
      if(length(train.indices.day)!=0){
        train.idx[day,1:length(train.indices.day)] <- train.indices.day
      }
      
      # selecting training data from Landsat
      train.LST[[day]] <- LST.block[train.indices.day]
      rm(LST.block, LST.NA)
      
      # selecting training data from MODIS
      MOD.ras <- raster(paste0(modis.path, modis.file.names[day]))
      MOD <- t(matrix(values(MOD.ras),nc=ncols ,nr=nrows))
      MOD.block <- MOD[i:(i+bl.size),j:(j+bl.size)]
      train.MOD[[day]] <- MOD.block[train.indices.day]
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
      
    }
    
    train.da.nNA <- unlist(lapply(train.LST,length))
    
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
                            da=rep(mod.da,times=train.da.nNA)
    )
    
    
    rm(train.LST, train.elev, train.slope, train.asp, train.MOD, train.windspeed, train.nw.barrierdist, train.sw.barrierdist, train.w.barrierdist, train.sw.waterdist, train.w.waterdist)
    rm(asp.block,elev.block, slope.block,windspeed.block, nw.barrierdist.block, sw.barrierdist.block, w.barrierdist.block, w.waterdist.block)
    
    # creating a categorical variable as zero, between and hundred based on Landsat
    vec <- rep("btwn",dim(train.dat)[1])
    zero_indices <- which(train.dat[,1]==0)
    hundred_indices <- which(train.dat[,1]==100)
    
    vec[zero_indices] <- "zero"
    vec[hundred_indices] <-  "hundred"
    train.dat$class <- as.factor(vec)
    
    rm(vec,zero_indices,hundred_indices)
    
    # If there is less than 1000 datapoints to train on then we do not train that block
    if ( dim(train.dat[train.dat$class=="btwn",])[1]<1000 ){
      next
    }
    
    Model.Start.time <- Sys.time()
    ranger.regression <- ranger(lst~ elev+ slope + mod + asp + windspeed + nw.barrierdist+ sw.barrierdist+ w.barrierdist+ sw.waterdist+ w.waterdist + da, data=train.dat[train.dat$class=="btwn",],num.trees=100,probability=TRUE, importance = 'impurity')
    Model.End.time <- Sys.time()
    
    var.importance[k,3:13] <- ranger.regression$variable.importance
    var.importance[k,1] <- model.time[k,1] <- i
    var.importance[k,2] <- model.time[k,2] <- j
    
    model.time[k,5] <- dim(train.dat[train.dat$class=="btwn",])[1]
    model.time[k,6] <- format(object.size(ranger.regression),units="GB")
    
    save(ranger.regression,file=paste0(Rdata_path,"forest/","ranger.regression.",i,".",j,".Rda"))
    rm(ranger.regression, train.dat, train.da.nNA)
    
    Block.End.time <- Sys.time()
    model.time[k,3] <- difftime(Block.End.time, Block.Start.time, units = "mins")
    model.time[k,4] <- difftime(Model.End.time, Model.Start.time, units = "mins")
    
    k <- k+1
    
    save(var.importance, file = paste0(Rdata_path,"variable.importance.regression.RData"))
    save(model.time, file = paste0(Rdata_path,"model.time.regression.RData"))
    save(train.idx,file=paste0(Rdata_path,"train_indices_locations/","train.indices.regression.",i,".",j,".RData"))
    rm(train.idx)
    
    print(paste0("Finished model training block: ",i," ",j," at: ",Sys.time()))
    
  }
}

