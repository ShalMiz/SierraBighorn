## Combining Shalini's classification with Will's regression locally fitted models from 2bRegressionForests.R/4bRegressionDownscale.R

library(rgdal)
library(fields)
library(raster)
library(ranger)

#local classification downscaled
local.classification.path <- "/pl/active/rittger_esp/SierraBighorn/downscaledv5_localForests/downscaled_classification/downscaled/"
local.classification.names <- list.files(local.classification.path)

#local regression downscaled # NOTE: here's the difference
local.regression.path <- "/pl/active/rittger_esp/SierraBighorn/downscaledv5_localForests/downscaled_regression_750/downscaled/"
local.regression.names <- list.files(local.regression.path)

#local downscales combined 
local.combined.path <- "/pl/active/rittger_esp/SierraBighorn/downscaledv5_localForests/downscaled_combined_750/downscaled/"

#getting date
lst.date <- NULL
for(i in 1:length(local.classification.names)){
  lst.date[i] <- regmatches(local.classification.names[i], regexpr("\\d{4}\\d{2}\\d{2}",local.classification.names[i]))
}

lst.date <- as.integer(lst.date)
lst.date <- strptime(x=as.character(lst.date),format="%Y%m%d")


ncols <- 14752
nrows <- 9712

n <- length(local.regression.names)

for (day in 1:n) {
  
  local.classfication <- raster(paste0(local.classification.path,local.classification.names[day]))
  local.classfication.mat <- t(matrix(values(local.classfication),nc=ncols,nr=nrows))
  local.classfication.mat[which(local.classfication.mat==2)] <- 100
  idx <- which(local.classfication.mat==1)
  
  local.regression <- raster(paste0(local.regression.path,local.regression.names[day]))
  local.regression.mat <- t(matrix(values(local.regression),nc=ncols,nr=nrows))
  
  local.classfication.mat[idx] <- local.regression.mat[idx]
  
  values(local.classfication) <- c(t(matrix(local.classfication.mat,nr=ncols,nc=nrows)))
  
  writeRaster(local.classfication, filename=paste0(local.combined.path,"downscaled.combined.750", format(lst.date[day],"%Y%m%d"),".tif"), format="GTiff",option="COMPRESS=LZW",datatype="INT1U",overwrite=TRUE)
  print(paste0("Done with day: ",day))
  rm(local.classfication, local.classfication.mat, local.regression, local.regression.mat)
}
