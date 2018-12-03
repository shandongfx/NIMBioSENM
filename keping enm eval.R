###############################################################
library("raster")
library("dismo")

# prepare spatial occ data
if(!file.exists("data/occ_raw.rdata")){
  occ_raw <- gbif(genus="Dasypus",species="novemcinctus",download=TRUE) 
  save(occ_raw,file = "data/occ_raw.rdata")
}else{
  load("data/occ_raw.rdata")
}
occ_clean <- subset(occ_raw,(!is.na(lat))&(!is.na(lon))) 
occ_unique <- occ_clean[!duplicated( occ_clean[c("lat","lon")]  ),]
occ_unique_specimen <- subset(occ_unique, basisOfRecord=="PRESERVED_SPECIMEN")
occ_final <- subset(occ_unique_specimen, year>=1950 & year <=2000)
coordinates(occ_final) <- ~ lon + lat
myCRS1 <- CRS("+init=epsg:4326") # WGS 84
crs(occ_final) <- myCRS1

# prepare raster data
if( !file.exists( paste0("data/bioclim/bio_10m_bil.zip")   )){
  utils::download.file(url="http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_10m_bil.zip",
                       destfile="data/bioclim/bio_10m_bil.zip"   ) 
  utils::unzip("data/bioclim/bio_10m_bil.zip",exdir="data/bioclim") 
}

# load rasters
clim_list <- list.files("data/bioclim/",pattern=".bil$",full.names = T)
clim <- raster::stack(clim_list) 

occ_buffer <- buffer(occ_final,width=4*10^5) #unit is meter
clim_mask <- mask(clim, occ_buffer)

set.seed(1) 
bg <- sampleRandom(x=clim_mask,
                   size=10000,
                   na.rm=T, #removes the 'Not Applicable' points  
                   sp=T) # return spatial points 

temp1 <- extract(clim_mask[[1]],occ_final)
occ_final <- occ_final[!is.na(temp1),]
#############################################################################################################################################################################################
###########################################################
#### Run analyses for My data study ####
###########################################################

library(ENMeval)


###loading occurrence and climate data, make sure that occurrence and background data only have two columns (i.e., long and lat), these are present-only data; For background data, you can use R or GIS software to generate these background data(e.g. random sampling 5000,10000),these background data should be within the accessible area (i.e. M in BAM diagram) of species in question###.


env <- clim_mask[[c("bio1","bio5","bio6","bio12")]]
#env <- stack("bio1.tif","bio5.tif","bio6.tif")

occ <- occ_final@coords
#occ <- read.csv("occ.csv", head=TRUE)

bg <- bg@coords
#bg<- read.csv("bc.csv", head=TRUE)


######start ENMeval; "RMvalues" is used to set a range of RM values, here we set RM ranged from 0.5 to 4 at at the interval of 0.5; "method" is used to spatial parting occurrence data, there are mainly two approaches available in ENM eval, i.e. block and checkerboard methods, the former method is used when your model in a transferred manner/need to be transferred (i.e. in the application of biological invasions, climate change), the latter is used in a none transfer manner (i.e. setting priority area for conservation); "fc" is used to set feature combination (for example, fc = c('L', 'LQ', 'H')), here we use default, which has six feature combinations; "overlap" is asking whether you are going to perform overlap measurements of Maxent prediction during the iterative running; "bin.output" is asking whether you are going to reserved the iterative prediction; "categoricals" is aksing wether you are allow categorical variables to fit model##### 

#   fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT")
coc.results <- ENMevaluate(occ = occ,  # add comments...............
                           env = env, 
                           bg.coords = bg, 
                           RMvalues=seq(0.5,4,0.5),
                           fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                           method='block', 
                           overlap=FALSE, 
                           bin.output=FALSE)



##################################################
#### Exploration the results ####
##################################################

# View ENMevaluation:
res <- coc.results


# Look at results table AND save it in working directory for later checking.
res@results
write.csv (res@results, file = "MyResult.csv")

# Which settings gave delta.AICc < 2?
aicmods <- which(res@results$delta.AICc < 2)
res@results[aicmods,]

# Visualize how data were partitioned
#1 Background points (error !!!!):
r <- raster("bio1.tif")
plot(r, xlim=c(-180,-20))
points(res@bg.pts, col= res@bg.grp, cex=.75)

#2 Occurrence localities:
plot(res@predictions[[which(res@results$delta.AICc == 0)]], xlim=c(-180,-20))
points(res@occ.pts, pch=16, col= res@occ.grp, cex=.75)
dev.print(tiff, "OCC_partition.tiff", res=600, units="in")


# View predictions in geographic space for these models
plot(res@predictions[[aicmods]])

# Plot delta.AICc for different settings
par(mfrow=c(2,2))
eval.plot(res@results)
eval.plot(res@results,      'avg.test.AUC')
eval.plot(res@results,      'avg.diff.AUC')
eval.plot(res@results,   'avg.test.or10pct')
eval.plot(res@resultDifs, 'avg.test.orMTP')



plot(res@results$avg.test.AUC, 
     res@results$delta.AICc, 
     bg=res@results$features, pch=21, cex= res@results$rm/2, ylim=c(0,30))
legend("topright", legend=unique(res@results$features), pt.bg=res@results$features, pch=21)
mtext("Circle size proportional to regularization multiplier value")


##################################################
#### Home work: use "block" and "checkerboard" methods to spatial parting occurrence records in ENMeval to find the best Maxent setting, and run Maxent model to compare their predictions based these setting, here are the codes####
##################################################


##preparing data###

env <- stack("bio1.tif","bio5.tif","bio6.tif")

occ <- read.csv("occ.csv", head=TRUE)

bg<- read.csv("bc.csv", head=TRUE)


###run ENMeval using two method to spatial parting method to 
coc.results <- ENMevaluate(occ, env, bg, RMvalues=seq(0.5,4,0.5),method='block', overlap=FALSE, bin.output=FALSE, categoricals=1)

coc.results1 <- ENMevaluate(occ, env, bg, RMvalues=seq(0.5,4,0.5),method='checkerboard1', overlap=FALSE, bin.output=FALSE, categoricals=1)



###Writing the results of block method#####

res <- coc.results
res@results
write.csv (res@results, file = "MyResult.csv")

###Writing the results of checkerboard method#####
res <- coc.results1
res@results1
write.csv (res@results1, file = "MyResult.csv")



# Selecting settings gave delta.AICc < 2 in block method####
aicmods <- which(res@results$delta.AICc < 2)
res@results[aicmods,]

# Selecting settings gave delta.AICc < 2 in checkerboard####
aicmods <- which(res@results$delta.AICc < 2)
res@results[aicmods,]
