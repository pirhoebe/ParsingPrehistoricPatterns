###############################################################################
# Supplementary Information for Hoebe Peeters & Arnoldussen 2022 
###############################################################################
# 
###############################################################################
#### 5 Analysis #### 
###############################################################################

# install and read relevant libraries
source(file="read libraries.R")

# run preparation script if necessary
if(exists("data_full")==FALSE | exists("data_full.cal")==FALSE){
  source(file="1 Dataset preparation.R")}

# Prepare bins
h_value <- 100
if(!"bins" %in% colnames(data_full)){
data_full$bins = binPrep(sites=data_full$site_id,
                         ages=data_full$c14age, 
                         h=h_value)
}

###############################################################################
#### 5.1 Prepare SPDs ####
###############################################################################

### Full dataset SPD aggregation ###
# parameters
runningmean <- 50
lowerlim <- 16000
upperlim <- 7500
normalised <- FALSE
if(exists("data_full.spd")==FALSE){
  bins = data_full$bins
  set.cal <- data_full.cal
  set.spd= spd(set.cal,timeRange=c(lowerlim,upperlim), bins = bins)
  data_full.spd <- set.spd}


### SPD aggregation with 200 year margins for model fitting ###
#parameters
if(exists("data_full.spd_m")==FALSE){
  runningmean <- 50
  lowerlim <- 16200
  upperlim <- 7300
  normalised <- FALSE
  bins = data_full$bins
  set.cal <- data_full.cal
  set.spd= spd(set.cal,timeRange=c(lowerlim,upperlim), bins = bins)
  data_full.spd_m <- set.spd
  lowerlim <- 16000
  upperlim <- 7500
  }


###############################################################################
#### 5.2 Prepare modeltested SPDs ####
###############################################################################

simulations <- 500
### exponential modeltest ###
if(exists("data_full.expmod")==FALSE){
set.cal <- data_full.cal
bins <- data_full$bins
myModel <- "exponential" #options are: "exponential", "linear", "uniform"
source("Modeltest.R")
data_full.expmod <- set.mod
}

### custom climate modeltest ###
simulations <- 500 
set <- data_full
set.cal <- data_full.cal
set.spd <- data_full.spd_m
bins <- data_full$bins
# open and run this script seperately for SI-7 output
source("Modeltest climate.R")
data_full.climod <- set.climod


###############################################################################
#### 5.3 Plotting ####
###############################################################################

###############################################################################
#### 5.3.1 Exponential vs climate modeltest ####
###############################################################################
par(mfrow=c(2,1))
par(oma=c(4,3,3,2))
###############################################################################
## a. Exponential modeltest ##
par(mar=c(1.5,0,0,0))
plot(data_full.expmod,
     xlim=c(lowerlim, upperlim),
     ylim=c(0, 1.3))
par(xpd=NA)
phaserect <- 5
phaselabels <- -0.2
source("phaserect.R")
par(xpd=FALSE)
phaserect <- 6
rm(phaselabels)
source("phaserect.R")
text(x=16000, y=(par("usr")[4]-(0.15*(par("usr")[4]))), 
     labels=c("a. Exponential modeltest"), 
     col="black", cex=1,font=2, pos=4)
axis(3, at = seq(8000,16000, by=1000), labels = FALSE)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

axis(1, at = seq(8000,16000, by=1000), labels = FALSE)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

###############################################################################
## b. climate modeltest ##
par(mar=c(0,0,1.5,0))
plot(data_full.climod,
     xlim=c(lowerlim, upperlim),
     ylim=c(0, 1.3))
text(x=16000, y=(par("usr")[4]-(0.15*(par("usr")[4]))), 
     labels=c("b. Climate modeltest"), 
     col="black", cex=1,font=2, pos=4)
par(xpd=NA)
phaserect <- 5
source("phaserect.R")
phaserect=6
source("phaserect.R")
axis(3, at = seq(8000,16000, by=1000), labels = FALSE)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

axis(1, at = seq(8000,16000, by=1000), labels = FALSE)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)


###############################################################################
#### 5.3.2 Permtest by landscape zone ####
###############################################################################

# categorise dataset
data_full$landscape <- NA_character_
data_full$landscape[data_full$elevation<=10] <- "Low"
data_full$landscape[data_full$elevation<=50 &
                      data_full$elevation>10] <- "Mid"
data_full$landscape[data_full$elevation>50 ] <- "High"

data_full %>% 
  group_by(landscape) %>% 
  summarise(n = n())

elevation10_50 <-data_full
elevation10_50 <- elevation10_50%>%filter(!is.na(landscape))
elevation10_50 %>% 
  group_by(landscape) %>% 
  summarise(n = n())

# recalibrate set with categories
set <- elevation10_50
source("calibrate data.R")
elevation10_50.cal <- set.cal
elevation10_50.cal$metadata$landscape <- elevation10_50$landscape
simulations=500

# permutation test by landscape zone 
perm.elevation10_50=permTest(x=elevation10_50.cal,
                             marks=elevation10_50$landscape,
                             bins = elevation10_50$bins,
                             timeRange=c(lowerlim, upperlim),
                             nsim=simulations,
                             runm=runningmean)

data <- filter(data_full,
               data_full$calBP_end <= lowerlim &
                 data_full$calBP_start >= upperlim)

# calculate right number of dates within lowlands and uplands

ndates$e <- list(
  "LL" = nrow(filter(elevation10_50,
                     elevation10_50$calBP_end<=lowerlim &
                       elevation10_50$calBP_start>=upperlim &
                       elevation10_50$landscape=="Low")),
  "ML" = nrow(filter(elevation10_50,
                      elevation10_50$calBP_end<=lowerlim &
                        elevation10_50$calBP_start>=upperlim &
                        elevation10_50$landscape=="Mid")),
  "HL" = nrow(filter(elevation10_50,
                      elevation10_50$calBP_end<=lowerlim &
                        elevation10_50$calBP_start>=upperlim &
                        elevation10_50$landscape=="High")))

## Plot landscape zone permutation test ##
###############################################################################
par(mfrow=c(3,1))
par(oma=c(2,3,2,2))
###############################################################################
## a. Low ## 
par(mar=c(1.5,0,1.5,0))
plot(perm.elevation10_50,focalm = "Low",
     ylim=c(0,0.5))
text(x=16000, y=(par("usr")[4]-(0.15*(par("usr")[4]))), 
     labels=paste("a. Low: 10 m and lower (n=",
                  ndates$e$LL, ")", sep = ""), 
     col="black", cex=1.5, pos=4, font=2)
axis(3, at = seq(8000,16000, by=1000), labels = FALSE)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

axis(1, at = seq(8000,16000, by=1000), labels = FALSE)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)
par(xpd=NA)
phaserect <- 5
phaselabels <- -0.2
source("phaserect.R")
par(xpd=FALSE)
phaserect <- 6
rm(phaselabels)
source("phaserect.R")

###############################################################################
## b. Mid ##
par(mar=c(1.5,0,1.5,0))
plot(perm.elevation10_50,focalm = "Mid",
     ylim=c(0,0.5))
text(x=16000, y=(par("usr")[4]-(0.15*(par("usr")[4]))), 
     labels=paste("b. Middle: between 10 and 50 (n=",
                  ndates$e$ML, ")", sep = ""), 
     col="black", cex=1.5, pos=4, font=2)
axis(3, at = seq(8000,16000, by=1000), labels = FALSE)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

axis(1, at = seq(8000,16000, by=1000), labels = FALSE)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)
par(xpd=NA)
phaserect <- 5
source("phaserect.R")
par(xpd=FALSE)
phaserect <- 6
rm(phaselabels)
source("phaserect.R")

###############################################################################
## c. High ##
par(mar=c(1.5,0,1.5,0))
plot(perm.elevation10_50,focalm = "High",
     ylim=c(0,0.5))
text(x=16000, y=(par("usr")[4]-(0.15*(par("usr")[4]))), 
     labels=paste("c. High: above 50 m (n=",
                  ndates$e$HL, ")", sep = ""), 
     col="black", cex=1.5, pos=4, font=2)
axis(3, at = seq(8000,16000, by=1000), labels = FALSE)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

axis(1, at = seq(8000,16000, by=1000), labels = FALSE)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)
par(xpd=NA)
phaserect <- 5
source("phaserect.R")
par(xpd=FALSE)
phaserect <- 6
rm(phaselabels)
source("phaserect.R")


###############################################################################
#### 5.3.3 Permtest by country ####
###############################################################################
# Preparation
simulations <- 500
runningmean <- 50

# permutation test by country
perm.country=permTest(x=data_full.cal,
                          marks=data_full$country,
                          timeRange=c(lowerlim, upperlim),
                          bins=data_full$bins,
                          nsim=simulations,
                          runm=runningmean)
summary(perm.country)


###############################################################################
## Plot country permutation test ##
par(mfrow=c(2,3))
par(oma=c(2,2,2,2))

## SPD for northwest Europe ##
settitle <- "Northwest Europe"
rm(phaselabels)
set.spd <- data_full.spd
set.ndates=ndates$full
par(mar=c(2,2,2,2))
plot(set.spd, main = paste( 
  "Full dataset NW Europe (n dates = ", set.ndates, ")", sep=""))
rm(phaselabels)
source("phaserect.R")
plot(set.spd, runm=runningmean, add=TRUE, type="simple", 
     col="black", lwd=1.5, lty=2)

## Netherlands ##
plot(perm.country,focalm = "Netherlands",
     main=paste( "Netherlands (n=", ndates$c$NL, ")", sep = ""))
source("phaserect.R")

## Denmark ##
plot(perm.country,focalm = "Denmark",
     main=paste("Denmark (n=", ndates$c$DK, ")", sep = ""))
source("phaserect.R")

## Britain ##
plot(perm.country,focalm = "United Kingdom",
     main=paste( "Britain (n=", ndates$c$UK, ")", sep = ""))
source("phaserect.R")

## Belgium ##
plot(perm.country,focalm = "Belgium",
     main=paste("Belgium (n=", ndates$c$BE, ")", sep = ""))
source("phaserect.R")

## Germany ##
plot(perm.country,focalm = "Germany",
     main=paste("Germany (n=", ndates$c$DE, ")", sep = ""))
source("phaserect.R")


###############################################################################
#### 5.3.4 Material category permtest ####
###############################################################################
# Preparation
# simplify material categories
data_full$MatCat <- NA
data_full$MatCat[data_full$material == "charcoal"] <- "charcoal"

animal <- c("antler","bone", "bone / antler", "bone collagen", "tooth / ivory", 
            "amino acid")
for (CatNr in 1:length(animal)) {
  xMat <- animal[CatNr]
  message(paste( "Material:", CatNr, xMat  ))
  data_full$MatCat[
    str_detect(data_full$material, 
               fixed(xMat, ignore_case = TRUE))] <- "animal remains"
}

plant <- c("botanical remains", "wood", "adhesives")
for (CatNr in 1:length(plant)) {
  xMat <- plant[CatNr]
  message(paste( "Material:", CatNr, xMat  ))
  data_full$MatCat[
    str_detect(data_full$material, 
    fixed(xMat, ignore_case = TRUE))] <- "plant remains"
}

other <- c("other", "bulk", "pottery", "same sample", "sediment", "shell")
for (CatNr in 1:length(other)) {
  xMat <- other[CatNr]
  message(paste( "Material:", CatNr, xMat  ))
  data_full$MatCat[
    str_detect(data_full$material, 
    fixed(xMat, ignore_case = TRUE))] <- "other"
}

# determine number of dates for each category
ndates$m <- list(
"animal"= nrow(filter(data_full,
                      data_full$calBP_end<=lowerlim &
                      data_full$calBP_start>=upperlim &
                      data_full$MatCat=="animal remains")),

"plant"= nrow(filter(data_full,
                     data_full$calBP_end<=lowerlim &
                     data_full$calBP_start>=upperlim &
                     data_full$MatCat=="plant remains")),

"charcoal"= nrow(filter(data_full,
                        data_full$calBP_end<=lowerlim &
                        data_full$calBP_start>=upperlim &
                        data_full$MatCat=="charcoal")),

"other"= nrow(filter(data_full,
                     data_full$calBP_end<=lowerlim &
                     data_full$calBP_start>=upperlim &
                     data_full$MatCat=="other"))
)

# recalibrate to include these categories
set <- data_full
source("calibrate data.R")
set.cal$metadata$MatCat <- set$MatCat
data_full_mat.cal <- set.cal

# parameters
simulations <- 500
runningmean <- 50

# permutation test by material
perm.material=permTest(x=data_full_mat.cal,marks=data_full$MatCat,
                           timeRange=c(lowerlim, upperlim),nsim=simulations,
                           runm=runningmean)

###############################################################################
## Plot material permutation test ##
par(mfrow=c(4,1))
par(oma=c(4,3,3,2))
par(mar=c(1.5,0,1,0))
## charcoal ##
plot(perm.material,focalm = "charcoal",
     ylim=c(0,1))
text(x=16000, y=(par("usr")[4]-(0.15*(par("usr")[4]))), 
     labels=paste("a. charcoal (n=", ndates$m$charcoal, ")", sep = ""), 
     col="black", cex=1.5, font=2, pos=4)
axis(3, at = seq(8000,16000, by=1000), labels = FALSE)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

axis(1, at = seq(8000,16000, by=1000), labels = FALSE)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)
phaserect <- 5
phaselabels <- -0.2
source("phaserect.R")
rm(phaselabels)
phaserect <- 6
source("phaserect.R")

## animal remains ##
par(mar=c(1.5,0,1.5,0))
plot(perm.material,focalm = "animal remains",
     ylim=c(0,0.6))
text(x=16000, y=(par("usr")[4]-(0.15*(par("usr")[4]))), 
     labels=paste("b. animal remains (n=", ndates$m$animal, ")", sep = ""), 
     col="black", cex=1.5, font=2, pos=4)
axis(3, at = seq(8000,16000, by=1000), labels = FALSE)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

axis(1, at = seq(8000,16000, by=1000), labels = FALSE)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

## plant remains ##
source("phaserect.R")
par(mar=c(1.5,0,1.5,0))
plot(perm.material,focalm = "plant remains",
     ylim=c(0,0.6))
text(x=16000, y=(par("usr")[4]-(0.15*(par("usr")[4]))), 
     labels=paste("c. plant remains (n=", ndates$m$plant, ")", sep = ""), 
     col="black", cex=1.5, font=2, pos=4)
axis(3, at = seq(8000,16000, by=1000), labels = FALSE)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

axis(1, at = seq(8000,16000, by=1000), labels = FALSE)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)
source("phaserect.R")

# other remains ##
par(mar=c(0,0,1.5,0))
plot(perm.material,focalm = "other",
     ylim=c(0,0.5))
text(x=16000, y=(par("usr")[4]-(0.15*(par("usr")[4]))), 
     labels=paste("d. other (n=", ndates$m$other, ")", sep = ""), 
     col="black", cex=1.5, font=2, pos=4)
axis(3, at = seq(8000,16000, by=1000), labels = FALSE)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

axis(1, at = seq(8000,16000, by=1000), labels = FALSE)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)
source("phaserect.R")


###############################################################################
save.image(file='14cEnvironment.RData')
message("Environment saved")

data_tf <- filter(data_full, data_full$calBP_end<=lowerlim &
                            data_full$calBP_start>=upperlim)
path <- "- Supplementary data post analysis.xlsx"
write.xlsx(data_tf, file=path)

file.edit("6 Calculations in text.R") 
