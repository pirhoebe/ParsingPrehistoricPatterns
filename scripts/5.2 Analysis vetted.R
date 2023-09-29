###############################################################################
# Supplementary Information for Hoebe Peeters & Arnoldussen 2023b 
###############################################################################
# 
###############################################################################
#### 5.2 Analysis vetted #### 
###############################################################################

# install and read relevant libraries
source(file="read libraries.R")

# run preparation script if necessary
if(exists("data_vet")==FALSE | exists("data_vet.cal")==FALSE){
  source(file="1 Dataset preparation.R")}

# Prepare bins
h_value <- 100
if(!"bins" %in% colnames(data_vet)){
  data_vet$bins = binPrep(sites=data_vet$site_id,
                           ages=data_vet$c14age, 
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
if(exists("data_vet.spd")==FALSE){
  bins = data_vet$bins
  set.cal <- data_vet.cal
  set.spd= spd(set.cal,timeRange=c(lowerlim,upperlim), bins = bins)
  data_vet.spd <- set.spd}


### SPD aggregation with 200 year margins for model fitting ###
#parameters
if(exists("data_vet.spd_m")==FALSE){
  runningmean <- 50
  lowerlim <- 16200
  upperlim <- 7300
  normalised <- FALSE
  bins = data_vet$bins
  set.cal <- data_vet.cal
  set.spd= spd(set.cal,timeRange=c(lowerlim,upperlim), bins = bins)
  data_vet.spd_m <- set.spd
  lowerlim <- 16000
  upperlim <- 7500
}


###############################################################################
#### 5.2 Prepare modeltested SPDs ####
###############################################################################

simulations <- 500
### exponential modeltest ###
if(exists("data_vet.expmod")==FALSE){
  set.cal <- data_vet.cal
  bins <- data_vet$bins
  myModel <- "exponential" #options are: "exponential", "linear", "uniform"
  source("Modeltest.R")
  data_vet.expmod <- set.mod
}

### custom climate modeltest ###
simulations <- 500 
set <- data_vet
set.cal <- data_vet.cal
set.spd <- data_vet.spd_m
bins <- data_vet$bins
# open and run this script seperately for SI-7 output
source("Modeltest climate.R")
data_vet.climod <- set.climod


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
plot(data_vet.expmod,
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
plot(data_vet.climod,
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
save.image(file='14cEnvironment.RData')
message("Environment saved")

data_tf <- filter(data_vet, data_vet$calBP_end<=lowerlim &
                    data_vet$calBP_start>=upperlim)
path <- "- Supplementary data post analysis.xlsx"
write.xlsx(data_tf, file=path)

