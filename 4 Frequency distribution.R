###############################################################################
# Supplementary Information for Hoebe Peeters & Arnoldussen 2023 
###############################################################################
#
###############################################################################
#### 4 Frequency distribution ####
###############################################################################

# install and read relevant libraries
source(file="read libraries.R")

# run preparation script if necessary
if(exists("data_full")==FALSE | exists("data_full.cal")==FALSE){
source(file="1 Dataset preparation.R")}


###############################################################################
#### 4.1 Binsense ####
###############################################################################

# Unfortunately rcarbon does not seem to support adding binsense to other plots
# save plot to add the plot with other software. Figure margins here work well 
# if the multiplot below is saved at 1000 x 1000.

png(file="binsense.png",
    width=1230, height=300)

binsense(x=set.cal,
         y=set.cal$metadata$site_id,
         h=c(50,100,200,500),
         timeRange=c(lowerlim, upperlim),
         calendar = "BP")
dev.off()

h_value <- 100 # h of 100 assigns dates from the same site that are within 100y
# of each other into the same bin 

#full dataset
data_full$bins = binPrep(sites=data_full$site_id,
                         ages=data_full$c14age, 
                         h=h_value)
data_full.bins.med=binMed(x = data_full.cal,bins=data_full$bins)
#vetted dataset
data_vet$bins = binPrep(sites=data_vet$site_id,
                        ages=data_vet$c14age,
                        h=h_value)


###############################################################################
#### 4.2 Prepare SPDs and KDE ####
###############################################################################

# PARAMETERS 
runningmean <- 50
lowerlim <- 16000
upperlim <- 7500
normalised <- FALSE

### Full dataset SPD aggregation ###
if(exists("data_full.spd")==FALSE){
  bins = data_full$bins
  set.cal <- data_full.cal
  set.spd= spd(set.cal,timeRange=c(lowerlim,upperlim), bins = bins)
  data_full.spd <- set.spd
  }

### Vetted dataset SPD aggregation ###
if(exists("data_vet.spd")==FALSE){
  bins = data_vet$bins
  set.cal <- data_vet.cal
  set.spd= spd(set.cal,timeRange=c(lowerlim,upperlim), bins = bins)
  data_vet.spd <- set.spd
  }

### Full dataset SPD normalised ###
if(exists("data_full_n.cal")==FALSE){
normalised <- TRUE
set <- data_full
source("calibrate data.R")
data_full_n.cal <- set.cal
normalised <- FALSE}

if(exists("data_full_n.spd")==FALSE){
bins = data_full$bins
set.cal <- data_full_n.cal
set.spd= spd(set.cal,timeRange=c(lowerlim,upperlim), bins = bins)
data_full_n.spd <- set.spd
}

### Full dataset KDE ###
# PARAMETERS 
simulations=500

# 50
bandwidth=50
bins = data_full$bins
set.cal <- data_full.cal
source("KDE analysis.R")
data_full.kde50 <- set.ckde

# 100
bandwidth=100
bins = data_full$bins
set.cal <- data_full.cal
source("KDE analysis.R")
data_full.kde100 <- set.ckde

# 200
bandwidth=200
bins = data_full$bins
set.cal <- data_full.cal
source("KDE analysis.R")
data_full.kde200 <- set.ckde

###############################################################################
#### 4.4 prepare other datasets ####
###############################################################################

### load intcal20 data ####
intcal20 <- read.csv('http://intcal.org/curves/intcal20.14c',
                     encoding="UTF-8",skip=11,header=F)
colnames(intcal20) <- c("BP","CRA","Error","D14C","Sigma")

### load NGRIP data ###
if(!require('pangaear')) install.packages('pangaear')
if(!"pangaear" %in% (.packages())){library(pangaear)} else 
  {message(paste("pangear is attached"))}
NGRIP<- pg_data(doi='10.1594/PANGAEA.824889')
NGRIPdata <- NGRIP[[1]]$data 
NGRIPdata <- rename(NGRIPdata, age = 2)
NGRIPdata <- rename(NGRIPdata, d18O = 3)
NGRIPdata$calBP <-NGRIPdata$age*1000 
NGRIPdata$calBP <- round(NGRIPdata$calBP, 0)
NGRIPdata <- NGRIPdata %>% filter(calBP>=7000 & calBP<=16500)
NGRIPdata$age <- NULL

### climate phases ###
# prepare cold climate phases as following Rasmussen et al. 2014 
# and 10.3ka (Bond et al 2001)
event <- data.frame(name=c("GS-2", "GI\n1d", "GI\n1c2","GI\n1b","GS-1", 
                           "11.4ka\nevent", "10.3ka\nevent", 
                           "9.3ka\nevent", "8.2ka\nevent"), 
                    start=c(16000, 14025, 13610, 13261, 12846, 
                            11470, 10375,  9300, 8250), 
                    end=c(14642, 13904, 13550, 13049, 11653, 11350, 
                          10225,  9190, 8090))
event$centre <- event$start-((event$start-event$end)/2)
chronozones$col <- c(rgb(0,0,0,alpha=0.2), "plum", "palegreen3",
                     rgb(0,0,0,alpha=0.2), "palegreen3", "palegreen4", 
                     "palegreen4", "green4")


###############################################################################
#### 4.5 Plotting ####
###############################################################################

###############################################################################
### Figure 3: NGRIP, Binsense, SPD, KDE ###

par(mfrow=c(4,1))
par(yaxs="i")
par(xaxs="i")
par(oma=c(4,0,4,0))
###############################################################################
## a. NGRIP plot with climate chronozones ##
par(mar=c(0,6,0,2)) # margins b, l, t, r
plot(x=NGRIPdata$calBP[NGRIPdata$calBP<=16000 &
                                      NGRIPdata$calBP>=7500],
                  y=NGRIPdata$d18O[NGRIPdata$calBP<=16000 &
                                     NGRIPdata$calBP>=7500],
                  xlim=rev(range(NGRIPdata$calBP[NGRIPdata$calBP<=16000 &
                                                   NGRIPdata$calBP>=7500])),
                  xaxt="n",
                  type="l",
                  ylab = "d18O",
                  ylim=c(min(NGRIPdata$d18O)-1, max(NGRIPdata$d18O)+4),
                  col="firebrick",lwd=2)
# chronozones
rect(x=chronozones$llim,
     xright = chronozones$ulim,
     ybottom = max(NGRIPdata$d18O)+1,
     ytop = max(NGRIPdata$d18O)+4,
     col=chronozones$col, # add chronozone colour to chronozones
     border="black")
text(x=chronozones$llim-((chronozones$llim-chronozones$ulim)/2), 
     y=max(NGRIPdata$d18O)+2.5, 
     labels=chronozones$name, col="black", cex=1.5)
# climate event phases
rect(x= event$start,
     xright = event$end,
     ybottom = min(NGRIPdata$d18O)-1,
     ytop = max(NGRIPdata$d18O)+1,
     col= rgb(0,0,0,alpha=0.2),
     border=NA)
phaselabels=0.2
text(x=event$centre, y=c(-38,-42.5,-42.5,-42.5,-38,-42.5,-42.5,-42.5,-42.5), 
labels=event$name, col="black", cex=1.5)
text(x= 16000, y=-34,
     labels="a. NGRIP, chronozones and events", 
     col="black", cex=1.5, font=2, pos=4)
# ticks on top axis
axis(3, at = seq(8000,16000, by=1000))
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)
# ticks on middle axis
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.01)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=0.01)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.02)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=0.02)
axis(1, at = seq(8000,16000, by=1000), labels=FALSE, tck=0.025)
axis(1, at = seq(8000,16000, by=1000), labels=FALSE, tck=-0.025)

###############################################################################
## b. empty plot for binsense plot ##
par(mar=c(0,6,0,2))
multi.plot <-plot(data_full.spd,
                 type="simple",
                 col=NA,
                 xaxt="n",
                 yaxt = "n",
                 lwd=0,lty=1)
rm(phaselabels)
phaserect=1
source("phaserect.R")
# visualise barcodes
barCodes(data_full.bins.med, yrng = c(0,0.07))
text(x= 16000, y=par("usr")[4]-(0.1*par("usr")[4]),
     labels="b. Binsense", col="black", cex=1.5, font=2, pos=4)

###############################################################################
## c. comparison of normal SPD and SPD with redundant normalisation ##
par(mar=c(0,6,0,2))
multi.plot <-plot(data_full_n.spd,
                 type= "simple",
                 col=NA,
                 xaxt="n",
                 lwd=2,lty=1)
rm(phaselabels)
phaserect=1
source("phaserect.R")
multi.plot <-plot(data_full_n.spd,
                 fill="tan1",
                 add=TRUE,
                 xaxt="n",
                 lwd=1.5,lty=1)

multi.plot <-plot(data_full.spd, 
                 add=TRUE,
                 type="simple",
                 col="black",
                 xlim=c(16000,7500),
                 xaxt="n")
# add Intcal20 to contextualise artificial spiking
lines(intcal20[intcal20$BP<=16000 & intcal20$BP>=7500,"BP"],
                   reScale(intcal20[intcal20$BP<=16000 & 
                                      intcal20$BP>=7500,"CRA"])*par('usr')[4],
                   col="royalblue",lwd=2)
text(x= 16000, y=par("usr")[4]-(0.1*par("usr")[4]),
     labels="c. SPD", col="black", cex=1.5, font=2, pos=4)
# legend
legend(x= 16000, 
       y=par("usr")[4]-(0.2*par("usr")[4]), 
       bty="n", 
       col=c("tan1", "black","royalblue"),
       legend=c("normalised", "non-normalised", "Intcal20"),lwd=c(6,2,2), cex=1)
# ticks on middle axis top
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.01)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=0.01)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=0.02)
axis(3, at = seq(8000,16000, by=1000), labels=FALSE, tck=0.025)
axis(3, at = seq(8000,16000, by=1000), labels=FALSE, tck=-0.025)
# ticks on middle axis bottom
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.01)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=0.01)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.015)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=0.015)
axis(1, at = seq(8000,16000, by=1000), labels=FALSE, tck=0.025)
axis(1, at = seq(8000,16000, by=1000), labels=FALSE, tck=-0.025)

###############################################################################
## d. KDE plot ##
par(mar=c(0,6,0,2))
plot(data_full.kde50,
     xlim=c(16000,7500))
source("phaserect.R")
axis(1, at = seq(9000,15000, by=2000))
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)
text(x= 16000, y=par("usr")[4]-(0.1*par("usr")[4]),
     labels="d. KDE", col="black", cex=1.5, font=2, pos=4)
legend(x=16000,y=par("usr")[4]-(0.2*par("usr")[4]),
       bty="n", col=c("darkgrey", "black"), legend = c("envelope", "median"), 
       lwd=c(6,2), lty=c(1, 3), cex=1)
rm(l, i, phaserect)

###############################################################################
### SI-6 KDE 50, 100, 200 y bw ###

par(mfrow=c(3,1))
par(yaxs="i")
par(xaxs="i")
par(oma=c(4,0,4,0))
###############################################################################
par(mar=c(1.5,6,1.5,2))
plot(data_full.kde50,
     xaxt="n",
     xlim=c(16000,7500))
source("phaserect.R")
text(x= 16000, y=par("usr")[4]-(0.1*par("usr")[4]),
     labels="a. 50 y bandwidth", col="black", cex=1.5, font=2, pos=4)
axis(1, at = seq(9000,15000, by=2000))
# ticks on top axis
axis(3, at = seq(8000,16000, by=1000))
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)
# ticks on middle axis
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.01)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=0.01)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.02)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=0.02)
axis(1, at = seq(8000,16000, by=1000), labels=FALSE, tck=0.025)
axis(1, at = seq(8000,16000, by=1000), labels=FALSE, tck=-0.025)
###############################################################################
par(mar=c(1.5,6,1.5,2))
plot(data_full.kde100,
     xlim=c(16000,7500))
source("phaserect.R")
text(x= 16000, y=par("usr")[4]-(0.1*par("usr")[4]),
     labels="b. 100 y bandwidth", col="black", cex=1.5, font=2, pos=4)
axis(1, at = seq(9000,15000, by=2000))
# ticks on middle axis top
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.01)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=0.01)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=0.02)
axis(3, at = seq(8000,16000, by=1000), labels=FALSE, tck=0.025)
axis(3, at = seq(8000,16000, by=1000), labels=FALSE, tck=-0.025)
# ticks on middle axis bottom
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.01)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=0.01)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.015)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=0.015)
axis(1, at = seq(8000,16000, by=1000), labels=FALSE, tck=0.025)
axis(1, at = seq(8000,16000, by=1000), labels=FALSE, tck=-0.025)
###############################################################################
par(mar=c(1.5,6,1.5,2))
plot(data_full.kde200,
     xaxt="n",
     xlim=c(16000,7500))
source("phaserect.R")
text(x= 16000, y=par("usr")[4]-(0.1*par("usr")[4]),
     labels="c. 200 y bandwidth", col="black", cex=1.5, font=2, pos=4)
axis(1, at = seq(9000,15000, by=2000))
# ticks on middle axis top
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.01)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=0.01)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.02)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=0.02)
axis(3, at = seq(8000,16000, by=1000), labels=FALSE, tck=0.025)
axis(3, at = seq(8000,16000, by=1000), labels=FALSE, tck=-0.025)

# ticks on middle axis bottom
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.01)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=0.01)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.015)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=0.015)
axis(1, at = seq(8000,16000, by=1000), labels=FALSE, tck=0.025)
axis(1, at = seq(8000,16000, by=1000), labels=FALSE, tck=-0.025)

###############################################################################
save.image(file='14cEnvironment.RData')
message("Environment saved for use in SPD workflow and other scripts")

file.edit("5 Analysis.R")
