###############################################################################
# Supplementary Information for Hoebe Peeters & Arnoldussen 2023 
###############################################################################
# 
###############################################################################
#### 2 Aoristic weight #### 
###############################################################################

# install and read relevant libraries
source(file="read libraries.R")

# run preparation script if necessary
if (exists("data_full")==FALSE){
  source("1 Dataset preparation.R")
}

# prepare and load data
#rewrite chronozone dataframe to include overflow categories . 
chronozones <- data.frame(row.names=1:8)
chronozones$name <- c("Older Dryas", "Bølling", "Allerød", 
                      "Younger Dryas", "Preboreal", "Early Boreal", 
                      "Late Boreal", "Atlantic")
chronozones$llim <- c(16000, 14650, 13904, 12846, 
                      11653, 10800, 9190, 8090)
chronozones$ulim <- c(14650, 13904, 12846, 11653, 
                      10800, 9190, 8090, 7500)

#chronozone duration in years and in centuries
chronozones$duration <-   chronozones$llim-chronozones$ulim
chronozones$centuries <- chronozones$duration/100


###############################################################################
#### 2.1 load the function aoristic()  #### 
###############################################################################
# aoristic() requires 6 parameters originating from two datasets
# 1 the data that you want weighed
# an identifier (data_id), and the start and end dates of each record
# 2 the periods you want the data distirbuted over
# an identifier (bin_id), and the start and end dates of each period.
# note that only dates that fall in three timeframes or less can be weighed
# using this function

source("function aoristic().R")


###############################################################################
#### 2.2 calculate aoristic weights  #### 
###############################################################################
#### a. aoristic weight per chronozone for the full dataset #####

# only dates in timeframe are used for aoristic calculation
lowerlim=16000
upperlim=7500
data <- filter(data_full,
               data_full$calBP_end<=lowerlim &
                 data_full$calBP_start>=upperlim) 

d_aor<- aoristic(
  data_id=data$labcode , 
  data_start = data$calBP_start,
  data_end = data$calBP_end,
  bin_id = chronozones$name,
  bin_start = chronozones$llim,
  bin_end = chronozones$ulim)

# add result to chronozone df
chronozones <- left_join(chronozones, d_aor$bins %>% 
                           select(bin_id,
                                  aor_n, 
                                  aor_w), 
                          by=c("name"="bin_id"))
rm(d_aor)

###############################################################################
#### b. aoristic weight per chronozone for site phases #####

phases <- data.frame(site_id=data_full$site_id,
                     phase=data_full$bins, 
                     site=data_full$site)
phases$calBP_start <- NA_real_
phases$calBP_end <- NA_real_

# calculate start and end dates per site phase
for (i in 1:nrow(phases)){
  # max age for each site phase
  phases$calBP_start[data_full$bins==phases$phase[i]] <- 
    max(data_full$calBP_start[data_full$bins==phases$phase[i]])
  # min age for each site phase
  phases$calBP_end[data_full$bins==phases$phase[i]] <- 
    min(data_full$calBP_end[data_full$bins==phases$phase[i]])
}
# remove duplicates and filter for timeframe
phases <- unique(phases)
phases <- filter(phases,phases$calBP_end<=lowerlim &
                   phases$calBP_start>=upperlim)
# aoristic function
p_aor <- aoristic(
  data_id=phases$phase , 
  data_start = phases$calBP_start,
  data_end = phases$calBP_end,
  bin_id = chronozones$name,
  bin_start = chronozones$llim,
  bin_end = chronozones$ulim)

# add result to chronozone df
chronozones <- left_join(chronozones, p_aor$bins %>% 
                           select(bin_id,
                           aor_p=aor_n, 
                           aor_pw=aor_w), 
                         by=c("name"="bin_id"))
rm(p_aor)

###############################################################################
#### c. aoristic weight per chronozone for the vetted dataset #####

data <- filter(data_vet,
               data_vet$calBP_end<=lowerlim &
                 data_vet$calBP_start>=upperlim)

# aoristic function
dv_aor<- aoristic(
  data_id=data$labcode , 
  data_start = data$calBP_start,
  data_end = data$calBP_end,
  bin_id = chronozones$name,
  bin_start = chronozones$llim,
  bin_end = chronozones$ulim)

# add result to chronozone df
chronozones <- left_join(chronozones, dv_aor$bins %>% 
                           select(bin_id, 
                           aor_n.v=aor_n, 
                           aor_w.v=aor_w), 
                         by=c("name"="bin_id"))
rm(dv_aor)

###############################################################################
#### d. aoristic weight per chronozone for site phases in vetted dataset #####

phases <- data.frame(site_id=data_vet$site_id,
                     phase=data_vet$bins, 
                     site=data_vet$site)
phases$calBP_start <- NA_real_
phases$calBP_end <- NA_real_

# calculate start and end dates per site phase
for (i in 1:nrow(phases)){
  phases$calBP_start[data_vet$bins==phases$phase[i]] <- 
    max(data_vet$calBP_start[data_vet$bins==phases$phase[i]])
  phases$calBP_end[data_vet$bins==phases$phase[i]] <- 
    min(data_vet$calBP_end[data_vet$bins==phases$phase[i]])
}
# filter for unique values and timeframe
phases <- unique(phases)
phases <- filter(phases,phases$calBP_end<=lowerlim &
                   phases$calBP_start>=upperlim)

# aoristic calculation
pv_aor <- aoristic(
  data_id=phases$phase , 
  data_start = phases$calBP_start,
  data_end = phases$calBP_end,
  bin_id = chronozones$name,
  bin_start = chronozones$llim,
  bin_end = chronozones$ulim)

# join result to dataframe
chronozones <- left_join(chronozones, pv_aor$bins %>% 
                           select(bin_id,
                           aor_p.v=aor_n, 
                           aor_pw.v=aor_w), 
                         by=c("name"="bin_id"))
rm(pv_aor, phases)
chronozones <- chronozones%>%filter(!name %in% c(">16000", "<7500"))
path <- "- chronozones aoristic distribution.xlsx"
write.xlsx(chronozones, file=path)


###############################################################################
#### 2.3 Plotting  #### 
###############################################################################
#### a. date density plot #####
par(mfrow=c(1,1))
par(oma=c(3,5,3,3))
par(mar=c(0,0,0,0))
par(xpd=NA)
plot(x=chronozones$duration, type="n", xlim=c(16000,7500), 
     ylim=c(0,max(chronozones$aor_w)+(0.1*max(chronozones$aor_w))),
     ylab = "date density per century",
     xaxs="i", yaxs="i")
phaserect <- 1
phaselabels <- -0.1
source("phaserect.R")
rect(x= chronozones$llim,
     xright = chronozones$ulim,
     ybottom = 0,
     ytop = chronozones$aor_w,
     col= "darkgrey",
     border="black")
text(x=chronozones$llim-(chronozones$duration/2), 
     y=chronozones$aor_w+(0.05*max(chronozones$aor_w)), 
     labels=chronozones$aor_w, col="black", cex=0.8)
text(x=chronozones$llim-(chronozones$duration/2), 
     y=(0.10*max(chronozones$aor_w)), 
     labels=chronozones$name, col="black", cex=0.8)
axis(1, at = seq(8000,16000, by=1000), labels = FALSE, tck=-0.04)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.015)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

###############################################################################
#### b. phase density plot ####
par(mfrow=c(1,1))
par(oma=c(3,5,3,3))
par(mar=c(0,0,0,0))
par(xpd=NA)
plot(x=chronozones$duration, type="n", xlim=c(16000,7500), 
     ylim=c(0,max(chronozones$aor_pw)+(0.20*max(chronozones$aor_pw))),
     ylab = "aoristic weight",
     xaxs="i", yaxs="i")
phaserect <- 1
phaselabels <- -0.1
source("phaserect.R")
rect(x= chronozones$llim,
     xright = chronozones$ulim,
     ybottom = 0,
     ytop = chronozones$aor_pw,
     col= "darkgrey",
     border="black")
text(x=chronozones$llim-(chronozones$duration/2), 
     y=chronozones$aor_pw+(0.05*max(chronozones$aor_pw)), 
     labels=chronozones$aor_pw, col="black", cex=0.8)
text(x=chronozones$llim-(chronozones$duration/2), 
     y=(0.10*max(chronozones$aor_pw)), 
     labels=chronozones$name, col="black", cex=0.8)
axis(1, at = seq(8000,16000, by=1000), labels = FALSE, tck=-0.04)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.015)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

###############################################################################
#### c. multiplot dates and phases ####
par(mfrow=c(2,1))
par(oma=c(3,5,3,3))
par(mar=c(1.5,0,0,0))
par(xpd=NA)
plot(x=chronozones$duration, type="n", xlim=c(16000,7500), 
     ylim=c(0,max(chronozones$aor_w)+(0.30*max(chronozones$aor_w))),
     xlab="",
     ylab = "aoristic weight",
     xaxs="i", xaxt="n",yaxs="i")
phaserect <- 1
phaselabels <- -0.15
source("phaserect.R")
rect(x= chronozones$llim,
     xright = chronozones$ulim,
     ybottom = 0,
     ytop = chronozones$aor_w,
     col= "darkgrey",
     border="black")
text(x=chronozones$llim-(chronozones$duration/2), 
     y=chronozones$aor_w+(0.1*max(chronozones$aor_w)), 
     labels=chronozones$aor_w, col="black", cex=0.8)
text(x=chronozones$llim-(chronozones$duration/2), 
     y=(0.10*max(chronozones$aor_w)), 
     labels=chronozones$name, col="black", cex=0.8)
text(x= 16000, y=par("usr")[4]-(0.15*max(chronozones$aor_w)),
     labels="a. Dates", col="black", cex=1, font=2, pos=4)
axis(3, at = seq(8000,16000, by=1000), labels = FALSE, tck=-0.04)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.015)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

axis(1, at = seq(8000,16000, by=1000), labels = TRUE, tck=-0.04)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.015)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

par(mar=c(3,0,1.5,0))
par(xpd=NA)
plot(x=chronozones$duration, type="n", xlim=c(16000,7500), 
     xlab="years cal BP",
     ylim=c(0,max(chronozones$aor_pw)+(0.30*max(chronozones$aor_pw))),
     ylab = "aoristic weight",
     xaxs="i", xaxt="n",yaxs="i")
phaserect <- 1
rm(phaselabels)
source("phaserect.R")
rect(x= chronozones$llim,
     xright = chronozones$ulim,
     ybottom = 0,
     ytop = chronozones$aor_pw,
     col= "darkgrey",
     border="black")
text(x=chronozones$llim-(chronozones$duration/2), 
     y=chronozones$aor_pw+(0.1*max(chronozones$aor_w)), 
     labels=chronozones$aor_pw, col="black", cex=0.8)
text(x=chronozones$llim-(chronozones$duration/2), 
     y=(0.10*max(chronozones$aor_pw)), 
     labels=chronozones$name, col="black",cex=0.8)
text(x= 16000, y=par("usr")[4]-(0.15*max(chronozones$aor_w)),
     labels="b. Site phases", col="black", cex=1, font=2, pos=4)
axis(3, at = seq(8000,16000, by=1000), labels = FALSE, tck=-0.04)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.015)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

axis(1, at = seq(8000,16000, by=1000), labels = TRUE, tck=-0.04)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.015)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

###############################################################################
#### d. multiplot vetting comparison dates ####
par(mfrow=c(2,1))
par(oma=c(3,5,3,3))
par(mar=c(1.5,0,0,0))
par(xpd=NA)
plot(x=chronozones$duration, type="n", xlim=c(16000,7500), 
     xlab="",
     ylim=c(0,max(chronozones$aor_w)+(0.30*max(chronozones$aor_w))),
     ylab = "aoristic weight",
     xaxs="i", xaxt="n",yaxs="i")
phaserect <- 1
phaselabels <- -0.15
source("phaserect.R")
rect(x= chronozones$llim,
     xright = chronozones$ulim,
     ybottom = 0,
     ytop = chronozones$aor_w,
     col= "darkgrey",
     border="black")
text(x=chronozones$llim-(chronozones$duration/2), 
     y=chronozones$aor_w+(0.1*max(chronozones$aor_w)), 
     labels=chronozones$aor_w, col="black", cex=0.8)
text(x=chronozones$llim-(chronozones$duration/2), 
     y=(0.10*max(chronozones$aor_w)), 
     labels=chronozones$name, col="black", cex=0.8)
text(x= 16000, y=par("usr")[4]-(0.15*max(chronozones$aor_w)),
     labels="a. Full dataset", col="black", cex=1, font=2, pos=4)
axis(3, at = seq(8000,16000, by=1000), labels = FALSE, tck=-0.04)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.015)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

axis(1, at = seq(8000,16000, by=1000), labels = TRUE, tck=-0.04)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.015)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

par(mar=c(3,0,1.5,0))
par(xpd=NA)
plot(x=chronozones$duration, type="n", xlim=c(16000,7500), 
     xlab="years cal BP",
     ylim=c(0,max(chronozones$aor_pw)+(0.30*max(chronozones$aor_pw))),
     #the max value for phase weight is used here instead of vetted n dates
     # because it ensures uniformity in the y axes between the two plots. 
     ylab = "aoristic weight",
     xaxs="i", xaxt="n",yaxs="i")
phaserect <- 1
rm(phaselabels)
source("phaserect.R")
rect(x= chronozones$llim,
     xright = chronozones$ulim,
     ybottom = 0,
     ytop = chronozones$aor_w.v,
     col= "darkgrey",
     border="black")
text(x=chronozones$llim-(chronozones$duration/2), 
     y=chronozones$aor_w.v+(0.1*max(chronozones$aor_w)), 
     labels=chronozones$aor_w.v, col="black", cex=0.8)
text(x=chronozones$llim-(chronozones$duration/2), 
     y=(0.10*max(chronozones$aor_w.v)), 
     labels=chronozones$name, col="black", cex=0.8)
text(x= 16000, y=par("usr")[4]-(0.15*max(chronozones$aor_w)),
     labels="b. Vetted dataset", col="black", cex=1, font=2, pos=4)
axis(3, at = seq(8000,16000, by=1000), labels = FALSE, tck=-0.04)
axis(3, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.015)
axis(3, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

axis(1, at = seq(8000,16000, by=1000), labels = TRUE, tck=-0.04)
axis(1, at = seq(7500,16000, by=100),labels=FALSE, tck=-0.015)
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)

rm(l, phaserect)

###############################################################################
file.edit("3 Vetting impact.R")