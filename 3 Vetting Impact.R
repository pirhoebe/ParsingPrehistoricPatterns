###############################################################################
# Supplementary Information for Hoebe Peeters & Arnoldussen 2023 
###############################################################################
# 
###############################################################################
#### 3 Vetting impact #### 
###############################################################################

# install and read relevant libraries
source(file="read libraries.R")

# run preparation script if necessary
if (exists("data_full")==FALSE){
  source("1 Dataset preparation.R")
}


###############################################################################
#### 3.1 Vetting impact by country ####
###############################################################################
# Comparison of the distribution of dates among the vetted and full dataset
# by country

###############################################################################
# a. number of dates

data <- filter(data_full,
               data_full$calBP_end<=lowerlim &
                 data_full$calBP_start>=upperlim)
freqCountry <- #summarise per country
  data %>% 
  group_by(country, vet) %>% 
  summarise(n = n())

# pivot table
pivot_d <- spread(freqCountry, vet, n)
pivot_d$country <- factor(pivot_d$country,
                          levels=c("United Kingdom", "Belgium", "Netherlands", 
                                   "Germany", "Denmark"))
pivot_d <- pivot_d %>% arrange(country)
pivot_d$dates_full <- pivot_d$`0` + pivot_d$`1`
pivot_d$`0` <- NULL
pivot_d<- rename(pivot_d, dates_vetted= `1` )
pivot_d$`%dates_lost` <- 
  round(100-(pivot_d$dates_vetted*100/pivot_d$dates_full),1)
pivot_d <-pivot_d%>% select(country, dates_full, dates_vetted, `%dates_lost`)
rm(freqCountry)

###############################################################################
# b. number of phases

phases <- data.frame(site_id=data$site_id,
                    phase_id=data$bins,
                    country=data$country)
phases$vet <- 0
phases$vet[phases$phase_id %in% data_vet$bins] <- 1

freqCountry <- 
  phases %>% 
  group_by(country, vet) %>% 
  summarise(n = n())
# pivot table
pivot_p <- spread(freqCountry, vet, n)
pivot_p$country <- factor(pivot_p$country,
                          levels=c("United Kingdom", "Belgium", "Netherlands", 
                                   "Germany", "Denmark"))
pivot_p <- pivot_p %>% arrange(country)
pivot_p$phases_full <- pivot_p$`0` + pivot_p$`1`
pivot_p$`0` <- NULL
pivot_p<- rename(pivot_p, phases_vetted= `1` )
pivot_p$`%phases_lost` <- 
  round(100-(pivot_p$phases_vetted*100/pivot_p$phases_full),1)
pivot_p <-pivot_p%>% select(country, phases_full, phases_vetted, `%phases_lost`)
rm(freqCountry, phases)

###############################################################################
# c. number of sites

sites <- data.frame(site_id=data$site_id,
                      site=data$site,
                    country=data$country)
sites$vet <- 0
sites$vet[sites$site_id %in% data_vet$site_id] <- 1
freqCountry <- #summarise per country
  sites %>% 
  group_by(country, vet) %>% 
  summarise(n = n())

pivot_s <- spread(freqCountry, vet, n)
pivot_s$country <- factor(pivot_s$country,
                          levels=c("United Kingdom", "Belgium", "Netherlands", 
                                   "Germany", "Denmark"))
pivot_s <- pivot_s %>% arrange(country)
pivot_s$sites_full <- pivot_s$`0` + pivot_s$`1`
pivot_s$`0` <- NULL
pivot_s<- rename(pivot_s, sites_vetted= `1` )
pivot_s$`%sites_lost` <- 
  round(100-(pivot_s$sites_vetted*100/pivot_s$sites_full),1)
pivot_s <-pivot_s%>% select(country, sites_full, sites_vetted, `%sites_lost`)
rm(freqCountry, sites)

###############################################################################
# d. combined pivot table 
pivot <- cbind(pivot_d, 
               pivot_p[,c("phases_full", "phases_vetted", "%phases_lost")], 
               pivot_s[,c("sites_full", "sites_vetted", "%sites_lost")])

# save as xlsx
path <- "- vetting comparison by country.xlsx"
write.xlsx(pivot, file=path)
rm(pivot_d, pivot_p, pivot_s)


###############################################################################
#### 3.2 Vetting impact on SPD ####
###############################################################################
# prepare SPDs
# parameters
ncores <- 3
runningmean <- 50
lowerlim <- 16000
upperlim <- 7500

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
# calculate the difference between full and vetted
spd.vetcompare <- data_full.spd$grid
spd.vetcompare$PrDens_v <- data_vet.spd$grid$PrDens
spd.vetcompare$diff <- spd.vetcompare$PrDens - spd.vetcompare$PrDens_v


###############################################################################
#### 3.3 Plotting ####
###############################################################################
# plot the difference between full and vetted dataset through time. ####

par(mfrow=c(2,1))
par(yaxs="i")
par(xaxs="i")
par(mar=c(0,4,3,2))
par(oma=c(0,0,0,0))
# plot full SPD
plot(x=spd.vetcompare$calBP, 
     y=spd.vetcompare$PrDens, 
     type="l",
     lwd=2,
     ylim=c(0,1.1),
     xlim = c(16000,7500),
     xaxt="n",
     ylab = "Probability Density")
# add vetted SPD
lines(x=spd.vetcompare$calBP, 
      y=spd.vetcompare$PrDens_v, 
      col="darkgrey",
      type="l",
      lwd=2,
      xlim = c(16000,7500),
      xaxt="n")
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.015)
axis(1, at = seq(8000,16000, by=1000), labels=FALSE, tck=-0.025)
legend("topleft", bty="n", col = c("black", "darkgrey"),
       legend=c("full dataset", "vetted dataset"), lwd=c(2, 2), cex=1)

# plot difference between full and vetted
par(mar=c(8,4,1.5,2))
plot(x=spd.vetcompare$calBP, 
     y=spd.vetcompare$diff, 
     type="l",
     col= "red",
     xlim = c(16000,7500),
     ylim=c(0,0.4),
     xaxt="n",
     ylab="Density difference",
     xlab ="years cal BP")
# indicating low and high impact areas (>1 sd)
abline(h=mean(spd.vetcompare$diff), col="grey40")
abline(h=mean(spd.vetcompare$diff)-sd(spd.vetcompare$diff), col="grey")
abline(h=mean(spd.vetcompare$diff)+sd(spd.vetcompare$diff), col="grey")
lines(x=spd.vetcompare$calBP, 
      y=spd.vetcompare$diff, 
      type="l",
      lwd=2,
      col= "black")
lines(x=spd.vetcompare$calBP[spd.vetcompare$diff>
      (mean(spd.vetcompare$diff)+sd(spd.vetcompare$diff))], 
      y=spd.vetcompare$diff[spd.vetcompare$diff>
      (mean(spd.vetcompare$diff)+sd(spd.vetcompare$diff))], 
      type="l",
      lwd=2,
      col= "red")
lines(x=spd.vetcompare$calBP[spd.vetcompare$diff<
      (mean(spd.vetcompare$diff)-sd(spd.vetcompare$diff))], 
      y=spd.vetcompare$diff[spd.vetcompare$diff<
       (mean(spd.vetcompare$diff)-sd(spd.vetcompare$diff))], 
      type="l",
      lwd=2,
      col= "blue")
axis(1, at = seq(7500,16000, by=500),labels=FALSE, tck=-0.03)
axis(1, at = seq(8000,16000, by=1000), tck=-0.06)
legend("topleft", bty="n", col = c("black"),
       legend=c("vetting impact"), lwd=c(2), cex=1)
# values of mean impact and sd
mean(spd.vetcompare$diff)
sd(spd.vetcompare$diff)

rm(set.cal, set.spd, spd.vetcompare)

###############################################################################
file.edit("4 Frequency distribution.R")
