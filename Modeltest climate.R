###########################################################################
#### Modeltest climate #### 
###########################################################################

# prepare NGRIPdata
if(!require('pangaear')) install.packages('pangaear')
if(!"pangaear" %in% (.packages())){library(pangaear)} else 
  {message(paste("pangear is attached"))}
if(exists("NGRIPdata")==FALSE){
  NGRIP<- pg_data(doi='10.1594/PANGAEA.824889')
  NGRIPdata <- NGRIP[[1]]$data #download NGRIP data
  NGRIPdata <- rename(NGRIPdata, age = 2)
  NGRIPdata <- rename(NGRIPdata, d18O = 3)
  NGRIPdata$calBP <-NGRIPdata$age*1000 
  NGRIPdata$calBP <- round(NGRIPdata$calBP, 0)
  NGRIPdata <- NGRIPdata %>% filter(calBP<=20000)
  NGRIPdata$age <- NULL
  }

# model fitting
lowerlim=16200
upperlim=7300
message("this script fits the range of NGRIPdata to that of set.SPD")
Sys.sleep(3)
par(mfrow=c(2,1))
par(oma=c(2,2,2,2))
par(mar=c(1,0,0,0))
# past climate
plot(x=NGRIPdata$calBP,
     y=NGRIPdata$d18O,
     xlim=c(16000,7500),
     xaxs="i",
     type="l", col="firebrick", lwd=2)
# spd
par(mar=c(0,0,1,0))
plot(x=set.spd$grid$calBP,
     y=set.spd$grid$PrDens,
     xlim=c(16000,7500),
     xaxs="i",
     type="l", lwd=2)

###########################################################################
# calculate and plot correlation coefficient
message(
"is there a significant correlation between d18O and the SPD?")
Sys.sleep(3)
#corr.coeff
c <- left_join(NGRIPdata, set.spd$grid, by="calBP")
c <- filter(c, c$calBP <=16000 & c$calBP>=7500)
cor.test(c$d18O, c$PrDens, 
         method = 'spearman', exact=FALSE)
corr <- cor.test(c$d18O, c$PrDens, 
                 method = 'spearman', exact=FALSE)

par(mfrow=c(1,1))
par(oma=c(1,1,1,1))
par(mar=c(5,5,3,1))
plot(y=c$d18O, x=c$PrDens,
     main=paste("Spearman corr. coeff r =", round(corr$estimate,3)),
     pch=19,
     ylab = "NGRIP d18O\n temperature proxy",
     xlab = "14C summed probability density\n population proxy")
abline(lm(c$d18O~c$PrDens), col= "blue")


NGRIPdata.int <- approx(NGRIPdata$calBP, 
                        NGRIPdata$d18O,
                        method="linear",
                        xout=c((lowerlim):(upperlim)))

cliModel=data.frame(calBP= NGRIPdata.int$x, 
                    rTemp= NGRIPdata.int$y)

cliModel <- cliModel %>% filter(calBP<=lowerlim & calBP>=upperlim)


# # calculate mutual information
# install.packages("entropy")
# library(entropy)
# ?mi.plugin()
# 
# install.packages("infotheo")
# library(infotheo)
# 
# min(c$d18O)
# c$Tpos=c$d18O+45
# min(c$Tpos)
# freqs=matrix(sapply(seq(max(c$Tpos)*max(c$PrDens)), function(x) length(which(((c$Tpos-1)*max(c$PrDens)+c$PrDens)==x))),ncol=max(c$Tpos))
# mutinformation(c$Tpos, c$PrDens, method="emp")
# 
# min(standard$rTemp)
# min(standard$PrDens)
# mi$T <- standard$rTemp+2.5
# mi$P <- standard$PrDens+2.5  
# mi$Tr <- round(mi$T,2)
# mi$Pr <- round(mi$P,2)
# mi$TrA <- mi$Tr*100
# mi$PrA <- mi$Pr*100
# freqs=matrix(sapply(seq(max(mi$TrA)*max(mi$PrA)), function(x) length(which(((mi$TrA-1)*max(mi$PrA)+mi$PrA)==x))),ncol=max(mi$TrA))
# mi.plugin(freqs)

###########################################################################
# standardize data
message(
"standardizing transforms both datasets to the same scale, with mean = 0")
Sys.sleep(3)
standard <- left_join(cliModel, set.spd$grid, by = "calBP")
standard <- standard %>% 
  mutate_at(c("rTemp", "PrDens"),~(scale(.) %>% as.vector))

# plot histogram distribution of standardized climate and spd data
par(mfrow=c(2,1))
par(oma=c(2,2,2,1))
par(mar=c(2,2,1.5,0))
hist(x=standard$rTemp,
     xlim=c(-3,3),
     ylim=c(0,2500),
     breaks=c(seq(from=-3, to=3, by=0.25)),
     xlab=NA,
     main="relative temperature level frequency")
par(mar=c(1,2,1.5,0))
hist(standard$PrDens,
     xlim=c(-3,3),
     ylim=c(0,2500),
     breaks=c(seq(from=-3, to=3, by=0.25)),
     xlab=NA,
     main="relative probability density frequency")
Sys.sleep(3)

###########################################################################
# transform the standardized data to the scale of set.spd
standardT<- standard
r1 <- range(set.spd$grid$PrDens)[2]-range(set.spd$grid$PrDens)[1]
r0 <- range(standardT$PrDens)[2]-range(standardT$PrDens)[1]
# Transforming data to the same range as set.spd using rT
rT <- r0/r1 
r0/rT
standardT$PrDens <- standardT$PrDens/rT
standardT$rTemp <- standardT$rTemp/rT
range(standardT$PrDens)
range(standardT$rTemp)
# add the difference between the minimums
Diff <- min(set.spd$grid$PrDens) - min(standardT$PrDens)
Diff
standardT$PrDens <- standardT$PrDens+Diff
standardT$rTemp <- standardT$rTemp+Diff
# rTemp is standardized at the same range as the SPD
summary(set.spd$grid$PrDens)
summary(standardT$PrDens)
summary(standardT$rTemp)

###########################################################################
# both datasets are bimodal distributions, because of the higher frequency 
# of population and temperature at a certain level during different climate 
# regimes
# obtain modes
getmode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]}

cliset1 <- round(standardT$rTemp[standardT$rTemp<= 0.5],2)
cliset2 <- round(standardT$rTemp[standardT$rTemp> 0.5],2)
c_m1 <- getmode(x=cliset1)
c_m2 <- getmode(x=cliset2)
popset1 <- round(standardT$PrDens[standardT$PrDens<= 0.5],2)
popset2 <- round(standardT$PrDens[standardT$PrDens> 0.5],2)
p_m1 <- getmode(x=popset1)
p_m2 <- getmode(x=popset2)

# plot with modes 
par(mfrow=c(2,1))
par(oma=c(2,2,2,3))
par(mar=c(4,4,3,0))
hist(x=standardT$rTemp,
     xlim=c(-0.15,1.15),
     ylim=c(0,1400),
     breaks=c(seq(from=-0.2, to=1.5, by=0.05)),
     main=NA,
     xlab="relative Temperature",
     ylab="Frequency",
     col="darkgreen")
abline(v=c_m1, col="blue", lty=3, lwd=2)
abline(v=c_m2, col="red", lty=3, lwd=2)
par(xpd=NA)
text(x=c_m1, y=1600, labels="stadial mode")
text(x=c_m2, y=1600, labels="interstadial mode")
par(xpd=FALSE)
par(mar=c(4,4,4,0))
hist(standardT$PrDens,
     xlim=c(rev(range(-0.15,1.15))),
     ylim=c(0,1400),
     breaks=c(seq(from=-0.2, to=1.5, by=0.05)),
     main=NA,
     xlab="Summed Probability Density",
     ylab="Frequency",
     col="grey10")
abline(v=p_m1, col="blue", lwd=2)
abline(v=p_m2, col="red", lwd=2)
par(xpd=NA)
text(x=p_m1, y=1600, labels="stadial mode")
text(x=p_m2, y=1600, labels="interstadial mode")
par(xpd=FALSE)

###########################################################################
# plot the model fit
cliModSt <- data.frame(calBP=standardT$calBP, PrDens=standardT$rTemp)
par(mfrow=c(1,1))
par(oma=c(2,2,2,2))
par(mar=c(0,0,0,0))
plot(x=cliModSt$calBP, y=cliModSt$PrDens, type="l", col="darkgreen",
     ylim=c(min(cliModSt$PrDens), max(standardT$PrDens)),
     xlim=c(16000,7500),
     xaxs="i")
lines(x=standardT$calBP, y=standardT$PrDens, type="l", col="grey10",
      xlim = c(16000,7500), xaxs="i")
abline(h=p_m1, col="blue", lwd=2)
abline(h=c_m1, col="blue", lwd=2, lty=3)
text(x=7500, y= p_m1+0.03, 
     label=paste("spd stadial mode =", p_m1), pos=2)
text(x=7500, y= c_m1-0.03, 
     label=paste("cliMod stadial mode =", c_m1), pos=2)
abline(h=p_m2, col= "red", lwd=2)
abline(h=c_m2, col="red", lwd=2, lty=3)
text(x=16000, y= p_m2+0.03, 
     label=paste("spd interstadial mode =", p_m2), pos=4)
text(x=16000, y= c_m2-0.03, 
     label=paste("cliMod interstadial mode =", round(c_m2,2)), pos=4)
abline(h=mean(cliModSt$PrDens))
text(x=7500, 
     y= round(mean(cliModSt$PrDens),2)+0.03, 
     label=paste("mean =", round(mean(cliModSt$PrDens),2)), pos=2)

###########################################################################
# perform modeltest
Sys.sleep(5)
message(paste("Chosen parameters:",
              "\n simulations : ", simulations,
              "\n running mean: ", runningmean,
              "\n time range  : ", lowerlim, " - ", upperlim,
              "\n model       : ", "custom",
              "\n redundant normalisation: ", normalised,
              "\n processor cores: ", ncores
))
set.climod <- modelTest(set.cal, 
                        errors=set$c14sd, 
                        bins=bins,
                        nsim=simulations,
                        timeRange=c(lowerlim,upperlim), 
                        model="custom",
                        predgrid=cliModSt, 
                        runm=runningmean, 
                        raw=TRUE)

###########################################################################
# reset temporal frame
lowerlim <- 16000
upperlim <- 7500
# 
rm(c, cliModel, corr, NGRIPdata.int, NGRIP, set, set.cal, standard, 
   standardT, TempDensModel, c_m1, c_m2, p_m1, p_m2, cliset1, cliset2, 
   Diff, popset1, popset2, r0, r1, rT, getmode, rFromWilcox)
message("complete!")
