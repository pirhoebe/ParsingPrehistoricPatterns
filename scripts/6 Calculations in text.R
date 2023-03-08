###############################################################################
# Supplementary Information for Hoebe Peeters & Arnoldussen 2022 
###############################################################################
#
###############################################################################
#### 6 Calculations in text ####
###############################################################################

# install and read relevant libraries
source(file="read libraries.R")

# to use this script, the main analysis must be completed first. 
if(!"landscape" %in% colnames(data_full)){
  file.edit("5 Analysis.R")  
}

lowerlim=16000
upperlim=7500
# calculations made within the timeframe
data_tf <- data_full %>% filter(data_full$calBP_end<=lowerlim &
                                  data_full$calBP_start>=upperlim)

###############################################################################
#### Research bias in landscape zones ####
###############################################################################

# n dates in the middle landscape zone during the Younger Dryas and Preboreal?
a <- nrow(data_tf %>% filter(
  landscape=="Mid" &
    (chronozone_start %in% c("Younger Dryas", "Preboreal")|
       chronozone_end %in% c("Younger Dryas", "Preboreal"))))

# portion of Dutch dates
b <-nrow(data_tf %>% filter(
  landscape=="Mid" & country == "Netherlands" &
    (chronozone_start %in% c("Younger Dryas", "Preboreal")|
       chronozone_end %in% c("Younger Dryas", "Preboreal"))))
a
b
x <- round(b/a*100,1)
message(paste("Of all dates in this zone and period",x,
              "% are from the Netherlands"))

###############################################################################
# n dates in the low landscape zones from the 8.2ka event onwards
a <- nrow(data_tf %>% filter(
  landscape=="Low" & calBP_end <8250
))

# portion of Dutch dates
b <- nrow(data_tf %>% filter(
  landscape=="Low" & calBP_end <8250 & country == "Netherlands"))
a
b
x <- round(b/a*100,1)
message(paste("Of all dates in this zone and period",x,
              "% are from the Netherlands"))

###############################################################################
#### Boreal peaks and sampling bias ####
###############################################################################
# hazelnut
# n Belgian dates in the Early Boreal
be <- nrow(data_tf %>% filter(
  country== "Belgium" & calBP_end <10700 & calBP_start>9190))

# n Dutch dates in the Early Boreal
nl <- nrow(data_tf %>% filter(
  country== "Netherlands" & calBP_end <10700 & calBP_start>9190))

# n corylus dates in the Early Boreal
a <- nrow(filter(data_tf, data_tf$species =="Corylus"
                 & calBP_end <10700 & calBP_start>9190))

# portion of Belgian Hazelnut
b <- nrow(filter(data_tf, data_tf$species =="Corylus"
                 & country=="Belgium"
                 & calBP_end <10700 & calBP_start>9190))
# portion of Dutch Hazelnut
c <- nrow(filter(data_tf, data_tf$species =="Corylus"
                 & country=="Netherlands"
                 & calBP_end <10700 & calBP_start>9190))
be
nl
a
b
c
x <- round(b/a*100,1)
y <- round(c/a*100,1)
message(paste(
  "Of all of the hazelnut dates in this period,\n",
  x,"% are from Belgium and",y,"% are from the Netherlands"))

x <- round(b/be*100,1)
y <- round(c/nl*100,1)
message(paste(
  "Of all of the Belgian dates in this period,\n",
  x,"% are on hazelnut"))
message(paste(
  "Of all of the Dutch dates in this period,\n",
  y,"% are on hazelnut"))

###############################################################################
# pit hearths

# n dates in the Netherlands during the Late Boreal
a <- nrow(data_tf %>% filter(country =="Netherlands"
                             & chronozone_end == "Late Boreal" & 
                               chronozone_start== "Late Boreal"))
# portion of pit hearths
b <- nrow(data_tf %>% filter(country =="Netherlands" & site_type=="pit hearths"
                             & chronozone_end == "Late Boreal" & 
                               chronozone_start== "Late Boreal"))

a
b
x <- round(b/a*100,1)
message(paste("Of all of the Dutch dates in this period",x,
              "% are from pit hearths"))

# n dates in Belgium during the Late Boreal
a <- nrow(data_tf %>% filter(country =="Belgium"
                             & chronozone_end == "Late Boreal" & 
                               chronozone_start== "Late Boreal"))
# portion of pit hearths
b <- nrow(data_tf %>% filter(country =="Belgium" & site_type=="pit hearths"
                             & chronozone_end == "Late Boreal" & 
                               chronozone_start== "Late Boreal"))

a
b
x <- round(b/a*100,1)
message(paste("Of all the Belgian dates in this period",x,
              "% are from pit hearths"))


###############################################################################
# plot the distribution of hazelnut and pit hearth dates (SI-8)
corylus <- filter(data_full, data_full$species =="Corylus")
set <- corylus
source("calibrate data.R")
corylus.cal <- set.cal
set.spd= spd(set.cal,timeRange=c(lowerlim,upperlim), bins = FALSE)
corylus.spd<- set.spd


pits <- filter(data_full, data_full$site_type =="pit hearths")
set <- pits
source("calibrate data.R")
pits.cal <- set.cal
set.spd= spd(set.cal,timeRange=c(lowerlim,upperlim), bins = FALSE)
pits.spd <- set.spd

# plot
par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
plot(data_full.spd, runm=50, xlim=c(11653,7500), xaxt= "n")
plot(corylus.spd, runm = 50, add=TRUE, type="simple", lwd=2,col="brown",
     )
plot(pits.spd, runm = 50, add=TRUE, type="simple",lwd=2,col="black")
axis(1, at = seq(8000,11000, by=1000))
axis(1, at = seq(7500,11600, by=100),labels=FALSE, tck=-0.02)
axis(1, at = seq(7500,11500, by=500),labels=FALSE, tck=-0.03)
rect(x= event$start,
     xright = event$end,
     ybottom = 0,
     ytop = par("usr")[3]+(l*(par("usr")[4])),
     col= rgb(0,0,0,alpha=0.3),
     border=NA)
abline(v=10800)
abline(v=9190)
abline(v=8090)
par(xpd=NA)
text(x=11653-((11653-10800)/2), y=1.2, labels= "Preboreal")
text(x=10800-((10800-9190)/2), y=1.2, labels= "Early Boreal")
text(x=9190-((9190-8090)/2), y=1.2, labels= "Late Boreal")
par(xpd=FALSE)
legend(x= 11653, 
       y=par("usr")[4]-(0.05*par("usr")[4]), 
       bty="n", 
       col=c("brown", "black"),
       legend=c("hazel", "pit hearths"),lwd=c(2,2), cex=1)

###############################################################################
###############################################################################
