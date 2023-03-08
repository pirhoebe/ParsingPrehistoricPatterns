###############################################################################
# Supplementary Information for Hoebe Peeters & Arnoldussen 2023
###############################################################################
#
# This SI consists of 6 main scripts that:
# 1 prepare the dataset
# 2 demonstrate the general temporal distribution in the dataset with aoristics
# 3 explore the impact of vetting on the spatiotemporal distribution
# 4 show the general frequency distribution through SPD and KDE
# 5 conduct the main analysis: modeltesting and permutation testing the SPDs
# 6 Conduct some additional calculations
#
# There are several supporting scripts that are called upon in the main script.
# These facilitate running certain functions and serve to declutter the scripts.
# The execution of these supporting scripts requires input parameters that are
# displayed with a message.
# these parameters have to be set manually if one wants to use the supporting 
# scripts seperately.
#
# Pir Hoebe (p.w.hoebe@rug.nl), p.w.hoebe@rug.nl

###############################################################################
#### 1 Dataset preparation #### 
###############################################################################

# install and read relevant libraries
source(file="read libraries.R")

# load data
data <- read_excel("SI-2 Supplementary data.xlsx") 


###############################################################################
##### 1.1 Vetting procedure #####
###############################################################################

data$vet <- 1 # adding vetting column

# remove dates on unreliable materials
data$vet[data$material %in% c("bulk","other", "same sample",
                                        "sediment", "shell")] <- 0

# remove dates on aquatic species, reptiles, omnivores and carnivores.
data$vet[data$species %in% c("Homo","Ursus", "Vipera", 
                                       "Lynx", "Panthera", "Canis lupus", 
                                       "Canis familiaris", "Frogs", 
                                       "Esox", "Anguis", "Toads", "Turtles",
                                       "Vulpes", "Meles")] <- 0
# remove dates from bulk samples.
data$vet[like (data$sample, "bulk")] <- 0

# seperate vetted and full datasets
data_vet <- filter(data, data$vet==1)
data_full <- data


###############################################################################
#### 1.2 Calibrating data and subsets ####
###############################################################################

# PARAMETERS 
ncores <- 3 # set the number of processor cores that you want to set to work.
curve <- "intcal20"
normalised <- FALSE # turn of redundant normalisation

# calibration
set <- data_full
source("calibrate data.R")
data_full.cal <- set.cal
rm(set, set.cal)

# seperate the vetted dataset
data_vet.cal <- data_full.cal[data_full.cal$metadata$vet==1]


###############################################################################
#### 1.3 appending calibrated date ranges to dataset ####
###############################################################################

# add calBP start and end to the main dataset for further selection purposes
# create dataframe with the calibrated from - to data. 
calsum <- summary(data_full.cal)

# extracting the full two sigma calibrated range from calsum
calsum <- select(calsum, "DateID", "MedianBP",starts_with("TwoSigma"))

# splitting the date ranges of all the individual peaks in each date 
# distribution into a start date and an end date
calsum <- calsum %>% 
separate(TwoSigma_BP_1, c("Two_sigma_1_start", "Two_sigma_1_end"), " to ")%>% 
separate(TwoSigma_BP_2, c("Two_sigma_2_start", "Two_sigma_2_end"), " to ")%>% 
separate(TwoSigma_BP_3, c("Two_sigma_3_start", "Two_sigma_3_end"), " to ")%>% 
separate(TwoSigma_BP_4, c("Two_sigma_4_start", "Two_sigma_4_end"), " to ")%>% 
separate(TwoSigma_BP_5, c("Two_sigma_5_start", "Two_sigma_5_end"), " to ")%>% 
separate(TwoSigma_BP_6, c("Two_sigma_6_start", "Two_sigma_6_end"), " to ")%>% 
separate(TwoSigma_BP_7, c("Two_sigma_7_start", "Two_sigma_7_end"), " to ")%>% 
separate(TwoSigma_BP_8, c("Two_sigma_8_start", "Two_sigma_8_end"), " to ")%>% 
separate(TwoSigma_BP_9, c("Two_sigma_9_start", "Two_sigma_9_end"), " to ")

calsum[calsum=="NA"] <- NA

# earliest date in column 1
calsum$calBP_first <- calsum$Two_sigma_1_start
# latest date in columns 2 to 9
calsum$calBP_last <- calsum$Two_sigma_9_end
calsum$calBP_last[is.na(calsum$calBP_last)] <- 
  calsum$Two_sigma_8_end[is.na(calsum$calBP_last)]
calsum$calBP_last[is.na(calsum$calBP_last)] <- 
  calsum$Two_sigma_7_end[is.na(calsum$calBP_last)]
calsum$calBP_last[is.na(calsum$calBP_last)] <- 
  calsum$Two_sigma_6_end[is.na(calsum$calBP_last)]
calsum$calBP_last[is.na(calsum$calBP_last)] <- 
  calsum$Two_sigma_5_end[is.na(calsum$calBP_last)]
calsum$calBP_last[is.na(calsum$calBP_last)] <- 
  calsum$Two_sigma_4_end[is.na(calsum$calBP_last)]
calsum$calBP_last[is.na(calsum$calBP_last)] <- 
  calsum$Two_sigma_3_end[is.na(calsum$calBP_last)]
calsum$calBP_last[is.na(calsum$calBP_last)] <- 
  calsum$Two_sigma_2_end[is.na(calsum$calBP_last)]
calsum$calBP_last[is.na(calsum$calBP_last)] <- 
  calsum$Two_sigma_1_end[is.na(calsum$calBP_last)]

caldates <- data_full.cal$metadata
caldates <- left_join(x=caldates, 
                      y=calsum[ , c("DateID","calBP_first", "calBP_last")], 
                      by = c("DateID"="DateID"))
caldates$calBP_start <- as.numeric(caldates$calBP_first)
caldates$calBP_end <- as.numeric(caldates$calBP_last)

# append dataset with calBP start and end
data_full <- left_join(x=data_full, 
                       y=caldates[ , c("labcode","calBP_start", "calBP_end")], 
                       by = c("labcode"="labcode"))
rm(caldates, calsum)
data_full$dlength=data_full$calBP_start-data_full$calBP_end

# remake the vetted dataset to include new columns
data_vet <- filter(data_full, data_full$vet==1)


###############################################################################
#### 1.4 add chronozone column to dataset ####
###############################################################################

# create aoristic chronozone dataframe. 
chronozones <- data.frame(row.names=1:8)
chronozones$name <- c("Older Dryas", "Bølling", "Allerød", "Younger Dryas", 
                      "Preboreal", "Early Boreal", "Late Boreal", "Atlantic")
chronozones$llim <- c(25000, 14650, 13904, 12846, 11653, 10800, 9190, 8090)
chronozones$ulim <- c(14650, 13904, 12846, 11653, 10800, 9190, 8090, 6000)


data_full$chronozone_start <- NA
data_full$chronozone_end <- NA
# check for start dates that are within the limits of a chronozone
for (i in 1:length(chronozones$name)){
  data_full$chronozone_start[
    data_full$calBP_start<=chronozones$llim[i]& 
      data_full$calBP_start>chronozones$ulim[i]] <- chronozones$name[i]
}

# check for end dates that are within the limits of a chronozone
for (i in 1:length(chronozones$name)){
  data_full$chronozone_end[
    data_full$calBP_end<chronozones$llim[i]& 
      data_full$calBP_end>=chronozones$ulim[i]] <- chronozones$name[i]
}

rm(i)

# collapse data into one column
data_full$chronozone_end2 <- data_full$chronozone_end
data_full$chronozone_end2[
  data_full$chronozone_start==data_full$chronozone_end2] <- NA
data_full$chronozones <- NA
data_full$chronozones <- apply(cbind(data_full$chronozone_start,
                                data_full$chronozone_end2), 1,
                                   function(x) paste(x[!is.na(x)], 
                                                     collapse = " - "))
data_full$chronozone_end2 <- NULL

# factorize the chronozone data
data_full$chronozone_start <- factor(data_full$chronozone_start,
                                     levels=c("Older Dryas", 
                                              "Bølling", "Allerød", 
                                              "Younger Dryas", 
                                              "Preboreal", 
                                              "Early Boreal",
                                              "Late Boreal", 
                                              "Atlantic"))
data_full$chronozone_end <- factor(data_full$chronozone_end,
                                   levels=c("Older Dryas", 
                                            "Bølling", "Allerød", 
                                            "Younger Dryas", 
                                            "Preboreal", 
                                            "Early Boreal",
                                            "Late Boreal", 
                                            "Atlantic"))

data_vet <- filter(data_full, data_full$vet==1)


###############################################################################
#### 1.6 Binning ####
###############################################################################

# timeframe:
lowerlim <- 16000
upperlim <- 7500

binsense(x=data_full.cal,
         y=data_full.cal$metadata$site_id,          
         h=c(50,100,200,500),          
         timeRange=c(lowerlim, upperlim),          
         calendar = "BP") 

h_value <- 100 # h of 100 assigns dates from the same site that are within 100y
# of each other into the same bin 

# determine the number of phases (number depends on set binwidth)
# full dataset
data_full$bins = binPrep(sites=data_full$site_id,
                         ages=data_full$c14age, 
                         h=h_value)


# vetted dataset
data_vet$bins = binPrep(sites=data_vet$site_id,
                        ages=data_vet$c14age,
                        h=h_value)


###############################################################################
#### 1.7 Data overview ####
###############################################################################
# As SPD aggregation includes data outside the timeframe to minimize edge 
# effects the number of dates within the timeframe is calculated here and used 
# for reference. The numbers are based on dates whose calibrated ranges
# overlap at least a year with the time frame. 
# Note: These numbers may vary based on the used calibration curve!



# determine number of dates within timeframe for the full and vetted dataset
count <- list(full=list(dates=list(total=c(nrow(filter(data_full,
                                       data_full$calBP_end<=lowerlim &
                                         data_full$calBP_start>=upperlim))))), 
              vetted=list(dates=list(total=c(nrow(filter(data_full,
                                         data_full$calBP_end<=lowerlim &
                                           data_full$calBP_start>=upperlim &
                                           data_full$vet==1))))))
###############################################################################
# Counting the distribution in the full dataset
# determine the number of dates by country
count$full$dates$country <- list(
  "UK" = nrow(filter(data_full,
                     data_full$calBP_end<=lowerlim &
                       data_full$calBP_start>=upperlim &
                       data_full$country=="United Kingdom")),
  "BE" = nrow(filter(data_full,
                     data_full$calBP_end<=lowerlim &
                       data_full$calBP_start>=upperlim &
                       data_full$country=="Belgium")),
  "NL" = nrow(filter(data_full,
                     data_full$calBP_end<=lowerlim &
                       data_full$calBP_start>=upperlim &
                       data_full$country=="Netherlands")),
  "DE" = nrow(filter(data_full,
                     data_full$calBP_end<=lowerlim &
                       data_full$calBP_start>=upperlim &
                       data_full$country=="Germany")),
  "DK" = nrow(filter(data_full,
                     data_full$calBP_end<=lowerlim &
                       data_full$calBP_start>=upperlim &
                       data_full$country=="Denmark"))
)
# determine the number of phases within the timeframe...
count$full$phases$total <- 
  length(unique(data_full$bins[data_full$calBP_end<=lowerlim &
                               data_full$calBP_start>=upperlim]))
# ... and phases by country
count$full$phases$country <- list(
  "UK" = length(unique(data_full$bins[data_full$calBP_end<=lowerlim &
                                        data_full$calBP_start>=upperlim&
                                        data_full$country=="United Kingdom"])),
  "BE" = length(unique(data_full$bins[data_full$calBP_end<=lowerlim &
                                        data_full$calBP_start>=upperlim&
                                        data_full$country=="Belgium"])),
  "NL" = length(unique(data_full$bins[data_full$calBP_end<=lowerlim &
                                        data_full$calBP_start>=upperlim&
                                        data_full$country=="Netherlands"])),
  "DE" =length(unique(data_full$bins[data_full$calBP_end<=lowerlim &
                                       data_full$calBP_start>=upperlim &
                                       data_full$country=="Germany"])),
  "DK" = length(unique(data_full$bins[data_full$calBP_end<=lowerlim &
                                        data_full$calBP_start>=upperlim &
                                        data_full$country=="Denmark"]))
)
# determine the number of sites within the timeframe...
count$full$sites$total <- length(unique(
  data_full$site[data_full$calBP_end<=lowerlim &
                 data_full$calBP_start>=upperlim]))
# ...and sites by country
count$full$sites$country <- list(
  "UK" = length(unique(data_full$site[data_full$calBP_end<=lowerlim &
                                        data_full$calBP_start>=upperlim&
                                        data_full$country=="United Kingdom"])),
  "BE" = length(unique(data_full$site[data_full$calBP_end<=lowerlim &
                                        data_full$calBP_start>=upperlim&
                                        data_full$country=="Belgium"])),
  "NL" = length(unique(data_full$site[data_full$calBP_end<=lowerlim &
                                        data_full$calBP_start>=upperlim&
                                        data_full$country=="Netherlands"])),
  "DE" =length(unique(data_full$site[data_full$calBP_end<=lowerlim &
                                       data_full$calBP_start>=upperlim &
                                       data_full$country=="Germany"])),
  "DK" = length(unique(data_full$site[data_full$calBP_end<=lowerlim &
                                        data_full$calBP_start>=upperlim &
                                        data_full$country=="Denmark"]))
)
###############################################################################
# Counting the distribution in the vetted dataset
# determine the number of dates by country

count$vetted$dates$country <- list(
  "UK" = nrow(filter(data_vet,
                     data_vet$calBP_end<=lowerlim &
                       data_vet$calBP_start>=upperlim &
                       data_vet$country=="United Kingdom")),
  "BE" = nrow(filter(data_vet,
                     data_vet$calBP_end<=lowerlim &
                       data_vet$calBP_start>=upperlim &
                       data_vet$country=="Belgium")),
  "NL" = nrow(filter(data_vet,
                     data_vet$calBP_end<=lowerlim &
                       data_vet$calBP_start>=upperlim &
                       data_vet$country=="Netherlands")),
  "DE" = nrow(filter(data_vet,
                     data_vet$calBP_end<=lowerlim &
                       data_vet$calBP_start>=upperlim &
                       data_vet$country=="Germany")),
  "DK" = nrow(filter(data_vet,
                     data_vet$calBP_end<=lowerlim &
                       data_vet$calBP_start>=upperlim &
                       data_vet$country=="Denmark"))
)
# determine the number of phases within the timeframe...
count$vetted$phases$total <- length(
  unique(data_vet$bins[data_vet$calBP_end<=lowerlim &
                       data_vet$calBP_start>=upperlim]))
# ... and phases by country
count$vetted$phases$country <- list(
  "UK" = length(unique(data_vet$bins[data_vet$calBP_end<=lowerlim &
                                        data_vet$calBP_start>=upperlim&
                                        data_vet$country=="United Kingdom"])),
  "BE" = length(unique(data_vet$bins[data_vet$calBP_end<=lowerlim &
                                        data_vet$calBP_start>=upperlim&
                                        data_vet$country=="Belgium"])),
  "NL" = length(unique(data_vet$bins[data_vet$calBP_end<=lowerlim &
                                        data_vet$calBP_start>=upperlim&
                                        data_vet$country=="Netherlands"])),
  "DE" =length(unique(data_vet$bins[data_vet$calBP_end<=lowerlim &
                                       data_vet$calBP_start>=upperlim &
                                       data_vet$country=="Germany"])),
  "DK" = length(unique(data_vet$bins[data_vet$calBP_end<=lowerlim &
                                        data_vet$calBP_start>=upperlim &
                                        data_vet$country=="Denmark"]))
)
# determine the number of sites within the timeframe...
count$vetted$sites$total <- length(
  unique(data_vet$site[data_vet$calBP_end<=lowerlim &
                       data_vet$calBP_start>=upperlim]))
# ...and sites by country
count$vetted$sites$country <- list(
  "UK" = length(unique(data_vet$site[data_vet$calBP_end<=lowerlim &
                                        data_vet$calBP_start>=upperlim&
                                        data_vet$country=="United Kingdom"])),
  "BE" = length(unique(data_vet$site[data_vet$calBP_end<=lowerlim &
                                        data_vet$calBP_start>=upperlim&
                                        data_vet$country=="Belgium"])),
  "NL" = length(unique(data_vet$site[data_vet$calBP_end<=lowerlim &
                                        data_vet$calBP_start>=upperlim&
                                        data_vet$country=="Netherlands"])),
  "DE" =length(unique(data_vet$site[data_vet$calBP_end<=lowerlim &
                                       data_vet$calBP_start>=upperlim &
                                       data_vet$country=="Germany"])),
  "DK" = length(unique(data_vet$site[data_vet$calBP_end<=lowerlim &
                                        data_vet$calBP_start>=upperlim &
                                        data_vet$country=="Denmark"]))
)
###############################################################################
# transforming this overview to an exportable format
countdf <- data.frame(row.names=1:12)
countdf$selection <- NA_character_
countdf$country <- NA_character_
countdf$selection[1:6] <- "full"
countdf$selection[7:12] <- "vetted"
countdf$country [1:6] <- c("United Kingdom",  "Belgium", "Netherlands",
                           "Germany", "Denmark", "Total")
countdf$country [7:12] <- c("United Kingdom",  "Belgium", "Netherlands",
                            "Germany", "Denmark", "Total")
countdf$n.dates <- NA_integer_
countdf$n.dates[1:5] <- unlist(count$full$dates$country, use.names=FALSE)
countdf$n.dates[6] <- sum(countdf$n.dates[1:5])
countdf$n.dates[7:11] <- unlist(count$vetted$dates$country, use.names=FALSE)
countdf$n.dates[12] <- sum(countdf$n.dates[7:11])

countdf$n.phases <- NA_integer_
countdf$n.phases[1:5] <- unlist(count$full$phases$country, use.names=FALSE)
countdf$n.phases[6] <- sum(countdf$n.phases[1:5])
countdf$n.phases[7:11] <- unlist(count$vetted$phases$country, use.names=FALSE)
countdf$n.phases[12] <- sum(countdf$n.phases[7:11])

countdf$n.sites <- NA_integer_
countdf$n.sites[1:5] <- unlist(count$full$sites$country, use.names=FALSE)
countdf$n.sites[6] <- sum(countdf$n.sites[1:5])
countdf$n.sites[7:11] <- unlist(count$vetted$sites$country, use.names=FALSE)
countdf$n.sites[12] <- sum(countdf$n.sites[7:11])


# the following numbers are called on for plotting in subsequent scripts
ndates <- list(full=count$full$dates$total, 
               vetted=count$vetted$dates$total, 
               c=count$full$dates$country)


###############################################################################
data <- data_full %>% filter(calBP_end<=lowerlim &
                               calBP_start>=upperlim)

path <- "- data in timeframe.xlsx"
write.xlsx(data, file=path)

path <- "- data count overview.xlsx"
write.xlsx(countdf,  file=path)

save.image(file='14cEnvironment.RData')
message("Environment saved for use in SPD workflow and other scripts")
file.edit("2 Aoristic weight.R")
