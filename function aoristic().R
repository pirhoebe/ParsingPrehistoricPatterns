###############################################################################
#### Function: aoristic() #### 
###############################################################################

# define function:
aoristic <- function (data_id, data_start, data_end, 
                      bin_id, bin_start, bin_end){

# create dataframes
  # data
aor_data <- data.frame(data_id, data_start, data_end) 
aor_bin <- data.frame(bin_id=factor(bin_id, levels=bin_id, 
                                    ordered=TRUE), bin_start, bin_end)
aor_bin$duration <- aor_bin$bin_start-aor_bin$bin_end
aor_bin$centuries <- aor_bin$duration/100
  # bins
aor_data$bin_start <- NA_character_
aor_data$bin_end <- NA_character_

# assign bin name based on data_start
for (i in 1:length(aor_bin$bin_id)){
  aor_data$bin_start[
    aor_data$data_start<=aor_bin$bin_start[i]& 
      aor_data$data_start>aor_bin$bin_end[i]] <- aor_bin$bin_id[i]}

# assign bin name based on data_end
for (i in 1:length(aor_bin$bin_id)){
  aor_data$bin_end[
    aor_data$data_end<aor_bin$bin_start[i]& 
      aor_data$data_end>=aor_bin$bin_end[i]] <- aor_bin$bin_id[i]}

aor_data$n_bins <- (as.numeric(aor_data$bin_end) - 
                      as.numeric(aor_data$bin_start))+1

# function can only be used when dateranges do not exceed 3 bins
if(max(aor_data$n_bins[!is.na(aor_data$n_bins)]>3)){
  message(
"Error: margins are too high,
data that exceed 3 bins are present")}else{
Sys.sleep(1)

# calculating aoristic count for the differnt bins
message("\nassigning weight of data confined to 1 bin...\n")
for(i in 1:nrow(aor_bin)){ 
  chron <- as.numeric(aor_bin$bin_id[i])
  set <- filter(aor_data, as.numeric(aor_data$bin_start)==chron & 
                       as.numeric(aor_data$bin_end)==chron)
  message(paste(chron,":", nrow(set)))
  aor_bin$aor_n[i] <- nrow(set)
  
} 
Sys.sleep(1)
message("\ndividing aoristic weight of data spanning 2 bins...\n")
for(i in 1:nrow(aor_bin)){
  chron <- as.numeric(aor_bin$bin_id[i])
  chron2 <- as.numeric(aor_bin$bin_id[i+1])
  set <- filter(aor_data, as.numeric(aor_data$bin_start)==chron & 
                       as.numeric(aor_data$bin_end)==chron2)
  
  set$aor_1 <- (as.numeric(set$data_start)-aor_bin$bin_end[i])/
    (as.numeric(set$data_start)-as.numeric(set$data_end))
  set$aor_2 <- (aor_bin$bin_start[i+1]-as.numeric(set$data_end))/
    (as.numeric(set$data_start)-as.numeric(set$data_end))
  
  if((sum(set$aor_1)!=0)==TRUE){
    message(paste(aor_bin$bin_id[i],":", aor_bin$aor_n[i]," + ", 
                  round(sum(set$aor_1),2)," = ", 
                  aor_bin$aor_n[i] + round(sum(set$aor_1),2))) 
    aor_bin$aor_n[i] <- aor_bin$aor_n[i] + round(sum(set$aor_1),2)}
  
  if((sum(set$aor_2)!=0)==TRUE ){if(!is.na(aor_bin$aor_n[i+1])){
    message(paste(aor_bin$bin_id[i+1],":", aor_bin$aor_n[i+1]," + ", 
                  round(sum(set$aor_2),2)," = ", 
                  aor_bin$aor_n[i+1] + round(sum(set$aor_2),2)))     
    aor_bin$aor_n[i+1] <- aor_bin$aor_n[i+1] + 
      round(sum(set$aor_2),2)}}
}
Sys.sleep(1)
message("\ndividing aoristic weight of data spanning 3 bins...\n")
for(i in 1:nrow(aor_bin)){
  chron <- as.numeric(aor_bin$bin_id[i])
  chron2 <- as.numeric(aor_bin$bin_id[i+2])
  set <- filter(aor_data, as.numeric(aor_data$bin_start)==chron & 
                       as.numeric(aor_data$bin_end)==chron2)
  
  set$aor_1 <- (as.numeric(set$data_start)-aor_bin$bin_end[i])/
    (as.numeric(set$data_start)-as.numeric(set$data_end))
  set$aor_2 <- aor_bin$duration[i+1]/
    (as.numeric(set$data_start)-as.numeric(set$data_end))
  set$aor_3 <- (aor_bin$bin_start[i+2]-as.numeric(set$data_end))/
    (as.numeric(set$data_start)-as.numeric(set$data_end))
  
  if((sum(set$aor_1)!=0)==TRUE){
    message(paste(aor_bin$bin_id[i],":", aor_bin$aor_n[i]," + ", 
                  round(sum(set$aor_1),2)," = ", 
                  aor_bin$aor_n[i] + round(sum(set$aor_1),2))) 
    aor_bin$aor_n[i] <- aor_bin$aor_n[i] + round(sum(set$aor_1),2)}
  
  if((sum(set$aor_2)!=0)==TRUE ){  if(!is.na(aor_bin$aor_n[i+1])){
    message(paste(aor_bin$bin_id[i+1],":", aor_bin$aor_n[i+1]," + ", 
                  round(sum(set$aor_2),2)," = ", 
                  aor_bin$aor_n[i+1] + round(sum(set$aor_2),2)))     
    aor_bin$aor_n[i+1] <- aor_bin$aor_n[i+1] + 
      round(sum(set$aor_2),2)}}
  
  if((sum(set$aor_3)!=0)==TRUE ){  if(!is.na(aor_bin$aor_n[i+2])){
    message(paste(aor_bin$bin_id[i+2],":", aor_bin$aor_n[i+2]," + ", 
                  round(sum(set$aor_3),2)," = ", 
                  aor_bin$aor_n[i+2] + round(sum(set$aor_3),2)))    
    aor_bin$aor_n[i+2] <- aor_bin$aor_n[i+2] + 
      round(sum(set$aor_3),2)}}
} 
# calculate aoristic weight 
aor_bin$aor_w <- round(aor_bin$aor_n/aor_bin$duration,2)
# message
message(
"\naoristic count and weight are calculated for each bin:
see $bins$aor_n and $bins$aor_w

data was categorised by bin:
see $data$bin_start and $bin_end")
Sys.sleep(1)
print(list(bins=aor_bin,data=aor_data))
}}