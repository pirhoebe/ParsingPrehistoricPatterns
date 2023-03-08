message("reading libraries...")

if (!require(readr)) install.packages("readr")
if (!require(readxl)) install.packages("readxl")
if (!require(openxlsx)) install.packages("openxlsx")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(data.table)) install.packages("data.table")
if (!require(rcarbon)) install.packages("rcarbon")
if (!require(stringr)) install.packages("stringr")
if (!require(beepr)) install.packages("beepr")

if(!"readr" %in% (.packages())){library(readr)} else 
  {message(paste("readr is attached"))}
if(!"readxl" %in% (.packages())){library(readxl)} else 
  {message(paste("readxl is attached"))}
if(!"openxlsx" %in% (.packages())){library(openxlsx)} else 
  {message(paste("openxlsx is attached"))}
if(!"dplyr" %in% (.packages())){library(dplyr)} else 
  {message(paste("dplyr is attached"))}
if(!"tidyr" %in% (.packages())){library(tidyr)} else 
  {message(paste("tidyr is attached"))}
if(!"data.table" %in% (.packages())){library(data.table)} else 
  {message(paste("data.table is attached"))}
if(!"rcarbon" %in% (.packages())){library(rcarbon)} else 
  {message(paste("rcarbon is attached"))}
if(!"stringr" %in% (.packages())){library(stringr)} else 
  {message(paste("stringr is attached"))}

message("complete")