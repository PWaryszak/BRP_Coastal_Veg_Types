#Functions to processs raw BRP_Coastal_Veg_Types CN-data======
#Install 2 R-packages we need if not on your comp by now:
if (!require(readxl)) library(readxl) # load that package if not done already
if (!require(tidyverse)) library(tidyverse) # load that package if not done already

#Load these packages into your working environment
library("readxl")
library("tidyverse")

#IMPORTANT:
#YOUR SAMPLE_ID (I named it CNCode) should be consistent across all your spreadsheets.
setwd("~/00DeakinUni/R/BCL_R/BCL/BRP/CN_BRP_COASTAL_VEG")
#Specify path to folder where you store youre CN spreadsheets and list them all:
files <- list.files(path = "C:/Users/BlueCarbon/Documents/00DeakinUni/R/BCL_R/BCL/BRP/CN_BRP_COASTAL_VEG",
                    pattern = "*.xls", full.names = T)

#Create function to export data from sheet = 1 (Sample Table)
ReadCN_smp <- function(x) read_xlsx (path = x,
                                     sheet = 1,
                                     skip = 7,
                                     range = cell_cols("C:E"))

#Export "Sample Table" data from all files in your folder:
tbl1 <- sapply(files, ReadCN_smp, simplify=FALSE) %>%
  bind_rows(.id = "id1")

tbl1
#Create function to export data from sheet = 2 (Element%)
ReadCN <- function(x) read_xlsx (path = x,sheet = 2,skip = 7, range = cell_cols("C:D"))

#Export "Element%" data from all files in your folder using sapply:
tbl2 <- sapply(files, 
               ReadCN, simplify=FALSE) %>% 
  bind_rows(.id = "id2")

tbl2

#Bind sample (tbl1) and CN (tbl2) data together
CN_DATA <- cbind(tbl1,tbl2)#bind smp and CN data

#Double check if data alligns well:
all.equal(tbl1$id1,tbl2$id2) #should be TRUE!

#Clean up the file to keep data you need in a form you/R likes (no special signs):
CN_DATA_clean <- CN_DATA %>%
  filter(Type == "Smp") %>%
  select("id1", "Sample Name","Weight", "(N) %", "(C) %" ) %>%
  rename(file = "id1", CNCode = "Sample Name", Weight.mg = "Weight",
         N.percent = "(N) %", C.percent = "(C) %")

#Check for duplicates:
anyDuplicated(CN_DATA_clean$CNCode)#Should be 0!!!
#If not 0, some samples are dupliacted,
#You have to decide what to do with duplicates (e.g., average them, remove them)
#Extract duplicate elements:
CN_DATA_clean[duplicated(CN_DATA_clean)]#Nothing!


#Merge new CN data with your MASTER file:
MASTER_DATA <- read.csv("~/00DeakinUni/R/BCL_R/BCL/BRP/CoastalVeg.csv")
NewDATA <- left_join(MASTER_DATA, CN_DATA_clean, by = "CNCode")

dim(MASTER_DATA)#2233 rows   44 cols
dim(NewDATA)#2233 rows  48 cols
#write.csv(NewDATA, file = "CN_BRP_CoastalVeg_NewDATA.csv") #SEEK THIS FILE

