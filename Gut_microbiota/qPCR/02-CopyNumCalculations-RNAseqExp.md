RNA-sequencing experiment: qPCR analyses of microbial loads in the gut
of gnotobiotic honey bees
================
Joanito Liberti and Lucie Kešnerová, University of Lausanne

## The following code processes raw qPCR data to calculate copies of the 16S rRNA gene per gut

``` r
setwd("/Users/joanitoliberti/The-gut-microbiota-affects-the-social-network-of-honeybees/Gut_microbiota/qPCR/")

# load data
cts = read.csv("RNAseqExp_guts_qPCR_rawdata.csv")
cts$X <- NULL
colnames(cts) <- c("Bee_ID", "Target", "Ct", "SD") # rename columns for better work

# rename levels of targets
cts$Target <- factor(cts$Target)
levels(cts$Target)
levels(cts$Target)[levels(cts$Target)=="actin_0113/0114"] <- "actin"
levels(cts$Target)[levels(cts$Target)=="UV_0356/0772"] <- "univ"

# extract the Hive and Box numbers based on Bee_ID format
control <- cts[grepl("[N]{1}[[:digit:]]{1}", cts$Bee_ID),]
cts <- cts[!grepl("[N]{1}[[:digit:]]{1}", cts$Bee_ID),] # delete negative control samples from the main dataset to be added later
```

``` r
#install.packages("tidyr")
library(dplyr)
library(tidyr)
# extract Bee number
cts <- separate(cts, col = 1, into = c("HiveBox","Bee_num"), sep = "[^[:alnum:]]+", remove = F) # gives number of bee (Bee_num) in the Box
# extract Hive number
ids <- cts$HiveBox
a <- strsplit(ids, "[0-9]")
b <- strsplit(ids, "[A-Z]")
cts$Hive <- paste0(lapply(a, "[", c(1)),lapply(b, "[", c(2))) # give Hive
# extract Box number
cts$Box <- sub("[H]{1}[[:digit:]]{1,}", "\\1", cts$HiveBox) # gives Box 

# merge together with info about treatment
library(readxl)
treatments <- read_excel("RNAseqExp-Treatments.xlsx", sheet=1) # loads the sheet with info
dt <- merge(cts,treatments) # adds Treatment and exp
```

## Check data distribution

``` r
#install.packages("ggplot2")
library(ggplot2)
ac <- ggplot(dt[dt$Target=="actin",], aes(x = Bee_ID, y = Ct)) + geom_bar(stat = "identity") + ggtitle("actin")
print(ac)
hist(dt$Ct[dt$Target=="actin"], breaks = 100) 
hist(dt$Ct[dt$Target=="univ"], breaks = 100)

# turn SD values to NA if Ct = NA
dt$SD[is.na(dt$Ct)] <- NA

# check proportion of not quality data
poor <- dt[dt$SD > 0.7 & dt$Ct < 30,]
poor <- poor[!is.na(poor$Ct),]
a = nrow(poor)
b = nrow(dt)
(a/b)*100 # 0 % of data can be doubtful
```

## Calculate copy numbers

``` r
# change order of columns
dt <- dt[c("Bee_ID","Hive","Box","HiveBox","Bee_num","Treatment","Target","Ct","SD")]

## add actin Ct value to each row with the same Bee_ID, regardless the target
actin = dt[dt$Target=="actin",c(1,8)] # generate new variable with only actin Cts
names(actin)[2] <- "actin_Ct" # rename this column, otherwise merge() wouldn't work as it should
n <- merge(dt[dt$Target!="actin",],actin) # merge 

## Calculate CopyNum for each target
StdCurves <- read_excel("StdCurves_RNAseqExp.xlsx", sheet=1) # export intercept, slope and threshold of detection values for each target
StdCurves <- StdCurves[c(1:2),c(2:5)] # select only data needed
StdCurves <- as.data.frame(StdCurves) 
n <- merge(n, StdCurves) # add them to the dataframe based on Target which is in common

# calculate CopyNums of actin
n$actin_CopyNum <- (StdCurves[1,4]^(StdCurves[1,2] - n$actin_Ct))*100 # n = E^(intercept - Cq), *gut sample volume: we used  half of the gut homogenate, and DNA was eluted in 50 ul
n$actin_CopyNum <- round(n$actin_CopyNum) # round it

# calculate 16S rRNA gene copies
n$CopyNum <- (n$Efficiency^(n$Intercept - n$Ct))*100 # n = E^(intercept - Cq), *gut sample volume: we used  half of the gut homogenate, and DNA was eluted in 50 ul
n$CopyNum <- round(n$CopyNum) # round it

# calculate median actin value
actin_mean <- round(mean(n$actin_CopyNum))
actin_median <- round(median(n$actin_CopyNum))

# calculate CopyNum_norm - normalized with median actin
n$CopyNum_norm <- round((n$CopyNum / n$actin_CopyNum) * actin_median)

n$Intercept <- NULL # get rid of extra columns
n$Slope <- NULL
n$Efficiency <- NULL
n$Threshold <- NULL
n$`LOD #copy` <- NULL

# EXPORT
write.csv(n, file="RNAseqExp_gut_qPCRdata.csv")
```
