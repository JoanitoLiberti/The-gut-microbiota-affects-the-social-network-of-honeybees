Automated behavioral tracking experiment: qPCR analyses of microbial
loads in the gut of gnotobiotic honey bees
================
Joanito Liberti, University of Lausanne

## The following code processes raw qPCR data to calculate copies of the 16S rRNA gene per gut

``` r
setwd("/Users/joanitoliberti/The-gut-microbiota-affects-the-social-network-of-honeybees/Gut_microbiota/qPCR/")

# load raw data
cts = read.csv("BeeTracking_guts_qPCR_rawdata.csv")
cts$X <- NULL
colnames(cts) <- c("Bee_ID", "Target", "Ct", "SD") # rename columns for better work

# rename levels of targets
cts$Target <- as.factor(cts$Target)
levels(cts$Target)
levels(cts$Target)[levels(cts$Target)=="actin_0113/0114"] <- "actin"
levels(cts$Target)[levels(cts$Target)=="UV_0356/0772"] <- "univ"

# merge metadata
library(readxl)
treatments <- read_excel("Metadata_BeeTracking_Guts.xlsx", sheet=1) # loads the sheet with info
selected<-c("MD","CL")
treatments <- treatments[treatments$Treatment %in% selected,]
dt <- merge(cts,treatments, by.x="Bee_ID", by.y="Sample_ID") # adds Treatment and exp#
```

## Check data distribution

``` r
#install.packages("ggplot2")
library(ggplot2)
ac <- ggplot(dt[dt$Target=="actin",], aes(x = Bee_ID, y = Ct)) + geom_bar(stat = "identity") + ggtitle("actin")
print(ac)
hist(dt$Ct[dt$Target=="actin"], breaks = 100)
hist(dt$Ct[dt$Target=="univ"], breaks = 100)

# turn remained SD values to NA if Ct = NA
dt$SD[is.na(dt$Ct)] <- NA

# check proportion of not quality data
poor <- dt[dt$SD > 0.7 & dt$Ct < 30,]
poor <- poor[!is.na(poor$Ct),]
a = nrow(poor)
b = nrow(dt)
(a/b)*100 
```

## Calculate copy numbers

``` r
actin = dt[dt$Target=="actin",c(1,3)] # generate new variable with only actin Cts
names(actin)[2] <- "actin_Ct" # rename this column, otherwise merge() wouldn't work as it should
n <- merge(dt[dt$Target!="actin",],actin) # merge 

## Calculate CopyNum for each target
StdCurves <- read_excel("StdCurves_BeeTracking.xlsx", sheet=1) # export intercept, slope and threshold of detection values for each target
StdCurves <- StdCurves[c(1:2),c(2:5)] # select only data needed
StdCurves <- as.data.frame(StdCurves) # table from excel is loaded as 'tibble' -> needs to be fixed
n <- merge(n, StdCurves) # add them to the dataframe based on Target which is in common

# calculate CopyNums of actin
n$actin_CopyNum <- (StdCurves[1,4]^(StdCurves[1,2] - n$actin_Ct))*400 # n = E^(intercept - Cq), *gut sample volume: we used  half of the gut homogenate, and DNA was eluted in 200 ul
n$actin_CopyNum <- round(n$actin_CopyNum) # round it

# calculate 16S rRNA gene copies
n$CopyNum <- (n$Efficiency^(n$Intercept - n$Ct))*400 # n = E^(intercept - Cq), *gut sample volume: we used  half of the gut homogenate, and DNA was eluted in 200 ul
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
write.csv(n, file="/Users/joanitoliberti/The-gut-microbiota-affects-the-social-network-of-honeybees/Gut_microbiota/qPCR/BeeTracking_gut_qPCRdata.csv")
```
