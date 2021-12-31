03-PlottingTimeCourses
================
Tomas Kay
12/21/2021

## This script plots average changes in behavior across the experiment overall and for each subcolony

### Set directories and load packages

``` r
MAINDIR <- "/Volumes/Lacie/BM"
INDIR   <- paste(MAINDIR, "RawData", sep = "/")
OUTDIR  <- paste(MAINDIR, "ProcessedData", sep = "/")
FIGSDIR <- paste(MAINDIR, "Figures", sep = "/")
library('jpeg'); library('gapminder'); library('ggfortify'); library('stringr');
library('igraph');library('pals');library('tidyverse'); library('lubridate')
```

### Import data

``` r
setwd(INDIR)
files <- list.files(pattern=".rds")
for (file in files){
  assign(sub("\\..*", "", file),readRDS(file))
}
```

### Trim data as in script 01

``` r
RawContacts_ALL_1C <- RawContacts_ALL_1C[RawContacts_ALL_1C$start > "2020-05-13 24:00:00",]
RawContacts_ALL_1D <- RawContacts_ALL_1D[RawContacts_ALL_1D$start > "2020-05-13 24:00:00",]
RawContacts_ALL_2C <- RawContacts_ALL_2C[RawContacts_ALL_2C$start > "2020-05-27 24:00:00",]
RawContacts_ALL_2D <- RawContacts_ALL_2D[RawContacts_ALL_2D$start > "2020-05-27 24:00:00",]
RawContacts_ALL_3C <- RawContacts_ALL_3C[RawContacts_ALL_3C$start > "2020-06-03 24:00:00",]
RawContacts_ALL_3D <- RawContacts_ALL_3D[RawContacts_ALL_3D$start > "2020-06-03 24:00:00",]
RawContacts_ALL_4C <- RawContacts_ALL_4C[RawContacts_ALL_4C$start > "2020-06-10 24:00:00",]
RawContacts_ALL_4D <- RawContacts_ALL_4D[RawContacts_ALL_4D$start > "2020-06-10 24:00:00",]
RawContacts_ALL_5C <- RawContacts_ALL_5C[RawContacts_ALL_5C$start > "2020-06-24 24:00:00",]
RawContacts_ALL_5D <- RawContacts_ALL_5D[RawContacts_ALL_5D$start > "2020-06-24 24:00:00",]
RawContacts_ALL_6C <- RawContacts_ALL_6C[RawContacts_ALL_6C$start > "2020-07-01 24:00:00",]
RawContacts_ALL_6D <- RawContacts_ALL_6D[RawContacts_ALL_6D$start > "2020-07-01 24:00:00",]
RawContacts_ALL_7C <- RawContacts_ALL_7C[RawContacts_ALL_7C$start > "2020-07-15 24:00:00",]
RawContacts_ALL_7D <- RawContacts_ALL_7D[RawContacts_ALL_7D$start > "2020-07-15 24:00:00",]
RawContacts_ALL_8C <- RawContacts_ALL_8C[RawContacts_ALL_8C$start > "2020-07-22 24:00:00",]
RawContacts_ALL_8D <- RawContacts_ALL_8D[RawContacts_ALL_8D$start > "2020-07-22 24:00:00",]
RawContacts_ALL_9C <- RawContacts_ALL_9C[RawContacts_ALL_9C$start > "2020-09-16 24:00:00",]
RawContacts_ALL_9D <- RawContacts_ALL_9D[RawContacts_ALL_9D$start > "2020-09-16 24:00:00",]
RawContacts_ALL_1C <- RawContacts_ALL_1C[RawContacts_ALL_1C$start < "2020-05-20 08:00:00",]
RawContacts_ALL_1D <- RawContacts_ALL_1D[RawContacts_ALL_1D$start < "2020-05-20 08:00:00",]
RawContacts_ALL_2C <- RawContacts_ALL_2C[RawContacts_ALL_2C$start < "2020-06-03 08:00:00",]
RawContacts_ALL_2D <- RawContacts_ALL_2D[RawContacts_ALL_2D$start < "2020-06-03 08:00:00",]
RawContacts_ALL_3C <- RawContacts_ALL_3C[RawContacts_ALL_3C$start < "2020-06-10 08:00:00",]
RawContacts_ALL_3D <- RawContacts_ALL_3D[RawContacts_ALL_3D$start < "2020-06-10 08:00:00",]
RawContacts_ALL_4C <- RawContacts_ALL_4C[RawContacts_ALL_4C$start < "2020-06-17 08:00:00",]
RawContacts_ALL_4D <- RawContacts_ALL_4D[RawContacts_ALL_4D$start < "2020-06-17 08:00:00",]
RawContacts_ALL_5C <- RawContacts_ALL_5C[RawContacts_ALL_5C$start < "2020-07-01 08:00:00",]
RawContacts_ALL_5D <- RawContacts_ALL_5D[RawContacts_ALL_5D$start < "2020-07-01 08:00:00",]
RawContacts_ALL_6C <- RawContacts_ALL_6C[RawContacts_ALL_6C$start < "2020-07-08 08:00:00",]
RawContacts_ALL_6D <- RawContacts_ALL_6D[RawContacts_ALL_6D$start < "2020-07-08 08:00:00",]
RawContacts_ALL_7C <- RawContacts_ALL_7C[RawContacts_ALL_7C$start < "2020-07-22 08:00:00",]
RawContacts_ALL_7D <- RawContacts_ALL_7D[RawContacts_ALL_7D$start < "2020-07-22 08:00:00",]
RawContacts_ALL_8C <- RawContacts_ALL_8C[RawContacts_ALL_8C$start < "2020-07-29 08:00:00",]
RawContacts_ALL_8D <- RawContacts_ALL_8D[RawContacts_ALL_8D$start < "2020-07-29 08:00:00",]
RawContacts_ALL_9C <- RawContacts_ALL_9C[RawContacts_ALL_9C$start < "2020-09-23 08:00:00",]
RawContacts_ALL_9D <- RawContacts_ALL_9D[RawContacts_ALL_9D$start < "2020-09-23 08:00:00",]
RawContacts_ALL_3C <- RawContacts_ALL_3C[!("2020-06-04 21:01:20" < RawContacts_ALL_3C$start & RawContacts_ALL_3C$start < "2020-06-04 21:01:20"),]
RawContacts_ALL_3C <- RawContacts_ALL_3C[!("2020-06-06 09:01:04" < RawContacts_ALL_3C$start & RawContacts_ALL_3C$start < "2020-06-06 11:56:55"),]
RawContacts_ALL_3C <- RawContacts_ALL_3C[!("2020-06-07 07:50:31" < RawContacts_ALL_3C$start & RawContacts_ALL_3C$start < "2020-06-07 09:29:01"),]
RawContacts_ALL_3C <- RawContacts_ALL_3C[!("2020-06-07 14:55:39" < RawContacts_ALL_3C$start & RawContacts_ALL_3C$start < "2020-06-07 15:38:22"),]
RawContacts_ALL_3D <- RawContacts_ALL_3D[!("2020-06-04 21:01:20" < RawContacts_ALL_3D$start & RawContacts_ALL_3D$start < "2020-06-04 21:01:20"),]
RawContacts_ALL_3D <- RawContacts_ALL_3D[!("2020-06-06 09:01:04" < RawContacts_ALL_3D$start & RawContacts_ALL_3D$start < "2020-06-06 11:56:55"),]
RawContacts_ALL_3D <- RawContacts_ALL_3D[!("2020-06-07 07:50:31" < RawContacts_ALL_3D$start & RawContacts_ALL_3D$start < "2020-06-07 09:29:01"),]
RawContacts_ALL_3D <- RawContacts_ALL_3D[!("2020-06-07 14:55:39" < RawContacts_ALL_3D$start & RawContacts_ALL_3D$start < "2020-06-07 15:38:22"),]
RawContacts_ALL_4D <- RawContacts_ALL_4D[!("2020-06-14 13:46:30" < RawContacts_ALL_4D$start & RawContacts_ALL_4D$start < "2020-06-14 20:08:30"),]
RawContacts_ALL_4D <- RawContacts_ALL_4D[!("2020-06-15 22:59:18" < RawContacts_ALL_4D$start & RawContacts_ALL_4D$start < "2020-06-16 06:35:05"),]
RawContacts_ALL_4C <- RawContacts_ALL_4C[!("2020-06-14 13:46:30" < RawContacts_ALL_4C$start & RawContacts_ALL_4C$start < "2020-06-14 20:08:30"),]
RawContacts_ALL_4C <- RawContacts_ALL_4C[!("2020-06-15 22:59:18" < RawContacts_ALL_4C$start & RawContacts_ALL_4C$start < "2020-06-16 06:35:05"),]
RawContacts_ALL_5C <- RawContacts_ALL_5C[!("2020-06-26 20:38:51" < RawContacts_ALL_5C$start & RawContacts_ALL_5C$start < "2020-06-27 07:58:18"),]
RawContacts_ALL_5C <- RawContacts_ALL_5C[!("2020-06-27 13:48:20" < RawContacts_ALL_5C$start & RawContacts_ALL_5C$start < "2020-06-27 16:39:02"),]
RawContacts_ALL_5C <- RawContacts_ALL_5C[!("2020-06-30 28:58:42" < RawContacts_ALL_5C$start),]
RawContacts_ALL_5D <- RawContacts_ALL_5D[!("2020-06-26 20:38:51" < RawContacts_ALL_5D$start & RawContacts_ALL_5D$start < "2020-06-27 07:58:18"),]
RawContacts_ALL_5D <- RawContacts_ALL_5D[!("2020-06-27 13:48:20" < RawContacts_ALL_5D$start & RawContacts_ALL_5D$start < "2020-06-27 16:39:02"),]
RawContacts_ALL_5D <- RawContacts_ALL_5D[!("2020-06-30 28:58:42" < RawContacts_ALL_5D$start),]
RawContacts_ALL_7C <- RawContacts_ALL_7C[!("2020-07-20 17:27:42" < RawContacts_ALL_7C$start & RawContacts_ALL_7C$start < "2020-07-20 18:24:29"),]
RawContacts_ALL_7C <- RawContacts_ALL_7C[!("2020-07-20 23:47:14" < RawContacts_ALL_7C$start & RawContacts_ALL_7C$start < "2020-07-21 13:04:28"),]
RawContacts_ALL_7D <- RawContacts_ALL_7D[!("2020-07-20 17:27:42" < RawContacts_ALL_7D$start & RawContacts_ALL_7D$start < "2020-07-20 18:24:29"),]
RawContacts_ALL_7D <- RawContacts_ALL_7D[!("2020-07-20 23:47:14" < RawContacts_ALL_7D$start & RawContacts_ALL_7D$start < "2020-07-21 13:04:28"),]

setwd(OUTDIR)
SocOut1C <- read.csv("SocOut1C.csv")$x
SocOut2C <- read.csv("SocOut2C.csv")$x
SocOut3C <- read.csv("SocOut3C.csv")$x
SocOut4C <- read.csv("SocOut4C.csv")$x
SocOut5C <- read.csv("SocOut5C.csv")$x
SocOut6C <- read.csv("SocOut6C.csv")$x
SocOut7C <- read.csv("SocOut7C.csv")$x
SocOut8C <- read.csv("SocOut8C.csv")$x
SocOut9C <- read.csv("SocOut9C.csv")$x
SocOut1D <- read.csv("SocOut1D.csv")$x
SocOut2D <- read.csv("SocOut2D.csv")$x
SocOut3D <- read.csv("SocOut3D.csv")$x
SocOut4D <- read.csv("SocOut4D.csv")$x
SocOut5D <- read.csv("SocOut5D.csv")$x
SocOut6D <- read.csv("SocOut6D.csv")$x
SocOut7D <- read.csv("SocOut7D.csv")$x
SocOut8D <- read.csv("SocOut8D.csv")$x
SocOut9D <- read.csv("SocOut9D.csv")$x
RawContacts_ALL_1C <- RawContacts_ALL_1C[!(RawContacts_ALL_1C$ant1 %in% SocOut1C) & !(RawContacts_ALL_1C$ant2 %in% SocOut1C),]
RawContacts_ALL_2C <- RawContacts_ALL_2C[!(RawContacts_ALL_2C$ant1 %in% SocOut2C) & !(RawContacts_ALL_2C$ant2 %in% SocOut2C),]
RawContacts_ALL_3C <- RawContacts_ALL_3C[!(RawContacts_ALL_3C$ant1 %in% SocOut3C) & !(RawContacts_ALL_3C$ant2 %in% SocOut3C),]
RawContacts_ALL_4C <- RawContacts_ALL_4C[!(RawContacts_ALL_4C$ant1 %in% SocOut4C) & !(RawContacts_ALL_4C$ant2 %in% SocOut4C),]
RawContacts_ALL_5C <- RawContacts_ALL_5C[!(RawContacts_ALL_5C$ant1 %in% SocOut5C) & !(RawContacts_ALL_5C$ant2 %in% SocOut5C),]
RawContacts_ALL_6C <- RawContacts_ALL_6C[!(RawContacts_ALL_6C$ant1 %in% SocOut6C) & !(RawContacts_ALL_6C$ant2 %in% SocOut6C),]
RawContacts_ALL_7C <- RawContacts_ALL_7C[!(RawContacts_ALL_7C$ant1 %in% SocOut7C) & !(RawContacts_ALL_7C$ant2 %in% SocOut7C),]
RawContacts_ALL_8C <- RawContacts_ALL_8C[!(RawContacts_ALL_8C$ant1 %in% SocOut8C) & !(RawContacts_ALL_8C$ant2 %in% SocOut8C),]
RawContacts_ALL_9C <- RawContacts_ALL_9C[!(RawContacts_ALL_9C$ant1 %in% SocOut9C) & !(RawContacts_ALL_9C$ant2 %in% SocOut9C),]
RawContacts_ALL_1D <- RawContacts_ALL_1D[!(RawContacts_ALL_1D$ant1 %in% SocOut1D) & !(RawContacts_ALL_1D$ant2 %in% SocOut1D),]
RawContacts_ALL_2D <- RawContacts_ALL_2D[!(RawContacts_ALL_2D$ant1 %in% SocOut2D) & !(RawContacts_ALL_2D$ant2 %in% SocOut2D),]
RawContacts_ALL_3D <- RawContacts_ALL_3D[!(RawContacts_ALL_3D$ant1 %in% SocOut3D) & !(RawContacts_ALL_3D$ant2 %in% SocOut3D),]
RawContacts_ALL_4D <- RawContacts_ALL_4D[!(RawContacts_ALL_4D$ant1 %in% SocOut4D) & !(RawContacts_ALL_4D$ant2 %in% SocOut4D),]
RawContacts_ALL_5D <- RawContacts_ALL_5D[!(RawContacts_ALL_5D$ant1 %in% SocOut5D) & !(RawContacts_ALL_5D$ant2 %in% SocOut5D),]
RawContacts_ALL_6D <- RawContacts_ALL_6D[!(RawContacts_ALL_6D$ant1 %in% SocOut6D) & !(RawContacts_ALL_6D$ant2 %in% SocOut6D),]
RawContacts_ALL_7D <- RawContacts_ALL_7D[!(RawContacts_ALL_7D$ant1 %in% SocOut7D) & !(RawContacts_ALL_7D$ant2 %in% SocOut7D),]
RawContacts_ALL_8D <- RawContacts_ALL_8D[!(RawContacts_ALL_8D$ant1 %in% SocOut8D) & !(RawContacts_ALL_8D$ant2 %in% SocOut8D),]
RawContacts_ALL_9D <- RawContacts_ALL_9D[!(RawContacts_ALL_9D$ant1 %in% SocOut9D) & !(RawContacts_ALL_9D$ant2 %in% SocOut9D),]
```

### Select head-head interactions (mostly encoded as 1-1 but inverted for some replicates)

``` r
RawContacts_11_1C <- RawContacts_ALL_1C[RawContacts_ALL_1C$types == "1-1",]
RawContacts_11_2C <- RawContacts_ALL_2C[RawContacts_ALL_2C$types == "1-1",]
RawContacts_11_3C <- RawContacts_ALL_3C[RawContacts_ALL_3C$types == "1-1",]
RawContacts_11_4C <- RawContacts_ALL_4C[RawContacts_ALL_4C$types == "1-1",]
RawContacts_11_5C <- RawContacts_ALL_5C[RawContacts_ALL_5C$types == "1-1",]
RawContacts_11_6C <- RawContacts_ALL_6C[RawContacts_ALL_6C$types == "2-2",]
RawContacts_11_7C <- RawContacts_ALL_7C[RawContacts_ALL_7C$types == "2-2",]
RawContacts_11_8C <- RawContacts_ALL_8C[RawContacts_ALL_8C$types == "2-2",]
RawContacts_11_9C <- RawContacts_ALL_9C[RawContacts_ALL_9C$types == "2-2",]
RawContacts_11_1D <- RawContacts_ALL_1D[RawContacts_ALL_1D$types == "1-1",]
RawContacts_11_2D <- RawContacts_ALL_2D[RawContacts_ALL_2D$types == "1-1",]
RawContacts_11_3D <- RawContacts_ALL_3D[RawContacts_ALL_3D$types == "1-1",]
RawContacts_11_4D <- RawContacts_ALL_4D[RawContacts_ALL_4D$types == "1-1",]
RawContacts_11_5D <- RawContacts_ALL_5D[RawContacts_ALL_5D$types == "1-1",]
RawContacts_11_6D <- RawContacts_ALL_6D[RawContacts_ALL_6D$types == "1-1",]
RawContacts_11_7D <- RawContacts_ALL_7D[RawContacts_ALL_7D$types == "2-2",]
RawContacts_11_8D <- RawContacts_ALL_8D[RawContacts_ALL_8D$types == "2-2",]
RawContacts_11_9D <- RawContacts_ALL_9D[RawContacts_ALL_9D$types == "2-2",]
```

# Separate data by arena

``` r
RawContacts_11_N_1C <- RawContacts_11_1C[RawContacts_11_1C$space == 1,]
RawContacts_11_F_1C <- RawContacts_11_1C[RawContacts_11_1C$space == 2,]
RawContacts_11_N_2C <- RawContacts_11_2C[RawContacts_11_2C$space == 2,]
RawContacts_11_F_2C <- RawContacts_11_2C[RawContacts_11_2C$space == 1,]
RawContacts_11_N_3C <- RawContacts_11_3C[RawContacts_11_3C$space == 2,]
RawContacts_11_F_3C <- RawContacts_11_3C[RawContacts_11_3C$space == 1,]
RawContacts_11_N_4C <- RawContacts_11_4C[RawContacts_11_4C$space == 2,]
RawContacts_11_F_4C <- RawContacts_11_4C[RawContacts_11_4C$space == 1,]
RawContacts_11_N_5C <- RawContacts_11_5C[RawContacts_11_5C$space == 2,]
RawContacts_11_F_5C <- RawContacts_11_5C[RawContacts_11_5C$space == 1,]
RawContacts_11_N_6C <- RawContacts_11_6C[RawContacts_11_6C$space == 2,]
RawContacts_11_F_6C <- RawContacts_11_6C[RawContacts_11_6C$space == 1,]
RawContacts_11_N_7C <- RawContacts_11_7C[RawContacts_11_7C$space == 2,]
RawContacts_11_F_7C <- RawContacts_11_7C[RawContacts_11_7C$space == 1,]
RawContacts_11_N_8C <- RawContacts_11_8C[RawContacts_11_8C$space == 2,]
RawContacts_11_F_8C <- RawContacts_11_8C[RawContacts_11_8C$space == 1,]
RawContacts_11_N_9C <- RawContacts_11_9C[RawContacts_11_9C$space == 2,]
RawContacts_11_F_9C <- RawContacts_11_9C[RawContacts_11_9C$space == 1,]
RawContacts_11_N_1D <- RawContacts_11_1D[RawContacts_11_1D$space == 2,]
RawContacts_11_F_1D <- RawContacts_11_1D[RawContacts_11_1D$space == 1,]
RawContacts_11_N_2D <- RawContacts_11_2D[RawContacts_11_2D$space == 2,]
RawContacts_11_F_2D <- RawContacts_11_2D[RawContacts_11_2D$space == 1,]
RawContacts_11_N_3D <- RawContacts_11_3D[RawContacts_11_3D$space == 2,]
RawContacts_11_F_3D <- RawContacts_11_3D[RawContacts_11_3D$space == 1,]
RawContacts_11_N_4D <- RawContacts_11_4D[RawContacts_11_4D$space == 2,]
RawContacts_11_F_4D <- RawContacts_11_4D[RawContacts_11_4D$space == 1,]
RawContacts_11_N_5D <- RawContacts_11_5D[RawContacts_11_5D$space == 2,]
RawContacts_11_F_5D <- RawContacts_11_5D[RawContacts_11_5D$space == 1,]
RawContacts_11_N_6D <- RawContacts_11_6D[RawContacts_11_6D$space == 2,]
RawContacts_11_F_6D <- RawContacts_11_6D[RawContacts_11_6D$space == 1,]
RawContacts_11_N_7D <- RawContacts_11_7D[RawContacts_11_7D$space == 2,]
RawContacts_11_F_7D <- RawContacts_11_7D[RawContacts_11_7D$space == 1,]
RawContacts_11_N_8D <- RawContacts_11_8D[RawContacts_11_8D$space == 2,]
RawContacts_11_F_8D <- RawContacts_11_8D[RawContacts_11_8D$space == 1,]
RawContacts_11_N_9D <- RawContacts_11_9D[RawContacts_11_9D$space == 2,]
RawContacts_11_F_9D <- RawContacts_11_9D[RawContacts_11_9D$space == 1,]
```

### Calculate the number of interactions per hour in each of the arena

``` r
TimeCourse1C_F <- table(RawContacts_11_F_1C$hour)/length(unique(c(RawContacts_11_1C$ant1, RawContacts_11_1C$ant2)))
names(TimeCourse1C_F) <- as.numeric(names(TimeCourse1C_F)) - as.numeric(names(TimeCourse1C_F)[1])
TimeCourse2C_F <- table(RawContacts_11_F_2C$hour)/length(unique(c(RawContacts_11_2C$ant1, RawContacts_11_2C$ant2)))
names(TimeCourse2C_F) <- as.numeric(names(TimeCourse2C_F)) - as.numeric(names(TimeCourse2C_F)[1])
TimeCourse3C_F <- table(RawContacts_11_F_3C$hour)/length(unique(c(RawContacts_11_3C$ant1, RawContacts_11_3C$ant2)))
names(TimeCourse3C_F) <- as.numeric(names(TimeCourse3C_F)) - as.numeric(names(TimeCourse3C_F)[1])
TimeCourse4C_F <- table(RawContacts_11_F_4C$hour)/length(unique(c(RawContacts_11_4C$ant1, RawContacts_11_4C$ant2)))
names(TimeCourse4C_F) <- as.numeric(names(TimeCourse4C_F)) - as.numeric(names(TimeCourse4C_F)[1])
TimeCourse5C_F <- table(RawContacts_11_F_5C$hour)/length(unique(c(RawContacts_11_5C$ant1, RawContacts_11_5C$ant2)))
names(TimeCourse5C_F) <- as.numeric(names(TimeCourse5C_F)) - as.numeric(names(TimeCourse5C_F)[1])
TimeCourse6C_F <- table(RawContacts_11_F_6C$hour)/length(unique(c(RawContacts_11_6C$ant1, RawContacts_11_6C$ant2)))
names(TimeCourse6C_F) <- as.numeric(names(TimeCourse6C_F)) - as.numeric(names(TimeCourse6C_F)[1])
TimeCourse7C_F <- table(RawContacts_11_F_7C$hour)/length(unique(c(RawContacts_11_7C$ant1, RawContacts_11_7C$ant2)))
names(TimeCourse7C_F) <- as.numeric(names(TimeCourse7C_F)) - as.numeric(names(TimeCourse7C_F)[1])
TimeCourse8C_F <- table(RawContacts_11_F_8C$hour)/length(unique(c(RawContacts_11_8C$ant1, RawContacts_11_8C$ant2)))
names(TimeCourse8C_F) <- as.numeric(names(TimeCourse8C_F)) - as.numeric(names(TimeCourse8C_F)[1])
TimeCourse9C_F <- table(RawContacts_11_F_9C$hour)/length(unique(c(RawContacts_11_9C$ant1, RawContacts_11_9C$ant2)))
names(TimeCourse9C_F) <- as.numeric(names(TimeCourse9C_F)) - as.numeric(names(TimeCourse9C_F)[1])
TimeCourse1D_F <- table(RawContacts_11_F_1D$hour)/length(unique(c(RawContacts_11_1D$ant1, RawContacts_11_1D$ant2)))
names(TimeCourse1D_F) <- as.numeric(names(TimeCourse1D_F)) - as.numeric(names(TimeCourse1D_F)[1])
TimeCourse2D_F <- table(RawContacts_11_F_2D$hour)/length(unique(c(RawContacts_11_2D$ant1, RawContacts_11_2D$ant2)))
names(TimeCourse2D_F) <- as.numeric(names(TimeCourse2D_F)) - as.numeric(names(TimeCourse2D_F)[1])
TimeCourse3D_F <- table(RawContacts_11_F_3D$hour)/length(unique(c(RawContacts_11_3D$ant1, RawContacts_11_3D$ant2)))
names(TimeCourse3D_F) <- as.numeric(names(TimeCourse3D_F)) - as.numeric(names(TimeCourse3D_F)[1])
TimeCourse4D_F <- table(RawContacts_11_F_4D$hour)/length(unique(c(RawContacts_11_4D$ant1, RawContacts_11_4D$ant2)))
names(TimeCourse4D_F) <- as.numeric(names(TimeCourse4D_F)) - as.numeric(names(TimeCourse4D_F)[1])
TimeCourse5D_F <- table(RawContacts_11_F_5D$hour)/length(unique(c(RawContacts_11_5D$ant1, RawContacts_11_5D$ant2)))
names(TimeCourse5D_F) <- as.numeric(names(TimeCourse5D_F)) - as.numeric(names(TimeCourse5D_F)[1])
TimeCourse6D_F <- table(RawContacts_11_F_6D$hour)/length(unique(c(RawContacts_11_6D$ant1, RawContacts_11_6D$ant2)))
names(TimeCourse6D_F) <- as.numeric(names(TimeCourse6D_F)) - as.numeric(names(TimeCourse6D_F)[1])
TimeCourse7D_F <- table(RawContacts_11_F_7D$hour)/length(unique(c(RawContacts_11_7D$ant1, RawContacts_11_7D$ant2)))
names(TimeCourse7D_F) <- as.numeric(names(TimeCourse7D_F)) - as.numeric(names(TimeCourse7D_F)[1])
TimeCourse8D_F <- table(RawContacts_11_F_8D$hour)/length(unique(c(RawContacts_11_8D$ant1, RawContacts_11_8D$ant2)))
names(TimeCourse8D_F) <- as.numeric(names(TimeCourse8D_F)) - as.numeric(names(TimeCourse8D_F)[1])
TimeCourse9D_F <- table(RawContacts_11_F_9D$hour)/length(unique(c(RawContacts_11_9D$ant1, RawContacts_11_9D$ant2)))

TimeCourse1C_N <- table(RawContacts_11_N_1C$hour)/length(unique(c(RawContacts_11_1C$ant1, RawContacts_11_1C$ant2)))
names(TimeCourse1C_N) <- as.numeric(names(TimeCourse1C_N)) - as.numeric(names(TimeCourse1C_F)[1])
TimeCourse2C_N <- table(RawContacts_11_N_2C$hour)/length(unique(c(RawContacts_11_2C$ant1, RawContacts_11_2C$ant2)))
names(TimeCourse2C_N) <- as.numeric(names(TimeCourse2C_N)) - as.numeric(names(TimeCourse2C_F)[1])
TimeCourse3C_N <- table(RawContacts_11_N_3C$hour)/length(unique(c(RawContacts_11_3C$ant1, RawContacts_11_3C$ant2)))
names(TimeCourse3C_N) <- as.numeric(names(TimeCourse3C_N)) - as.numeric(names(TimeCourse3C_F)[1])
TimeCourse4C_N <- table(RawContacts_11_N_4C$hour)/length(unique(c(RawContacts_11_4C$ant1, RawContacts_11_4C$ant2)))
names(TimeCourse4C_N) <- as.numeric(names(TimeCourse4C_N)) - as.numeric(names(TimeCourse4C_F)[1])
TimeCourse5C_N <- table(RawContacts_11_N_5C$hour)/length(unique(c(RawContacts_11_5C$ant1, RawContacts_11_5C$ant2)))
names(TimeCourse5C_N) <- as.numeric(names(TimeCourse5C_N)) - as.numeric(names(TimeCourse5C_F)[1])
TimeCourse6C_N <- table(RawContacts_11_N_6C$hour)/length(unique(c(RawContacts_11_6C$ant1, RawContacts_11_6C$ant2)))
names(TimeCourse6C_N) <- as.numeric(names(TimeCourse6C_N)) - as.numeric(names(TimeCourse6C_F)[1])
TimeCourse7C_N <- table(RawContacts_11_N_7C$hour)/length(unique(c(RawContacts_11_7C$ant1, RawContacts_11_7C$ant2)))
names(TimeCourse7C_N) <- as.numeric(names(TimeCourse7C_N)) - as.numeric(names(TimeCourse7C_F)[1])
TimeCourse8C_N <- table(RawContacts_11_N_8C$hour)/length(unique(c(RawContacts_11_8C$ant1, RawContacts_11_8C$ant2)))
names(TimeCourse8C_N) <- as.numeric(names(TimeCourse8C_N)) - as.numeric(names(TimeCourse8C_F)[1])
TimeCourse9C_N <- table(RawContacts_11_N_9C$hour)/length(unique(c(RawContacts_11_9C$ant1, RawContacts_11_9C$ant2)))
names(TimeCourse9C_N) <- as.numeric(names(TimeCourse9C_N)) - as.numeric(names(TimeCourse9C_F)[1])
TimeCourse1D_N <- table(RawContacts_11_N_1D$hour)/length(unique(c(RawContacts_11_1D$ant1, RawContacts_11_1D$ant2)))
names(TimeCourse1D_N) <- as.numeric(names(TimeCourse1D_N)) - as.numeric(names(TimeCourse1D_F)[1])
TimeCourse2D_N <- table(RawContacts_11_N_2D$hour)/length(unique(c(RawContacts_11_2D$ant1, RawContacts_11_2D$ant2)))
names(TimeCourse2D_N) <- as.numeric(names(TimeCourse2D_N)) - as.numeric(names(TimeCourse2D_F)[1])
TimeCourse3D_N <- table(RawContacts_11_N_3D$hour)/length(unique(c(RawContacts_11_3D$ant1, RawContacts_11_3D$ant2)))
names(TimeCourse3D_N) <- as.numeric(names(TimeCourse3D_N)) - as.numeric(names(TimeCourse3D_F)[1])
TimeCourse4D_N <- table(RawContacts_11_N_4D$hour)/length(unique(c(RawContacts_11_4D$ant1, RawContacts_11_4D$ant2)))
names(TimeCourse4D_N) <- as.numeric(names(TimeCourse4D_N)) - as.numeric(names(TimeCourse4D_F)[1])
TimeCourse5D_N <- table(RawContacts_11_N_5D$hour)/length(unique(c(RawContacts_11_5D$ant1, RawContacts_11_5D$ant2)))
names(TimeCourse5D_N) <- as.numeric(names(TimeCourse5D_N)) - as.numeric(names(TimeCourse5D_F)[1])
TimeCourse6D_N <- table(RawContacts_11_N_6D$hour)/length(unique(c(RawContacts_11_6D$ant1, RawContacts_11_6D$ant2)))
names(TimeCourse6D_N) <- as.numeric(names(TimeCourse6D_N)) - as.numeric(names(TimeCourse6D_F)[1])
TimeCourse7D_N <- table(RawContacts_11_N_7D$hour)/length(unique(c(RawContacts_11_7D$ant1, RawContacts_11_7D$ant2)))
names(TimeCourse7D_N) <- as.numeric(names(TimeCourse7D_N)) - as.numeric(names(TimeCourse7D_F)[1])
TimeCourse8D_N <- table(RawContacts_11_N_8D$hour)/length(unique(c(RawContacts_11_8D$ant1, RawContacts_11_8D$ant2)))
names(TimeCourse8D_N) <- as.numeric(names(TimeCourse8D_N)) - as.numeric(names(TimeCourse8D_F)[1])
TimeCourse9D_N <- table(RawContacts_11_N_9D$hour)/length(unique(c(RawContacts_11_9D$ant1, RawContacts_11_9D$ant2)))
names(TimeCourse9D_N) <- as.numeric(names(TimeCourse9D_N)) - as.numeric(names(TimeCourse9D_F)[1])
```

### Plot time-series

``` r
setwd(FIGSDIR)
jpeg('Timecourse.jpg', width=9000, height=2000, unit='px')
par(mfrow = c(2,9), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
plot(y = TimeCourse1C_F, x = as.numeric(names(TimeCourse1C_F)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,90), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,90), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse1C_F~ as.numeric(names(TimeCourse1C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse1D_F~ as.numeric(names(TimeCourse1D_F)), col="blue" , lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse2C_F, x = as.numeric(names(TimeCourse2C_F)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,90), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,90), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse2C_F~ as.numeric(names(TimeCourse2C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse2D_F~ as.numeric(names(TimeCourse2D_F)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse3C_F, x = as.numeric(names(TimeCourse3C_F)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,90), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,90), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse3C_F~ as.numeric(names(TimeCourse3C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse3D_F~ as.numeric(names(TimeCourse3D_F)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse4C_F, x = as.numeric(names(TimeCourse4C_F)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,90), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,90), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse4C_F~ as.numeric(names(TimeCourse4C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse4D_F~ as.numeric(names(TimeCourse4D_F)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse5C_F, x = as.numeric(names(TimeCourse5C_F)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,90), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,90), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse5C_F~ as.numeric(names(TimeCourse5C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse5D_F~ as.numeric(names(TimeCourse5D_F)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse6C_F, x = as.numeric(names(TimeCourse6C_F)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,90), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,90), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse6C_F~ as.numeric(names(TimeCourse6C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse6D_F~ as.numeric(names(TimeCourse6D_F)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse7C_F, x = as.numeric(names(TimeCourse7C_F)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,90), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,90), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse7C_F~ as.numeric(names(TimeCourse7C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse7D_F~ as.numeric(names(TimeCourse7D_F)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse8C_F, x = as.numeric(names(TimeCourse8C_F)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,90), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,90), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse8D_F~ as.numeric(names(TimeCourse8D_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse8C_F~ as.numeric(names(TimeCourse8C_F)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse9C_F, x = as.numeric(names(TimeCourse9C_F)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,90), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,90), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse9C_F~ as.numeric(names(TimeCourse9C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse9D_F~ as.numeric(names(TimeCourse9D_F)), col="blue", lwd=8 , pch=19 , type="b" )

plot(y = TimeCourse1C_N, x = as.numeric(names(TimeCourse1C_N)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,4), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,4), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse1C_N~ as.numeric(names(TimeCourse1C_N)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse1D_N~ as.numeric(names(TimeCourse1D_N)), col="blue" , lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse2C_N, x = as.numeric(names(TimeCourse2C_N)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,4), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,4), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse2C_N~ as.numeric(names(TimeCourse2C_N)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse2D_N~ as.numeric(names(TimeCourse2D_N)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse3C_N, x = as.numeric(names(TimeCourse3C_N)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,4), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,4), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse3C_N~ as.numeric(names(TimeCourse3C_N)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse3D_N~ as.numeric(names(TimeCourse3D_N)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse4C_N, x = as.numeric(names(TimeCourse4C_N)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,4), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,4), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse4C_N~ as.numeric(names(TimeCourse4C_N)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse4D_N~ as.numeric(names(TimeCourse4D_N)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse5C_N, x = as.numeric(names(TimeCourse5C_N)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,4), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,4), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse5C_N~ as.numeric(names(TimeCourse5C_N)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse5D_N~ as.numeric(names(TimeCourse5D_N)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse6C_N, x = as.numeric(names(TimeCourse6C_N)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,4), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,4), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse6C_N~ as.numeric(names(TimeCourse6C_N)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse6D_N~ as.numeric(names(TimeCourse6D_N)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse7C_N, x = as.numeric(names(TimeCourse7C_N)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,4), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,4), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse7C_N~ as.numeric(names(TimeCourse7C_N)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse7D_N~ as.numeric(names(TimeCourse7D_N)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse8C_N, x = as.numeric(names(TimeCourse8C_N)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,4), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,4), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse8D_N~ as.numeric(names(TimeCourse8D_N)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse8C_N~ as.numeric(names(TimeCourse8C_N)), col="blue", lwd=8 , pch=19 , type="b" )
plot(y = TimeCourse9C_N, x = as.numeric(names(TimeCourse9C_N)), pch = 16, xlab = "Hour", ylab = "Interactions",
     ylim = c(0,4), xlim = c(0,150), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, 0, 8, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, 0, 32, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, 0, 56, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, 0, 80, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, 0, 104, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, 0, 128, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, 0, 158, 90,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(0,4), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines( TimeCourse9C_N~ as.numeric(names(TimeCourse9C_N)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines( TimeCourse9D_N~ as.numeric(names(TimeCourse9D_N)), col="blue", lwd=8 , pch=19 , type="b" )
dev.off()
```

### Repeat using residuals with respect to subcolony size rather than absolute values

``` r
TimeCourse1C_Ft <- rbind(TimeCourse1C_F, length(unique(c(RawContacts_11_1C$ant1, RawContacts_11_1C$ant2))), "1C")
TimeCourse2C_Ft <- rbind(TimeCourse2C_F, length(unique(c(RawContacts_11_2C$ant1, RawContacts_11_2C$ant2))), "2C")
TimeCourse3C_Ft <- rbind(TimeCourse3C_F, length(unique(c(RawContacts_11_3C$ant1, RawContacts_11_3C$ant2))), "3C")
TimeCourse4C_Ft <- rbind(TimeCourse4C_F, length(unique(c(RawContacts_11_4C$ant1, RawContacts_11_4C$ant2))), "4C")
TimeCourse5C_Ft <- rbind(TimeCourse5C_F, length(unique(c(RawContacts_11_5C$ant1, RawContacts_11_5C$ant2))), "5C")
TimeCourse6C_Ft <- rbind(TimeCourse6C_F, length(unique(c(RawContacts_11_6C$ant1, RawContacts_11_6C$ant2))), "6C")
TimeCourse7C_Ft <- rbind(TimeCourse7C_F, length(unique(c(RawContacts_11_7C$ant1, RawContacts_11_7C$ant2))), "7C")
TimeCourse8C_Ft <- rbind(TimeCourse8C_F, length(unique(c(RawContacts_11_8C$ant1, RawContacts_11_8C$ant2))), "8C")
TimeCourse9C_Ft <- rbind(TimeCourse9C_F, length(unique(c(RawContacts_11_9C$ant1, RawContacts_11_9C$ant2))), "9C")
TimeCourse1D_Ft <- rbind(TimeCourse1D_F, length(unique(c(RawContacts_11_1D$ant1, RawContacts_11_1D$ant2))), "1D")
TimeCourse2D_Ft <- rbind(TimeCourse2D_F, length(unique(c(RawContacts_11_2D$ant1, RawContacts_11_2D$ant2))), "2D")
TimeCourse3D_Ft <- rbind(TimeCourse3D_F, length(unique(c(RawContacts_11_3D$ant1, RawContacts_11_3D$ant2))), "3D")
TimeCourse4D_Ft <- rbind(TimeCourse4D_F, length(unique(c(RawContacts_11_4D$ant1, RawContacts_11_4D$ant2))), "4D")
TimeCourse5D_Ft <- rbind(TimeCourse5D_F, length(unique(c(RawContacts_11_5D$ant1, RawContacts_11_5D$ant2))), "5D")
TimeCourse6D_Ft <- rbind(TimeCourse6D_F, length(unique(c(RawContacts_11_6D$ant1, RawContacts_11_6D$ant2))), "6D")
TimeCourse7D_Ft <- rbind(TimeCourse7D_F, length(unique(c(RawContacts_11_7D$ant1, RawContacts_11_7D$ant2))), "7D")
TimeCourse8D_Ft <- rbind(TimeCourse8D_F, length(unique(c(RawContacts_11_8D$ant1, RawContacts_11_8D$ant2))), "8D")
TimeCourse9D_Ft <- rbind(TimeCourse9D_F, length(unique(c(RawContacts_11_9D$ant1, RawContacts_11_9D$ant2))), "9D")

ALL <-rbind(t(TimeCourse1C_Ft),
            t(TimeCourse2C_Ft),
            t(TimeCourse3C_Ft),
            t(TimeCourse4C_Ft),
            t(TimeCourse5C_Ft),
            t(TimeCourse6C_Ft),
            t(TimeCourse7C_Ft),
            t(TimeCourse8C_Ft),
            t(TimeCourse9C_Ft),
            t(TimeCourse1D_Ft),
            t(TimeCourse2D_Ft),
            t(TimeCourse3D_Ft),
            t(TimeCourse4D_Ft),
            t(TimeCourse5D_Ft),
            t(TimeCourse6D_Ft),
            t(TimeCourse7D_Ft),
            t(TimeCourse8D_Ft),
            t(TimeCourse9D_Ft))

ALL <- as.data.frame(ALL, row.names = 1)
ALL$TimeCourse1C_F <- as.numeric(ALL$TimeCourse1C_F)
ALL$V2 <- as.numeric(ALL$V2)
colnames(ALL) <- c("HHPerHour", "N_Bee", "Rep")
ALL_AV <- aggregate(.~Rep, data=ALL, mean)

Averages <- c()
for (row in 1:nrow(ALL)){
  rep       <- ALL[row,3]
  Averages  <- c(Averages, ALL_AV[ALL_AV$Rep == rep,2])
}
ALL$Averages <- Averages
ALL$residuals <- ALL$HHPerHour - ALL$Averages
```

### Plot the same plot as previously but using residuals

``` r
jpeg('TimecourseResiduals.jpg', width=9000, height=2000, unit='px')
par(mfrow = c(1,9), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
plot(y = ALL[ALL$Rep == "1C",]$residuals, x = as.numeric(names(TimeCourse1C_F)), pch = 16, xlab = "Hour", ylab = "Interaction Residuals",
     ylim = c(-25,35), xlim = c(0,160), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, -30, 8, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, -30, 32, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, -30, 56, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, -30, 80, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, -30, 104, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, -30, 128, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, -30, 158, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(-25,25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(ALL[ALL$Rep == "1C",]$residuals ~ as.numeric(names(TimeCourse1C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines(ALL[ALL$Rep == "1D",]$residuals ~ as.numeric(names(TimeCourse1D_F)), col = "blue", lwd=8 , pch=19 , type="b" )

plot(y = ALL[ALL$Rep == "2C",]$residuals, x = as.numeric(names(TimeCourse2C_F)), pch = 16, xlab = "Hour", ylab = "Interaction Residuals",
     ylim = c(-25,35), xlim = c(0,160), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, -30, 8, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, -30, 32, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, -30, 56, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, -30, 80, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, -30, 104, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, -30, 128, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, -30, 158, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(-25,25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(ALL[ALL$Rep == "2C",]$residuals ~ as.numeric(names(TimeCourse2C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines(ALL[ALL$Rep == "2D",]$residuals ~ as.numeric(names(TimeCourse2D_F)), col = "blue", lwd=8 , pch=19 , type="b" )

plot(y = ALL[ALL$Rep == "3C",]$residuals, x = as.numeric(names(TimeCourse3C_F)), pch = 16, xlab = "Hour", ylab = "Interaction Residuals",
     ylim = c(-25,35), xlim = c(0,160), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, -30, 8, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, -30, 32, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, -30, 56, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, -30, 80, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, -30, 104, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, -30, 128, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, -30, 158, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(-25,25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(ALL[ALL$Rep == "3C",]$residuals ~ as.numeric(names(TimeCourse3C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines(ALL[ALL$Rep == "3D",]$residuals ~ as.numeric(names(TimeCourse3D_F)), col = "blue", lwd=8 , pch=19 , type="b" )

plot(y = ALL[ALL$Rep == "4C",]$residuals, x = as.numeric(names(TimeCourse4C_F)), pch = 16, xlab = "Hour", ylab = "Interaction Residuals",
     ylim = c(-25,35), xlim = c(0,160), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, -30, 8, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, -30, 32, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, -30, 56, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, -30, 80, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, -30, 104, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, -30, 128, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, -30, 158, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(-25,25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(ALL[ALL$Rep == "4C",]$residuals ~ as.numeric(names(TimeCourse4C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines(ALL[ALL$Rep == "4D",]$residuals ~ as.numeric(names(TimeCourse4D_F)), col = "blue", lwd=8 , pch=19 , type="b" )

plot(y = ALL[ALL$Rep == "5C",]$residuals, x = as.numeric(names(TimeCourse5C_F)), pch = 16, xlab = "Hour", ylab = "Interaction Residuals",
     ylim = c(-25,35), xlim = c(0,160), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, -30, 8, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, -30, 32, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, -30, 56, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, -30, 80, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, -30, 104, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, -30, 128, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, -30, 158, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(-25,25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(ALL[ALL$Rep == "5C",]$residuals ~ as.numeric(names(TimeCourse5C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines(ALL[ALL$Rep == "5D",]$residuals ~ as.numeric(names(TimeCourse5D_F)), col = "blue", lwd=8 , pch=19 , type="b" )

plot(y = ALL[ALL$Rep == "6C",]$residuals, x = as.numeric(names(TimeCourse6C_F)), pch = 16, xlab = "Hour", ylab = "Interaction Residuals",
     ylim = c(-25,35), xlim = c(0,160), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, -30, 8, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, -30, 32, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, -30, 56, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, -30, 80, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, -30, 104, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, -30, 128, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, -30, 158, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(-25,25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(ALL[ALL$Rep == "6C",]$residuals ~ as.numeric(names(TimeCourse6C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines(ALL[ALL$Rep == "6D",]$residuals ~ as.numeric(names(TimeCourse6D_F)), col = "blue", lwd=8 , pch=19 , type="b" )

plot(y = ALL[ALL$Rep == "7C",]$residuals, x = as.numeric(names(TimeCourse7C_F)), pch = 16, xlab = "Hour", ylab = "Interaction Residuals",
     ylim = c(-25,35), xlim = c(0,160), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, -30, 8, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, -30, 32, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, -30, 56, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, -30, 80, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, -30, 104, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, -30, 128, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, -30, 158, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(-25,25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(ALL[ALL$Rep == "7C",]$residuals ~ as.numeric(names(TimeCourse7C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines(ALL[ALL$Rep == "7D",]$residuals ~ as.numeric(names(TimeCourse7D_F)), col = "blue", lwd=8 , pch=19 , type="b" )

plot(y = ALL[ALL$Rep == "8C",]$residuals, x = as.numeric(names(TimeCourse8C_F)), pch = 16, xlab = "Hour", ylab = "Interaction Residuals",
     ylim = c(-25,35), xlim = c(0,160), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, -30, 8, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, -30, 32, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, -30, 56, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, -30, 80, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, -30, 104, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, -30, 128, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, -30, 158, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(-25,25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(ALL[ALL$Rep == "8C",]$residuals ~ as.numeric(names(TimeCourse8C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines(ALL[ALL$Rep == "8D",]$residuals ~ as.numeric(names(TimeCourse8D_F)), col = "blue", lwd=8 , pch=19 , type="b" )

plot(y = ALL[ALL$Rep == "9C",]$residuals, x = as.numeric(names(TimeCourse9C_F)), pch = 16, xlab = "Hour", ylab = "Interaction Residuals",
     ylim = c(-25,35), xlim = c(0,160), cex.axis = 10, yaxt="n", xaxt="n", cex.lab = 10)
rect(0, -30, 8, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(20, -30, 32, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(44, -30, 56, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(68, -30, 80, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(92, -30, 104, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(116, -30, 128, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
rect(140, -30, 158, 35,col = "gray", border = NULL, lty = par("lty"), lwd = par("lwd"))
axis(2, at = c(-25,25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(ALL[ALL$Rep == "9C",]$residuals ~ as.numeric(names(TimeCourse9C_F)), col=rgb(0.8,0.4,0.1,0.7) , lwd=8 , pch=19 , type="b" )
lines(ALL[ALL$Rep == "9D",]$residuals ~ as.numeric(names(TimeCourse9D_F)), col = "blue", lwd=8 , pch=19 , type="b" )
dev.off()
```

### Construct and plot time-course data for individual speed and social contacts across colonies

``` r
# Reset start dates so all begin on the 1st
RawContacts_11_1C$start <- RawContacts_11_1C$start - days(13)
RawContacts_11_2C$start <- RawContacts_11_2C$start - days(27)
RawContacts_11_3C$start <- RawContacts_11_3C$start - days(3)
RawContacts_11_4C$start <- RawContacts_11_4C$start - days(10)
RawContacts_11_5C$start <- RawContacts_11_5C$start - days(24)
RawContacts_11_6C$start <- RawContacts_11_6C$start - days(1)
RawContacts_11_7C$start <- RawContacts_11_7C$start - days(15)
RawContacts_11_8C$start <- RawContacts_11_8C$start - days(22)
RawContacts_11_9C$start <- RawContacts_11_9C$start - days(16)
RawContacts_11_1D$start <- RawContacts_11_1D$start - days(13)
RawContacts_11_2D$start <- RawContacts_11_2D$start - days(27)
RawContacts_11_3D$start <- RawContacts_11_3D$start - days(3)
RawContacts_11_4D$start <- RawContacts_11_4D$start - days(10)
RawContacts_11_5D$start <- RawContacts_11_5D$start - days(24)
RawContacts_11_6D$start <- RawContacts_11_6D$start - days(1)
RawContacts_11_7D$start <- RawContacts_11_7D$start - days(15)
RawContacts_11_8D$start <- RawContacts_11_8D$start - days(22)
RawContacts_11_9D$start <- RawContacts_11_9D$start - days(16)
```

### Reset down the start hour and sum per hour

``` r
RawContacts_11_1C$start <- format(floor_date(RawContacts_11_1C$start, unit="hour"), format="%d %H:%M")
RawContacts_11_2C$start <- format(floor_date(RawContacts_11_2C$start, unit="hour"), format="%d %H:%M")
RawContacts_11_3C$start <- format(floor_date(RawContacts_11_3C$start, unit="hour"), format="%d %H:%M")
RawContacts_11_4C$start <- format(floor_date(RawContacts_11_4C$start, unit="hour"), format="%d %H:%M")
RawContacts_11_5C$start <- format(floor_date(RawContacts_11_5C$start, unit="hour"), format="%d %H:%M")
RawContacts_11_6C$start <- format(floor_date(RawContacts_11_6C$start, unit="hour"), format="%d %H:%M")
RawContacts_11_7C$start <- format(floor_date(RawContacts_11_7C$start, unit="hour"), format="%d %H:%M")
RawContacts_11_8C$start <- format(floor_date(RawContacts_11_8C$start, unit="hour"), format="%d %H:%M")
RawContacts_11_9C$start <- format(floor_date(RawContacts_11_9C$start, unit="hour"), format="%d %H:%M")
RawContacts_11_1D$start <- format(floor_date(RawContacts_11_1D$start, unit="hour"), format="%d %H:%M")
RawContacts_11_2D$start <- format(floor_date(RawContacts_11_2D$start, unit="hour"), format="%d %H:%M")
RawContacts_11_3D$start <- format(floor_date(RawContacts_11_3D$start, unit="hour"), format="%d %H:%M")
RawContacts_11_4D$start <- format(floor_date(RawContacts_11_4D$start, unit="hour"), format="%d %H:%M")
RawContacts_11_5D$start <- format(floor_date(RawContacts_11_5D$start, unit="hour"), format="%d %H:%M")
RawContacts_11_6D$start <- format(floor_date(RawContacts_11_6D$start, unit="hour"), format="%d %H:%M")
RawContacts_11_7D$start <- format(floor_date(RawContacts_11_7D$start, unit="hour"), format="%d %H:%M")
RawContacts_11_8D$start <- format(floor_date(RawContacts_11_8D$start, unit="hour"), format="%d %H:%M")
RawContacts_11_9D$start <- format(floor_date(RawContacts_11_9D$start, unit="hour"), format="%d %H:%M")
TimeCourse1C <- table(RawContacts_11_1C$start)/length(unique(c(RawContacts_11_1C$ant1, RawContacts_11_1C$ant2)))
TimeCourse2C <- table(RawContacts_11_2C$start)/length(unique(c(RawContacts_11_2C$ant1, RawContacts_11_2C$ant2)))
TimeCourse3C <- table(RawContacts_11_3C$start)/length(unique(c(RawContacts_11_3C$ant1, RawContacts_11_3C$ant2)))
TimeCourse4C <- table(RawContacts_11_4C$start)/length(unique(c(RawContacts_11_4C$ant1, RawContacts_11_4C$ant2)))
TimeCourse5C <- table(RawContacts_11_5C$start)/length(unique(c(RawContacts_11_5C$ant1, RawContacts_11_5C$ant2)))
TimeCourse6C <- table(RawContacts_11_6C$start)/length(unique(c(RawContacts_11_6C$ant1, RawContacts_11_6C$ant2)))
TimeCourse7C <- table(RawContacts_11_7C$start)/length(unique(c(RawContacts_11_7C$ant1, RawContacts_11_7C$ant2)))
TimeCourse8C <- table(RawContacts_11_8C$start)/length(unique(c(RawContacts_11_8C$ant1, RawContacts_11_8C$ant2)))
TimeCourse9C <- table(RawContacts_11_9C$start)/length(unique(c(RawContacts_11_9C$ant1, RawContacts_11_9C$ant2)))
TimeCourse1D <- table(RawContacts_11_1D$start)/length(unique(c(RawContacts_11_1D$ant1, RawContacts_11_1D$ant2)))
TimeCourse2D <- table(RawContacts_11_2D$start)/length(unique(c(RawContacts_11_2D$ant1, RawContacts_11_2D$ant2)))
TimeCourse3D <- table(RawContacts_11_3D$start)/length(unique(c(RawContacts_11_3D$ant1, RawContacts_11_3D$ant2)))
TimeCourse4D <- table(RawContacts_11_4D$start)/length(unique(c(RawContacts_11_4D$ant1, RawContacts_11_4D$ant2)))
TimeCourse5D <- table(RawContacts_11_5D$start)/length(unique(c(RawContacts_11_5D$ant1, RawContacts_11_5D$ant2)))
TimeCourse6D <- table(RawContacts_11_6D$start)/length(unique(c(RawContacts_11_6D$ant1, RawContacts_11_6D$ant2)))
TimeCourse7D <- table(RawContacts_11_7D$start)/length(unique(c(RawContacts_11_7D$ant1, RawContacts_11_7D$ant2)))
TimeCourse8D <- table(RawContacts_11_8D$start)/length(unique(c(RawContacts_11_8D$ant1, RawContacts_11_8D$ant2)))
TimeCourse9D <- table(RawContacts_11_9D$start)/length(unique(c(RawContacts_11_9D$ant1, RawContacts_11_9D$ant2)))
```

### Compile into one dataframe per treatment and average

``` r
StartHours_DF_C <- data.frame(matrix(0, ncol = 9, nrow = 152))
rownames(StartHours_DF_C) <- names(TimeCourse1C)
StartHours_DF_D <- data.frame(matrix(0, ncol = 9, nrow = 152))
rownames(StartHours_DF_D) <- names(TimeCourse1C)
for (t in 1:length(TimeCourse1C)){
  timepoint <- names(TimeCourse1C)[t]
  StartHours_DF_C[which(rownames(StartHours_DF_C) == timepoint),1] <- as.numeric(TimeCourse1C)[t]
}
for (t in 1:length(TimeCourse2C)){
  timepoint <- names(TimeCourse2C)[t]
  StartHours_DF_C[which(rownames(StartHours_DF_C) == timepoint),2] <- as.numeric(TimeCourse2C)[t]
}
for (t in 1:length(TimeCourse3C)){
  timepoint <- names(TimeCourse3C)[t]
  StartHours_DF_C[which(rownames(StartHours_DF_C) == timepoint),3] <- as.numeric(TimeCourse3C)[t]
}
for (t in 1:length(TimeCourse4C)){
  timepoint <- names(TimeCourse4C)[t]
  StartHours_DF_C[which(rownames(StartHours_DF_C) == timepoint),4] <- as.numeric(TimeCourse4C)[t]
}
for (t in 1:length(TimeCourse5C)){
  timepoint <- names(TimeCourse5C)[t]
  StartHours_DF_C[which(rownames(StartHours_DF_C) == timepoint),5] <- as.numeric(TimeCourse5C)[t]
}
for (t in 1:length(TimeCourse6C)){
  timepoint <- names(TimeCourse6C)[t]
  StartHours_DF_C[which(rownames(StartHours_DF_C) == timepoint),6] <- as.numeric(TimeCourse6C)[t]
}
for (t in 1:length(TimeCourse7C)){
  timepoint <- names(TimeCourse7C)[t]
  StartHours_DF_C[which(rownames(StartHours_DF_C) == timepoint),7] <- as.numeric(TimeCourse7C)[t]
}
for (t in 1:length(TimeCourse8C)){
  timepoint <- names(TimeCourse8C)[t]
  StartHours_DF_C[which(rownames(StartHours_DF_C) == timepoint),8] <- as.numeric(TimeCourse8C)[t]
}
for (t in 1:length(TimeCourse9C)){
  timepoint <- names(TimeCourse9C)[t]
  StartHours_DF_C[which(rownames(StartHours_DF_C) == timepoint),9] <- as.numeric(TimeCourse9C)[t]
}
for (t in 1:length(TimeCourse1D)){
  timepoint <- names(TimeCourse1D)[t]
  StartHours_DF_D[which(rownames(StartHours_DF_D) == timepoint),1] <- as.numeric(TimeCourse1D)[t]
}
for (t in 1:length(TimeCourse2D)){
  timepoint <- names(TimeCourse2D)[t]
  StartHours_DF_D[which(rownames(StartHours_DF_D) == timepoint),2] <- as.numeric(TimeCourse2D)[t]
}
for (t in 1:length(TimeCourse3D)){
  timepoint <- names(TimeCourse3D)[t]
  StartHours_DF_D[which(rownames(StartHours_DF_D) == timepoint),3] <- as.numeric(TimeCourse3D)[t]
}
for (t in 1:length(TimeCourse4D)){
  timepoint <- names(TimeCourse4D)[t]
  StartHours_DF_D[which(rownames(StartHours_DF_D) == timepoint),4] <- as.numeric(TimeCourse4D)[t]
}
for (t in 1:length(TimeCourse5D)){
  timepoint <- names(TimeCourse5D)[t]
  StartHours_DF_D[which(rownames(StartHours_DF_D) == timepoint),5] <- as.numeric(TimeCourse5D)[t]
}
for (t in 1:length(TimeCourse6D)){
  timepoint <- names(TimeCourse6D)[t]
  StartHours_DF_D[which(rownames(StartHours_DF_D) == timepoint),6] <- as.numeric(TimeCourse6D)[t]
}
for (t in 1:length(TimeCourse7D)){
  timepoint <- names(TimeCourse7D)[t]
  StartHours_DF_D[which(rownames(StartHours_DF_D) == timepoint),7] <- as.numeric(TimeCourse7D)[t]
}
for (t in 1:length(TimeCourse8D)){
  timepoint <- names(TimeCourse8D)[t]
  StartHours_DF_D[which(rownames(StartHours_DF_D) == timepoint),8] <- as.numeric(TimeCourse8D)[t]
}
for (t in 1:length(TimeCourse9D)){
  timepoint <- names(TimeCourse9D)[t]
  StartHours_DF_D[which(rownames(StartHours_DF_D) == timepoint),9] <- as.numeric(TimeCourse9D)[t]
}
StartHours_DF_C[StartHours_DF_C == 0] <- NA
StartHours_DF_D[StartHours_DF_D == 0] <- NA
Hourly_total_C <- as.numeric(rowMeans(StartHours_DF_C, na.rm = TRUE))
Hourly_total_D <- as.numeric(rowMeans(StartHours_DF_D, na.rm = TRUE))
```

### Import speed data

``` r
setwd(OUTDIR)
Hourly_Speed_C <- read.csv("HourlySpeed_C.csv")[-153,]
Hourly_Speed_D <- read.csv("HourlySpeed_D.csv")[-153,]
```

### Plot Summary timecourses

``` r
setwd(FIGSDIR)
pdf('OverallTimecourse.pdf', width=60, height=25)
par(mfrow = c(2,1), mar = c(10,16,5,0), oma = c(0,0,0,0), bty="n", mgp = c(0,9,0), family = "sans")
plot(seq(1,152,1), Hourly_total_C, type = "l", col = "#00BFC4",
     xlab = "", ylab = "", lwd = 15, yaxt="n", xaxt="n",
     cex.main = 5, cex.axis = 5, cex.lab = 1, xlim = c(1,152))
rect(0, -5, 8, 38, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))
rect(20, -5, 32, 38, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))
rect(44, -5, 56, 38, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))
rect(68, -5, 80, 38, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))
rect(92, -5, 104, 38, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))
rect(116, -5, 128, 38, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))
rect(140, -5, 151, 38, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))

lines(seq(1,152,1), Hourly_total_D, col = "#C77CFF", lwd = 15)
lines(seq(1,152,1), Hourly_total_C, col = "#00BFC4", lwd = 15)
axis(2, at = c(12,34), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(1,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
mtext("Contacts", side=2, line=2.2, cex=10)

plot(Hourly_Speed_D$speed ~ seq(1,length(Hourly_Speed_D$speed),1), type = "l", lwd = 15, col = "#00BFC4",
     ylab = "", xlab = "", yaxt="n", xaxt="n",
     cex.main = 5, cex.axis = 5, cex.lab = 1)
rect(0, -50, 8, 180, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))
rect(20, -50, 32, 180, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))
rect(44, -50, 56, 180, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))
rect(68, -50, 80, 180, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))
rect(92, -50, 104, 180, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))
rect(116, -50, 128, 180, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))
rect(140, -50, 151, 180, col = "lightgray", border = "lightgray", lty = par("lty"), lwd = par("lwd"))

lines(Hourly_Speed_D$speed ~ seq(1,length(Hourly_Speed_D$speed),1), lwd = 15, col = "#C77CFF")
lines(Hourly_Speed_C$speed ~ seq(1,length(Hourly_Speed_C$speed),1), lwd = 15, col = "#00BFC4")
axis(2, at = c(40,160), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(1,150), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
mtext("Hour", side=1, line=6, cex=10)
mtext("Speed", side=2, line=2.2, cex=10)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] lubridate_1.8.0  forcats_0.5.1    dplyr_1.0.7      purrr_0.3.4     
    ##  [5] readr_2.1.1      tidyr_1.1.4      tibble_3.1.6     tidyverse_1.3.1 
    ##  [9] pals_1.7         igraph_1.2.9     stringr_1.4.0    ggfortify_0.4.13
    ## [13] ggplot2_3.3.5    gapminder_0.3.0  jpeg_0.1-9      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.7       assertthat_0.2.1 digest_0.6.29    utf8_1.2.2      
    ##  [5] R6_2.5.1         cellranger_1.1.0 backports_1.4.1  reprex_2.0.1    
    ##  [9] evaluate_0.14    httr_1.4.2       pillar_1.6.4     rlang_0.4.12    
    ## [13] readxl_1.3.1     rstudioapi_0.13  rmarkdown_2.11   munsell_0.5.0   
    ## [17] broom_0.7.10     compiler_4.1.2   modelr_0.1.8     xfun_0.28       
    ## [21] pkgconfig_2.0.3  htmltools_0.5.2  tidyselect_1.1.1 gridExtra_2.3   
    ## [25] fansi_0.5.0      crayon_1.4.2     tzdb_0.2.0       dbplyr_2.1.1    
    ## [29] withr_2.4.3      grid_4.1.2       jsonlite_1.7.2   gtable_0.3.0    
    ## [33] lifecycle_1.0.1  DBI_1.1.1        magrittr_2.0.1   scales_1.1.1    
    ## [37] cli_3.1.0        stringi_1.7.6    mapproj_1.2.7    fs_1.5.1        
    ## [41] xml2_1.3.3       ellipsis_0.3.2   generics_0.1.1   vctrs_0.3.8     
    ## [45] tools_4.1.2      dichromat_2.0-0  glue_1.5.1       maps_3.4.0      
    ## [49] hms_1.1.1        fastmap_1.1.0    yaml_2.2.1       colorspace_2.0-2
    ## [53] rvest_1.0.2      knitr_1.36       haven_2.4.3
