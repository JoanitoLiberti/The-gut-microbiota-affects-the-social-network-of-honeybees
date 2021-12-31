Main analysis of automated behavioral tracking data
================
Tomas Kay

## Comparison of the behavior of colonized (C) and microbiota-depleted (D) bees

## Study design is paired: One C and one D subcolony for each of 9 hives

## This first script processed raw data saves summary stats for each subcolony

### Set directories and load packages

``` r
MAINDIR <- "/Volumes/Lacie/BM"
INDIR   <- paste(MAINDIR, "RawData", sep = "/")
OUTDIR  <- paste(MAINDIR, "ProcessedData", sep = "/")
FIGSDIR <- paste(MAINDIR, "Figures", sep = "/")


library('jpeg'); library('gapminder'); library('ggfortify'); library('stringr'); library('igraph')
library('pals'); library('OneR'); library('lubridate')
```

### Import all tracking data

``` r
### For each sub-colony there are two input .RDS files: one with contact data and the other with trajectory data
setwd(INDIR)
files <- list.files(pattern=".rds")
for (file in files){
  assign(sub("\\..*", "", file),readRDS(file))
}
```

### Set start and stop times

``` r
# All of these data files are trimmed to begin at 12 pm the night tracking was initiated (always on a Wednesday)
# and to end 152 hours later (8am the following Wednesday).
# Start and stop times for each replicate are as follows:
# 1C = 2020-05-13 21:54:00; 1D = 2020-05-13 21:07:00: END = 2020-05-20
# 2C = 2020-05-27 19:13:00; 2D = 2020-05-27 18:52:00: END = 2020-06-03
# 3C = 2020-06-03 16:50:00; 3D = 2020-06-03 16:48:00: END = 2020-06-10
# 4C = 2020-06-10 17:53:00; 4D = 2020-06-10 17:36:00: END = 2020-06-17
# 5C = 2020-06-24 18:33:00; 5D = 2020-06-24 18:31:00: END = 2020-07-01
# 6C = 2020-07-01 18:04:00; 6D = 2020-07-01 17:52:00: END = 2020-07-08
# 7C = 2020-07-15 17:56:00; 7D = 2020-07-15 17:44:00: END = 2020-07-22
# 8C = 2020-07-22 20:20:00; 8D = 2020-07-22 20:18:00: END = 2020-07-29
# 9C = 2020-09-16 18:38:00; 9D = 2020-09-16 18:42:00: END = 2020-09-23

# Fix start times
AntMinuteMeans_ALL_1C <- AntMinuteMeans_ALL_1C[AntMinuteMeans_ALL_1C$time > "2020-05-13 24:00:00",]
AntMinuteMeans_ALL_1D <- AntMinuteMeans_ALL_1D[AntMinuteMeans_ALL_1D$time > "2020-05-13 24:00:00",]
AntMinuteMeans_ALL_2C <- AntMinuteMeans_ALL_2C[AntMinuteMeans_ALL_2C$time > "2020-05-27 24:00:00",]
AntMinuteMeans_ALL_2D <- AntMinuteMeans_ALL_2D[AntMinuteMeans_ALL_2D$time > "2020-05-27 24:00:00",]
AntMinuteMeans_ALL_3C <- AntMinuteMeans_ALL_3C[AntMinuteMeans_ALL_3C$time > "2020-06-03 24:00:00",]
AntMinuteMeans_ALL_3D <- AntMinuteMeans_ALL_3D[AntMinuteMeans_ALL_3D$time > "2020-06-03 24:00:00",]
AntMinuteMeans_ALL_4C <- AntMinuteMeans_ALL_4C[AntMinuteMeans_ALL_4C$time > "2020-06-10 24:00:00",]
AntMinuteMeans_ALL_4D <- AntMinuteMeans_ALL_4D[AntMinuteMeans_ALL_4D$time > "2020-06-10 24:00:00",]
AntMinuteMeans_ALL_5C <- AntMinuteMeans_ALL_5C[AntMinuteMeans_ALL_5C$time > "2020-06-24 24:00:00",]
AntMinuteMeans_ALL_5D <- AntMinuteMeans_ALL_5D[AntMinuteMeans_ALL_5D$time > "2020-06-24 24:00:00",]
AntMinuteMeans_ALL_6C <- AntMinuteMeans_ALL_6C[AntMinuteMeans_ALL_6C$time > "2020-07-01 24:00:00",]
AntMinuteMeans_ALL_6D <- AntMinuteMeans_ALL_6D[AntMinuteMeans_ALL_6D$time > "2020-07-01 24:00:00",]
AntMinuteMeans_ALL_7C <- AntMinuteMeans_ALL_7C[AntMinuteMeans_ALL_7C$time > "2020-07-15 24:00:00",]
AntMinuteMeans_ALL_7D <- AntMinuteMeans_ALL_7D[AntMinuteMeans_ALL_7D$time > "2020-07-15 24:00:00",]
AntMinuteMeans_ALL_8C <- AntMinuteMeans_ALL_8C[AntMinuteMeans_ALL_8C$time > "2020-07-22 24:00:00",]
AntMinuteMeans_ALL_8D <- AntMinuteMeans_ALL_8D[AntMinuteMeans_ALL_8D$time > "2020-07-22 24:00:00",]
AntMinuteMeans_ALL_9C <- AntMinuteMeans_ALL_9C[AntMinuteMeans_ALL_9C$time > "2020-09-16 24:00:00",]
AntMinuteMeans_ALL_9D <- AntMinuteMeans_ALL_9D[AntMinuteMeans_ALL_9D$time > "2020-09-16 24:00:00",]
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

# Fix stop times
AntMinuteMeans_ALL_1C <- AntMinuteMeans_ALL_1C[AntMinuteMeans_ALL_1C$time < "2020-05-20 08:00:00",]
AntMinuteMeans_ALL_1D <- AntMinuteMeans_ALL_1D[AntMinuteMeans_ALL_1D$time < "2020-05-20 08:00:00",]
AntMinuteMeans_ALL_2C <- AntMinuteMeans_ALL_2C[AntMinuteMeans_ALL_2C$time < "2020-06-03 08:00:00",]
AntMinuteMeans_ALL_2D <- AntMinuteMeans_ALL_2D[AntMinuteMeans_ALL_2D$time < "2020-06-03 08:00:00",]
AntMinuteMeans_ALL_3C <- AntMinuteMeans_ALL_3C[AntMinuteMeans_ALL_3C$time < "2020-06-10 08:00:00",]
AntMinuteMeans_ALL_3D <- AntMinuteMeans_ALL_3D[AntMinuteMeans_ALL_3D$time < "2020-06-10 08:00:00",]
AntMinuteMeans_ALL_4C <- AntMinuteMeans_ALL_4C[AntMinuteMeans_ALL_4C$time < "2020-06-17 08:00:00",]
AntMinuteMeans_ALL_4D <- AntMinuteMeans_ALL_4D[AntMinuteMeans_ALL_4D$time < "2020-06-17 08:00:00",]
AntMinuteMeans_ALL_5C <- AntMinuteMeans_ALL_5C[AntMinuteMeans_ALL_5C$time < "2020-07-01 08:00:00",]
AntMinuteMeans_ALL_5D <- AntMinuteMeans_ALL_5D[AntMinuteMeans_ALL_5D$time < "2020-07-01 08:00:00",]
AntMinuteMeans_ALL_6C <- AntMinuteMeans_ALL_6C[AntMinuteMeans_ALL_6C$time < "2020-07-08 08:00:00",]
AntMinuteMeans_ALL_6D <- AntMinuteMeans_ALL_6D[AntMinuteMeans_ALL_6D$time < "2020-07-08 08:00:00",]
AntMinuteMeans_ALL_7C <- AntMinuteMeans_ALL_7C[AntMinuteMeans_ALL_7C$time < "2020-07-22 08:00:00",]
AntMinuteMeans_ALL_7D <- AntMinuteMeans_ALL_7D[AntMinuteMeans_ALL_7D$time < "2020-07-22 08:00:00",]
AntMinuteMeans_ALL_8C <- AntMinuteMeans_ALL_8C[AntMinuteMeans_ALL_8C$time < "2020-07-29 08:00:00",]
AntMinuteMeans_ALL_8D <- AntMinuteMeans_ALL_8D[AntMinuteMeans_ALL_8D$time < "2020-07-29 08:00:00",]
AntMinuteMeans_ALL_9C <- AntMinuteMeans_ALL_9C[AntMinuteMeans_ALL_9C$time < "2020-09-23 08:00:00",]
AntMinuteMeans_ALL_9D <- AntMinuteMeans_ALL_9D[AntMinuteMeans_ALL_9D$time < "2020-09-23 08:00:00",]
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
```

### In the tracking of a few subcolonies the tracking system crashed and we have no data for some short time periods. We remove the data from those times from the corresponding sub-colony in the replicate

``` r
# True of 3D (4x), 4C (2x), 5D (3x), and 7D (2x)

# 3D gap 1: "2020-06-04 17:05:02" -> "2020-06-04 21:01:20"
# 3D gap 2: "2020-06-06 09:01:04" -> "2020-06-06 11:56:55"
# 3D gap 3: "2020-06-07 07:50:31" -> "2020-06-07 09:29:01"
# 3D gap 4: "2020-06-07 14:55:39" -> "2020-06-07 15:38:22"
AntMinuteMeans_ALL_3C <- AntMinuteMeans_ALL_3C[!("2020-06-04 17:05:02" < AntMinuteMeans_ALL_3C$time & AntMinuteMeans_ALL_3C$time < "2020-06-04 21:01:20"),]
RawContacts_ALL_3C <- RawContacts_ALL_3C[!("2020-06-04 21:01:20" < RawContacts_ALL_3C$start & RawContacts_ALL_3C$start < "2020-06-04 21:01:20"),]
AntMinuteMeans_ALL_3C <- AntMinuteMeans_ALL_3C[!("2020-06-06 09:01:04" < AntMinuteMeans_ALL_3C$time & AntMinuteMeans_ALL_3C$time < "2020-06-06 11:56:55"),]
RawContacts_ALL_3C <- RawContacts_ALL_3C[!("2020-06-06 09:01:04" < RawContacts_ALL_3C$start & RawContacts_ALL_3C$start < "2020-06-06 11:56:55"),]
AntMinuteMeans_ALL_3C <- AntMinuteMeans_ALL_3C[!("2020-06-07 07:50:31" < AntMinuteMeans_ALL_3C$time & AntMinuteMeans_ALL_3C$time < "2020-06-07 09:29:01"),]
RawContacts_ALL_3C <- RawContacts_ALL_3C[!("2020-06-07 07:50:31" < RawContacts_ALL_3C$start & RawContacts_ALL_3C$start < "2020-06-07 09:29:01"),]
AntMinuteMeans_ALL_3C <- AntMinuteMeans_ALL_3C[!("2020-06-07 14:55:39" < AntMinuteMeans_ALL_3C$time & AntMinuteMeans_ALL_3C$time < "2020-06-07 15:38:22"),]
RawContacts_ALL_3C <- RawContacts_ALL_3C[!("2020-06-07 14:55:39" < RawContacts_ALL_3C$start & RawContacts_ALL_3C$start < "2020-06-07 15:38:22"),]
AntMinuteMeans_ALL_3D <- AntMinuteMeans_ALL_3D[!("2020-06-04 17:05:02" < AntMinuteMeans_ALL_3D$time & AntMinuteMeans_ALL_3D$time < "2020-06-04 21:01:20"),]
RawContacts_ALL_3D <- RawContacts_ALL_3D[!("2020-06-04 21:01:20" < RawContacts_ALL_3D$start & RawContacts_ALL_3D$start < "2020-06-04 21:01:20"),]
AntMinuteMeans_ALL_3D <- AntMinuteMeans_ALL_3D[!("2020-06-06 09:01:04" < AntMinuteMeans_ALL_3D$time & AntMinuteMeans_ALL_3D$time < "2020-06-06 11:56:55"),]
RawContacts_ALL_3D <- RawContacts_ALL_3D[!("2020-06-06 09:01:04" < RawContacts_ALL_3D$start & RawContacts_ALL_3D$start < "2020-06-06 11:56:55"),]
AntMinuteMeans_ALL_3D <- AntMinuteMeans_ALL_3D[!("2020-06-07 07:50:31" < AntMinuteMeans_ALL_3D$time & AntMinuteMeans_ALL_3D$time < "2020-06-07 09:29:01"),]
RawContacts_ALL_3D <- RawContacts_ALL_3D[!("2020-06-07 07:50:31" < RawContacts_ALL_3D$start & RawContacts_ALL_3D$start < "2020-06-07 09:29:01"),]
AntMinuteMeans_ALL_3D <- AntMinuteMeans_ALL_3D[!("2020-06-07 14:55:39" < AntMinuteMeans_ALL_3D$time & AntMinuteMeans_ALL_3D$time < "2020-06-07 15:38:22"),]
RawContacts_ALL_3D <- RawContacts_ALL_3D[!("2020-06-07 14:55:39" < RawContacts_ALL_3D$start & RawContacts_ALL_3D$start < "2020-06-07 15:38:22"),]

# 4C gap 1: "2020-06-14 13:46:30" -> "2020-06-14 20:08:30"
# 4C gap 2: "2020-06-15 22:59:18" -> "2020-06-16 06:35:05"
AntMinuteMeans_ALL_4D <- AntMinuteMeans_ALL_4D[!("2020-06-14 13:46:30" < AntMinuteMeans_ALL_4D$time & AntMinuteMeans_ALL_4D$time < "2020-06-14 20:08:30"),]
RawContacts_ALL_4D <- RawContacts_ALL_4D[!("2020-06-14 13:46:30" < RawContacts_ALL_4D$start & RawContacts_ALL_4D$start < "2020-06-14 20:08:30"),]
AntMinuteMeans_ALL_4D <- AntMinuteMeans_ALL_4D[!("2020-06-15 22:59:18" < AntMinuteMeans_ALL_4D$time & AntMinuteMeans_ALL_4D$time < "2020-06-16 06:35:05"),]
RawContacts_ALL_4D <- RawContacts_ALL_4D[!("2020-06-15 22:59:18" < RawContacts_ALL_4D$start & RawContacts_ALL_4D$start < "2020-06-16 06:35:05"),]
AntMinuteMeans_ALL_4C <- AntMinuteMeans_ALL_4C[!("2020-06-14 13:46:30" < AntMinuteMeans_ALL_4C$time & AntMinuteMeans_ALL_4C$time < "2020-06-14 20:08:30"),]
RawContacts_ALL_4C <- RawContacts_ALL_4C[!("2020-06-14 13:46:30" < RawContacts_ALL_4C$start & RawContacts_ALL_4C$start < "2020-06-14 20:08:30"),]
AntMinuteMeans_ALL_4C <- AntMinuteMeans_ALL_4C[!("2020-06-15 22:59:18" < AntMinuteMeans_ALL_4C$time & AntMinuteMeans_ALL_4C$time < "2020-06-16 06:35:05"),]
RawContacts_ALL_4C <- RawContacts_ALL_4C[!("2020-06-15 22:59:18" < RawContacts_ALL_4C$start & RawContacts_ALL_4C$start < "2020-06-16 06:35:05"),]


# 5D gap 1: "2020-06-26 20:38:51" -> "2020-06-27 07:58:18"
# 5D gap 2: "2020-06-27 13:48:20" -> "2020-06-27 16:39:02"
# 5D gap 3: "2020-06-30 28:58:42" -> Not restarted
AntMinuteMeans_ALL_5C <- AntMinuteMeans_ALL_5C[!("2020-06-26 20:38:51" < AntMinuteMeans_ALL_5C$time & AntMinuteMeans_ALL_5C$time < "2020-06-27 07:58:18"),]
RawContacts_ALL_5C <- RawContacts_ALL_5C[!("2020-06-26 20:38:51" < RawContacts_ALL_5C$start & RawContacts_ALL_5C$start < "2020-06-27 07:58:18"),]
AntMinuteMeans_ALL_5C <- AntMinuteMeans_ALL_5C[!("2020-06-27 13:48:20" < AntMinuteMeans_ALL_5C$time & AntMinuteMeans_ALL_5C$time < "2020-06-27 16:39:02"),]
RawContacts_ALL_5C <- RawContacts_ALL_5C[!("2020-06-27 13:48:20" < RawContacts_ALL_5C$start & RawContacts_ALL_5C$start < "2020-06-27 16:39:02"),]
AntMinuteMeans_ALL_5C <- AntMinuteMeans_ALL_5C[!("2020-06-30 28:58:42" < AntMinuteMeans_ALL_5C$time),]
RawContacts_ALL_5C <- RawContacts_ALL_5C[!("2020-06-30 28:58:42" < RawContacts_ALL_5C$start),]
AntMinuteMeans_ALL_5D <- AntMinuteMeans_ALL_5D[!("2020-06-26 20:38:51" < AntMinuteMeans_ALL_5D$time & AntMinuteMeans_ALL_5D$time < "2020-06-27 07:58:18"),]
RawContacts_ALL_5D <- RawContacts_ALL_5D[!("2020-06-26 20:38:51" < RawContacts_ALL_5D$start & RawContacts_ALL_5D$start < "2020-06-27 07:58:18"),]
AntMinuteMeans_ALL_5D <- AntMinuteMeans_ALL_5D[!("2020-06-27 13:48:20" < AntMinuteMeans_ALL_5D$time & AntMinuteMeans_ALL_5D$time < "2020-06-27 16:39:02"),]
RawContacts_ALL_5D <- RawContacts_ALL_5D[!("2020-06-27 13:48:20" < RawContacts_ALL_5D$start & RawContacts_ALL_5D$start < "2020-06-27 16:39:02"),]
AntMinuteMeans_ALL_5D <- AntMinuteMeans_ALL_5D[!("2020-06-30 28:58:42" < AntMinuteMeans_ALL_5D$time),]
RawContacts_ALL_5D <- RawContacts_ALL_5D[!("2020-06-30 28:58:42" < RawContacts_ALL_5D$start),]

# 7D gap 1: "2020-07-20 17:27:42" -> "2020-07-20 18:24:29"
# 7D gap 2: "2020-07-20 23:47:14" -> "2020-07-21 13:04:28"

AntMinuteMeans_ALL_7C <- AntMinuteMeans_ALL_7C[!("2020-07-20 17:27:42" < AntMinuteMeans_ALL_7C$time & AntMinuteMeans_ALL_7C$time < "2020-07-20 18:24:29"),]
RawContacts_ALL_7C <- RawContacts_ALL_7C[!("2020-07-20 17:27:42" < RawContacts_ALL_7C$start & RawContacts_ALL_7C$start < "2020-07-20 18:24:29"),]
AntMinuteMeans_ALL_7C <- AntMinuteMeans_ALL_7C[!("2020-07-20 23:47:14" < AntMinuteMeans_ALL_7C$time & AntMinuteMeans_ALL_7C$time < "2020-07-21 13:04:28"),]
RawContacts_ALL_7C <- RawContacts_ALL_7C[!("2020-07-20 23:47:14" < RawContacts_ALL_7C$start & RawContacts_ALL_7C$start < "2020-07-21 13:04:28"),]
AntMinuteMeans_ALL_7D <- AntMinuteMeans_ALL_7D[!("2020-07-20 17:27:42" < AntMinuteMeans_ALL_7D$time & AntMinuteMeans_ALL_7D$time < "2020-07-20 18:24:29"),]
RawContacts_ALL_7D <- RawContacts_ALL_7D[!("2020-07-20 17:27:42" < RawContacts_ALL_7D$start & RawContacts_ALL_7D$start < "2020-07-20 18:24:29"),]
AntMinuteMeans_ALL_7D <- AntMinuteMeans_ALL_7D[!("2020-07-20 23:47:14" < AntMinuteMeans_ALL_7D$time & AntMinuteMeans_ALL_7D$time < "2020-07-21 13:04:28"),]
RawContacts_ALL_7D <- RawContacts_ALL_7D[!("2020-07-20 23:47:14" < RawContacts_ALL_7D$start & RawContacts_ALL_7D$start < "2020-07-21 13:04:28"),]
```

### Identify and exclude outlier bees from all data sets

``` r
# All replicates began with 100 bees but there were usually ~ 10 that were dead at the end of the replicate
# Others were possibly dying
# We want to keep bees that behaved normally for most of the experiment
# The distribution of total contacts is negatively skewed
# We filter out bees that had fewer interactions than 2 SDs below the mean interactions for their subcolony
# No bees has more interactions than 2 SDs above mean of their subcolony
# This filtration removes 0 to 8 bees per subcolony
SocOut1C <- names(table(c(RawContacts_ALL_1C$ant1, RawContacts_ALL_1C$ant2)))[table(c(RawContacts_ALL_1C$ant1, RawContacts_ALL_1C$ant2)) < mean(table(c(RawContacts_ALL_1C$ant1, RawContacts_ALL_1C$ant2))) - (2*sd(table(c(RawContacts_ALL_1C$ant1, RawContacts_ALL_1C$ant2))))]
SocOut2C <- names(table(c(RawContacts_ALL_2C$ant1, RawContacts_ALL_2C$ant2)))[table(c(RawContacts_ALL_2C$ant1, RawContacts_ALL_2C$ant2)) < mean(table(c(RawContacts_ALL_2C$ant1, RawContacts_ALL_2C$ant2))) - (2*sd(table(c(RawContacts_ALL_2C$ant1, RawContacts_ALL_2C$ant2))))]
SocOut3C <- names(table(c(RawContacts_ALL_3C$ant1, RawContacts_ALL_3C$ant2)))[table(c(RawContacts_ALL_3C$ant1, RawContacts_ALL_3C$ant2)) < mean(table(c(RawContacts_ALL_3C$ant1, RawContacts_ALL_3C$ant2))) - (2*sd(table(c(RawContacts_ALL_3C$ant1, RawContacts_ALL_3C$ant2))))]
SocOut4C <- names(table(c(RawContacts_ALL_4C$ant1, RawContacts_ALL_4C$ant2)))[table(c(RawContacts_ALL_4C$ant1, RawContacts_ALL_4C$ant2)) < mean(table(c(RawContacts_ALL_4C$ant1, RawContacts_ALL_4C$ant2))) - (2*sd(table(c(RawContacts_ALL_4C$ant1, RawContacts_ALL_4C$ant2))))]
SocOut5C <- names(table(c(RawContacts_ALL_5C$ant1, RawContacts_ALL_5C$ant2)))[table(c(RawContacts_ALL_5C$ant1, RawContacts_ALL_5C$ant2)) < mean(table(c(RawContacts_ALL_5C$ant1, RawContacts_ALL_5C$ant2))) - (2*sd(table(c(RawContacts_ALL_5C$ant1, RawContacts_ALL_5C$ant2))))]
SocOut6C <- names(table(c(RawContacts_ALL_6C$ant1, RawContacts_ALL_6C$ant2)))[table(c(RawContacts_ALL_6C$ant1, RawContacts_ALL_6C$ant2)) < mean(table(c(RawContacts_ALL_6C$ant1, RawContacts_ALL_6C$ant2))) - (2*sd(table(c(RawContacts_ALL_6C$ant1, RawContacts_ALL_6C$ant2))))]
SocOut7C <- names(table(c(RawContacts_ALL_7C$ant1, RawContacts_ALL_7C$ant2)))[table(c(RawContacts_ALL_7C$ant1, RawContacts_ALL_7C$ant2)) < mean(table(c(RawContacts_ALL_7C$ant1, RawContacts_ALL_7C$ant2))) - (2*sd(table(c(RawContacts_ALL_7C$ant1, RawContacts_ALL_7C$ant2))))]
SocOut8C <- names(table(c(RawContacts_ALL_8C$ant1, RawContacts_ALL_8C$ant2)))[table(c(RawContacts_ALL_8C$ant1, RawContacts_ALL_8C$ant2)) < mean(table(c(RawContacts_ALL_8C$ant1, RawContacts_ALL_8C$ant2))) - (2*sd(table(c(RawContacts_ALL_8C$ant1, RawContacts_ALL_8C$ant2))))]
SocOut9C <- names(table(c(RawContacts_ALL_9C$ant1, RawContacts_ALL_9C$ant2)))[table(c(RawContacts_ALL_9C$ant1, RawContacts_ALL_9C$ant2)) < mean(table(c(RawContacts_ALL_9C$ant1, RawContacts_ALL_9C$ant2))) - (2*sd(table(c(RawContacts_ALL_9C$ant1, RawContacts_ALL_9C$ant2))))]
SocOut1D <- names(table(c(RawContacts_ALL_1D$ant1, RawContacts_ALL_1D$ant2)))[table(c(RawContacts_ALL_1D$ant1, RawContacts_ALL_1D$ant2)) < mean(table(c(RawContacts_ALL_1D$ant1, RawContacts_ALL_1D$ant2))) - (2*sd(table(c(RawContacts_ALL_1D$ant1, RawContacts_ALL_1D$ant2))))]
SocOut2D <- names(table(c(RawContacts_ALL_2D$ant1, RawContacts_ALL_2D$ant2)))[table(c(RawContacts_ALL_2D$ant1, RawContacts_ALL_2D$ant2)) < mean(table(c(RawContacts_ALL_2D$ant1, RawContacts_ALL_2D$ant2))) - (2*sd(table(c(RawContacts_ALL_2D$ant1, RawContacts_ALL_2D$ant2))))]
SocOut3D <- names(table(c(RawContacts_ALL_3D$ant1, RawContacts_ALL_3D$ant2)))[table(c(RawContacts_ALL_3D$ant1, RawContacts_ALL_3D$ant2)) < mean(table(c(RawContacts_ALL_3D$ant1, RawContacts_ALL_3D$ant2))) - (2*sd(table(c(RawContacts_ALL_3D$ant1, RawContacts_ALL_3D$ant2))))]
SocOut4D <- names(table(c(RawContacts_ALL_4D$ant1, RawContacts_ALL_4D$ant2)))[table(c(RawContacts_ALL_4D$ant1, RawContacts_ALL_4D$ant2)) < mean(table(c(RawContacts_ALL_4D$ant1, RawContacts_ALL_4D$ant2))) - (2*sd(table(c(RawContacts_ALL_4D$ant1, RawContacts_ALL_4D$ant2))))]
SocOut5D <- names(table(c(RawContacts_ALL_5D$ant1, RawContacts_ALL_5D$ant2)))[table(c(RawContacts_ALL_5D$ant1, RawContacts_ALL_5D$ant2)) < mean(table(c(RawContacts_ALL_5D$ant1, RawContacts_ALL_5D$ant2))) - (2*sd(table(c(RawContacts_ALL_5D$ant1, RawContacts_ALL_5D$ant2))))]
SocOut6D <- names(table(c(RawContacts_ALL_6D$ant1, RawContacts_ALL_6D$ant2)))[table(c(RawContacts_ALL_6D$ant1, RawContacts_ALL_6D$ant2)) < mean(table(c(RawContacts_ALL_6D$ant1, RawContacts_ALL_6D$ant2))) - (2*sd(table(c(RawContacts_ALL_6D$ant1, RawContacts_ALL_6D$ant2))))]
SocOut7D <- names(table(c(RawContacts_ALL_7D$ant1, RawContacts_ALL_7D$ant2)))[table(c(RawContacts_ALL_7D$ant1, RawContacts_ALL_7D$ant2)) < mean(table(c(RawContacts_ALL_7D$ant1, RawContacts_ALL_7D$ant2))) - (2*sd(table(c(RawContacts_ALL_7D$ant1, RawContacts_ALL_7D$ant2))))]
SocOut8D <- names(table(c(RawContacts_ALL_8D$ant1, RawContacts_ALL_8D$ant2)))[table(c(RawContacts_ALL_8D$ant1, RawContacts_ALL_8D$ant2)) < mean(table(c(RawContacts_ALL_8D$ant1, RawContacts_ALL_8D$ant2))) - (2*sd(table(c(RawContacts_ALL_8D$ant1, RawContacts_ALL_8D$ant2))))]
SocOut9D <- names(table(c(RawContacts_ALL_9D$ant1, RawContacts_ALL_9D$ant2)))[table(c(RawContacts_ALL_9D$ant1, RawContacts_ALL_9D$ant2)) < mean(table(c(RawContacts_ALL_9D$ant1, RawContacts_ALL_9D$ant2))) - (2*sd(table(c(RawContacts_ALL_9D$ant1, RawContacts_ALL_9D$ant2))))]

# These outliers are saved so the same bees can be filtered out elsewhere
setwd(OUTDIR)
#write.csv(SocOut1C, "SocOut1C.csv", row.names = FALSE)
#write.csv(SocOut2C, "SocOut2C.csv", row.names = FALSE)
#write.csv(SocOut3C, "SocOut3C.csv", row.names = FALSE)
#write.csv(SocOut4C, "SocOut4C.csv", row.names = FALSE)
#write.csv(SocOut5C, "SocOut5C.csv", row.names = FALSE)
#write.csv(SocOut6C, "SocOut6C.csv", row.names = FALSE)
#write.csv(SocOut7C, "SocOut7C.csv", row.names = FALSE)
#write.csv(SocOut8C, "SocOut8C.csv", row.names = FALSE)
#write.csv(SocOut9C, "SocOut9C.csv", row.names = FALSE)
#write.csv(SocOut1D, "SocOut1D.csv", row.names = FALSE)
#write.csv(SocOut2D, "SocOut2D.csv", row.names = FALSE)
#write.csv(SocOut3D, "SocOut3D.csv", row.names = FALSE)
#write.csv(SocOut4D, "SocOut4D.csv", row.names = FALSE)
#write.csv(SocOut5D, "SocOut5D.csv", row.names = FALSE)
#write.csv(SocOut6D, "SocOut6D.csv", row.names = FALSE)
#write.csv(SocOut7D, "SocOut7D.csv", row.names = FALSE)
#write.csv(SocOut8D, "SocOut8D.csv", row.names = FALSE)
#write.csv(SocOut9D, "SocOut9D.csv", row.names = FALSE)

# Remove Outliers
AntMinuteMeans_ALL_1C <- AntMinuteMeans_ALL_1C[!(AntMinuteMeans_ALL_1C$ID %in% SocOut1C),]
AntMinuteMeans_ALL_2C <- AntMinuteMeans_ALL_2C[!(AntMinuteMeans_ALL_2C$ID %in% SocOut2C),]
AntMinuteMeans_ALL_3C <- AntMinuteMeans_ALL_3C[!(AntMinuteMeans_ALL_3C$ID %in% SocOut3C),]
AntMinuteMeans_ALL_4C <- AntMinuteMeans_ALL_4C[!(AntMinuteMeans_ALL_4C$ID %in% SocOut4C),]
AntMinuteMeans_ALL_5C <- AntMinuteMeans_ALL_5C[!(AntMinuteMeans_ALL_5C$ID %in% SocOut5C),]
AntMinuteMeans_ALL_6C <- AntMinuteMeans_ALL_6C[!(AntMinuteMeans_ALL_6C$ID %in% SocOut6C),]
AntMinuteMeans_ALL_7C <- AntMinuteMeans_ALL_7C[!(AntMinuteMeans_ALL_7C$ID %in% SocOut7C),]
AntMinuteMeans_ALL_8C <- AntMinuteMeans_ALL_8C[!(AntMinuteMeans_ALL_8C$ID %in% SocOut8C),]
AntMinuteMeans_ALL_9C <- AntMinuteMeans_ALL_9C[!(AntMinuteMeans_ALL_9C$ID %in% SocOut9C),]
RawContacts_ALL_1C <- RawContacts_ALL_1C[!(RawContacts_ALL_1C$ant1 %in% SocOut1C) & !(RawContacts_ALL_1C$ant2 %in% SocOut1C),]
RawContacts_ALL_2C <- RawContacts_ALL_2C[!(RawContacts_ALL_2C$ant1 %in% SocOut2C) & !(RawContacts_ALL_2C$ant2 %in% SocOut2C),]
RawContacts_ALL_3C <- RawContacts_ALL_3C[!(RawContacts_ALL_3C$ant1 %in% SocOut3C) & !(RawContacts_ALL_3C$ant2 %in% SocOut3C),]
RawContacts_ALL_4C <- RawContacts_ALL_4C[!(RawContacts_ALL_4C$ant1 %in% SocOut4C) & !(RawContacts_ALL_4C$ant2 %in% SocOut4C),]
RawContacts_ALL_5C <- RawContacts_ALL_5C[!(RawContacts_ALL_5C$ant1 %in% SocOut5C) & !(RawContacts_ALL_5C$ant2 %in% SocOut5C),]
RawContacts_ALL_6C <- RawContacts_ALL_6C[!(RawContacts_ALL_6C$ant1 %in% SocOut6C) & !(RawContacts_ALL_6C$ant2 %in% SocOut6C),]
RawContacts_ALL_7C <- RawContacts_ALL_7C[!(RawContacts_ALL_7C$ant1 %in% SocOut7C) & !(RawContacts_ALL_7C$ant2 %in% SocOut7C),]
RawContacts_ALL_8C <- RawContacts_ALL_8C[!(RawContacts_ALL_8C$ant1 %in% SocOut8C) & !(RawContacts_ALL_8C$ant2 %in% SocOut8C),]
RawContacts_ALL_9C <- RawContacts_ALL_9C[!(RawContacts_ALL_9C$ant1 %in% SocOut9C) & !(RawContacts_ALL_9C$ant2 %in% SocOut9C),]
AntMinuteMeans_ALL_1D <- AntMinuteMeans_ALL_1D[!(AntMinuteMeans_ALL_1D$ID %in% SocOut1D),]
AntMinuteMeans_ALL_2D <- AntMinuteMeans_ALL_2D[!(AntMinuteMeans_ALL_2D$ID %in% SocOut2D),]
AntMinuteMeans_ALL_3D <- AntMinuteMeans_ALL_3D[!(AntMinuteMeans_ALL_3D$ID %in% SocOut3D),]
AntMinuteMeans_ALL_4D <- AntMinuteMeans_ALL_4D[!(AntMinuteMeans_ALL_4D$ID %in% SocOut4D),]
AntMinuteMeans_ALL_5D <- AntMinuteMeans_ALL_5D[!(AntMinuteMeans_ALL_5D$ID %in% SocOut5D),]
AntMinuteMeans_ALL_6D <- AntMinuteMeans_ALL_6D[!(AntMinuteMeans_ALL_6D$ID %in% SocOut6D),]
AntMinuteMeans_ALL_7D <- AntMinuteMeans_ALL_7D[!(AntMinuteMeans_ALL_7D$ID %in% SocOut7D),]
AntMinuteMeans_ALL_8D <- AntMinuteMeans_ALL_8D[!(AntMinuteMeans_ALL_8D$ID %in% SocOut8D),]
AntMinuteMeans_ALL_9D <- AntMinuteMeans_ALL_9D[!(AntMinuteMeans_ALL_9D$ID %in% SocOut9D),]
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

### Summarise the speed and interaction data for each subcolony

``` r
# We next construct a single dataframe of summary statistics for each subcolony
# Considering both social interaction and trajectory data

# Number of bees and the total number of contacts for each subcolony
ResDF <- data.frame(condition = rep(c("C","D"), 9),
                    replicate = rep(1:9, each = 2),
                    Bees = c(length(unique(c(AntMinuteMeans_ALL_1C$ID))), length(unique(c(AntMinuteMeans_ALL_1D$ID))),
                             length(unique(c(AntMinuteMeans_ALL_2C$ID))), length(unique(c(AntMinuteMeans_ALL_2D$ID))),
                             length(unique(c(AntMinuteMeans_ALL_3C$ID))), length(unique(c(AntMinuteMeans_ALL_3D$ID))),
                             length(unique(c(AntMinuteMeans_ALL_4C$ID))), length(unique(c(AntMinuteMeans_ALL_4D$ID))),
                             length(unique(c(AntMinuteMeans_ALL_5C$ID))), length(unique(c(AntMinuteMeans_ALL_5D$ID))),
                             length(unique(c(AntMinuteMeans_ALL_6C$ID))), length(unique(c(AntMinuteMeans_ALL_6D$ID))),
                             length(unique(c(AntMinuteMeans_ALL_7C$ID))), length(unique(c(AntMinuteMeans_ALL_7D$ID))),
                             length(unique(c(AntMinuteMeans_ALL_8C$ID))), length(unique(c(AntMinuteMeans_ALL_8D$ID))),
                             length(unique(c(AntMinuteMeans_ALL_9C$ID))), length(unique(c(AntMinuteMeans_ALL_9D$ID)))),
                    TotContacts = c(nrow(RawContacts_ALL_1C), nrow(RawContacts_ALL_1D), nrow(RawContacts_ALL_2C), nrow(RawContacts_ALL_2D),
                                    nrow(RawContacts_ALL_3C), nrow(RawContacts_ALL_3D), nrow(RawContacts_ALL_4C), nrow(RawContacts_ALL_4D),
                                    nrow(RawContacts_ALL_5C), nrow(RawContacts_ALL_5D), nrow(RawContacts_ALL_6C), nrow(RawContacts_ALL_6D),
                                    nrow(RawContacts_ALL_7C), nrow(RawContacts_ALL_7D), nrow(RawContacts_ALL_8C), nrow(RawContacts_ALL_8D),
                                    nrow(RawContacts_ALL_9C), nrow(RawContacts_ALL_9D)))
# The average number of contacts per bee in the subcolony

ResDF$ContPerBee <- ResDF$TotContacts/ResDF$Bees

# The average duration of contact
ResDF$MedianDur <- c(as.numeric(median(RawContacts_ALL_1C$duration)), as.numeric(median(RawContacts_ALL_1D$duration)),
                     as.numeric(median(RawContacts_ALL_2C$duration)), as.numeric(median(RawContacts_ALL_2D$duration)),
                     as.numeric(median(RawContacts_ALL_3C$duration)), as.numeric(median(RawContacts_ALL_3D$duration)),
                     as.numeric(median(RawContacts_ALL_4C$duration)), as.numeric(median(RawContacts_ALL_4D$duration)),
                     as.numeric(median(RawContacts_ALL_5C$duration)), as.numeric(median(RawContacts_ALL_5D$duration)),
                     as.numeric(median(RawContacts_ALL_6C$duration)), as.numeric(median(RawContacts_ALL_6D$duration)),
                     as.numeric(median(RawContacts_ALL_7C$duration)), as.numeric(median(RawContacts_ALL_7D$duration)),
                     as.numeric(median(RawContacts_ALL_8C$duration)), as.numeric(median(RawContacts_ALL_8D$duration)),
                     as.numeric(median(RawContacts_ALL_9C$duration)), as.numeric(median(RawContacts_ALL_9D$duration)))

# In fort-studio two body parts were annotated onto each bee: head and body
# We define social interacts as head-head contact
# And use the comparison of head-head contact and body-body contact to look whether certain pairs interact deliberately
# Or are just bumping into whichever bee is nearby

RawContacts_11_1C <- RawContacts_ALL_1C[RawContacts_ALL_1C$types == "1-1",]
RawContacts_11_2C <- RawContacts_ALL_2C[RawContacts_ALL_2C$types == "1-1",]
RawContacts_11_3C <- RawContacts_ALL_3C[RawContacts_ALL_3C$types == "1-1",]
RawContacts_11_4C <- RawContacts_ALL_4C[RawContacts_ALL_4C$types == "1-1",]
RawContacts_11_5C <- RawContacts_ALL_5C[RawContacts_ALL_5C$types == "1-1",]
RawContacts_11_6C <- RawContacts_ALL_6C[RawContacts_ALL_6C$types == "1-1",]
RawContacts_11_7C <- RawContacts_ALL_7C[RawContacts_ALL_7C$types == "1-1",]
RawContacts_11_8C <- RawContacts_ALL_8C[RawContacts_ALL_8C$types == "1-1",]
RawContacts_11_9C <- RawContacts_ALL_9C[RawContacts_ALL_9C$types == "1-1",]
RawContacts_11_1D <- RawContacts_ALL_1D[RawContacts_ALL_1D$types == "1-1",]
RawContacts_11_2D <- RawContacts_ALL_2D[RawContacts_ALL_2D$types == "1-1",]
RawContacts_11_3D <- RawContacts_ALL_3D[RawContacts_ALL_3D$types == "1-1",]
RawContacts_11_4D <- RawContacts_ALL_4D[RawContacts_ALL_4D$types == "1-1",]
RawContacts_11_5D <- RawContacts_ALL_5D[RawContacts_ALL_5D$types == "1-1",]
RawContacts_11_6D <- RawContacts_ALL_6D[RawContacts_ALL_6D$types == "1-1",]
RawContacts_11_7D <- RawContacts_ALL_7D[RawContacts_ALL_7D$types == "1-1",]
RawContacts_11_8D <- RawContacts_ALL_8D[RawContacts_ALL_8D$types == "1-1",]
RawContacts_11_9D <- RawContacts_ALL_9D[RawContacts_ALL_9D$types == "1-1",]

ResDF$Conts11 <- c(nrow(RawContacts_11_1C), nrow(RawContacts_11_1D), nrow(RawContacts_11_2C), nrow(RawContacts_11_2D),
                   nrow(RawContacts_11_3C), nrow(RawContacts_11_3D), nrow(RawContacts_11_4C), nrow(RawContacts_11_4D),
                   nrow(RawContacts_11_5C), nrow(RawContacts_11_5D), nrow(RawContacts_11_6C), nrow(RawContacts_11_6D),
                   nrow(RawContacts_11_7C), nrow(RawContacts_11_7D), nrow(RawContacts_11_8C), nrow(RawContacts_11_8D),
                   nrow(RawContacts_11_9C), nrow(RawContacts_11_9D))
ResDF$Prop11 <- ResDF$Conts11/ResDF$TotContacts

RawContacts_22_1C <- RawContacts_ALL_1C[RawContacts_ALL_1C$types == "2-2",]
RawContacts_22_2C <- RawContacts_ALL_2C[RawContacts_ALL_2C$types == "2-2",]
RawContacts_22_3C <- RawContacts_ALL_3C[RawContacts_ALL_3C$types == "2-2",]
RawContacts_22_4C <- RawContacts_ALL_4C[RawContacts_ALL_4C$types == "2-2",]
RawContacts_22_5C <- RawContacts_ALL_5C[RawContacts_ALL_5C$types == "2-2",]
RawContacts_22_6C <- RawContacts_ALL_6C[RawContacts_ALL_6C$types == "2-2",]
RawContacts_22_7C <- RawContacts_ALL_7C[RawContacts_ALL_7C$types == "2-2",]
RawContacts_22_8C <- RawContacts_ALL_8C[RawContacts_ALL_8C$types == "2-2",]
RawContacts_22_9C <- RawContacts_ALL_9C[RawContacts_ALL_9C$types == "2-2",]
RawContacts_22_1D <- RawContacts_ALL_1D[RawContacts_ALL_1D$types == "2-2",]
RawContacts_22_2D <- RawContacts_ALL_2D[RawContacts_ALL_2D$types == "2-2",]
RawContacts_22_3D <- RawContacts_ALL_3D[RawContacts_ALL_3D$types == "2-2",]
RawContacts_22_4D <- RawContacts_ALL_4D[RawContacts_ALL_4D$types == "2-2",]
RawContacts_22_5D <- RawContacts_ALL_5D[RawContacts_ALL_5D$types == "2-2",]
RawContacts_22_6D <- RawContacts_ALL_6D[RawContacts_ALL_6D$types == "2-2",]
RawContacts_22_7D <- RawContacts_ALL_7D[RawContacts_ALL_7D$types == "2-2",]
RawContacts_22_8D <- RawContacts_ALL_8D[RawContacts_ALL_8D$types == "2-2",]
RawContacts_22_9D <- RawContacts_ALL_9D[RawContacts_ALL_9D$types == "2-2",]

ResDF$Conts22 <- c(nrow(RawContacts_22_1C), nrow(RawContacts_22_1D), nrow(RawContacts_22_2C), nrow(RawContacts_22_2D),
                   nrow(RawContacts_22_3C), nrow(RawContacts_22_3D), nrow(RawContacts_22_4C), nrow(RawContacts_22_4D),
                   nrow(RawContacts_22_5C), nrow(RawContacts_22_5D), nrow(RawContacts_22_6C), nrow(RawContacts_22_6D),
                   nrow(RawContacts_22_7C), nrow(RawContacts_22_7D), nrow(RawContacts_22_8C), nrow(RawContacts_22_8D),
                   nrow(RawContacts_22_9C), nrow(RawContacts_22_9D))

ResDF$Prop22 <- ResDF$Conts22/ResDF$TotContacts

# In fort-studio in some replicates the head was labelled first (and encoded as 1) and the body second (2)
# This was inverted for other replicates
# Since there were always more body-body contacts than head-head contacts, we relabel these as below

ResDF$HH <- pmin(ResDF$Conts11, ResDF$Conts22)
ResDF$BB <- pmax(ResDF$Conts11, ResDF$Conts22)

# Calculate the number of HH and BB contacts per bee 
ResDF$HHPerBee <- ResDF$HH/ResDF$Bees
ResDF$BBPerBee <- ResDF$BB/ResDF$Bees

# NB: The number of contacts per bee increases with the number of bees
par(mfrow=c(1,1))
plot(ResDF$Bees, ResDF$HHPerBee, col = as.numeric(ResDF$condition), pch = 16, cex = 2,
     xlab = "Number of bees", ylab = "Contacts per bee")
abline(lm(ResDF$HHPerBee ~ ResDF$Bees), col = "dark gray", lwd = 3)
```

![](01-TrackingOverview_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# Many of the trajectories are very short; it does not make so much sense to use these to calculate movement dynamics
# The trajectories are filtered to keep only those with more than 100 observations

AntMinuteMeans_abbrv_1C <- AntMinuteMeans_ALL_1C[AntMinuteMeans_ALL_1C$N_observations > 100,]
AntMinuteMeans_abbrv_2C <- AntMinuteMeans_ALL_2C[AntMinuteMeans_ALL_2C$N_observations > 100,]
AntMinuteMeans_abbrv_3C <- AntMinuteMeans_ALL_3C[AntMinuteMeans_ALL_3C$N_observations > 100,]
AntMinuteMeans_abbrv_4C <- AntMinuteMeans_ALL_4C[AntMinuteMeans_ALL_4C$N_observations > 100,]
AntMinuteMeans_abbrv_5C <- AntMinuteMeans_ALL_5C[AntMinuteMeans_ALL_5C$N_observations > 100,]
AntMinuteMeans_abbrv_6C <- AntMinuteMeans_ALL_6C[AntMinuteMeans_ALL_6C$N_observations > 100,]
AntMinuteMeans_abbrv_7C <- AntMinuteMeans_ALL_7C[AntMinuteMeans_ALL_7C$N_observations > 100,]
AntMinuteMeans_abbrv_8C <- AntMinuteMeans_ALL_8C[AntMinuteMeans_ALL_8C$N_observations > 100,]
AntMinuteMeans_abbrv_9C <- AntMinuteMeans_ALL_9C[AntMinuteMeans_ALL_9C$N_observations > 100,]
AntMinuteMeans_abbrv_1D <- AntMinuteMeans_ALL_1D[AntMinuteMeans_ALL_1D$N_observations > 100,]
AntMinuteMeans_abbrv_2D <- AntMinuteMeans_ALL_2D[AntMinuteMeans_ALL_2D$N_observations > 100,]
AntMinuteMeans_abbrv_3D <- AntMinuteMeans_ALL_3D[AntMinuteMeans_ALL_3D$N_observations > 100,]
AntMinuteMeans_abbrv_4D <- AntMinuteMeans_ALL_4D[AntMinuteMeans_ALL_4D$N_observations > 100,]
AntMinuteMeans_abbrv_5D <- AntMinuteMeans_ALL_5D[AntMinuteMeans_ALL_5D$N_observations > 100,]
AntMinuteMeans_abbrv_6D <- AntMinuteMeans_ALL_6D[AntMinuteMeans_ALL_6D$N_observations > 100,]
AntMinuteMeans_abbrv_7D <- AntMinuteMeans_ALL_7D[AntMinuteMeans_ALL_7D$N_observations > 100,]
AntMinuteMeans_abbrv_8D <- AntMinuteMeans_ALL_8D[AntMinuteMeans_ALL_8D$N_observations > 100,]
AntMinuteMeans_abbrv_9D <- AntMinuteMeans_ALL_9D[AntMinuteMeans_ALL_9D$N_observations > 100,]

# A few trajectories give infinite speed values - these we remove
AntMinuteMeans_abbrv_7C <- AntMinuteMeans_abbrv_7C[-(which(is.infinite(AntMinuteMeans_abbrv_7C$speed))),]
AntMinuteMeans_abbrv_8C <- AntMinuteMeans_abbrv_8C[-(which(is.infinite(AntMinuteMeans_abbrv_8C$speed))),]
AntMinuteMeans_abbrv_9C <- AntMinuteMeans_abbrv_9C[-(which(is.infinite(AntMinuteMeans_abbrv_9C$speed))),]
AntMinuteMeans_abbrv_7D <- AntMinuteMeans_abbrv_7D[-(which(is.infinite(AntMinuteMeans_abbrv_7D$speed))),]
AntMinuteMeans_abbrv_8D <- AntMinuteMeans_abbrv_8D[-(which(is.infinite(AntMinuteMeans_abbrv_8D$speed))),]
AntMinuteMeans_abbrv_9D <- AntMinuteMeans_abbrv_9D[-(which(is.infinite(AntMinuteMeans_abbrv_9D$speed))),]

# Add the average speed and SD of speed to the summary dataframe
ResDF$speed <- c(mean(AntMinuteMeans_abbrv_1C$speed), mean(AntMinuteMeans_abbrv_1D$speed),
                 mean(AntMinuteMeans_abbrv_2C$speed), mean(AntMinuteMeans_abbrv_2D$speed),
                 mean(AntMinuteMeans_abbrv_3C$speed), mean(AntMinuteMeans_abbrv_3D$speed),
                 mean(AntMinuteMeans_abbrv_4C$speed), mean(AntMinuteMeans_abbrv_4D$speed),
                 mean(AntMinuteMeans_abbrv_5C$speed), mean(AntMinuteMeans_abbrv_5D$speed),
                 mean(AntMinuteMeans_abbrv_6C$speed), mean(AntMinuteMeans_abbrv_6D$speed),
                 mean(AntMinuteMeans_abbrv_7C$speed), mean(AntMinuteMeans_abbrv_7D$speed),
                 mean(AntMinuteMeans_abbrv_8C$speed), mean(AntMinuteMeans_abbrv_8D$speed),
                 mean(AntMinuteMeans_abbrv_9C$speed), mean(AntMinuteMeans_abbrv_9D$speed))

ResDF$speed_sd <- c(mean(AntMinuteMeans_abbrv_1C$speed_sd), mean(AntMinuteMeans_abbrv_1D$speed_sd),
                    mean(AntMinuteMeans_abbrv_2C$speed_sd), mean(AntMinuteMeans_abbrv_2D$speed_sd),
                    mean(AntMinuteMeans_abbrv_3C$speed_sd), mean(AntMinuteMeans_abbrv_3D$speed_sd),
                    mean(AntMinuteMeans_abbrv_4C$speed_sd), mean(AntMinuteMeans_abbrv_4D$speed_sd),
                    mean(AntMinuteMeans_abbrv_5C$speed_sd), mean(AntMinuteMeans_abbrv_5D$speed_sd),
                    mean(AntMinuteMeans_abbrv_6C$speed_sd), mean(AntMinuteMeans_abbrv_6D$speed_sd),
                    mean(AntMinuteMeans_abbrv_7C$speed_sd), mean(AntMinuteMeans_abbrv_7D$speed_sd),
                    mean(AntMinuteMeans_abbrv_8C$speed_sd), mean(AntMinuteMeans_abbrv_8D$speed_sd),
                    mean(AntMinuteMeans_abbrv_9C$speed_sd), mean(AntMinuteMeans_abbrv_9D$speed_sd))


# Construct a social network for each subcolony
AggContacts1C       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_11_1C)
AggContacts2C       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_11_2C)
AggContacts3C       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_11_3C)
AggContacts4C       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_11_4C)
AggContacts5C       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_11_5C)
AggContacts6C       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_22_6C)
AggContacts7C       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_22_7C)
AggContacts8C       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_22_8C)
AggContacts9C       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_22_9C)
AggContacts1D       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_11_1D)
AggContacts2D       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_11_2D)
AggContacts3D       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_11_3D)
AggContacts4D       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_11_4D)
AggContacts5D       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_11_5D)
AggContacts6D       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_11_6D)
AggContacts7D       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_22_7D)
AggContacts8D       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_22_8D)
AggContacts9D       <- aggregate(end ~ ant1 + ant2, FUN=length,RawContacts_22_9D)

colnames(AggContacts1C)[match("end",colnames(AggContacts1C))] <- "weight"
colnames(AggContacts2C)[match("end",colnames(AggContacts2C))] <- "weight"
colnames(AggContacts3C)[match("end",colnames(AggContacts3C))] <- "weight"
colnames(AggContacts4C)[match("end",colnames(AggContacts4C))] <- "weight"
colnames(AggContacts5C)[match("end",colnames(AggContacts5C))] <- "weight"
colnames(AggContacts6C)[match("end",colnames(AggContacts6C))] <- "weight"
colnames(AggContacts7C)[match("end",colnames(AggContacts7C))] <- "weight"
colnames(AggContacts8C)[match("end",colnames(AggContacts8C))] <- "weight"
colnames(AggContacts9C)[match("end",colnames(AggContacts9C))] <- "weight"
colnames(AggContacts1D)[match("end",colnames(AggContacts1D))] <- "weight"
colnames(AggContacts2D)[match("end",colnames(AggContacts2D))] <- "weight"
colnames(AggContacts3D)[match("end",colnames(AggContacts3D))] <- "weight"
colnames(AggContacts4D)[match("end",colnames(AggContacts4D))] <- "weight"
colnames(AggContacts5D)[match("end",colnames(AggContacts5D))] <- "weight"
colnames(AggContacts6D)[match("end",colnames(AggContacts6D))] <- "weight"
colnames(AggContacts7D)[match("end",colnames(AggContacts7D))] <- "weight"
colnames(AggContacts8D)[match("end",colnames(AggContacts8D))] <- "weight"
colnames(AggContacts9D)[match("end",colnames(AggContacts9D))] <- "weight"

# Save edgelists for other analyses
#write.csv(AggContacts1C, "Network_1C.csv", row.names = FALSE)
#write.csv(AggContacts1D, "Network_1D.csv", row.names = FALSE)
#write.csv(AggContacts2C, "Network_2C.csv", row.names = FALSE)
#write.csv(AggContacts2D, "Network_2D.csv", row.names = FALSE)
#write.csv(AggContacts3C, "Network_3C.csv", row.names = FALSE)
#write.csv(AggContacts3D, "Network_3D.csv", row.names = FALSE)
#write.csv(AggContacts4C, "Network_4C.csv", row.names = FALSE)
#write.csv(AggContacts4D, "Network_4D.csv", row.names = FALSE)
#write.csv(AggContacts5C, "Network_5C.csv", row.names = FALSE)
#write.csv(AggContacts5D, "Network_5D.csv", row.names = FALSE)
#write.csv(AggContacts6C, "Network_6C.csv", row.names = FALSE)
#write.csv(AggContacts6D, "Network_6D.csv", row.names = FALSE)
#write.csv(AggContacts7C, "Network_7C.csv", row.names = FALSE)
#write.csv(AggContacts7D, "Network_7D.csv", row.names = FALSE)
#write.csv(AggContacts8C, "Network_8C.csv", row.names = FALSE)
#write.csv(AggContacts8D, "Network_8D.csv", row.names = FALSE)
#write.csv(AggContacts9C, "Network_9C.csv", row.names = FALSE)
#write.csv(AggContacts9D, "Network_9D.csv", row.names = FALSE)

# Calculate a specialization index based on the individual-level variation in edge-weight
# For each node We take a vector of its edge-wights (including 0 values for nodes it is not connected to)
# We then take the SD of this vector, and average within sub-colonies

AntIDs1C <- unique(c(AggContacts1C$ant1, AggContacts1C$ant2))
SpecInd1C <- c()
for (ID in AntIDs1C){
  weights <- AggContacts1C[AggContacts1C$ant1 == ID | AggContacts1C$ant2 == ID,]$weight
  missing <- length(AntIDs1C) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd1C <- c(SpecInd1C, var(weights))
}
AntIDs2C <- unique(c(AggContacts2C$ant1, AggContacts2C$ant2))
SpecInd2C <- c()
for (ID in AntIDs2C){
  weights <- AggContacts2C[AggContacts2C$ant1 == ID | AggContacts2C$ant2 == ID,]$weight
  missing <- length(AntIDs2C) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd2C <- c(SpecInd2C, var(weights))
}
AntIDs3C <- unique(c(AggContacts3C$ant1, AggContacts3C$ant2))
SpecInd3C <- c()
for (ID in AntIDs3C){
  weights <- AggContacts3C[AggContacts3C$ant1 == ID | AggContacts3C$ant2 == ID,]$weight
  missing <- length(AntIDs3C) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd3C <- c(SpecInd3C, var(weights))
}
AntIDs4C <- unique(c(AggContacts4C$ant1, AggContacts4C$ant2))
SpecInd4C <- c()
for (ID in AntIDs4C){
  weights <- AggContacts4C[AggContacts4C$ant1 == ID | AggContacts4C$ant2 == ID,]$weight
  missing <- length(AntIDs4C) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd4C <- c(SpecInd4C, var(weights))
}
AntIDs5C <- unique(c(AggContacts5C$ant1, AggContacts5C$ant2))
SpecInd5C <- c()
for (ID in AntIDs5C){
  weights <- AggContacts5C[AggContacts5C$ant1 == ID | AggContacts5C$ant2 == ID,]$weight
  missing <- length(AntIDs5C) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd5C <- c(SpecInd5C, var(weights))
}
AntIDs6C <- unique(c(AggContacts6C$ant1, AggContacts6C$ant2))
SpecInd6C <- c()
for (ID in AntIDs6C){
  weights <- AggContacts6C[AggContacts6C$ant1 == ID | AggContacts6C$ant2 == ID,]$weight
  missing <- length(AntIDs6C) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd6C <- c(SpecInd6C, var(weights))
}
AntIDs7C <- unique(c(AggContacts7C$ant1, AggContacts7C$ant2))
SpecInd7C <- c()
for (ID in AntIDs7C){
  weights <- AggContacts7C[AggContacts7C$ant1 == ID | AggContacts7C$ant2 == ID,]$weight
  missing <- length(AntIDs7C) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd7C <- c(SpecInd7C, var(weights))
}
AntIDs8C <- unique(c(AggContacts8C$ant1, AggContacts8C$ant2))
SpecInd8C <- c()
for (ID in AntIDs8C){
  weights <- AggContacts8C[AggContacts8C$ant1 == ID | AggContacts8C$ant2 == ID,]$weight
  missing <- length(AntIDs8C) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd8C <- c(SpecInd8C, var(weights))
}
AntIDs9C <- unique(c(AggContacts9C$ant1, AggContacts9C$ant2))
SpecInd9C <- c()
for (ID in AntIDs9C){
  weights <- AggContacts9C[AggContacts9C$ant1 == ID | AggContacts9C$ant2 == ID,]$weight
  missing <- length(AntIDs9C) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd9C <- c(SpecInd9C, var(weights))
}
AntIDs1D <- unique(c(AggContacts1D$ant1, AggContacts1D$ant2))
SpecInd1D <- c()
for (ID in AntIDs1D){
  weights <- AggContacts1D[AggContacts1D$ant1 == ID | AggContacts1D$ant2 == ID,]$weight
  missing <- length(AntIDs1D) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd1D <- c(SpecInd1D, var(weights))
}
AntIDs2D <- unique(c(AggContacts2D$ant1, AggContacts2D$ant2))
SpecInd2D <- c()
for (ID in AntIDs2D){
  weights <- AggContacts2D[AggContacts2D$ant1 == ID | AggContacts2D$ant2 == ID,]$weight
  missing <- length(AntIDs2D) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd2D <- c(SpecInd2D, var(weights))
}
AntIDs3D <- unique(c(AggContacts3D$ant1, AggContacts3D$ant2))
SpecInd3D <- c()
for (ID in AntIDs3D){
  weights <- AggContacts3D[AggContacts3D$ant1 == ID | AggContacts3D$ant2 == ID,]$weight
  missing <- length(AntIDs3D) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd3D <- c(SpecInd3D, var(weights))
}
AntIDs4D <- unique(c(AggContacts4D$ant1, AggContacts4D$ant2))
SpecInd4D <- c()
for (ID in AntIDs4D){
  weights <- AggContacts4D[AggContacts4D$ant1 == ID | AggContacts4D$ant2 == ID,]$weight
  missing <- length(AntIDs4D) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd4D <- c(SpecInd4D, var(weights))
}
AntIDs5D <- unique(c(AggContacts5D$ant1, AggContacts5D$ant2))
SpecInd5D <- c()
for (ID in AntIDs5D){
  weights <- AggContacts5D[AggContacts5D$ant1 == ID | AggContacts5D$ant2 == ID,]$weight
  missing <- length(AntIDs5D) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd5D <- c(SpecInd5D, var(weights))
}
AntIDs6D <- unique(c(AggContacts6D$ant1, AggContacts6D$ant2))
SpecInd6D <- c()
for (ID in AntIDs6D){
  weights <- AggContacts6D[AggContacts6D$ant1 == ID | AggContacts6D$ant2 == ID,]$weight
  missing <- length(AntIDs6D) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd6D <- c(SpecInd6D, var(weights))
}
AntIDs7D <- unique(c(AggContacts7D$ant1, AggContacts7D$ant2))
SpecInd7D <- c()
for (ID in AntIDs7D){
  weights <- AggContacts7D[AggContacts7D$ant1 == ID | AggContacts7D$ant2 == ID,]$weight
  missing <- length(AntIDs7D) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd7D <- c(SpecInd7D, var(weights))
}
AntIDs8D <- unique(c(AggContacts8D$ant1, AggContacts8D$ant2))
SpecInd8D <- c()
for (ID in AntIDs8D){
  weights <- AggContacts8D[AggContacts8D$ant1 == ID | AggContacts8D$ant2 == ID,]$weight
  missing <- length(AntIDs8D) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd8D <- c(SpecInd8D, var(weights))
}
AntIDs9D <- unique(c(AggContacts9D$ant1, AggContacts9D$ant2))
SpecInd9D <- c()
for (ID in AntIDs9D){
  weights <- AggContacts9D[AggContacts9D$ant1 == ID | AggContacts9D$ant2 == ID,]$weight
  missing <- length(AntIDs9D) - length(weights)
  weights <- c(weights, rep(0, missing))
  SpecInd9D <- c(SpecInd9D, var(weights))
}

# Add specialization indexes to the summary dataframe
ResDF$Spec_Ind <- c(mean(SpecInd1C), mean(SpecInd1D),mean(SpecInd2C), mean(SpecInd2D),
                    mean(SpecInd3C), mean(SpecInd3D),mean(SpecInd4C), mean(SpecInd4D),
                    mean(SpecInd5C), mean(SpecInd5D),mean(SpecInd6C), mean(SpecInd6D),
                    mean(SpecInd7C), mean(SpecInd6D),mean(SpecInd8C), mean(SpecInd8D),
                    mean(SpecInd9C), mean(SpecInd9D))

# Calculate the proportion of time that each bee spent in each arena
Space1C <- c()
for (ID in AntIDs1C){
  Space1C <- c(Space1C, mean(RawContacts_11_1C[RawContacts_11_1C$ant1 == ID | RawContacts_11_1C$ant2 == ID,]$space))
}
Space2C <- c()
for (ID in AntIDs2C){
  Space2C <- c(Space2C, mean(RawContacts_11_2C[RawContacts_11_2C$ant1 == ID | RawContacts_11_2C$ant2 == ID,]$space))
}
Space3C <- c()
for (ID in AntIDs3C){
  Space3C <- c(Space3C, mean(RawContacts_11_3C[RawContacts_11_3C$ant1 == ID | RawContacts_11_3C$ant2 == ID,]$space))
}
Space4C <- c()
for (ID in AntIDs4C){
  Space4C <- c(Space4C, mean(RawContacts_11_4C[RawContacts_11_4C$ant1 == ID | RawContacts_11_4C$ant2 == ID,]$space))
}
Space5C <- c()
for (ID in AntIDs5C){
  Space5C <- c(Space5C, mean(RawContacts_11_5C[RawContacts_11_5C$ant1 == ID | RawContacts_11_5C$ant2 == ID,]$space))
}
Space6C <- c()
for (ID in AntIDs6C){
  Space6C <- c(Space6C, mean(RawContacts_11_6C[RawContacts_11_6C$ant1 == ID | RawContacts_11_6C$ant2 == ID,]$space))
}
Space7C <- c()
for (ID in AntIDs7C){
  Space7C <- c(Space7C, mean(RawContacts_11_7C[RawContacts_11_7C$ant1 == ID | RawContacts_11_7C$ant2 == ID,]$space))
}
Space8C <- c()
for (ID in AntIDs8C){
  Space8C <- c(Space8C, mean(RawContacts_11_8C[RawContacts_11_8C$ant1 == ID | RawContacts_11_8C$ant2 == ID,]$space))
}
Space9C <- c()
for (ID in AntIDs9C){
  Space9C <- c(Space9C, mean(RawContacts_11_9C[RawContacts_11_9C$ant1 == ID | RawContacts_11_9C$ant2 == ID,]$space))
}
Space1D <- c()
for (ID in AntIDs1D){
  Space1D <- c(Space1D, mean(RawContacts_11_1D[RawContacts_11_1D$ant1 == ID | RawContacts_11_1D$ant2 == ID,]$space))
}
Space2D <- c()
for (ID in AntIDs2D){
  Space2D <- c(Space2D, mean(RawContacts_11_2D[RawContacts_11_2D$ant1 == ID | RawContacts_11_2D$ant2 == ID,]$space))
}
Space3D <- c()
for (ID in AntIDs3D){
  Space3D <- c(Space3D, mean(RawContacts_11_3D[RawContacts_11_3D$ant1 == ID | RawContacts_11_3D$ant2 == ID,]$space))
}
Space4D <- c()
for (ID in AntIDs4D){
  Space4D <- c(Space4D, mean(RawContacts_11_4D[RawContacts_11_4D$ant1 == ID | RawContacts_11_4D$ant2 == ID,]$space))
}
Space5D <- c()
for (ID in AntIDs5D){
  Space5D <- c(Space5D, mean(RawContacts_11_5D[RawContacts_11_5D$ant1 == ID | RawContacts_11_5D$ant2 == ID,]$space))
}
Space6D <- c()
for (ID in AntIDs6D){
  Space6D <- c(Space6D, mean(RawContacts_11_6D[RawContacts_11_6D$ant1 == ID | RawContacts_11_6D$ant2 == ID,]$space))
}
Space7D <- c()
for (ID in AntIDs7D){
  Space7D <- c(Space7D, mean(RawContacts_11_7D[RawContacts_11_7D$ant1 == ID | RawContacts_11_7D$ant2 == ID,]$space))
}
Space8D <- c()
for (ID in AntIDs8D){
  Space8D <- c(Space8D, mean(RawContacts_11_8D[RawContacts_11_8D$ant1 == ID | RawContacts_11_8D$ant2 == ID,]$space))
}
Space9D <- c()
for (ID in AntIDs9D){
  Space9D <- c(Space9D, mean(RawContacts_11_9D[RawContacts_11_9D$ant1 == ID | RawContacts_11_9D$ant2 == ID,]$space))
}

# Add space-use date to summary dataframe
ResDF$Space <- c(mean(Space1C), mean(Space1D),mean(Space2C), mean(Space2D),
                 mean(Space3C), mean(Space3D),mean(Space4C), mean(Space4D),
                 mean(Space5C), mean(Space5D),mean(Space6C), mean(Space6D),
                 mean(Space7C), mean(Space6D),mean(Space8C), mean(Space8D),
                 mean(Space9C), mean(Space9D))

# Correct for the one replicate where the foraging box and nest box numerical identities were inverted
ResDF$Space[1] <- 3 - ResDF$Space[1]

# Save data
#write.csv(ResDF, "ResDF.csv")

# Save speed and time data for timeseries analsis
AntMinuteMeans_abbrv_1C$time <- round(as.numeric(AntMinuteMeans_abbrv_1C$time - AntMinuteMeans_abbrv_1C$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_2C$time <- round(as.numeric(AntMinuteMeans_abbrv_2C$time - AntMinuteMeans_abbrv_2C$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_3C$time <- round(as.numeric(AntMinuteMeans_abbrv_3C$time - AntMinuteMeans_abbrv_3C$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_4C$time <- round(as.numeric(AntMinuteMeans_abbrv_4C$time - AntMinuteMeans_abbrv_4C$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_5C$time <- round(as.numeric(AntMinuteMeans_abbrv_5C$time - AntMinuteMeans_abbrv_5C$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_6C$time <- round(as.numeric(AntMinuteMeans_abbrv_6C$time - AntMinuteMeans_abbrv_6C$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_7C$time <- round(as.numeric(AntMinuteMeans_abbrv_7C$time - AntMinuteMeans_abbrv_7C$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_8C$time <- round(as.numeric(AntMinuteMeans_abbrv_8C$time - AntMinuteMeans_abbrv_8C$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_9C$time <- round(as.numeric(AntMinuteMeans_abbrv_9C$time - AntMinuteMeans_abbrv_9C$time[1])/(60*60), 0)



AntMinuteMeans_abbrv_C <- rbind(AntMinuteMeans_abbrv_1C,
                                AntMinuteMeans_abbrv_2C,
                                AntMinuteMeans_abbrv_3C,
                                AntMinuteMeans_abbrv_4C,
                                AntMinuteMeans_abbrv_5C,
                                AntMinuteMeans_abbrv_6C,
                                AntMinuteMeans_abbrv_7C,
                                AntMinuteMeans_abbrv_8C,
                                AntMinuteMeans_abbrv_9C)

AntMinuteMeans_abbrv_C$time <- as.character(AntMinuteMeans_abbrv_C$time)

AntMinuteMeans_abbrv_C$time <- as.factor(AntMinuteMeans_abbrv_C$time)
Agg_C <- aggregate(.~time, data=AntMinuteMeans_abbrv_C, mean)
Agg_C$time <- as.numeric(as.character(Agg_C$time))
Agg_C <- Agg_C[order(Agg_C$time),]

AntMinuteMeans_abbrv_1D$time <- round(as.numeric(AntMinuteMeans_abbrv_1D$time - AntMinuteMeans_abbrv_1D$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_2D$time <- round(as.numeric(AntMinuteMeans_abbrv_2D$time - AntMinuteMeans_abbrv_2D$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_3D$time <- round(as.numeric(AntMinuteMeans_abbrv_3D$time - AntMinuteMeans_abbrv_3D$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_4D$time <- round(as.numeric(AntMinuteMeans_abbrv_4D$time - AntMinuteMeans_abbrv_4D$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_5D$time <- round(as.numeric(AntMinuteMeans_abbrv_5D$time - AntMinuteMeans_abbrv_5D$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_6D$time <- round(as.numeric(AntMinuteMeans_abbrv_6D$time - AntMinuteMeans_abbrv_6D$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_7D$time <- round(as.numeric(AntMinuteMeans_abbrv_7D$time - AntMinuteMeans_abbrv_7D$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_8D$time <- round(as.numeric(AntMinuteMeans_abbrv_8D$time - AntMinuteMeans_abbrv_8D$time[1])/(60*60), 0)
AntMinuteMeans_abbrv_9D$time <- round(as.numeric(AntMinuteMeans_abbrv_9D$time - AntMinuteMeans_abbrv_9D$time[1])/(60*60), 0)

AntMinuteMeans_abbrv_D <- rbind(AntMinuteMeans_abbrv_1D,
                                AntMinuteMeans_abbrv_2D,
                                AntMinuteMeans_abbrv_3D,
                                AntMinuteMeans_abbrv_4D,
                                AntMinuteMeans_abbrv_5D,
                                AntMinuteMeans_abbrv_6D,
                                AntMinuteMeans_abbrv_7D,
                                AntMinuteMeans_abbrv_8D,
                                AntMinuteMeans_abbrv_9D)

AntMinuteMeans_abbrv_D$time <- as.character(AntMinuteMeans_abbrv_D$time)

AntMinuteMeans_abbrv_D$time <- as.factor(AntMinuteMeans_abbrv_D$time)
Agg_D <- aggregate(.~time, data=AntMinuteMeans_abbrv_D, mean)
Agg_D$time <- as.numeric(as.character(Agg_D$time))
Agg_D <- Agg_D[order(Agg_D$time),]

setwd(OUTDIR)
#write.csv(Agg_C, "HourlySpeed_C.csv", row.names = FALSE)
#write.csv(Agg_D, "HourlySpeed_D.csv", row.names = FALSE)
```

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
    ## [1] lubridate_1.8.0  OneR_2.2         pals_1.7         igraph_1.2.9    
    ## [5] stringr_1.4.0    ggfortify_0.4.13 ggplot2_3.3.5    gapminder_0.3.0 
    ## [9] jpeg_0.1-9      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] highr_0.9        pillar_1.6.4     compiler_4.1.2   tools_4.1.2     
    ##  [5] digest_0.6.29    evaluate_0.14    lifecycle_1.0.1  tibble_3.1.6    
    ##  [9] gtable_0.3.0     pkgconfig_2.0.3  rlang_0.4.12     DBI_1.1.1       
    ## [13] mapproj_1.2.7    yaml_2.2.1       xfun_0.28        fastmap_1.1.0   
    ## [17] gridExtra_2.3    withr_2.4.3      dplyr_1.0.7      knitr_1.36      
    ## [21] maps_3.4.0       generics_0.1.1   vctrs_0.3.8      grid_4.1.2      
    ## [25] tidyselect_1.1.1 glue_1.5.1       R6_2.5.1         fansi_0.5.0     
    ## [29] rmarkdown_2.11   purrr_0.3.4      tidyr_1.1.4      magrittr_2.0.1  
    ## [33] scales_1.1.1     ellipsis_0.3.2   htmltools_0.5.2  dichromat_2.0-0 
    ## [37] colorspace_2.0-2 utf8_1.2.2       stringi_1.7.6    munsell_0.5.0   
    ## [41] crayon_1.4.2
