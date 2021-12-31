02-IndividualLevelAnalysis
================
Tomas Kay
12/21/2021

## This script calculates individual-level social interaction data

## For comparison with molecular data on individual brains and guts

### Set directories and load packages

``` r
MAINDIR <- "/Volumes/Lacie/BM"
INDIR   <- paste(MAINDIR, "RawData", sep = "/")
OUTDIR  <- paste(MAINDIR, "ProcessedData", sep = "/")
FIGSDIR <- paste(MAINDIR, "Figures", sep = "/")
CODEDIR <- paste(MAINDIR, "Code", sep = "/")
library('igraph')
```

### Establishing correspondence between the unique IDs written on the tubes for the molecular lab, and the ‘TagID’ and ’AntID’s in FortStudio

### Read in the IDs of the bees processed in the molecular lab

``` r
setwd(INDIR)
guts <- read.csv("ExtractedGuts.csv")
source(paste(CODEDIR, "to_FORT_tag36ARTag_ID.r", sep = "/"))
```

### Convert Tube IDs to Fort TagIDs

``` r
correctedID <- c()
for(ID in guts$ID){
  correctedID <- c(correctedID, to_FORT_tag36ARTag_ID(ID))
}
guts$corrected <- correctedID
```

# Read in tag-bee conversion files

``` r
setwd(INDIR)
conversion1C <- read.csv("Conversion_1C.csv")
conversion1D <- read.csv("Conversion_1D.csv")
conversion2C <- read.csv("Conversion_2C.csv")
conversion2D <- read.csv("Conversion_2D.csv")
conversion3C <- read.csv("Conversion_3C.csv")
conversion3D <- read.csv("Conversion_3D.csv")
conversion4C <- read.csv("Conversion_4C.csv")
conversion4D <- read.csv("Conversion_4D.csv")
conversion5C <- read.csv("Conversion_5C.csv")
conversion5D <- read.csv("Conversion_5D.csv")
conversion6C <- read.csv("Conversion_6C.csv")
conversion6D <- read.csv("Conversion_6D.csv")
conversion7C <- read.csv("Conversion_7C.csv")
conversion7D <- read.csv("Conversion_7D.csv")
conversion8C <- read.csv("Conversion_8C.csv")
conversion8D <- read.csv("Conversion_8D.csv")
conversion9C <- read.csv("Conversion_9C.csv")
conversion9D <- read.csv("Conversion_9D.csv")
```

# Separate by tube IDs by sub-colony and theck tube IDs are present in tracking data

``` r
extracted1C <- guts[guts$Rep == 1,]
extracted1D <- guts[guts$Rep == 2,]
extracted2C <- guts[guts$Rep == 3,]
extracted2D <- guts[guts$Rep == 4,]
extracted3C <- guts[guts$Rep == 6,]
extracted3D <- guts[guts$Rep == 5,]
extracted4C <- guts[guts$Rep == 7,]
extracted4D <- guts[guts$Rep == 8,]
extracted5C <- guts[guts$Rep == 10,]
extracted5D <- guts[guts$Rep == 9,]
extracted6C <- guts[guts$Rep == 12,]
extracted6D <- guts[guts$Rep == 11,]
extracted7C <- guts[guts$Rep == 14,]
extracted7D <- guts[guts$Rep == 13,]
extracted8C <- guts[guts$Rep == 16,]
extracted8D <- guts[guts$Rep == 15,]
extracted9C <- guts[guts$Rep == 18,]
extracted9D <- guts[guts$Rep == 17,]
extracted1C$corrected %in% conversion1C$tagDecimalValue
extracted1D$corrected %in% conversion1D$tagDecimalValue
extracted2C$corrected %in% conversion2C$tagDecimalValue
extracted2D$corrected %in% conversion2D$tagDecimalValue
extracted3C$corrected %in% conversion3C$tagDecimalValue
extracted3D$corrected %in% conversion3D$tagDecimalValue
extracted4C$corrected %in% conversion4C$tagDecimalValue # missing 1 in rep 4C
extracted4C <- extracted4C[extracted4C$corrected %in% conversion4C$tagDecimalValue,]
extracted4D$corrected %in% conversion4D$tagDecimalValue
extracted5C$corrected %in% conversion5C$tagDecimalValue
extracted5D$corrected %in% conversion5D$tagDecimalValue
extracted6C$corrected %in% conversion6C$tagDecimalValue
extracted6D$corrected %in% conversion6D$tagDecimalValue
extracted7C$corrected %in% conversion7C$tagDecimalValue
extracted7D$corrected %in% conversion7D$tagDecimalValue
extracted8C$corrected %in% conversion8C$tagDecimalValue
extracted8D$corrected %in% conversion8D$tagDecimalValue
extracted9C$corrected %in% conversion9C$tagDecimalValue
extracted9D$corrected %in% conversion9D$tagDecimalValue
```

# Create a data frame for each subcolony of tube ID -> bee ID

``` r
conversion1C_trimmed <- conversion1C[conversion1C$tagDecimalValue %in% extracted1C$corrected,]
extracted1C          <- extracted1C[order(extracted1C$corrected),]
DF_1C <- data.frame(Tagscan = extracted1C$ID, tag = extracted1C$corrected, bee = conversion1C_trimmed$antID)
conversion1D_trimmed <- conversion1D[conversion1D$tagDecimalValue %in% extracted1D$corrected,]
extracted1D          <- extracted1D[order(extracted1D$corrected),]
DF_1D <- data.frame(Tagscan = extracted1D$ID, tag = extracted1D$corrected, bee = conversion1D_trimmed$antID)
conversion2C_trimmed <- conversion2C[conversion2C$tagDecimalValue %in% extracted2C$corrected,]
extracted2C          <- extracted2C[order(extracted2C$corrected),]
DF_2C <- data.frame(Tagscan = extracted2C$ID, tag = extracted2C$corrected, bee = conversion2C_trimmed$antID)
conversion2D_trimmed <- conversion2D[conversion2D$tagDecimalValue %in% extracted2D$corrected,]
extracted2D          <- extracted2D[order(extracted2D$corrected),]
DF_2D <- data.frame(Tagscan = extracted2D$ID, tag = extracted2D$corrected, bee = conversion2D_trimmed$antID)
conversion3C_trimmed <- conversion3C[conversion3C$tagDecimalValue %in% extracted3C$corrected,]
extracted3C          <- extracted3C[order(extracted3C$corrected),]
DF_3C <- data.frame(Tagscan = extracted3C$ID, tag = extracted3C$corrected, bee = conversion3C_trimmed$antID)
conversion3D_trimmed <- conversion3D[conversion3D$tagDecimalValue %in% extracted3D$corrected,]
extracted3D          <- extracted3D[order(extracted3D$corrected),]
DF_3D <- data.frame(Tagscan = extracted3D$ID, tag = extracted3D$corrected, bee = conversion3D_trimmed$antID)
conversion4C_trimmed <- conversion4C[conversion4C$tagDecimalValue %in% extracted4C$corrected,]
extracted4C          <- extracted4C[order(extracted4C$corrected),]
extracted4C         <- extracted4C[extracted4C$corrected  %in% conversion4C$tagDecimalValue,]
DF_4C <- data.frame(Tagscan = extracted4C$ID, tag = extracted4C$corrected, bee = conversion4C_trimmed$antID)
conversion4D_trimmed <- conversion4D[conversion4D$tagDecimalValue %in% extracted4D$corrected,]
extracted4D          <- extracted4D[order(extracted4D$corrected),]
DF_4D <- data.frame(Tagscan = extracted4D$ID, tag = extracted4D$corrected, bee = conversion4D_trimmed$antID)
conversion5C_trimmed <- conversion5C[conversion5C$tagDecimalValue %in% extracted5C$corrected,]
extracted5C          <- extracted5C[order(extracted5C$corrected),]
DF_5C <- data.frame(Tagscan = extracted5C$ID, tag = extracted5C$corrected, bee = conversion5C_trimmed$antID)
conversion5D_trimmed <- conversion5D[conversion5D$tagDecimalValue %in% extracted5D$corrected,]
extracted5D          <- extracted5D[order(extracted5D$corrected),]
DF_5D <- data.frame(Tagscan = extracted5D$ID, tag = extracted5D$corrected, bee = conversion5D_trimmed$antID)
conversion6C_trimmed <- conversion6C[conversion6C$tagDecimalValue %in% extracted6C$corrected,]
extracted6C          <- extracted6C[order(extracted6C$corrected),]
DF_6C <- data.frame(Tagscan = extracted6C$ID, tag = extracted6C$corrected, bee = conversion6C_trimmed$antID)
conversion6D_trimmed <- conversion6D[conversion6D$tagDecimalValue %in% extracted6D$corrected,]
extracted6D          <- extracted6D[order(extracted6D$corrected),]
DF_6D <- data.frame(Tagscan = extracted6D$ID, tag = extracted6D$corrected, bee = conversion6D_trimmed$antID)
conversion7C_trimmed <- conversion7C[conversion7C$tagDecimalValue %in% extracted7C$corrected,]
extracted7C          <- extracted7C[order(extracted7C$corrected),]
DF_7C <- data.frame(Tagscan = extracted7C$ID, tag = extracted7C$corrected, bee = conversion7C_trimmed$antID)
conversion7D_trimmed <- conversion7D[conversion7D$tagDecimalValue %in% extracted7D$corrected,]
extracted7D          <- extracted7D[order(extracted7D$corrected),]
DF_7D <- data.frame(Tagscan = extracted7D$ID, tag = extracted7D$corrected, bee = conversion7D_trimmed$antID)
conversion8C_trimmed <- conversion8C[conversion8C$tagDecimalValue %in% extracted8C$corrected,]
extracted8C          <- extracted8C[order(extracted8C$corrected),]
DF_8C <- data.frame(Tagscan = extracted8C$ID, tag = extracted8C$corrected, bee = conversion8C_trimmed$antID)
conversion8D_trimmed <- conversion8D[conversion8D$tagDecimalValue %in% extracted8D$corrected,]
extracted8D          <- extracted8D[order(extracted8D$corrected),]
DF_8D <- data.frame(Tagscan = extracted8D$ID, tag = extracted8D$corrected, bee = conversion8D_trimmed$antID)
conversion9C_trimmed <- conversion9C[conversion9C$tagDecimalValue %in% extracted9C$corrected,]
extracted9C          <- extracted9C[order(extracted9C$corrected),]
DF_9C <- data.frame(Tagscan = extracted9C$ID, tag = extracted9C$corrected, bee = conversion9C_trimmed$antID)
conversion9D_trimmed <- conversion9D[conversion9D$tagDecimalValue %in% extracted9D$corrected,]
extracted9D          <- extracted9D[order(extracted9D$corrected),]
DF_9D <- data.frame(Tagscan = extracted9D$ID, tag = extracted9D$corrected, bee = conversion9D_trimmed$antID)
```

# Read in social interatction data

``` r
setwd(OUTDIR)
net_1C <- read.csv("Network_1C.csv")
net_1D <- read.csv("Network_1D.csv")
net_2C <- read.csv("Network_2C.csv")
net_2D <- read.csv("Network_2D.csv")
net_3C <- read.csv("Network_3C.csv")
net_3D <- read.csv("Network_3D.csv")
net_4C <- read.csv("Network_4C.csv")
net_4D <- read.csv("Network_4D.csv")
net_5C <- read.csv("Network_5C.csv")
net_5D <- read.csv("Network_5D.csv")
net_6C <- read.csv("Network_6C.csv")
net_6D <- read.csv("Network_6D.csv")
net_7C <- read.csv("Network_7C.csv")
net_7D <- read.csv("Network_7D.csv")
net_8C <- read.csv("Network_8C.csv")
net_8D <- read.csv("Network_8D.csv")
net_9C <- read.csv("Network_9C.csv")
net_9D <- read.csv("Network_9D.csv")
```

# Calculate strength for each bee

``` r
strength_1C <- colSums(as.matrix(get.adjacency(graph.data.frame(net_1C, directed = FALSE), attr='weight', type = 'both')))
strength_1D <- colSums(as.matrix(get.adjacency(graph.data.frame(net_1D, directed = FALSE), attr='weight', type = 'both')))
strength_1D <- strength_1D[order(as.numeric(names(strength_1D)))]
strength_2C <- colSums(as.matrix(get.adjacency(graph.data.frame(net_2C, directed = FALSE), attr='weight', type = 'both')))
strength_2D <- colSums(as.matrix(get.adjacency(graph.data.frame(net_2D, directed = FALSE), attr='weight', type = 'both')))
strength_3C <- colSums(as.matrix(get.adjacency(graph.data.frame(net_3C, directed = FALSE), attr='weight', type = 'both')))
strength_3D <- colSums(as.matrix(get.adjacency(graph.data.frame(net_3D, directed = FALSE), attr='weight', type = 'both')))
strength_4C <- colSums(as.matrix(get.adjacency(graph.data.frame(net_4C, directed = FALSE), attr='weight', type = 'both')))
strength_4D <- colSums(as.matrix(get.adjacency(graph.data.frame(net_4D, directed = FALSE), attr='weight', type = 'both')))
strength_5C <- colSums(as.matrix(get.adjacency(graph.data.frame(net_5C, directed = FALSE), attr='weight', type = 'both')))
strength_5D <- colSums(as.matrix(get.adjacency(graph.data.frame(net_5D, directed = FALSE), attr='weight', type = 'both')))
strength_6C <- colSums(as.matrix(get.adjacency(graph.data.frame(net_6C, directed = FALSE), attr='weight', type = 'both')))
strength_6D <- colSums(as.matrix(get.adjacency(graph.data.frame(net_6D, directed = FALSE), attr='weight', type = 'both')))
strength_7C <- colSums(as.matrix(get.adjacency(graph.data.frame(net_7C, directed = FALSE), attr='weight', type = 'both')))
strength_7D <- colSums(as.matrix(get.adjacency(graph.data.frame(net_7D, directed = FALSE), attr='weight', type = 'both')))
strength_7D <- strength_7D[order(as.numeric(names(strength_7D)))]
strength_8C <- colSums(as.matrix(get.adjacency(graph.data.frame(net_8C, directed = FALSE), attr='weight', type = 'both')))
strength_8D <- colSums(as.matrix(get.adjacency(graph.data.frame(net_8D, directed = FALSE), attr='weight', type = 'both')))
strength_9C <- colSums(as.matrix(get.adjacency(graph.data.frame(net_9C, directed = FALSE), attr='weight', type = 'both')))
strength_9D <- colSums(as.matrix(get.adjacency(graph.data.frame(net_9D, directed = FALSE), attr='weight', type = 'both')))
```

# Create a dataframe of all bees and trim to keep only bees that have gut data

# Five bees processed in the lab are not in the social network

``` r
strengthDF <- data.frame(Rep = c(rep("1C", length(strength_1C)), rep("2C", length(strength_2C)), rep("3C", length(strength_3C)),
                                 rep("4C", length(strength_4C)), rep("5C", length(strength_5C)), rep("6C", length(strength_6C)),
                                 rep("7C", length(strength_7C)), rep("8C", length(strength_8C)), rep("9C", length(strength_9C)),
                                 rep("1D", length(strength_1D)), rep("2D", length(strength_2D)), rep("3D", length(strength_3D)),
                                 rep("4D", length(strength_4D)), rep("5D", length(strength_5D)), rep("6D", length(strength_6D)),
                                 rep("7D", length(strength_7D)), rep("8D", length(strength_8D)), rep("9D", length(strength_9D))),
                         ID = names(c(strength_1C, strength_2C, strength_3C,
                                      strength_4C, strength_5C, strength_6C,
                                      strength_7C, strength_8C, strength_9C,
                                      strength_1D, strength_2D, strength_3D,
                                      strength_4D, strength_5D, strength_6D,
                                      strength_7D, strength_8D, strength_9D)),
                         strength = c(strength_1C, strength_2C, strength_3C,
                                      strength_4C, strength_5C, strength_6C,
                                      strength_7C, strength_8C, strength_9C,
                                      strength_1D, strength_2D, strength_3D,
                                      strength_4D, strength_5D, strength_6D,
                                      strength_7D, strength_8D, strength_9D))
strength_1C_trimmed <- strength_1C[names(strength_1C)%in% DF_1C$bee]
DF_1C$strength <- strength_1C_trimmed
DF_1C$rep <- 1
DF_1C$Rep <- "1C"
strength_1D_trimmed <- strength_1D[names(strength_1D)%in% DF_1D$bee]
DF_1D$strength <- strength_1D_trimmed
DF_1D$rep <- 2
DF_1D$Rep <- "1D"
strength_2C_trimmed <- strength_2C[names(strength_2C)%in% DF_2C$bee]
DF_2C <- DF_2C[DF_2C$bee %in% names(strength_2C_trimmed),] # 1 missing
DF_2C$strength <- strength_2C_trimmed
DF_2C$rep <- 3
DF_2C$Rep <- "2C"
strength_2D_trimmed <- strength_2D[names(strength_2D)%in% DF_2D$bee]
DF_2D <- DF_2D[DF_2D$bee %in% names(strength_2D_trimmed),] # 1 missing
DF_2D$strength <- strength_2D_trimmed
DF_2D$rep <- 4
DF_2D$Rep <- "2D"
strength_3C_trimmed <- strength_3C[names(strength_3C)%in% DF_3C$bee]
DF_3C$strength <- strength_3C_trimmed
DF_3C$rep <- 6
DF_3C$Rep <- "3C"
strength_3D_trimmed <- strength_3D[names(strength_3D)%in% DF_3D$bee]
DF_3D <- DF_3D[DF_3D$bee %in% names(strength_3D_trimmed),] # 1 missing
DF_3D$strength <- strength_3D_trimmed
DF_3D$rep <- 5
DF_3D$Rep <- "3D"
strength_4C_trimmed <- strength_4C[names(strength_4C)%in% DF_4C$bee]
DF_4C$strength <- strength_4C_trimmed
DF_4C$rep <- 7
DF_4C$Rep <- "4C"
strength_4D_trimmed <- strength_4D[names(strength_4D)%in% DF_4D$bee]
DF_4D <- DF_4D[DF_4D$bee %in% names(strength_4D_trimmed),] # 1 missing
DF_4D$strength <- strength_4D_trimmed
DF_4D$rep <- 8
DF_4D$Rep <- "4D"
strength_5C_trimmed <- strength_5C[names(strength_5C)%in% DF_5C$bee]
DF_5C$strength <- strength_5C_trimmed
DF_5C$rep <- 10
DF_5C$Rep <- "5C"
strength_5D_trimmed <- strength_5D[names(strength_5D)%in% DF_5D$bee]
DF_5D$strength <- strength_5D_trimmed
DF_5D$rep <- 9
DF_5D$Rep <- "5D"
strength_6C_trimmed <- strength_6C[names(strength_6C)%in% DF_6C$bee]
DF_6C$strength <- strength_6C_trimmed
DF_6C$rep <- 12
DF_6C$Rep <- "6C"
strength_6D_trimmed <- strength_6D[names(strength_6D)%in% DF_6D$bee]
DF_6D <- DF_6D[DF_6D$bee %in% names(strength_6D_trimmed),] # 1 missing
DF_6D$strength <- strength_6D_trimmed
DF_6D$rep <- 11
DF_6D$Rep <- "6D"
strength_7C_trimmed <- strength_7C[names(strength_7C)%in% DF_7C$bee]
DF_7C$strength <- strength_7C_trimmed
DF_7C$rep <- 14
DF_7C$Rep <- "7C"
strength_7D_trimmed <- strength_7D[names(strength_7D)%in% DF_7D$bee]
DF_7D$strength <- strength_7D_trimmed
DF_7D$rep <- 13
DF_7D$Rep <- "7D"
strength_8C_trimmed <- strength_8C[names(strength_8C)%in% DF_8C$bee]
DF_8C$strength <- strength_8C_trimmed
DF_8C$rep <- 16
DF_8C$Rep <- "8C"
strength_8D_trimmed <- strength_8D[names(strength_8D)%in% DF_8D$bee]
DF_8D$strength <- strength_8D_trimmed
DF_8D$rep <- 15
DF_8D$Rep <- "8D"
strength_9C_trimmed <- strength_9C[names(strength_9C)%in% DF_9C$bee]
DF_9C$strength <- strength_9C_trimmed
DF_9C$rep <- 18
DF_9C$Rep <- "9C"
strength_9D_trimmed <- strength_9D[names(strength_9D)%in% DF_9D$bee]
DF_9D$strength <- strength_9D_trimmed
DF_9D$rep <- 17
DF_9D$Rep <- "9D"
Strength_Guts <- rbind(DF_1C, DF_1D,
                       DF_2C, DF_2D,
                       DF_3C, DF_3D,
                       DF_4C, DF_4D,
                       DF_5C, DF_5D,
                       DF_6C, DF_6D,
                       DF_7C, DF_7D,
                       DF_8C, DF_8D,
                       DF_9C, DF_9D)
```

# Add the total number of bees in the replicate to the dataframe

``` r
Strength_Guts$Tot <- c(rep(length(unique(c(net_1C$ant1, net_1C$ant2))), nrow(DF_1C)),
                       rep(length(unique(c(net_1D$ant1, net_1D$ant2))), nrow(DF_1D)),
                       rep(length(unique(c(net_2C$ant1, net_2C$ant2))), nrow(DF_2C)),
                       rep(length(unique(c(net_2D$ant1, net_2D$ant2))), nrow(DF_2D)),
                       rep(length(unique(c(net_3C$ant1, net_3C$ant2))), nrow(DF_3C)),
                       rep(length(unique(c(net_3D$ant1, net_3D$ant2))), nrow(DF_3D)),
                       rep(length(unique(c(net_4C$ant1, net_4C$ant2))), nrow(DF_4C)),
                       rep(length(unique(c(net_4D$ant1, net_4D$ant2))), nrow(DF_4D)),
                       rep(length(unique(c(net_5C$ant1, net_5C$ant2))), nrow(DF_5C)),
                       rep(length(unique(c(net_5D$ant1, net_5D$ant2))), nrow(DF_5D)),
                       rep(length(unique(c(net_6C$ant1, net_6C$ant2))), nrow(DF_6C)),
                       rep(length(unique(c(net_6D$ant1, net_6D$ant2))), nrow(DF_6D)),
                       rep(length(unique(c(net_7C$ant1, net_7C$ant2))), nrow(DF_7C)),
                       rep(length(unique(c(net_7D$ant1, net_7D$ant2))), nrow(DF_7D)),
                       rep(length(unique(c(net_8C$ant1, net_8C$ant2))), nrow(DF_8C)),
                       rep(length(unique(c(net_8D$ant1, net_8D$ant2))), nrow(DF_8D)),
                       rep(length(unique(c(net_9C$ant1, net_9C$ant2))), nrow(DF_9C)),
                       rep(length(unique(c(net_9D$ant1, net_9D$ant2))), nrow(DF_9D)))
colnames(Strength_Guts)[3] <- "ID"
MetaData <- merge(strengthDF, Strength_Guts, by = c("Rep", "ID", "strength"), all = TRUE)
MetaData <- MetaData[,c(1,2,4,7,3)]
```

### Add total number of bees per subcolony

### And weight strength by number of bees

``` r
Totals <- c()
for(row in 1:nrow(MetaData)){
  rep <- MetaData[row,1]
  Totals <- c(Totals, as.numeric(table(MetaData$Rep)[names(table(MetaData$Rep)) == rep]))
}
MetaData$Tot <- Totals
MetaData$StrengthPerBee <- MetaData$strength / (MetaData$Tot - 1)
```

# calculate residuals and save dataframe

``` r
StrengthResidual       <- c()
StrengthPerBeeResidual <- c()
for(row in 1:nrow(MetaData)){
  rep      <- MetaData[row,1]
  repDATA  <- MetaData[MetaData$Rep == rep,]
  ave.strength       <- mean(repDATA$strength)
  ave.StrengthPerBee <- mean(repDATA$StrengthPerBee)
  StrengthResidual       <- c(StrengthResidual, (MetaData[row,5]- ave.strength))
  StrengthPerBeeResidual <- c(StrengthPerBeeResidual, (MetaData[row,6] - ave.StrengthPerBee))
}
MetaData$strengthResidual <- StrengthResidual 
MetaData$StrengthPerBeeResidual <- StrengthPerBeeResidual
#write.csv(MetaData, "BehaviouralMetadata.csv", row.names = FALSE)
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
    ## [1] igraph_1.2.9
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] lattice_0.20-45 digest_0.6.29   grid_4.1.2      magrittr_2.0.1 
    ##  [5] evaluate_0.14   rlang_0.4.12    stringi_1.7.6   Matrix_1.3-4   
    ##  [9] rmarkdown_2.11  tools_4.1.2     stringr_1.4.0   xfun_0.28      
    ## [13] yaml_2.2.1      fastmap_1.1.0   compiler_4.1.2  pkgconfig_2.0.3
    ## [17] htmltools_0.5.2 knitr_1.36
