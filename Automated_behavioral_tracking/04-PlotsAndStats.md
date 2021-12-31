Plots and statistical analyses of the automated behavioral tracking
experiment
================
Joanito Liberti

## Load data

``` r
setwd("/Volumes/gr_Engel/jliberti/BeeTracking/BM/Figures/")
dd <- read.table("/Volumes/gr_Engel/jliberti/BeeTracking/BM/ProcessedData/ResDF.csv", fill = TRUE ,sep="," , header=T)
dd$condition <- as.character(dd$condition)
dd$condition[dd$condition=="D"] <- "MD"
dd$condition[dd$condition=="C"] <- "CL"
dd$condition <- as.factor(dd$condition)
dd$replicate<-factor(dd$replicate)
```

## Plots main figure

``` r
library(ggplot2)
pl1 <- ggplot(dd, aes(x=relevel(condition, "MD"), y=HHPerBee, color=replicate, group=replicate)) +
  theme_light() +
  xlab(element_blank()) + # Remove the x axis title
  ylab("Head to head interactions per bee") +
  geom_line() +
  geom_boxplot(aes(x=relevel(condition, "MD"), group=condition), outlier.shape=NA, width=0.2, alpha = 0)

pl2 <- ggplot(dd, aes(x=relevel(condition, "MD"), y=Spec_Ind, color=replicate, group=replicate)) +
  theme_light() +
  xlab(element_blank()) + # Remove the x axis title
  ylab("Interaction bias") +
  geom_boxplot(aes(x=relevel(condition, "MD"), group=condition), outlier.shape=NA, width=0.2, alpha = 0) +
  geom_line() 

pl3 <- ggplot(dd, aes(x=relevel(condition, "MD"), y=BBPerBee, color=replicate, group=replicate)) +
  theme_light() +
  xlab(element_blank()) + # Remove the x axis title
  ylab("Body to body interactions per bee") +
  geom_boxplot(aes(x=relevel(condition, "MD"), group=condition), outlier.shape=NA, width=0.2, alpha = 0) +
  geom_line() 
  
pl4 <- ggplot(dd, aes(x=relevel(condition, "MD"), y=speed_sd, color=replicate, group=replicate)) +
  theme_light() +
  xlab(element_blank()) + # Remove the x axis title
  ylab("Average speed") +
  geom_boxplot(aes(x=relevel(condition, "MD"), group=condition), outlier.shape=NA, width=0.2, alpha = 0) +
  geom_line() 

library(ggpubr)
ggarrange(pl1,pl2,pl3,pl4, nrow=1, widths=c(1,1,1,1), common.legend = T)
```

![](AutomatedTracking_PlotsAndStats_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggsave(height=2,width=5,dpi=300, filename="MainTrackingFigures.pdf", useDingbats=FALSE)
```

## Statistical analyses

``` r
##############################################################################
## Paired t test for num of head to head interactions per bee
library(reshape2)
data_wide <- dcast(dd, replicate ~ condition, value.var="HHPerBee")

# Shapiro-Wilk normality test for the difference between paired values
d <- with(data_wide, 
          CL - MD)
shapiro.test(d)# W = 0.92415, p-value = 0.4277
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  d
    ## W = 0.92415, p-value = 0.4277

``` r
# Compute t-test
res <- t.test(data_wide$CL, data_wide$MD, paired = TRUE, alternative =  "two.sided")
res # t = 2.8235, df = 8, p-value = 0.02237
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  data_wide$CL and data_wide$MD
    ## t = 2.8235, df = 8, p-value = 0.02237
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##   94.14657 933.14395
    ## sample estimates:
    ## mean of the differences 
    ##                513.6453

``` r
##############################################################################
## Paired t test for interaction bias
data_wide <- dcast(dd, replicate ~ condition, value.var="Spec_Ind")

# Shapiro-Wilk normality test for the difference between paired values
d <- with(data_wide, 
          CL - MD)
shapiro.test(d)# W = 0.9142, p-value = 0.3464
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  d
    ## W = 0.9142, p-value = 0.3464

``` r
# Compute t-test
res <- t.test(data_wide$CL, data_wide$MD, paired = TRUE, alternative =  "two.sided")
res # t = 2.928, df = 8, p-value = 0.01906
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  data_wide$CL and data_wide$MD
    ## t = 2.928, df = 8, p-value = 0.01906
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##   124.2145 1045.3159
    ## sample estimates:
    ## mean of the differences 
    ##                584.7652

``` r
##############################################################################
## Paired t test for number of body to body interactions per bee
data_wide <- dcast(dd, replicate ~ condition, value.var="BBPerBee")

# Shapiro-Wilk normality test for the difference between paired values
d <- with(data_wide, 
          CL - MD)
shapiro.test(d)# W = 0.93682, p-value = 0.5488
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  d
    ## W = 0.93682, p-value = 0.5488

``` r
# Compute t-test
res <- t.test(data_wide$CL, data_wide$MD, paired = TRUE, alternative =  "two.sided")
res # t = 0.31677, df = 8, p-value = 0.7595
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  data_wide$CL and data_wide$MD
    ## t = 0.31677, df = 8, p-value = 0.7595
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -1269.376  1673.659
    ## sample estimates:
    ## mean of the differences 
    ##                202.1419

``` r
##############################################################################
## Paired t test for average speed
data_wide <- dcast(dd, replicate ~ condition, value.var="speed")

# Shapiro-Wilk normality test for the difference between paired values
d <- with(data_wide, 
          CL - MD)
shapiro.test(d)# W = 0.93655, p-value = 0.5461
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  d
    ## W = 0.93655, p-value = 0.5461

``` r
# Compute t-test
res <- t.test(data_wide$CL, data_wide$MD, paired = TRUE, alternative =  "two.sided")
res # t = 0.35305, df = 8, p-value = 0.7332
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  data_wide$CL and data_wide$MD
    ## t = 0.35305, df = 8, p-value = 0.7332
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -15.02610  20.45891
    ## sample estimates:
    ## mean of the differences 
    ##                2.716402

## Supplementary figures

``` r
S1 <- ggplot(dd, aes(x=relevel(condition, "MD"), y=speed_sd, color=replicate, group=replicate)) +
  theme_light() +
  xlab(element_blank()) + # Remove the x axis title
  ylab("Speed SD") +
  geom_line() +
geom_boxplot(aes(x=relevel(condition, "MD"), group=condition), outlier.shape=NA, width=0.2, alpha = 0)

S2 <- ggplot(dd, aes(x=relevel(condition, "MD"), y=Dead, color=replicate, group=replicate)) +
  theme_light() +
  xlab(element_blank()) + # Remove the x axis title
  ylab("Number of dead bees") +
  geom_line() +
  geom_boxplot(aes(x=relevel(condition, "MD"), group=condition), outlier.shape=NA, width=0.2, alpha = 0) 

ggarrange(S1,S2, nrow=1, widths=c(1,1), common.legend = T)
```

![](AutomatedTracking_PlotsAndStats_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggsave(height=2,width=2.5,dpi=300, filename="SupplementTrackingFigs.pdf", useDingbats=FALSE)
```

## Additional statistical tests

``` r
##############################################################################
## Paired t test for speed SD
data_wide <- dcast(dd, replicate ~ condition, value.var="speed_sd")

# Shapiro-Wilk normality test for the difference between paired values
d <- with(data_wide, 
          CL - MD)
shapiro.test(d)# W = 0.94911, p-value = 0.6803
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  d
    ## W = 0.94911, p-value = 0.6803

``` r
# Compute t-test
res <- t.test(data_wide$CL, data_wide$MD, paired = TRUE, alternative =  "two.sided")
res # t = 0.3349, df = 8, p-value = 0.7463
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  data_wide$CL and data_wide$MD
    ## t = 0.3349, df = 8, p-value = 0.7463
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -11.33426  15.18577
    ## sample estimates:
    ## mean of the differences 
    ##                1.925755

``` r
##############################################################################
## Paired t test for mortality
data_wide <- dcast(dd, replicate ~ condition, value.var="Dead")

# Shapiro-Wilk normality test for the difference between paired values
d <- with(data_wide, 
          CL - MD)
shapiro.test(d)# W = 0.98709, p-value = 0.9907
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  d
    ## W = 0.98709, p-value = 0.9907

``` r
# Compute t-test
res <- t.test(data_wide$CL, data_wide$MD, paired = TRUE, alternative =  "two.sided")
res # t = 1.1011, df = 8, p-value = 0.3029
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  data_wide$CL and data_wide$MD
    ## t = 1.1011, df = 8, p-value = 0.3029
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -4.13407 11.68963
    ## sample estimates:
    ## mean of the differences 
    ##                3.777778

``` r
sessionInfo()
```

    ## R version 4.1.0 (2021-05-18)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Mojave 10.14.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] reshape2_1.4.4 ggpubr_0.4.0   ggplot2_3.3.5 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.7       plyr_1.8.6       highr_0.9        pillar_1.6.4    
    ##  [5] compiler_4.1.0   tools_4.1.0      digest_0.6.28    evaluate_0.14   
    ##  [9] lifecycle_1.0.1  tibble_3.1.6     gtable_0.3.0     pkgconfig_2.0.3 
    ## [13] rlang_0.4.12     DBI_1.1.1        yaml_2.2.1       xfun_0.28       
    ## [17] fastmap_1.1.0    gridExtra_2.3    withr_2.4.2      stringr_1.4.0   
    ## [21] dplyr_1.0.7      knitr_1.36       generics_0.1.1   vctrs_0.3.8     
    ## [25] cowplot_1.1.1    grid_4.1.0       tidyselect_1.1.1 glue_1.5.0      
    ## [29] R6_2.5.1         rstatix_0.7.0    fansi_0.5.0      rmarkdown_2.11  
    ## [33] carData_3.0-4    farver_2.1.0     car_3.0-12       tidyr_1.1.4     
    ## [37] purrr_0.3.4      magrittr_2.0.1   backports_1.3.0  scales_1.1.1    
    ## [41] ellipsis_0.3.2   htmltools_0.5.2  abind_1.4-5      assertthat_0.2.1
    ## [45] colorspace_2.0-2 ggsignif_0.6.3   labeling_0.4.2   utf8_1.2.2      
    ## [49] stringi_1.7.5    munsell_0.5.0    broom_0.7.10     crayon_1.4.2
