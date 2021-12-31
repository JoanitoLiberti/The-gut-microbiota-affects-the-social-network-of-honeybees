Metabolite differences in the brains of gnotobiotic honeybees
================
Andrew Quinn, University of Lausanne

## Data Collection:

### MSTFA derivatized brain extracts analyzed via GC-MS (8890-5977B, Agilent) operated in scan mode

## Import Sample info and data from MassHunter Quantitative Analysis (Agilent)

### Set missing values to LOD = 1 and remove sample with failed derivatization

``` r
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(EnvStats)
library(ggpubr)
library(lubridate)
library(stats)
library(writexl)
library(scales)
wd <- "/Volumes/gr_Engel/aquinn/Joanito_collaboration/Publication_data"
setwd(wd)

brain_results_wide <- read_csv("Brain_data_trim.csv") %>%
  replace(is.na(.), 1)

sample_info <- read_csv("Sample_Info_brain.csv")

brain_results_wide <- full_join(brain_results_wide,sample_info)

brain_results <- brain_results_wide %>%
  pivot_longer(
   cols =  c("Lactic Acid":"Lactose"),
   names_to = "Metabolite",
   values_to = "Response",
   values_drop_na = FALSE
  ) 

brain_results$Type <- factor(brain_results$Type) 
brain_results$Colonization <- factor(brain_results$Colonization)
brain_results$Replicate <- factor(brain_results$Replicate)

brain_results <- brain_results %>%
  filter(Name != "B124")
```

## Check that Brain mass is not different between MD and CL bees.

### Use a linear mixed model to fit the data with colonization (fixed) and batch (random)

``` r
library(lme4)

lmm_b1 <-lmer(Mass ~ Colonization +(1|Replicate),data=subset(brain_results_wide, Type =="Sample"))
summary(lmm_b1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: Mass ~ Colonization + (1 | Replicate)
    ##    Data: subset(brain_results_wide, Type == "Sample")
    ## 
    ## REML criterion at convergence: -108.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -5.1480 -0.4770  0.1033  0.6153  2.4033 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance  Std.Dev.
    ##  Replicate (Intercept) 0.0003334 0.01826 
    ##  Residual              0.0299706 0.17312 
    ## Number of obs: 180, groups:  Replicate, 9
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error t value
    ## (Intercept)     1.39556    0.01924  72.546
    ## ColonizationMD -0.01667    0.02581  -0.646
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## ColoniztnMD -0.671

``` r
confint(lmm_b1)
```

    ##                      2.5 %     97.5 %
    ## .sig01          0.00000000 0.05801600
    ## .sigma          0.15581756 0.19245338
    ## (Intercept)     1.35776649 1.43334462
    ## ColonizationMD -0.06738439 0.03405106

``` r
ggplot(subset(brain_results_wide, Type %in% c("Sample")), mapping = aes(x = Colonization, y = Mass, color = Colonization)) +
  geom_boxplot(scale = "count", na.rm = TRUE, alpha = 0.5, trim = TRUE, aes(x = Colonization, y = Mass)) +
  geom_dotplot(aes(x = Colonization, y = Mass, color = Colonization, fill = Colonization),binaxis='y', stackdir='center', dotsize = 0.33, binwidth = 0.05, na.rm = TRUE) +
  stat_n_text(na.rm = TRUE, size = 2) +
  geom_hline(yintercept=0.8, col="grey50", linetype = "dashed") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        legend.position = "none", axis.text.x = element_text(color = "black", size=7)) +
  scale_y_continuous(name="Brain Mass (mg)", limits = c(0,2), breaks = seq(0,2,0.1)) +
  stat_compare_means(label = "p.format", label.x = 1.25, label.y = 2, size = 2)+
  facet_wrap(~Replicate, nrow=1)
```

![](Brain_MetaboliteAnalyses_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Remove outlier data.

### 1) Group by replicate and filter out samples \< mean - 2*sd \| \> mean + 2*sd

### 2) Remove brain samples with mass \< 0.8 mg

## Normalize Response

### 1) Normalize to ISTD (Norleucine) and brain mass

### 2) Convert normalized response to z-score

``` r
Internal_std <- filter(brain_results, Metabolite == "Norleucine") %>%
  select(Name, Type, Mass, Replicate, Response) %>%
  rename(Norleucine = Response) %>%
  #group_by(Replicate) %>%
  filter(Type == "Sample") %>%
  mutate(error = 2*sd(Norleucine),
         LB = median(Norleucine) - error,
         UB = median(Norleucine) + error,
         Decision = case_when(
           Norleucine > LB & Norleucine < UB & Mass > 0.8 ~ "ok",
           TRUE ~ "error"
         ))

brain_results_norm <-  merge(brain_results, Internal_std)%>%
  mutate(Resp_norm = Response / Norleucine,
         Resp_mnorm = Resp_norm / Mass) %>%
  filter(Decision == "ok")

brain_results_norm <- brain_results_norm %>%
  filter(Colonization != "None") %>%
  group_by(Metabolite) %>%
  mutate(z_score = (Resp_mnorm - mean(Resp_mnorm))/sd(Resp_mnorm))

brain_results_output <- brain_results_norm %>%
  filter(Metabolite != "Norleucine")
```

## Test for differentially abundant metabolites using mixed linear models

### Implement mixed linear model of lmm2met

#### Variables: Fixed (colonization and injection order), Random (batch)

``` r
metabdf <- brain_results_output %>%
  select(Name, Metabolite) %>%
  filter(Name == "B1")
metabdf$metabID<-1:nrow(metabdf)
metabdf$metabID = paste0('X', metabdf$metabID)
metabdf$metabID <- factor(metabdf$metabID)

metabdf <- metabdf %>%
  select(Metabolite, metabID)

brain_lmm <- brain_results_output %>%
  select(Name, Replicate, Colonization, Injection, Metabolite, z_score) %>%
  full_join(metabdf) %>%
  ungroup() %>%
  group_by(metabID) %>%
  select(Name, Replicate, Colonization, Injection, metabID, z_score) %>%
  pivot_wider(names_from = metabID, values_from = z_score)


brain_lmm$Colonization <- factor(brain_lmm$Colonization, levels = c("MD", "CL"))
brain_lmm$Replicate <- factor(brain_lmm$Replicate)
brain_lmm$Injection <- as.integer(brain_lmm$Injection)

library(lmm2met)
fitMet = fitLmm(fix= c('Colonization','Injection'), random='(1|Replicate)', data=brain_lmm, start=5, end=64)

plot(fitMet, type='coeff')

plot(fitMet, type='residual') 

plot(fitMet, type='randeff') 

plot(fitMet, type='chisq')

getFixCoeff(fitMet, save=TRUE) 

list2env(fitMet,envir = .GlobalEnv)

fittedDat_long <- fittedDat %>%
  pivot_longer(cols =  c("X1":"X60"),
               names_to = "Metabolite",
               values_to = "FittedDat",
               values_drop_na = FALSE)

lmm_output_csv <- read.csv("output.csv", check.names = F)
lmm_output <- lmm_output_csv %>% 
  dplyr::rename(p = "Pr(Chi).Colonization", metabID = "") %>%
  mutate(p_adj = p.adjust(p, method="BH",n = length(p)),
    Log10_p = -log10(p_adj)) %>%
  left_join(metabdf)


lmm_output$Significance <- "NS"

# if pvalue < 0.05, 0.01, 0.001 set as "*", "**", "***" 

lmm_output$Significance[lmm_output$p_adj < 0.05] <- "*"
lmm_output$Significance[lmm_output$p_adj < 0.01] <- "**"
lmm_output$Significance[lmm_output$p_adj < 0.001] <- "***"
lmm_output$Significant <- "No"
lmm_output$Significant[lmm_output$p_adj < 0.05] <- "Yes"
lmm_output$siglabel <- NA
lmm_output$Metabolite <- as.character(lmm_output$Metabolite)
lmm_output$siglabel[lmm_output$Significance == "*"] <- lmm_output$Metabolite[lmm_output$Significance == "*"]
lmm_output$siglabel[lmm_output$Significance == "**"] <- lmm_output$Metabolite[lmm_output$Significance == "**"]
lmm_output$siglabel[lmm_output$Significance == "***"] <- lmm_output$Metabolite[lmm_output$Significance == "***"]

Metabolite_info <- read_csv("Metabolite_Info.csv") %>%
  filter(Tissue == "Brain") %>%
  select(Metabolite, Class)

lmm_output <- lmm_output %>%
  left_join(Metabolite_info)
```

## Generate volcano plot of effect size vs p-value

``` r
library(ggrepel)

#dev.new()
#pdf("brain_lmm_effect_volcano-viridis.pdf",  width=4, height=5, useDingbats = F)
ggplot(data=lmm_output, aes(x=ColonizationCL, y=Log10_p, col=Significance, label = siglabel)) + 
  geom_point() + theme_classic() + 
  geom_text_repel(size = 3, max.overlaps = Inf)+
  geom_vline(xintercept=0, col="grey", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.01), col="grey", linetype = "dotdash") +
  geom_hline(yintercept=-log10(0.001), col="grey", linetype = "dotted") +
  scale_x_continuous(name="Effect Size (colonization)",limits=c(-0.6,0.6), breaks=seq(-0.6, 0.6, 0.2), labels = scales::comma)+
  scale_y_continuous(name="-Log10 (p-Value)",limits=c(0,6), breaks=seq(0, 6, 1)) +
  scale_color_viridis_d(option="plasma")
```

![](Brain_MetaboliteAnalyses_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
#dev.off()
```

## Generate barplot of significant metabolites by metabolite class

``` r
lmm_output$Class <- factor(lmm_output$Class, levels = c("Non-Essential Amino Acid", "Essential Amino Acid", "Amino Acid Metabolism",
"Carboxylic Acid", "Carbohydrate", "Nucleobase Metabolism", "Other"))

ggplot(data=lmm_output, aes(x=Significant))+
  geom_bar(aes(fill = Class),position = position_stack(reverse = TRUE))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        legend.position = "right", axis.text.x = element_text(color = "black", size=10), 
        axis.text.y = element_text(color = "black", size=10))+
  scale_fill_brewer(palette = "YlGnBu", direction = -1) +
  scale_y_continuous(name="Count",limits=c(0,40), breaks=seq(0, 40, 5))
```

![](Brain_MetaboliteAnalyses_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
#ggsave("Sig_metabs_bar_2.pdf")
```

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
    ##  [1] ggrepel_0.9.1      lmm2met_1.0        lme4_1.1-27.1      Matrix_1.3-4      
    ##  [5] scales_1.1.1       writexl_1.4.0      lubridate_1.8.0    ggpubr_0.4.0      
    ##  [9] EnvStats_2.4.0     RColorBrewer_1.1-2 readxl_1.3.1       forcats_0.5.1     
    ## [13] stringr_1.4.0      dplyr_1.0.7        purrr_0.3.4        readr_2.1.0       
    ## [17] tidyr_1.1.4        tibble_3.1.6       ggplot2_3.3.5      tidyverse_1.3.1   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-153       bitops_1.0-7       fs_1.5.0           bit64_4.0.5       
    ##  [5] webshot_0.5.2      httr_1.4.2         tools_4.1.0        backports_1.3.0   
    ##  [9] utf8_1.2.2         R6_2.5.1           KernSmooth_2.23-20 mgcv_1.8-38       
    ## [13] DBI_1.1.1          colorspace_2.0-2   withr_2.4.2        gridExtra_2.3     
    ## [17] tidyselect_1.1.1   bit_4.0.4          compiler_4.1.0     cli_3.1.0         
    ## [21] rvest_1.0.2        xml2_1.3.2         labeling_0.4.2     caTools_1.18.2    
    ## [25] systemfonts_1.0.2  digest_0.6.28      minqa_1.2.4        svglite_2.0.0     
    ## [29] rmarkdown_2.11     pkgconfig_2.0.3    htmltools_0.5.2    dbplyr_2.1.1      
    ## [33] fastmap_1.1.0      highr_0.9          rlang_0.4.12       rstudioapi_0.13   
    ## [37] generics_0.1.1     farver_2.1.0       jsonlite_1.7.2     gtools_3.9.2      
    ## [41] vroom_1.5.6        car_3.0-12         magrittr_2.0.1     kableExtra_1.3.4  
    ## [45] Rcpp_1.0.7         munsell_0.5.0      fansi_0.5.0        abind_1.4-5       
    ## [49] lifecycle_1.0.1    stringi_1.7.5      yaml_2.2.1         carData_3.0-4     
    ## [53] MASS_7.3-54        gplots_3.1.1       grid_4.1.0         parallel_4.1.0    
    ## [57] crayon_1.4.2       lattice_0.20-45    haven_2.4.3        splines_4.1.0     
    ## [61] hms_1.1.1          knitr_1.36         pillar_1.6.4       boot_1.3-28       
    ## [65] ggsignif_0.6.3     reprex_2.0.1       glue_1.5.0         evaluate_0.14     
    ## [69] modelr_0.1.8       vctrs_0.3.8        nloptr_1.2.2.3     tzdb_0.2.0        
    ## [73] cellranger_1.1.0   gtable_0.3.0       assertthat_0.2.1   xfun_0.28         
    ## [77] broom_0.7.10       rstatix_0.7.0      viridisLite_0.4.0  ellipsis_0.3.2
