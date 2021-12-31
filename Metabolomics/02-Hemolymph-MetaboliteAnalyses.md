Metabolite differences in the hemolymph of gnotobiotic honeybees
================
Andrew Quinn, University of Lausanne

# Data Collection:

## MSTFA derivatized brain extracts analyzed via GC-MS (8890-5977B, Agilent) operated in scan mode

# Import Sample info and data from MassHunter Quantitative Analysis (Agilent)

## Set missing values to LOD = 1 and remove sample with failed derivatization or little hemolyte recovery

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

hemo_results_wide <- read_csv("hemolymph_data_trim.csv") %>%
  replace(is.na(.), 1)

sample_info <- read_csv("Sample_Info_hemo.csv")

hemo_results_wide <- full_join(hemo_results_wide,sample_info)

hemo_results <- hemo_results_wide %>%
  pivot_longer(
   cols =  c("Propylene glycol":"Lactose"),
   names_to = "Metabolite",
   values_to = "Response",
   values_drop_na = FALSE
  ) 

hemo_results$Type <- factor(hemo_results$Type) 
hemo_results$Colonization <- factor(hemo_results$Colonization)
hemo_results$Replicate <- factor(hemo_results$Replicate)

hemo_results <- hemo_results %>%
  filter(Name != "H32" & Name != "H41" & Name != "H78" & Name != "H104" &
           Name != "H113" & Name != "H122" & Name != "H88")
```

# Remove outlier data.

## 1) Group by replicate and filter out samples \< mean - 2*sd \| \> mean + 2*sd

# Normalize Response

## 1) Normalize to ISTD (Norleucine)

## 2) Convert normalized response to z-score

``` r
Internal_std <- filter(hemo_results, Metabolite == "Norleucine") %>%
  select(Name, Type, Replicate, Response) %>%
  rename(Norleucine = Response) %>%
  filter(Type == "Sample") %>%
  mutate(error = 2*sd(Norleucine),
         LB = median(Norleucine) - error,
         UB = median(Norleucine) + error,
         Decision = case_when(
           Norleucine > LB & Norleucine < UB ~ "ok",
           TRUE ~ "error"
         ))

hemo_results_norm <-  merge(hemo_results, Internal_std)%>%
  mutate(Resp_norm = Response / Norleucine) %>%
  filter(Decision == "ok")

hemo_results_norm <- hemo_results_norm %>%
  filter(Colonization != "None") %>%
  group_by(Metabolite) %>%
  mutate(z_score = (Resp_norm - mean(Resp_norm))/sd(Resp_norm))

hemolymph_results_output <- hemo_results_norm %>%
  filter(Metabolite != "Norleucine")
```

# Test for differentially abundant metabolites using mixed linear models

## Implement mixed linear model of lmm2met

### Variables: Fixed (colonization and injection order), Random (batch)

``` r
metabdf <- hemolymph_results_output %>%
  select(Name, Metabolite) %>%
  filter(Name == "H1")
metabdf$metabID<-1:nrow(metabdf)
metabdf$metabID = paste0('X', metabdf$metabID)
metabdf$metabID <- factor(metabdf$metabID)

metabdf <- metabdf %>%
  select(Metabolite, metabID)

hemolymph_lmm <- hemolymph_results_output %>%
  select(Name, Replicate, Colonization, Injection, Metabolite, z_score) %>%
  full_join(metabdf) %>%
  ungroup() %>%
  group_by(metabID) %>%
  select(Name, Replicate, Colonization, Injection, metabID, z_score) %>%
  pivot_wider(names_from = metabID, values_from = z_score)


hemolymph_lmm$Colonization <- factor(hemolymph_lmm$Colonization, levels = c("MD", "CL"))
hemolymph_lmm$Replicate <- factor(hemolymph_lmm$Replicate)
hemolymph_lmm$Injection <- as.integer(hemolymph_lmm$Injection)

library(lmm2met)
fitMet = fitLmm(fix= c('Colonization','Injection'), random='(1|Replicate)', data=hemolymph_lmm, start=5, end=80)

plot(fitMet, type='coeff')

plot(fitMet, type='residual') 

plot(fitMet, type='randeff') 

plot(fitMet, type='chisq')

getFixCoeff(fitMet, save=TRUE) 

list2env(fitMet,envir = .GlobalEnv)

fittedDat_long <- fittedDat %>%
  pivot_longer(cols =  c("X1":"X76"),
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
```

# Generate volcano plot of effect size vs p-value

``` r
library(ggrepel)

#dev.new()
#pdf("hemolymph_lmm_effect_volcano-viridis.pdf",  width=4, height=5, useDingbats = F)
ggplot(data=lmm_output, aes(x=ColonizationCL, y=Log10_p, col=Significance, label = siglabel)) + 
  geom_point() + theme_classic() + 
  geom_text_repel(size = 3, max.overlaps = Inf)+
  geom_vline(xintercept=0, col="grey", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey", linetype = "dashed") +
  geom_hline(yintercept=-log10(0.01), col="grey", linetype = "dotdash") +
  geom_hline(yintercept=-log10(0.001), col="grey", linetype = "dotted") +
  scale_x_continuous(name="Effect Size (colonization)",limits=c(-0.8,0.8), breaks=seq(-0.8, 0.8, 0.2), labels = scales::comma)+
  scale_y_continuous(name="-Log10 (p-Value)",limits=c(0,6), breaks=seq(0, 6, 1)) +
  scale_color_viridis_d(option="plasma")
```

![](Hemolymph_MetaboliteAnalyses_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
#dev.off()
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
    ##  [1] ggrepel_0.9.1      lmm2met_1.0        scales_1.1.1       writexl_1.4.0     
    ##  [5] lubridate_1.8.0    ggpubr_0.4.0       EnvStats_2.4.0     RColorBrewer_1.1-2
    ##  [9] readxl_1.3.1       forcats_0.5.1      stringr_1.4.0      dplyr_1.0.7       
    ## [13] purrr_0.3.4        readr_2.1.0        tidyr_1.1.4        tibble_3.1.6      
    ## [17] ggplot2_3.3.5      tidyverse_1.3.1   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-153       bitops_1.0-7       fs_1.5.0           bit64_4.0.5       
    ##  [5] webshot_0.5.2      httr_1.4.2         tools_4.1.0        backports_1.3.0   
    ##  [9] utf8_1.2.2         R6_2.5.1           KernSmooth_2.23-20 mgcv_1.8-38       
    ## [13] DBI_1.1.1          colorspace_2.0-2   withr_2.4.2        tidyselect_1.1.1  
    ## [17] gridExtra_2.3      bit_4.0.4          compiler_4.1.0     cli_3.1.0         
    ## [21] rvest_1.0.2        xml2_1.3.2         labeling_0.4.2     caTools_1.18.2    
    ## [25] systemfonts_1.0.2  digest_0.6.28      minqa_1.2.4        rmarkdown_2.11    
    ## [29] svglite_2.0.0      pkgconfig_2.0.3    htmltools_0.5.2    lme4_1.1-27.1     
    ## [33] highr_0.9          dbplyr_2.1.1       fastmap_1.1.0      rlang_0.4.12      
    ## [37] rstudioapi_0.13    farver_2.1.0       generics_0.1.1     jsonlite_1.7.2    
    ## [41] gtools_3.9.2       vroom_1.5.6        car_3.0-12         magrittr_2.0.1    
    ## [45] kableExtra_1.3.4   Matrix_1.3-4       Rcpp_1.0.7         munsell_0.5.0     
    ## [49] fansi_0.5.0        abind_1.4-5        lifecycle_1.0.1    stringi_1.7.5     
    ## [53] yaml_2.2.1         carData_3.0-4      MASS_7.3-54        gplots_3.1.1      
    ## [57] grid_4.1.0         parallel_4.1.0     crayon_1.4.2       lattice_0.20-45   
    ## [61] haven_2.4.3        splines_4.1.0      hms_1.1.1          knitr_1.36        
    ## [65] pillar_1.6.4       boot_1.3-28        ggsignif_0.6.3     reprex_2.0.1      
    ## [69] glue_1.5.0         evaluate_0.14      modelr_0.1.8       nloptr_1.2.2.3    
    ## [73] vctrs_0.3.8        tzdb_0.2.0         cellranger_1.1.0   gtable_0.3.0      
    ## [77] assertthat_0.2.1   xfun_0.28          broom_0.7.10       rstatix_0.7.0     
    ## [81] viridisLite_0.4.0  ellipsis_0.3.2
