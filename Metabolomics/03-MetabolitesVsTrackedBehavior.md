Linear mixed models testing the effect of individual metabolites on
tracked social interactions of individual bees
================
Joanito Liberti, University of Lausanne

### Testing associations between brain and hemolymph metabolites and the behavioral outcomes under the tracking system

``` r
library(nlme)
library(ggplot2)
library(dplyr)

dt = read.table("/Users/joanitoliberti/The-gut-microbiota-affects-the-social-network-of-honeybees/Metabolomics/brain_zscores_wMetadata.txt", header = T, fill=TRUE, sep="\t", na.strings=c(""," ","NA")) # fill=TRUE allows to read a table with missing entries
rownames(dt) = dt$Sample_ID
dt$Colonization<-factor(dt$Colonization,levels = c("MD", "CL"))

dt <- subset(dt, !is.na(HH)) # we filter out samples for which we could not retrieve the number of head to head interactions (tag was missing)

# Functions to run multiple mixed effects models, correct p values and generate a table of the results
dep<-"HH"
indep<-names(dt)[6:ncol(dt)]

run_anova <- function(independent_variable, dependent_variable, data) {
  
  formula<-as.formula(paste0(dependent_variable, "~", independent_variable))
  test_results <- lme(formula, random=~1|Replicate, data=data)
  anova_results <- anova(test_results)
  output = list(f_stat = anova_results[2, 3], p_value = anova_results[2, 4])
  
}

generate_table_P <- function(independent_variables, dependent_variables, data) {
  
  output <- matrix(ncol = length(independent_variables), nrow = length(dependent_variables))
  colnames(output) <- independent_variables
  rownames(output) <- dependent_variables
  
  for (independent_variable in independent_variables) {
    
    for (dependent_variable in dependent_variables) {
      
      temp <- run_anova(independent_variable, dependent_variable, data)
      output[dependent_variable, independent_variable] <- temp$p_value
      
    }
  }
  
  return(output)
  
}


generate_table_F <- function(independent_variables, dependent_variables, data) {
  
  output <- matrix(ncol = length(independent_variables), nrow = length(dependent_variables))
  colnames(output) <- independent_variables
  rownames(output) <- dependent_variables
  
  for (independent_variable in independent_variables) {
    
    for (dependent_variable in dependent_variables) {
      
      temp <- run_anova(independent_variable, dependent_variable, data)
      output[dependent_variable, independent_variable] <- temp$f_stat
      
    }
  }
  
  return(output)
  
}

metabolitePvalue = generate_table_P(indep, dep, dt)
metaboliteFvalue = generate_table_F(indep, dep, dt)

padj <- matrix(p.adjust(as.vector(as.matrix(metabolitePvalue)), method='fdr'),ncol=length(colnames(metabolitePvalue)))
colnames(padj) <- colnames(metabolitePvalue)
rownames(padj) <- rownames(metabolitePvalue)
resTab<-as.data.frame(cbind(t(metaboliteFvalue),t(metabolitePvalue),t(padj)))
colnames(resTab)<-c("F-value", "p-value", "Adjusted p-value")
resTab<-resTab[order(resTab$`Adjusted p-value`),,drop=F] 
resTab
```

    ##                                F-value      p-value Adjusted p-value
    ## X3.Amino.2.piperidone     12.781812221 0.0004704575       0.01398933
    ## Phosphorylethanolamine    13.235208871 0.0003765033       0.01398933
    ## Tyrosine                  12.013312955 0.0006879996       0.01398933
    ## Ornithine_2               10.689491098 0.0013345324       0.01705442
    ## Serine                    10.597576393 0.0013979033       0.01705442
    ## Glycerol.3.phosphate_2     9.843451875 0.0020496297       0.02083790
    ## B.Alanine                  8.886223392 0.0033502639       0.02901025
    ## Myo.Inositol_2             8.640608754 0.0038046228       0.02901025
    ## Glycerol.3.phosphate_1     7.945010816 0.0054686238       0.03706512
    ## Aspartic.Acid              5.278763197 0.0229605618       0.14005943
    ## Glycine                    4.263785257 0.0406432092       0.16961552
    ## Erythritol                 4.225296566 0.0415481620       0.16961552
    ## Citric.Acid                4.228726969 0.0414666441       0.16961552
    ## Uracil                     4.264522967 0.0406260698       0.16961552
    ## Phosphate                  4.218560263 0.0417087347       0.16961552
    ## Isoleucine                 3.527106684 0.0623012149       0.23752338
    ## Malic.Acid                 3.270778232 0.0725138258       0.25876343
    ## Inosine                    3.184308217 0.0763564220       0.25876343
    ## Ornithine_1                3.060431823 0.0822520143       0.26082997
    ## Glutamic.Acid              2.995904628 0.0855180235       0.26082997
    ## Threitol                   2.723467948 0.1009614868       0.29326908
    ## X2.8.Dihydroxyquinoline    2.088870980 0.1504478208       0.37406157
    ## Valine                     1.832693388 0.1778325493       0.37406157
    ## Pipecolic.Acid             2.038379997 0.1554377742       0.37406157
    ## X5.Oxoproline              2.043949408 0.1548781690       0.37406157
    ## X4.Aminobutanoic.Acid      1.967416190 0.1627750525       0.37406157
    ## Asparagine                 1.921364933 0.1677497394       0.37406157
    ## Putrescine                 1.985175387 0.1609021911       0.37406157
    ## Glutamine                  1.836468504 0.1773891887       0.37406157
    ## Sorbitol                   1.685932688 0.1961175426       0.39877234
    ## Lactic.Acid                1.408828411 0.2371149561       0.46658104
    ## Proline                    1.224960802 0.2701493487       0.49936698
    ## Lactose                    1.261938732 0.2630672683       0.49936698
    ## Glucose                    1.119210717 0.2917779554       0.52348398
    ## Succinic.Acid              0.949970680 0.3312855768       0.57517375
    ## Hydroxyproline             0.867382424 0.3531676300       0.57517375
    ## Fructose                   0.910271231 0.3415673008       0.57517375
    ## Acetylglutamic.Acid        0.849000778 0.3583049620       0.57517375
    ## Phenylalanine              0.788977545 0.3758221920       0.58782445
    ## Glycerol                   0.660569281 0.4176386065       0.63643880
    ## Threonine                  0.632274493 0.4277703377       0.63643880
    ## Galactose                  0.445662748 0.5054208655       0.73406364
    ## Alanine                    0.388103728 0.5342376141       0.75787196
    ## Taurine                    0.319578440 0.5727000848       0.77632678
    ## Myo.Inositol_1             0.327553099 0.5679545937       0.77632678
    ## X2.4.Di.tert.butylphenoxy  0.238086923 0.6262990161       0.83052696
    ## Sucrose                    0.197249867 0.6575857162       0.85346231
    ## Glucose.6.phosphate        0.171901427 0.6790158708       0.86291600
    ## Timonacic.Acid             0.134828236 0.7139916565       0.87106982
    ## Unknown_6_94               0.140322980 0.7084856940       0.87106982
    ## Cysteine                   0.078160303 0.7801887462       0.88367800
    ## Norleucine                 0.104361481 0.7471040866       0.88367800
    ## Unknown_16_62              0.070137871 0.7914973115       0.88367800
    ## Adenosine                  0.063311357 0.8016788639       0.88367800
    ## Fumaric.Acid               0.076231609 0.7828482252       0.88367800
    ## Methylmaleic.Acid          0.057234869 0.8112453736       0.88367800
    ## Sarcosine                  0.039358705 0.8430065701       0.89933502
    ## Unknown_11_73              0.033460070 0.8551054250       0.89933502
    ## Glucitol                   0.024599843 0.8755778038       0.90525841
    ## Urea                       0.009317497 0.9232298793       0.93861704
    ## Gluconic.Acid              0.004367263 0.9473973503       0.94739735

### Plots

``` r
dt$Replicate<-factor(dt$Replicate)
plot_for_loop <- function(df, .x_var, .y_var, .z_var) {

  # convert strings to variable
  x_var <- sym(.x_var)
  y_var <- sym(.y_var)
  z_var <- sym(.z_var)

  # unquote variables using !! 
  ggplot(df, aes(x = !! x_var, y = !! y_var, color = !! z_var)) + 
  geom_point(aes(), shape = 16, size = 1, show.legend = T, alpha = 1) +
  theme_bw() +
  labs(y= "Head to head interactions") +
  geom_smooth(se=F,method=lm, size=1)
}


library(purrr)
plot_list <- rownames(resTab[resTab$`Adjusted p-value` < 0.05,,drop=F]) %>% 
  map( ~ plot_for_loop(dt, .x,"HH", "Replicate"))

# Combine all plots
library(ggpubr)
ggarrange(plotlist = plot_list,
          nrow=2,ncol = 5, common.legend = T)
```

![](03-MetabolitesVsTrackedBehavior_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# ggsave(height=5,width=12,filename="/Volumes/gr_Engel/jliberti/BeeTracking/Metabolomics/Figures/MetabolitesVsHH.pdf", useDingbats=FALSE)
```

### Analysis of the hemolymph metabolome

``` r
library(nlme)
library(ggplot2)
library(dplyr)

dt = read.table("/Users/joanitoliberti/The-gut-microbiota-affects-the-social-network-of-honeybees/Metabolomics/hemolymph_zscores_wMetadata.txt", header = T, fill=TRUE, sep="\t", na.strings=c(""," ","NA")) # fill=TRUE allows to read a table with missing entries
rownames(dt) = dt$Sample_ID

dt <- subset(dt, !is.na(HH)) # we filter out samples for which we could not retrieve the number of head to head interactions (tag was missing)

dep<-"HH"
indep<-names(dt)[6:ncol(dt)]


metabolitePvalue = generate_table_P(indep, dep, dt)

metaboliteFvalue = generate_table_F(indep, dep, dt)

padj <- matrix(p.adjust(as.vector(as.matrix(metabolitePvalue)), method='fdr'),ncol=length(colnames(metabolitePvalue)))
colnames(padj) <- colnames(metabolitePvalue)
rownames(padj) <- rownames(metabolitePvalue)
resTab2<-as.data.frame(cbind(t(metaboliteFvalue),t(metabolitePvalue),t(padj)))
colnames(resTab2)<-c("F-value", "p-value", "Adjusted p-value")
resTab2<-resTab2[order(resTab2$`Adjusted p-value`),,drop=F] 
resTab2
```

    ##                                       F-value     p-value Adjusted p-value
    ## Threitol                         9.8134901513 0.002088633        0.1587361
    ## Succinic.Acid                    5.0561452455 0.026015693        0.3291692
    ## B.Alanine                        4.5458485447 0.034649392        0.3291692
    ## Erythritol                       6.0791145919 0.014821558        0.3291692
    ## Tromethamine                     6.3525189828 0.012781805        0.3291692
    ## Asparagine                       4.6171050100 0.033281204        0.3291692
    ## Unknown.16.6246                  4.6349176071 0.032948182        0.3291692
    ## Glutamine                        5.2933461892 0.022803731        0.3291692
    ## Sarcosine                        4.0245037638 0.046665851        0.3414007
    ## Serine                           4.1954045788 0.042300985        0.3414007
    ## Pipecolic.Acid                   3.8881272485 0.050492874        0.3414007
    ## X5.Oxoproline                    3.7755284179 0.053905369        0.3414007
    ## Glucitol                         3.2081852227 0.075314377        0.4402994
    ## Leucine                          2.2606918397 0.134825671        0.6404219
    ## Fructose                         2.3781385715 0.125179545        0.6404219
    ## Lactose                          2.2657334872 0.134394787        0.6404219
    ## Valine                           2.1173059142 0.147760279        0.6605754
    ## Fumaric.Acid                     1.5679829009 0.212475390        0.7020926
    ## X3.Hydroxyproline                1.8816032557 0.172227906        0.7020926
    ## X4.Aminobutanoic.Acid            1.5849385402 0.210032865        0.7020926
    ## Methionine                       1.6487390925 0.201137120        0.7020926
    ## Adenine                          1.7560142951 0.187164591        0.7020926
    ## Inosine                          1.7623932834 0.186370352        0.7020926
    ## Glycerol                         1.4256863602 0.234379970        0.7044270
    ## Malic.Acid                       1.2810968524 0.259525751        0.7044270
    ## Aspartic.Acid                    1.3208296909 0.252297092        0.7044270
    ## Unknown.19.5708                  1.5003712183 0.222560199        0.7044270
    ## Tyrosine                         1.3195989918 0.252517216        0.7044270
    ## Unknown.20.7569                  1.2296393108 0.269275865        0.7056885
    ## Phenylalanine                    0.9913122820 0.321046912        0.7624864
    ## Unknown.17.1507                  1.0297886875 0.311864400        0.7624864
    ## Galacto.Hexodialdose             1.0173693362 0.314790033        0.7624864
    ## Urea                             0.7655662030 0.383009898        0.7652926
    ## Isoleucine                       0.7775919923 0.379308504        0.7652926
    ## Ornithine                        0.8810912673 0.349432323        0.7652926
    ## Ribitol                          0.7348053375 0.392715937        0.7652926
    ## Phosphorylethanolamine           0.8210223129 0.366354237        0.7652926
    ## Galactose.oxime                  0.7400128654 0.391048125        0.7652926
    ## Uric.Acid                        0.7737248513 0.380493194        0.7652926
    ## Proline                          0.6998924485 0.404168592        0.7679203
    ## Arabinofuranose                  0.4940302939 0.483239092        0.7739150
    ## Homocysteine                     0.5132506200 0.474864013        0.7739150
    ## Glucose                          0.6040897766 0.438262387        0.7739150
    ## Sorbitol                         0.6335088211 0.427345947        0.7739150
    ## Myo.inositol_1                   0.6156331373 0.433928775        0.7739150
    ## Maltose                          0.4815952778 0.488788429        0.7739150
    ## Myo.inositol_2                   0.5071187036 0.477509892        0.7739150
    ## Adenosine                        0.5770906602 0.448663293        0.7739150
    ## Glycine                          0.4421085250 0.507141273        0.7835534
    ## Taurine                          0.4249354213 0.515495627        0.7835534
    ## Putrescine                       0.3615759844 0.548552175        0.8081657
    ## Glycerol.3.phosphate             0.3536621472 0.552955446        0.8081657
    ## Cadaverine                       0.3220155036 0.571258680        0.8191634
    ## Glutamic.Acid                    0.2659778921 0.606812257        0.8540321
    ## X2.3.Butanediol                  0.0625888736 0.802796882        0.9406087
    ## Lactic.Acid                      0.0819209344 0.775111151        0.9406087
    ## Alanine                          0.1292804530 0.719691771        0.9406087
    ## Methylsuccinimide                0.0724611932 0.788161334        0.9406087
    ## Unknown.4.5703                   0.1484318716 0.700591589        0.9406087
    ## Diethylene.glycol..n.butyl.ether 0.0497550457 0.823797627        0.9406087
    ## Unknown.6.94                     0.0466917288 0.829220826        0.9406087
    ## X3.Aminoisobutyric.Acid          0.1192918928 0.730294169        0.9406087
    ## Timonacic.Acid                   0.0962648668 0.756795318        0.9406087
    ## Citric.Acid                      0.1531622738 0.696094379        0.9406087
    ## Altrose                          0.0516252680 0.820572496        0.9406087
    ## Gluconic.Acid                    0.0615448618 0.804414552        0.9406087
    ## Sucrose                          0.0569235713 0.811756184        0.9406087
    ## X2.Pyrrolidinone..1.methy        0.0346381962 0.852610686        0.9423392
    ## para.Isopropylbenzoic.Acid       0.0332574934 0.855544773        0.9423392
    ## Propylene.glycol                 0.0274432272 0.868650624        0.9431064
    ## X2.Oxoglutaric.Acid              0.0206738057 0.885866720        0.9482517
    ## X2.8.Dihydroxyquinoline          0.0160475757 0.899366387        0.9493312
    ## Phosphate                        0.0051920431 0.942654758        0.9698588
    ## Threonine                        0.0048915332 0.944336228        0.9698588
    ## X3.Pyridinol                     0.0001167140 0.991394832        0.9913948
    ## Gluconolactone                   0.0001505417 0.990227089        0.9913948

### Export results to Excel

``` r
library(xlsx)
write.xlsx(resTab, file="/Volumes/gr_Engel/jliberti/BeeTracking/Metabolomics/ProcessedData/MetabolitesVsHH.xlsx", sheetName="Brain", row.names=TRUE)
write.xlsx(resTab2, file="/Volumes/gr_Engel/jliberti/BeeTracking/Metabolomics/ProcessedData/MetabolitesVsHH.xlsx", sheetName="Hemolymph", append=TRUE, row.names=TRUE)
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
    ## [1] ggpubr_0.4.0  purrr_0.3.4   dplyr_1.0.7   ggplot2_3.3.5 nlme_3.1-153 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] highr_0.9        pillar_1.6.4     compiler_4.1.0   tools_4.1.0     
    ##  [5] digest_0.6.29    evaluate_0.14    lifecycle_1.0.1  tibble_3.1.6    
    ##  [9] gtable_0.3.0     lattice_0.20-45  mgcv_1.8-38      pkgconfig_2.0.3 
    ## [13] rlang_0.4.12     Matrix_1.3-4     DBI_1.1.1        yaml_2.2.1      
    ## [17] xfun_0.28        fastmap_1.1.0    gridExtra_2.3    withr_2.4.3     
    ## [21] stringr_1.4.0    knitr_1.36       generics_0.1.1   vctrs_0.3.8     
    ## [25] cowplot_1.1.1    grid_4.1.0       tidyselect_1.1.1 glue_1.6.0      
    ## [29] R6_2.5.1         rstatix_0.7.0    fansi_1.0.2      rmarkdown_2.11  
    ## [33] carData_3.0-4    farver_2.1.0     car_3.0-12       tidyr_1.1.4     
    ## [37] magrittr_2.0.1   splines_4.1.0    backports_1.3.0  scales_1.1.1    
    ## [41] ellipsis_0.3.2   htmltools_0.5.2  abind_1.4-5      assertthat_0.2.1
    ## [45] colorspace_2.0-2 ggsignif_0.6.3   labeling_0.4.2   utf8_1.2.2      
    ## [49] stringi_1.7.6    munsell_0.5.0    broom_0.7.10     crayon_1.4.2
