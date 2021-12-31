Linear mixed models testing the effect of individual metabolites on
tracked social interactions of individual bees
================
Joanito Liberti, University of Lausanne

### Testing associations between brain and hemolymph metabolites and the behavioral outcomes under the tracking system

``` r
library(nlme)
library(ggplot2)
library(dplyr)

dt = read.table("/Users/joanitoliberti/The-gut-microbiota-affects-the-social-network-of-honey-bees/Metabolomics/brain_zscores_wMetadata.txt", header = T, fill=TRUE, sep="\t", na.strings=c(""," ","NA")) # fill=TRUE allows to read a table with missing entries
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

![](MetabolitesVsTrackedBehavior_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# ggsave(height=5,width=12,filename="/Volumes/gr_Engel/jliberti/BeeTracking/Metabolomics/Figures/MetabolitesVsHH.pdf", useDingbats=FALSE)
```

### Analysis of the hemolymph metabolome

``` r
library(nlme)
library(ggplot2)
library(dplyr)

dt = read.table("/Users/joanitoliberti/The-gut-microbiota-affects-the-social-network-of-honey-bees/Metabolomics/hemolymph_zscores_wMetadata.txt", header = T, fill=TRUE, sep="\t", na.strings=c(""," ","NA")) # fill=TRUE allows to read a table with missing entries
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

    ##                                       F-value      p-value Adjusted p-value
    ## Threitol                         1.148644e+01 0.0009061623       0.06886834
    ## Erythritol                       7.014820e+00 0.0089921477       0.18414239
    ## Tromethamine                     6.461343e+00 0.0120898486       0.18414239
    ## Unknown.16.6246                  6.457536e+00 0.0121146312       0.18414239
    ## Glutamine                        7.313734e+00 0.0076738405       0.18414239
    ## Asparagine                       5.547119e+00 0.0198697895       0.25168400
    ## Succinic.Acid                    4.829297e+00 0.0295890604       0.27841324
    ## B.Alanine                        4.636810e+00 0.0329699895       0.27841324
    ## X5.Oxoproline                    5.056524e+00 0.0260622604       0.27841324
    ## Sarcosine                        4.050687e+00 0.0460307930       0.28814058
    ## Serine                           4.189292e+00 0.0425114905       0.28814058
    ## Pipecolic.Acid                   3.932144e+00 0.0492872048       0.28814058
    ## Lactose                          4.011969e+00 0.0470684013       0.28814058
    ## Sorbitol                         3.508186e+00 0.0631076150       0.34258420
    ## Valine                           3.017857e+00 0.0845053523       0.39433549
    ## Leucine                          2.899993e+00 0.0907524626       0.39433549
    ## X3.Hydroxyproline                3.143197e+00 0.0783732032       0.39433549
    ## Glucitol                         2.852782e+00 0.0933952486       0.39433549
    ## Glycerol                         1.683802e+00 0.1965093696       0.59738848
    ## Proline                          1.867619e+00 0.1738929457       0.59738848
    ## X4.Aminobutanoic.Acid            1.733417e+00 0.1900813578       0.59738848
    ## Methionine                       2.006602e+00 0.1587904681       0.59738848
    ## Unknown.19.5708                  1.748289e+00 0.1882030146       0.59738848
    ## Fructose                         1.904388e+00 0.1697417016       0.59738848
    ## Tyrosine                         1.861015e+00 0.1746508432       0.59738848
    ## Isoleucine                       1.484309e+00 0.2251078488       0.62671476
    ## Fumaric.Acid                     1.484464e+00 0.2250836329       0.62671476
    ## Unknown.20.7569                  1.397281e+00 0.2391411578       0.62671476
    ## Inosine                          1.443034e+00 0.2316348723       0.62671476
    ## Malic.Acid                       1.256768e+00 0.2641428489       0.65364120
    ## Phenylalanine                    1.207086e+00 0.2737564228       0.65364120
    ## Galacto.Hexodialdose             1.157396e+00 0.2838178894       0.65364120
    ## Adenine                          1.177156e+00 0.2797614923       0.65364120
    ## Unknown.17.1507                  1.104578e+00 0.2950363907       0.65949311
    ## Urea                             1.018668e+00 0.3145399385       0.66402876
    ## Aspartic.Acid                    1.055353e+00 0.3060119423       0.66402876
    ## Ribitol                          8.758406e-01 0.3509211932       0.70184239
    ## Phosphorylethanolamine           9.052214e-01 0.3429908372       0.70184239
    ## Ornithine                        7.948733e-01 0.3741282948       0.72858895
    ## Galactose.oxime                  7.642490e-01 0.3834678698       0.72858895
    ## Cadaverine                       7.306855e-01 0.3940909289       0.73051001
    ## Arabinofuranose                  5.535506e-01 0.4580913448       0.79124869
    ## Homocysteine                     5.754555e-01 0.4493476516       0.79124869
    ## Uric.Acid                        5.571868e-01 0.4566213642       0.79124869
    ## Glucose                          5.262403e-01 0.4693772913       0.79272609
    ## Lactic.Acid                      3.659292e-01 0.5461903038       0.80445093
    ## Glutamic.Acid                    3.406442e-01 0.5603770315       0.80445093
    ## Taurine                          4.003344e-01 0.5279270592       0.80445093
    ## Putrescine                       3.215218e-01 0.5715835544       0.80445093
    ## Glycerol.3.phosphate             4.601627e-01 0.4986441894       0.80445093
    ## Myo.inositol_1                   3.484020e-01 0.5559507331       0.80445093
    ## Maltose                          4.575587e-01 0.4998597500       0.80445093
    ## Myo.inositol_2                   3.369398e-01 0.5625145777       0.80445093
    ## Adenosine                        3.651342e-01 0.5466260722       0.80445093
    ## Timonacic.Acid                   2.646652e-01 0.6077272685       0.82608338
    ## Sucrose                          2.632424e-01 0.6086930165       0.82608338
    ## X2.3.Butanediol                  2.249761e-01 0.6359977297       0.84799697
    ## Unknown.4.5703                   1.238192e-01 0.7254465678       0.93447354
    ## Gluconic.Acid                    1.334663e-01 0.7154046286       0.93447354
    ## Propylene.glycol                 1.857326e-02 0.8917886395       0.93624070
    ## X2.Pyrrolidinone..1.methy        1.862044e-02 0.8916521336       0.93624070
    ## Alanine                          6.431923e-02 0.8001593175       0.93624070
    ## Methylsuccinimide                5.144061e-02 0.8208991018       0.93624070
    ## X3.Pyridinol                     5.578688e-02 0.8136212880       0.93624070
    ## Glycine                          9.085418e-02 0.7635319075       0.93624070
    ## Diethylene.glycol..n.butyl.ether 4.072087e-02 0.8403644574       0.93624070
    ## Unknown.6.94                     7.295495e-02 0.7874713903       0.93624070
    ## Threonine                        3.780357e-02 0.8461140407       0.93624070
    ## X3.Aminoisobutyric.Acid          2.106607e-02 0.8848034945       0.93624070
    ## X2.8.Dihydroxyquinoline          4.310797e-02 0.8358174877       0.93624070
    ## Citric.Acid                      2.337944e-02 0.8786900139       0.93624070
    ## Gluconolactone                   1.607599e-02 0.8992838310       0.93624070
    ## Altrose                          4.031526e-02 0.8411507376       0.93624070
    ## X2.Oxoglutaric.Acid              7.555929e-03 0.9308528468       0.95601103
    ## Phosphate                        2.491746e-03 0.9602578715       0.97306131
    ## para.Isopropylbenzoic.Acid       4.508047e-04 0.9830900551       0.98309006

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
    ##  [5] digest_0.6.28    evaluate_0.14    lifecycle_1.0.1  tibble_3.1.6    
    ##  [9] gtable_0.3.0     lattice_0.20-45  mgcv_1.8-38      pkgconfig_2.0.3 
    ## [13] rlang_0.4.12     Matrix_1.3-4     DBI_1.1.1        yaml_2.2.1      
    ## [17] xfun_0.28        fastmap_1.1.0    gridExtra_2.3    withr_2.4.2     
    ## [21] stringr_1.4.0    knitr_1.36       generics_0.1.1   vctrs_0.3.8     
    ## [25] cowplot_1.1.1    grid_4.1.0       tidyselect_1.1.1 glue_1.5.0      
    ## [29] R6_2.5.1         rstatix_0.7.0    fansi_0.5.0      rmarkdown_2.11  
    ## [33] carData_3.0-4    farver_2.1.0     car_3.0-12       tidyr_1.1.4     
    ## [37] magrittr_2.0.1   splines_4.1.0    backports_1.3.0  scales_1.1.1    
    ## [41] ellipsis_0.3.2   htmltools_0.5.2  abind_1.4-5      assertthat_0.2.1
    ## [45] colorspace_2.0-2 ggsignif_0.6.3   labeling_0.4.2   utf8_1.2.2      
    ## [49] stringi_1.7.5    munsell_0.5.0    broom_0.7.10     crayon_1.4.2
