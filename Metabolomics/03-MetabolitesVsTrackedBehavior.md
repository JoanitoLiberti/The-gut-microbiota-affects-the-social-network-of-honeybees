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

dt <- subset(dt, !is.na(HH_norm)) # we filter out samples for which we could not retrieve the number of head to head interactions (tag was missing)

# Functions to run multiple mixed effects models, correct p values and generate a table of the results
dep<-"HH_norm"
indep<-names(dt)[9:ncol(dt)]

run_anova <- function(independent_variable, dependent_variable, data) {
  
  formula<-as.formula(paste0(dependent_variable, "~", independent_variable))
  test_results <- lme(formula, random=~1|Subcolony, data=data)
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
  geom_point(aes(color = Colonization), shape = 16, size = 1, show.legend = T, alpha = 1) +
  theme_bw() +
  labs(y= "Normalised head-to-head interactions") +
  geom_smooth(aes(group=Subcolony, color = Colonization), se=F,method=lm, size=1) +
  scale_color_manual(values = c("#C77CFF", "#00BFC4")) +
  labs(colour = "Treatment")
}


library(purrr)
plot_list <- rownames(resTab[resTab$`Adjusted p-value` < 0.05,,drop=F]) %>% 
  map( ~ plot_for_loop(dt, .x,"HH_norm", "Subcolony"))

# Combine all plots
library(ggpubr)
ggarrange(plotlist = plot_list,
          nrow=1, common.legend = T)
# ggsave(height=3,width=14,filename="/Volumes/gr_Engel/jliberti/BeeTracking/Metabolomics/Figures/MetabolitesVsHH-2.pdf", useDingbats=FALSE)
```

### Analysis of the hemolymph metabolome

``` r
library(nlme)
library(ggplot2)
library(dplyr)

dt = read.table("/Users/joanitoliberti/The-gut-microbiota-affects-the-social-network-of-honeybees/Metabolomics/hemolymph_zscores_wMetadata.txt", header = T, fill=TRUE, sep="\t", na.strings=c(""," ","NA")) # fill=TRUE allows to read a table with missing entries
rownames(dt) = dt$Sample_ID

dt <- subset(dt, !is.na(HH_norm)) # we filter out samples for which we could not retrieve the number of head to head interactions (tag was missing)

dep<-"HH_norm"
indep<-names(dt)[9:ncol(dt)]


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

### Export results to Excel

``` r
library(xlsx)
write.xlsx(resTab, file="/Volumes/gr_Engel/jliberti/BeeTracking/Metabolomics/ProcessedData/MetabolitesVsHHnorm.xlsx", sheetName="Brain", row.names=TRUE)
write.xlsx(resTab2, file="/Volumes/gr_Engel/jliberti/BeeTracking/Metabolomics/ProcessedData/MetabolitesVsHHnorm.xlsx", sheetName="Hemolymph", append=TRUE, row.names=TRUE)
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
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.1.0  magrittr_2.0.1  fastmap_1.1.0   tools_4.1.0    
    ##  [5] htmltools_0.5.2 yaml_2.2.1      stringi_1.7.6   rmarkdown_2.11 
    ##  [9] knitr_1.36      stringr_1.4.0   xfun_0.28       digest_0.6.29  
    ## [13] rlang_0.4.12    evaluate_0.14
