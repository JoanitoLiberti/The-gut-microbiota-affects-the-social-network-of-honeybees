Automated behavioral tracking experiment: gut microbiota analyses
================
Joanito Liberti, University of Lausanne

## Number of reads per sample

``` bash
cd /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/RawFiles/
for f in *fastq.gz; do
echo $f
echo $(zless $f | wc -l)/4|bc
done
```

## Check data quality with FastQC

``` bash
mkdir /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis
mkdir /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/00_FASTQC/
  cd /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/RawFiles
for f in *.fastq.gz; do 
fastqc $f -o /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/00_FASTQC/; 
done
```

## Clean up sample names

``` bash
cd /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/RawFiles
for f in *fastq.gz; do mv $f ${f/_001_*/.fastq.gz}; done
for f in *fastq.gz; do mv $f ${f//_/.}; done
for f in H2O*; do mv $f ${f/./_}; done
for f in CLH*; do mv $f ${f/./_}; done
```

## Remove primers and anything before and after with Cutadapt

``` bash
mkdir /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/02_dada2
cd /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/RawFiles
  for f in *R1.fastq.gz; do 
/Users/joanitoliberti/Library/Python/2.7/bin//cutadapt -G GGACTACHVGGGTWTCTAAT -g GTGCCAGCMGCCGCGGTAA -o ${f/.fastq/.cut.fastq} -p ${f/R1.fastq/R2.cut.fastq} $f ${f/R1.fastq/R2.fastq} --pair-filter=any; done
```

## Move trimmed reads into a separate folder

``` bash
mkdir /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/01_trimmed
cd /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/RawFiles
for f in *cut.fastq.gz; do mv $f /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/01_trimmed/$f; done
```

## Number of reads per sample after cutadapt

``` bash
cd /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/RawFiles/

echo -e "Sample""\t""Raw""\t""Cutadapt"

for f in $(ls *R1.fastq.gz | cut -d '.' -f 1|sort|uniq); do
raw=`echo $(zless /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/RawFiles/$f.L1.R1.fastq.gz | wc -l)/4|bc`
cut=`echo $(zless /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/01_trimmed/$f.L1.R1.cut.fastq.gz | wc -l)/4|bc`

echo -e $f"\t"$raw"\t"$cut
done
```

## Now we can feed the quality-filtered reads to the DADA2 pipeline

## First, install required packages if needed (uncomment all installation lines)

``` r
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")
# BiocManager::install("dada2", version = "3.10")
# BiocManager::install("DECIPHER")
# BiocManager::install("phyloseq")
# BiocManager::install("biomformat")
# BiocManager::install("ShortRead")
# BiocManager::install("Biostrings")
# BiocManager::install("genefilter")
# BiocManager::install("decontam")
# BiocManager::install("DESeq2")

# install.packages("ggplot2")
# install.packages("phangorn")
# install.packages("vegan")
# install.packages("dplyr")
# install.packages("scales")
# install.packages("RColorBrewer")
# install.packages(“reshape2”)
# install.packages("cowplot")
# install.packages("tidyverse")
# install.packages("readxl")
# install.packages("ggbeeswarm")
# install.packages("magrittr")
# install.packages("ggpubr")
# install.packages('dendextend')
# install.packages("multcomp")
# install.packages("GUniFrac")

# Load libraries
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(biomformat)
```

## Set the working directory

``` r
path <- "/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/02_dada2"
knitr::opts_knit$set(root.dir = normalizePath(path)) 
```

``` r
library(dada2); packageVersion("dada2")
pathraw <-"/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/01_trimmed"

list.files(pathraw)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.trim.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(pathraw, pattern="R1.cut.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(pathraw, pattern="R2.cut.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)
sample.names
```

## Quality scores

``` r
path<-"/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/01_trimmed"
dev.new()
# Quality scores of R1 reads
plotQualityProfile(fnFs[1:20]) 
plotQualityProfile(fnFs[21:40])
plotQualityProfile(fnFs[41:60])
plotQualityProfile(fnFs[61:80])
plotQualityProfile(fnFs[81:100])
plotQualityProfile(fnFs[101:120])
plotQualityProfile(fnFs[121:140])
plotQualityProfile(fnFs[141:160])
plotQualityProfile(fnFs[161:180])
plotQualityProfile(fnFs[181:194])

# Quality scores of R2 reads
plotQualityProfile(fnRs[1:20]) 
plotQualityProfile(fnRs[21:40])
plotQualityProfile(fnRs[41:60])
plotQualityProfile(fnRs[61:80])
plotQualityProfile(fnRs[81:100])
plotQualityProfile(fnRs[101:120])
plotQualityProfile(fnRs[121:140])
plotQualityProfile(fnRs[141:160])
plotQualityProfile(fnRs[161:180])
plotQualityProfile(fnRs[181:194])
```

## Trim the data

``` r
# Place filtered files in filtered/ subdirectory
path<-"/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/02_dada2"
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

sys_str <- Sys.time()
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(232,231), # truncLen[[1]] + truncLen[[2]] > amplicon_length+25
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,      
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
sys_str[2] <- Sys.time()
sys_str
rm(sys_str)

derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

head(out)
out
```

## Learn the error rates

``` r
sys_str <- Sys.time()
set.seed(1) # Since we randomize what samples we use for learning the error rates, we use set.seed() to make the analysis reproducible
errF <- learnErrors(filtFs, randomize=TRUE, nbases=3e8, multithread=TRUE) # here we need to increase nbases to sample more than only a few samples, default is 1e8
set.seed(1)
errR <- learnErrors(filtRs, randomize=TRUE, nbases=3e8, multithread=TRUE) # here we need to increase nbases to sample more than only a few samples, default is 1e8
plotErrors(errF, nominalQ=TRUE)

# In the plots, the black line is the error model, the dots are the actual errors
sys_str[2] <- Sys.time()
sys_str
rm(sys_str)
```

## Sample inference

``` r
sys_str <- Sys.time()
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
sys_str[2] <- Sys.time()
sys_str
rm(sys_str)

dadaFs[[1]]
dadaRs[[1]]
```

## Merge paired reads

``` r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, trimOverhang=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

## Construct sequence table

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

## Remove sequences that appear too short or long to be real

``` r
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:255]
```

## Remove chimeras

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)
```

## Track reads through the pipeline

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)
```

## Assign taxonomy

``` r
#assignTaxonomy using Silva 
taxa <- assignTaxonomy(seqtab.nochim, "/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/silva_species_assignment_v138.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

``` bash
mkdir /Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/03_Taxonomy
```

``` r
pathtax <- "/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/03_Taxonomy/"

# Save taxonomy assignments from dada2 pipeline
write.csv2(file=paste(pathtax,"Taxtable_dada2",sep=""),taxa)
```

# Convert to fasta

``` r
uniquesToFasta(seqtab.nochim, fout='/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/03_Taxonomy/seqtab.nochim.fasta', ids=colnames(seqtab.nochim))
```

# Prepare a phylogenetic tree with phangorn

``` r
library(DECIPHER) # Make alignment with DECIPHER...
library(phangorn) # ...then build tree with phangorn
asv.seqs <- getSequences(seqtab.nochim )
names(asv.seqs) <- asv.seqs
asv.alignment <- AlignSeqs(DNAStringSet(asv.seqs), anchor=NA)

# convert to phangorn and run phylogenetic tree reconstruction
phang.align <- phyDat(as(asv.alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
```

## Evaluate accuracy

``` r
unqs.mock <- seqtab.nochim["MOCK",] 
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/mock_seq_plasmids.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

## Load and filter the data in phyloseq

``` r
library(ggplot2)
library(vegan) # ecological diversity analysis
library(dplyr)
library(scales) # scale functions for vizualizations
library(grid)
library(RColorBrewer)
library(reshape2) # data manipulation package
library(cowplot)
library(phyloseq)
library(tidyverse)
library(readxl)

# setting the working directory
setwd("/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/")

# Set plotting theme
theme_set(theme_bw())

#Data frame containing sample information
samdf = read.table("/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/GutMicrobiota_BeeTracking_metadata.txt", header = T, fill=TRUE, sep="\t") # fill=TRUE allows to read a table with missing entries
rownames(samdf) = samdf$Sample_ID

#Create a phyloseq object
ps.raw <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa), 
               phy_tree(fitGTR$tree))
sample_data(ps.raw) 

# Make phyloseq objects containing Mitochondria and Chloroplast ASVs so we can see how many there were
ps.mit <- subset_taxa(ps.raw, Family=="Mitochondria", prunesamples=TRUE) # Depending on the classifier used, double check if your Chloroplasts or Mitochondria are in the Order and Family levels...
tax_table(ps.mit)
ps.chl <- subset_taxa(ps.raw, Order=="Chloroplast", prunesamples=TRUE)
tail(tax_table(ps.chl))

dput(sample_names(ps)) # Check sample names
# As some sample names started with 0, it creates problems in some of the downstream functions, so we change the sample names within the ps object and modify the mapping file accordingly
samples = c("S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10", 
"S100", "S101", "S102", "S103", "S104", "S105", "S106", "S107", "S108", 
"S109", "S11", "S110", "S111", "S112", "S113", "S114", "S115", "S116", 
"S117", "S118", "S119", "S12", "S120", "S121", "S122", "S123", "S124", 
"S125", "S126", "S127", "S128", "S129", "S13", "S130", "S131", "S132", 
"S133", "S134", "S135", "S136", "S137", "S138", "S139", "S14", "S140", 
"S141", "S142", "S143", "S144", "S145", "S146", "S147", "S148", "S149", 
"S15", "S150", "S151", "S152", "S153", "S154", "S155", "S156", "S157", 
"S158", "S159", "S16", "S160", "S161", "S162", "S163", "S164", "S165", 
"S166", "S167", "S168", "S169", "S17", "S170", "S171", "S172", "S173", 
"S174", "S175", "S176", "S177", "S178", "S179", "S18", "S180", "S181", 
"S182", "S183", "S184", "S185", "S186", "S187", "S188", "S189", "S19", 
"S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27", "S28", "S29", "S30", 
"S31", "S32", "S33", "S34", "S35", "S36", "S37", "S38", "S39", "S40", "S41", 
"S42", "S43", "S44", "S45", "S46", "S47", "S48", "S49", "S50", "S51", "S52", 
"S53", "S54", "S55", "S56", "S57", "S58", "S59", "S60", "S61", "S62", "S63", 
"S64", "S65", "S66", "S67", "S68", "S69", "S70", "S71", "S72", "S73", "S74", 
"S75", "S76", "S77", "S78", "S79", "S80", "S81", "S82", "S83", "S84", "S85", 
"S86", "S87", "S88", "S89", "S90", "S91", "S92", "S93", "S94", "S95", "S96", 
"S97", "S98", "S99", "CLH_inoc", "H2O_1", "H2O_2", "H2O_3", "MOCK")

sample_names(ps) = samples
sample_names(ps)

samdf$Sample_ID<-c("S01","S02","S03","S04","S05","S06","S07","S08","S09","S10","S11",
                   "S12","S13","S14","S15","S16","S17","S18","S19","S20","S21","S22",
                   "S23","S24","S25","S26","S27","S28","S29","S30","S31","S32","S33",
                   "S34","S35","S36","S37","S38","S39","S40","S41","S42","S43","S44",
                   "S45","S46","S47","S48","S49","S50","S51","S52","S53","S54","S55",
                   "S56","S57","S58","S59","S60","S61","S62","S63","S64","S65","S66",
                   "S67","S68","S69","S70","S71","S72","S73","S74","S75","S76","S77",
                   "S78","S79","S80","S81","S82","S83","S84","S85","S86","S87","S88",
                   "S89","S90","S91","S92","S93","S94","S95","S96","S97","S98","S99",
                   "S100","S101","S102","S103","S104","S105","S106","S107","S108","S109","S110",
                   "S111","S112","S113","S114","S115","S116","S117","S118","S119","S120","S121",
                   "S122","S123","S124","S125","S126","S127","S128","S129","S130","S131","S132",
                   "S133","S134","S135","S136","S137","S138","S139","S140","S141","S142","S143",
                   "S144","S145","S146","S147","S148","S149","S150","S151","S152","S153","S154",
                   "S155","S156","S157","S158","S159","S160","S161","S162","S163","S164","S165",
                   "S166","S167","S168","S169","S170","S171","S172","S173","S174","S175","S176",
                   "S177","S178","S179","S180","S181","S182","S183","S184","S185","S186","S187",
                   "S188","S189","MOCK","H2O_1","H2O_3","CLH_inoc","H2O_2")
rownames(samdf) = samdf$Sample_ID
```

## Export raw ASV table

``` r
table = merge( tax_table(ps.raw),t(otu_table(ps.raw)), by="row.names")

write.table(table, "/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/RawASVtable_dada2taxa.txt", sep="\t")
```

## Filter the table and clean the classification

``` r
#Filter out Eukaryota, mitochondria and chloroplast reads
ps <- subset_taxa(ps.raw,
                      Kingdom !="Eukaryota" | is.na(Kingdom) &
                      Kingdom !="Unclassified" | is.na(Kingdom), 
                      prunesamples=TRUE)

# give new names to ASV's
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Instead of having NAs in the taxa names, at all levels of taxonomical assignment we replace them with the lowest level that the classifier could achieve
tax <- data.frame(tax_table(ps))

tax.clean <- data.frame(row.names = row.names(tax),
Kingdom = str_replace(tax[,1], "D_0__",""),
Phylum = str_replace(tax[,2], "D_1__",""),
Class = str_replace(tax[,3], "D_2__",""),
Order = str_replace(tax[,4], "D_3__",""),
Family = str_replace(tax[,5], "D_4__",""),
Genus = str_replace(tax[,6], "D_5__",""),
Species = str_replace(tax[,7], "D_6__",""), # uncheck this if using the Silva classification from DADA2, as that one also has the species column
stringsAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])} 
####### Fill holes in the tax table
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){

#Fill in missing taxonomy
if (tax.clean[i,2] == ""){
kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
tax.clean[i, 2:7] <- kingdom  
} else if (tax.clean[i,3] == ""){
phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
tax.clean[i, 3:7] <- phylum
} else if (tax.clean[i,4] == ""){
class <- paste("Class_", tax.clean[i,3], sep = "")
tax.clean[i, 4:7] <- class
} else if (tax.clean[i,5] == ""){
order <- paste("Order_", tax.clean[i,4], sep = "")
tax.clean[i, 5:7] <- order
} else if (tax.clean[i,6] == ""){
family <- paste("Family_", tax.clean[i,5], sep = "")
tax.clean[i, 6:7] <- family
} else if (tax.clean[i,6] == ""){
tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
}
}

tax_table(ps) <- as.matrix(tax.clean)
```

## Save the phyloseq object so next time you do not need to rerun the analysis from scratch

``` r
saveRDS(ps, "/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/ps_dada2taxa_BeeTracking.rds")

# Grab the XStringSet object with the OTU/ASV sequences
rs <- refseq(ps)
rs # inspect
# Get strings with the full taxonomy for each OTU/ASV
tax <- tax_table(ps)
head(tax) # inspect
tax_strings <- apply(tax, 1, paste, collapse=";")
head(tax_strings) # inspect
# Create new names for the sequences in the refseq object; these will become
# the sequence headers in the fasta file. Here, I set these to be the OTU/ASV
# name and the tax string separated by one space (determined by the sep
# parameter)
new_names <- paste(taxa_names(ps), tax_strings, sep = " ")
head(new_names) # inspect
# Update the refeq names and save as a fasta file
names(rs) <- new_names
rs # inspect
Biostrings::writeXStringSet(rs, "/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/03_Taxonomy/AllRawASVs.fasta")
```

## Blast all ASVs to complement the taxonomy

``` bash
cat AllRawASVs.fasta | parallel --block 2000k --recstart '>' --pipe /work/Joanito_new/Blast/ncbi-blast-2.11.0+/bin/blastn -query AllRawASVs.fasta -db /work/nucleotide/nt -outfmt \"6 std stitle\" -evalue 1e-3 -max_target_seqs 5 -num_threads 35 > AllRawASVs_Blast.txt
```

## Filter Blast results to retain only first hit of each ASV

``` r
bl <- read.table("/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/03_Taxonomy/AllRawASVs_Blast.txt", header=F, sep = "\t")
bl <- bl[!duplicated(bl$V1),]
write.table(bl, "/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/03_Taxonomy/AllRawASVs_Blast_OnlyFirstHit.txt", sep = "\t")
write.table(tax_table(ps), "/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/03_Taxonomy/TaxTable_ps.txt", sep = "\t")
```

## The analysis can be started from here if the ps object was saved

``` r
library(ggplot2)
library(vegan) # ecological diversity analysis
library(dplyr)
library(scales) # scale functions for vizualizations
library(grid)
library(RColorBrewer)
library(reshape2) # data manipulation package
library(cowplot)
library(phyloseq)
library(tidyverse)
library(readxl)
library(dada2)

ps <- readRDS("/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/ps_dada2taxa_BeeTracking.rds")
tax_tab <- read.table("/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/analysis/03_Taxonomy/TaxTable_ps_BeeTracking_Annotated.txt", header = T, sep = "\t") # read in the manually annotated taxonomy
row.names(tax_tab) <- tax_tab[,1]
tax_tab[,1] <- NULL
tax_table(ps) <- as.matrix(tax_tab)

#Data frame containing sample information
samdf = read.table("/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/GutMicrobiota_BeeTracking_metadata.txt", header = T, fill=TRUE,  sep="\t", na.strings=c(""," ","NA")) # fill=TRUE allows to read a table with missing entries

samdf$Sample_ID # Check the order is the same as the modified sample names in next lines
samdf$Sample_ID<-c("S01","S02","S03","S04","S05","S06","S07","S08","S09","S10","S11",
                   "S12","S13","S14","S15","S16","S17","S18","S19","S20","S21","S22",
                   "S23","S24","S25","S26","S27","S28","S29","S30","S31","S32","S33",
                   "S34","S35","S36","S37","S38","S39","S40","S41","S42","S43","S44",
                   "S45","S46","S47","S48","S49","S50","S51","S52","S53","S54","S55",
                   "S56","S57","S58","S59","S60","S61","S62","S63","S64","S65","S66",
                   "S67","S68","S69","S70","S71","S72","S73","S74","S75","S76","S77",
                   "S78","S79","S80","S81","S82","S83","S84","S85","S86","S87","S88",
                   "S89","S90","S91","S92","S93","S94","S95","S96","S97","S98","S99",
                   "S100","S101","S102","S103","S104","S105","S106","S107","S108","S109","S110",
                   "S111","S112","S113","S114","S115","S116","S117","S118","S119","S120","S121",
                   "S122","S123","S124","S125","S126","S127","S128","S129","S130","S131","S132",
                   "S133","S134","S135","S136","S137","S138","S139","S140","S141","S142","S143",
                   "S144","S145","S146","S147","S148","S149","S150","S151","S152","S153","S154",
                   "S155","S156","S157","S158","S159","S160","S161","S162","S163","S164","S165",
                   "S166","S167","S168","S169","S170","S171","S172","S173","S174","S175","S176",
                   "S177","S178","S179","S180","S181","S182","S183","S184","S185","S186","S187",
                   "S188","S189","MOCK","H2O_1","H2O_3","CLH_inoc","H2O_2")
rownames(samdf) = samdf$Sample_ID
```

## Or if you saved the whole session, you can load it like this

``` r
#load("/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/BeeTracking_Ampliseq.rdata")
```

## Plot distribution of reads per sample

``` r
# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(ps))

# Histogram of sample read counts
dev.new()
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# mean, max and min of sample read counts
smin <- min(sample_sums(ps))
smean <- mean(sample_sums(ps))
smax <- max(sample_sums(ps))

# printing the results
cat("The minimum sample read count is:",smin)
cat("The average sample read count is:",smean)
cat("The maximum sample read count is:",smax)

# setting the seed to one value in order to create reproducible results
set.seed(1)  
```

## Rarefaction curves

``` r
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 6
cols = gg_color_hue(n)

ps@sam_data$Sample_Cols<-cols[as.factor(ps@sam_data$Treatment)]

dev.new()
pdf(file="Rarefactions.pdf",width=5,height=4, pointsize=4)
rarecurve(otu_table(ps), ylab = "Number of ASVs", xlab = "Number of sequences", col= ps@sam_data$Sample_Cols, step= 50, cex=0.5)
dev.off()

rarecurve(otu_table(ps), ylab = "Number of ASVs", xlab = "Number of sequences", col= ps@sam_data$Sample_Cols, step= 50, cex=0.5)
```

## Retain only ASVs that have at least 1% relative abundance in minimum 1 sample

``` r
library(genefilter)
mothur_proportions = transform_sample_counts(ps, function(x) {x/sum(x)})
Filtered = filter_taxa(mothur_proportions, filterfun(kOverA(1, 0.01)), TRUE)
```

## Plot bacterial stacked bars

``` r
mdf = psmelt(Filtered)
mdf <- mdf[order(factor(mdf$Treatment, levels = c("MD", "CL", "inoculum", "blank", "H2O", "mock"), ordered=TRUE), mdf$Sample), ]  # Let's first sort by Treatment, then hive, then sample ID
mdf$Sample <- factor(mdf$Sample, levels = unique(mdf$Sample), ordered=TRUE)  # we need to do this to retain the order in plot.
tail(mdf, 4)

library(RColorBrewer)
cbPalette<-colorRampPalette(brewer.pal(12, "Paired"))

set.seed(7) # Here we use set.seed to make the colors reproducible, because in the scale_fill_manual we use "sample" to shuffle the order of colors in the palette
p1 = ggplot(mdf, aes_string(x = "Sample", y = "Abundance", fill = "Annotated")) +
        geom_bar(stat = "identity", position = "stack", color = "black") +
        ylab("Relative abundance") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position="bottom") +
        scale_fill_manual(values=sample(cbPalette(90))) 
  
print(p1, width = 1000, height = 200)
```

## Plot stacked bars of the various control samples (blanks, mock, negative PCR controls)

## Plotting only ASVs that have at least 0.1% relative abundance in 1 sample

``` r
ps.blanks <- subset_samples(ps, Treatment == "blank") 
ps.mock <- subset_samples(ps, Treatment == "mock")
ps.H2O <- subset_samples(ps, Treatment == "H2O")
ps.inoculum <- subset_samples(ps, Treatment == "inoculum")

library(genefilter)
proportions.blanks = transform_sample_counts(ps.blanks, function(x) {x/sum(x)})
Filtered.blanks = filter_taxa(proportions.blanks, filterfun(kOverA(1, 0.001)), TRUE)
proportions.mock = transform_sample_counts(ps.mock, function(x) {x/sum(x)})
Filtered.mock = filter_taxa(proportions.mock, filterfun(kOverA(1, 0.001)), TRUE)
proportions.H2O = transform_sample_counts(ps.H2O, function(x) {x/sum(x)})
Filtered.H2O = filter_taxa(proportions.H2O, filterfun(kOverA(1, 0.001)), TRUE)
proportions.inoculum = transform_sample_counts(ps.inoculum, function(x) {x/sum(x)})
Filtered.inoculum = filter_taxa(proportions.inoculum, filterfun(kOverA(1, 0.001)), TRUE)
```

``` r
mdf.blanks = psmelt(Filtered.blanks)
mdf.mock = psmelt(Filtered.mock)
mdf.H2O = psmelt(Filtered.H2O)
mdf.inoculum = psmelt(Filtered.inoculum)

# Make a color palette
library(RColorBrewer)
cbPalette<-colorRampPalette(brewer.pal(15, "Paired"))

p2 = ggplot(mdf.blanks, aes_string(x = "Sample", y = "Abundance", fill = "Annotated")) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
  ylab("Relative abundance in blank DNA") +
  xlab(element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position="bottom") +
  scale_fill_manual(values=cbPalette(100)) 

p3 = ggplot(mdf.mock, aes_string(x = "Sample", y = "Abundance", fill = "Annotated")) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
  ylab("Relative abundance in mock DNA") +
  xlab(element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position="bottom") +
  scale_fill_manual(values=cbPalette(9)) 

p4 = ggplot(mdf.H2O, aes_string(x = "Sample", y = "Abundance", fill = "Annotated")) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
  ylab("Relative abundance in negative PCR control") +
  xlab(element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position="bottom") +
  scale_fill_manual(values=cbPalette(44)) 

p5 = ggplot(mdf.inoculum, aes_string(x = "Sample", y = "Abundance", fill = "Annotated")) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
  ylab("Relative abundance in inocula") +
  xlab(element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position="bottom") +
  scale_fill_manual(values=cbPalette(23)) 


print(p2, width = 1000, height = 200)
print(p3, width = 1000, height = 200)
print(p4, width = 1000, height = 200)
print(p5, width = 1000, height = 200)
```

## Save plots to pdf

``` r
dev.new()
pdf(file="Barplots_allsample-BeeTrackings.pdf",width=12,height=6, pointsize=4)
p1
dev.off()

dev.new()
pdf(file="Barplots_ControlSamples-BeeTracking.pdf", width=14,height=25, pointsize=4)
plot_grid(p2, p3, p4, p5, labels = c("A", "B", "C", "D"), ncol = 1, nrow = 4, rel_heights = c(2,1,1.5,1))
dev.off()
```

## Distribution of library sizes

``` r
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Negative)) + geom_point() 
ggsave(file="DistributionLibrarySizes.pdf", width=8, height=6, useDingbats=FALSE)
```

## Check which sequences are likely to be contaminants (frequency method)

``` r
library(decontam); packageVersion("decontam")
contamdf.freq <- isContaminant(ps, method="frequency", conc="Picogreen")
head(contamdf.freq)
table(contamdf.freq$contaminant)

which(contamdf.freq$contaminant)
```

## Plot frequencies of identified contaminants vs. their concentrations after PCR

``` r
plot_frequency(ps, taxa_names(ps)[sample(which(contamdf.freq$contaminant), )], conc="Picogreen") + 
  xlab("DNA Concentration (PicoGreen ng/μl)")
ggsave(file="ContaminantFreqVsPicogreenConc.pdf", width=12, height=10, useDingbats=FALSE)
```

## Check which sequences are likely to be contaminants (prevalence method)

``` r
contamdf.prev <- isContaminant(ps, neg = "Negative", normalize = TRUE, threshold = 0.1, detailed = T)
table(contamdf.prev$contaminant)

head(which(contamdf.prev$contaminant))
```

## Plot prevalences of identified contaminants

``` r
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Negative == "TRUE", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Negative == "FALSE", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point(size=3, position=position_jitter(h=0.05,w=0.05)) +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ggsave(file="ContaminantPrevalences.pdf", width=8, height=6, useDingbats=FALSE)
```

## Both methods combined

``` r
contamdf.both <- isContaminant(ps, conc = "Picogreen", neg = "Negative", method = "either", normalize = TRUE, detailed = T)
table(contamdf.both$contaminant)

head(which(contamdf.both$contaminant))
```

## New phyloseq object after filtering out contaminants

``` r
ps.noncontam <- prune_taxa(!contamdf.both$contaminant, ps)
ps.noncontam
```

## Remove non experimental samples before replotting qPCR and MiSeq results

``` r
ps.noncontam.filt <- prune_samples(ps.noncontam@sam_data$Treatment == "MD" | ps.noncontam@sam_data$Treatment == "CL", ps.noncontam)
ps.noncontam.filt = filter_taxa(ps.noncontam.filt, function(x) sum(x) > 0, TRUE)
ps.noncontam.filt
ps.noncontam.filt.wInocula <- prune_samples(ps.noncontam@sam_data$Treatment == "MD" |  ps.noncontam@sam_data$Treatment == "CL" | ps.noncontam@sam_data$Treatment == "inoculum", ps.noncontam)
ps.noncontam.filt.wInocula = filter_taxa(ps.noncontam.filt.wInocula, function(x) sum(x) > 0, TRUE)
ps.noncontam.filt.wInocula
```

## Save filtered ASV table (not normalised by qPCR yet)

``` r
write.table(otu_table(ps.noncontam.filt), "BeeTracking_ASVtable_OnlyMiseq_Filtered.txt", sep="\t")
```

## Histogram of qPCR 16S rRNA copy numbers

``` r
mdf = sample_data(ps.noncontam.filt) 

mdf$Treatment <- factor(mdf$Treatment, levels = c("MD", "CL"))

p1 <- ggplot(mdf, aes(x = Sample_ID, y= CopyNum_norm)) + 
  geom_bar(stat="identity", color = "Black", fill = "indianred") +
  ylab("Normalised 16S rRNA gene copies") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_y_log10(limits=c(1,1e10), breaks = c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10)) + 
  coord_cartesian(ylim = c(1e5,1e10)) +
  facet_grid(. ~Treatment+Replicate, scale="free", space="fixed") +
  annotation_logticks(base = 10, sides = "l", scaled = TRUE,
                    short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"),
                    colour = "black", size = 0.5, linetype = 1, alpha = 1, color = NULL)
  
print(p1, width = 100, height = 20)
```

## Retain only OTUs that have at least 1% relative abundance in minimum 5 samples

``` r
library(genefilter)
mothur_proportions = transform_sample_counts(ps.noncontam.filt, function(x) {x/sum(x)})
Filtered = filter_taxa(mothur_proportions, filterfun(kOverA(5, 0.01)), TRUE)
```

## Plot bacterial stacked bars

``` r
mdf = psmelt(Filtered)

mdf$Treatment <- factor(mdf$Treatment, levels = c("MD", "CL"))

library(RColorBrewer)
cbPalette<-colorRampPalette(brewer.pal(15, "Paired"))

set.seed(5) # Change the number within parentheses to change the order of colors in the plot. Here we use set.seed to make the colors reproducible, because in scale_fill_manual we use "sample" to shuffle the color order.
p2 = ggplot(mdf, aes_string(x = "Sample_ID", y = "Abundance", fill = "Annotated")) +
        geom_bar(stat = "identity", position = "stack", color = "black") +
        ylab("Relative abundance") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position="bottom") +
        scale_fill_manual(values=sample(cbPalette(22)))

p3 = ggplot(mdf, aes_string(x = "Sample_ID", y = "Abundance", fill = "Annotated")) +
        geom_bar(stat = "identity", position = "stack", color = "black") +
        ylab("Relative abundance") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position="bottom", strip.background = element_blank()) + # ,  strip.text.x = element_blank()
        facet_grid(. ~ Treatment + Replicate, scale="free", space="fixed") +
        scale_fill_manual(values=sample(cbPalette(22)))  
print(p2, width = 1000, height = 200)
print(p3, width = 1000, height = 200)

# Get the ggplot grob
gt <- ggplotGrob(p3)

# Locate the tops of the plot panels
panels <- grep("panel", gt$layout$name)
top <- unique(gt$layout$t[panels])

# Remove the rows immediately above the plot panel
gt = gt[-(top-1), ]

# Draw it
grid.newpage()
grid.draw(gt)
```

## Plot qPCR and ampliseq bars together (and export to pdf)

``` r
dev.new()
pdf(file="Barplots_Experimentalsamples.pdf",width=12,height=7, pointsize=4)

plot_grid(p1, gt, ncol = 1, nrow = , rel_widths=c(2,2), rel_heights = c(2,6), align = 'v')

dev.off()

plot_grid(p1, gt, ncol = 1, nrow = , rel_widths=c(2,2), rel_heights = c(2,6), align = 'v')
```

## Let’s plot ordinations. First combine qPCR and MiSeq data

``` r
proportions = transform_sample_counts(ps.noncontam.filt, function(x) {x/sum(x)})

#calculate ratios
ratios = merge( tax_table(proportions),t(otu_table(proportions)), by="row.names")
copynum = data.frame(sample_data(proportions)[,"CopyNum_norm"])

row.names(ratios) <- ratios[,1]

ratios <- ratios[,11:ncol(ratios)]

copies = data.frame(row.names = (row.names(ratios)))
for (i in 1:ncol(ratios)){
  numbers = ratios[,i]*copynum$CopyNum[i]  # or copynum$CopyNum_norm[i]
  copies = cbind(copies,numbers)
}
colnames(copies) = colnames(ratios)

# Phyloseq object with the combined data 
ps_copies <- phyloseq(otu_table(copies, taxa_are_rows=T), 
               sample_data(samdf), 
               tax_table(proportions),
               phy_tree(proportions))
```

## Save qPCR normalised ASV table

``` r
write.table(otu_table(ps_copies), "BeeTracking_ASVtable_qPCRnormalised_Filtered.txt", sep="\t")
```

## Ordinations using qPCR-normalized counts for the whole dataset

``` r
set.seed(91193)
# Ordinate using Principal Coordinate analysis of Bray-Curtis dissimilarities (and make pdf file)
bee_pcoa <- ordinate(
  physeq = ps_copies, 
  method = "PCoA", 
  distance = "bray", 
  trymax = 100
)


# checking the ordination by making a scree plot
plot_ordination(
  physeq = ps_copies,
  ordination = bee_pcoa,
  type="scree")

# Plot PCoA
p1 <- plot_ordination(
  physeq = ps_copies,
  ordination = bee_pcoa,
  axes=c(1,2),   # this selects which axes to plot from the ordination
  color = "Treatment",
  title = "PCoA of Bray-Curtis dissimilarities"
)  +
 scale_color_manual(values = c("#00BFC4","#C77CFF")
 ) +
  theme_bw()+
  geom_point(aes(color = Treatment), alpha = 1, size = 3) 
ggsave(height=5,width=7,dpi=300, filename="Ordinations_qPCRnormalized_PCoA_Bray-BeeTracking.pdf", useDingbats=FALSE)
```

## Run a Permanova test with Adonis.

``` r
set.seed(8650)

sink("Adonis-Betadisper_WholeDataset.txt")
print("###############################################################")
print("Effect of microbiota treatment - Bray-Curtis dissimilarities")
# Calculate bray curtis distance matrix
bray <- phyloseq::distance(ps_copies, method = "bray", weighted=TRUE)

# make a data frame from the scaled sample_data
sampledf <- data.frame(sample_data(ps_copies))

# Adonis test
adonis(bray ~ Treatment, data = sampledf)

# test of Homegeneity of dispersion
beta <- betadisper(bray, sampledf$Treatment)

# run a permutation test to get a statistic and a significance score
permutest(beta)

sink()
file.show("Adonis-Betadisper_WholeDataset.txt")
```

## Plot number of head to head interactions of each bee

``` r
library(ggplot2)
library(ggbeeswarm)
library(magrittr)
library(ggpubr)
library(scales)

samdf1 <- subset(samdf, !is.na(HH_norm))
samdf1$Treatment <- factor(samdf1$Treatment, levels= c("MD","CL"))

give.n <- function(x){return(c(y = -0.65, label = length(x))) # experiment with the multiplier to find the perfect position
}

pq9 <- ggplot(samdf1, aes(x = Treatment, y = HH_norm))+ 
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(aes(colour = Treatment),cex=3) +
  scale_x_discrete()+ 
  scale_y_log10()+
  annotation_logticks(base = 10, sides = "l", scaled = TRUE,
                      short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.2, "cm"),
                      colour = "black", size = 0.5, linetype = 1, alpha = 1, color = NULL)+
  facet_grid(.~Replicate)+
  stat_summary(fun.data = give.n, geom = "text") + # add number of observations
  scale_color_manual(values = c("#C77CFF","#00BFC4")) +
  ylab("Normalised head-to-head interactions") +
  theme_bw()
print(pq9)

ggsave(height=4,width=8,dpi=300, filename="HHbyTreatmentbyReplicate-BeeTracking.pdf", useDingbats=FALSE)
```

## Stats

``` r
samdf1$Subcolony <- factor(paste0(samdf1$Treatment, samdf1$Replicate))

library(lmerTest)
resul <- lmer(samdf1$HH_norm ~ samdf1$Treatment  + (1|samdf1$Replicate) + (1|samdf1$Subcolony))
summary(resul)
anova(resul)
```

## To save the session

``` r
save.image("/Volumes/My\ Passport\ for\ Mac/MiSeqRun_Analyses/Bee_tracking/DADA2/BeeTracking_Ampliseq.rdata")
# Close R, Re-open R
#load("path/to/my_session.rdata")
```

## Session info

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
