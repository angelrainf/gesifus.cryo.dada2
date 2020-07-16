#Modified by @Angel,@GUI following https://benjjneb.github.io/dada2/tutorial.html: 
library(dada2)
library(phyloseq)
library(DECIPHER)
library(phangorn)

setwd("/Users/macbookpro/Dropbox/PhD Angel Rain/paper_cryopreservation/PrimerClipped") #Local path
path<-("/Users/macbookpro/Dropbox/PhD Angel Rain/paper_cryopreservation/PrimerClipped") #Local path

list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Quality assesment
#plotQualityProfile(fnFs)
#plotQualityProfile(fnRs)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#all other arguments are on their default levels
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,190),orient.fwd ="TACG",
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
#reads.in reads.out
#515FY-926Rjed-CANET_R1.fastq    29409     23290
#515FY-926Rjed-CM1_R1.fastq      30842     22712
#515FY-926Rjed-CM2_R1.fastq      53018     40342
#515FY-926Rjed-CM3_R1.fastq      51221     39331
#515FY-926Rjed-CS1_R1.fastq      55352     44292
#515FY-926Rjed-CS2_R1.fastq      42024     31896

#Learn the error rates of each base
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Remove replications
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE,pool=T)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE,pool=T)

dadaFs[[1]]
#dada-class: object describing DADA2 denoising results
#493 sequence variants were inferred from 8188 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
dadaRs[[1]]
#dada-class: object describing DADA2 denoising results
#523 sequence variants were inferred from 6944 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#[1]   25 9420
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#[1]  25 700
sum(seqtab.nochim)/sum(seqtab)
#[1] 0.813922

#Track the number of reads that made it through so far
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
###

#Taxonomy assignation using GTDB (https://gtdb.ecogenomic.org/)
taxa_gtdb <- assignTaxonomy(seqtab.nochim, "/Users/macbookpro/Dropbox/PhD Angel Rain/paper_cryopreservation/GTBD_database/GTDB_bac-arc_ssu_r86.fa",multithread=TRUE) #path local directory
taxa=taxa_gtdb

#Inspect taxonomic assignments
taxa.print <- taxa 
rownames(taxa.print) <- NULL # Removing sequence rownames for display only taxa
head(taxa.print)
taxa.print.names<-as.data.frame(taxa.print)

#Simple data.frame construction from the information encoded in the filenames
#Sample names are the 16S sequences -->ps0
samples.out <- rownames(seqtab.nochim)
community <-c("C", "C", "C", "C", "C", "C", "C","S","S","C","C","SC","SC", "SC", "SC", 
              "S", "S", "S", "S", "S", "S", "S","S", "S", "S") 
DOM<-c("C_initial", "M", "M", "M", "S", "S", "S","S_initial","S_initial","C_initial","C_initial", "M","M", "S", "S", 
       "SW", "SW", "SW", "M", "M", "M", "S_initial","S" ,"S", "S")
replicate<-c("1", "1" ,  "2"  , "3" ,  "1"  , "2" ,  "3" ,"2" ,"3","2","3", "1" , "2",  "2",  "3",
             "1",  "2",  "3",  "1",  "2"  , "3",   "1",  "1",   "2",   "3"  ) 
treatment<-c("C_initial", "CM", "CM", "CM", "CS", "CS", "CS","S_initial","S_initial","C_initial","C_initial" ,"SCM","SCM", "SCS", "SCS", 
             "SW", "SW", "SW", "SM", "SM", "SM", "S_initial","SS", "SS", "SS")

samdf <- data.frame(Community=community, DOM=DOM, Replicate=replicate, Treatment=treatment)
rownames(samdf) <- samples.out

# Phyloseq object from dada2 output
ps0 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(samdf), 
                tax_table(as.matrix(taxa)))

# Construction of a phylogenetic tree
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree.
names(seqs)<- paste0("SV_", seq(ntaxa(ps0)), "_", taxa.print.names$Order)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
seqs.tab<-as.data.frame(seqs)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
###

# New phyloseq object including phylogenetic tree
seqtab.nochim2<-seqtab.nochim
colnames(seqtab.nochim2)<-rownames(seqs.tab)
taxa2<-taxa
rownames(taxa2)<-rownames(seqs.tab)
ps <- phyloseq(otu_table(seqtab.nochim2, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(as.matrix(taxa2)),phy_tree(fitGTR$tree))

# Phyloseq object with relative abundance table
ps.rel = transform_sample_counts(ps, function(x) x/sum(x))
