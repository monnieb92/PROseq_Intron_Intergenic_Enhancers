# eRNA peak calling 
This is performed for each timepoint and replicate. 
## Call peaks using macs2 

```{bash}
$MACS2 callpeak -t ${BAM}/10767-MB-${ids[$i]}.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam -q 0.01 -f BAM -g mm --call-summits -n $OUTPUT/10767-MB-${times[$i]}-PS.q0.01.fBAM.gmm.callsummits
```
## Annotate Peaks with HOMER

```{bash}
annotatePeaks.pl $OUTPUT/10767-MB-${times[$i]}-PS.q0.01.fBAM.gmm.callsummits_peaks.narrowPeak $FASTA -size 200 -strand both -gtf $ANNO > $OUTPUT2/10767-MB-${times[$i]}-PS.q0.01.fBAM.gmm.callsummits_peaks.annotated.bothstrand.txt
```

## Filter for only peaks in introns and intergenic regions and not within 750bp of TSS using R
### added 150bp to each peak (this will merge peaks that might very close together)
```{r}
T0hA <- read_delim("/Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-macs2.q0.01/10767-PS-annotations/10767-MB-0hA-PS.q0.01.fBAM.gmm.callsummits_peaks.annotated.bothstrand.txt")
names(T0hA)[1] <- "PeakID" 
T0hA_macs <- read_delim("/Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-macs2.q0.01/10767-MB-0hA-PS.q0.01.fBAM.gmm.callsummits_peaks.narrowPeak", col_names = FALSE)
Macs2colnames <- c("Chrom","Start", "End", "PeakID", "score", "Strand", "signalValue", "pValue", "qValue", "peak")
colnames(T0hA_macs) <- Macs2colnames

T0hA_int <- subset(T0hA, grepl("(intron|intergenic)", tolower(Annotation))) %>%
  dplyr::filter(`Distance to TSS` > 750 | `Distance to TSS` < -750) %>%
  left_join(T0hA_macs, by = "PeakID", suffix = c("_anno", "_macs2")) %>% 
  dplyr::select(Chrom, Start_macs2, End_macs2,PeakID, "Peak Score", Strand_macs2) %>%
  dplyr::mutate(Start_macs2 = Start_macs2 - 150, 
                End_macs2 = End_macs2 + 150) %>%
  as.data.frame()
write_tsv(T0hA_int,"10767-MB-0hA-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.outsideTSS750.plusminus150.bothstrand.txt", col_names = FALSE)
```
## Concatenate peaks from all time points
```{bash}
##concatentate
cat /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-R_PROseq/10767-MB-PROseq-Routput/10767-MB-0hA-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.outsideTSS750.plusminus150.bothstrand.txt /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-R_PROseq/10767-MB-PROseq-Routput/10767-MB-0hB-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.outsideTSS750.plusminus150.bothstrand.txt /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-R_PROseq/10767-MB-PROseq-Routput/10767-MB-0hC-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.outsideTSS750.plusminus150.bothstrand.txt /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-R_PROseq/10767-MB-PROseq-Routput/10767-MB-15mA-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.outsideTSS750.plusminus150.bothstrand.txt /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-R_PROseq/10767-MB-PROseq-Routput/10767-MB-15mB-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.outsideTSS750.plusminus150.bothstrand.txt /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-R_PROseq/10767-MB-PROseq-Routput/10767-MB-15mC-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.outsideTSS750.plusminus150.bothstrand.txt /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-R_PROseq/10767-MB-PROseq-Routput/10767-MB-30mA-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.outsideTSS750.plusminus150.bothstrand.txt /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-R_PROseq/10767-MB-PROseq-Routput/10767-MB-30mB-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.outsideTSS750.plusminus150.bothstrand.txt /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-R_PROseq/10767-MB-PROseq-Routput/10767-MB-30mC-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.outsideTSS750.plusminus150.bothstrand.txt /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-R_PROseq/10767-MB-PROseq-Routput/10767-MB-1hA-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.outsideTSS750.plusminus150.bothstrand.txt /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-R_PROseq/10767-MB-PROseq-Routput/10767-MB-1hB-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.outsideTSS750.plusminus150.bothstrand.txt /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-R_PROseq/10767-MB-PROseq-Routput/10767-MB-1hC-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.outsideTSS750.plusminus150.bothstrand.txt > 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150.bed.txt
#sort bed files
sort -k1,1 -k2,2n 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150.bed.txt > 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150.sorted.bed.txt
#mergepeaks
bedtools merge -i 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150.sorted.bed.txt > 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150.sorted.merged.bed
#saf file for featurecounts
awk -v OFS="\t" 'BEGIN {print "GeneID","Chr","Start","End","Strand"} { print "Peak_"NR,$1,$2,$3,"."}' 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150.sorted.merged.bed > 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150.sorted.merged.saf
#replace positive strand
sed 's/\./+/g' 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150.sorted.merged.saf >10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150.sorted.merged.pos.saf
```

## Feature counts of the possible intron and intergenic regions 
```{bash}
#!/bin/bash 

featureCounts -F SAF -T 18 -s 2 -O -a 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150.sorted.merged.pos.saf  -o 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150_SubreadCounts.reversestrand.possaf.txt \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0001.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0002.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0003.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0004.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0005.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0006.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0007.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0008.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0009.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0010.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0011.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0012.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam

featureCounts -F SAF -T 18 -s 1 -O -a 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150.sorted.merged.pos.saf  -o 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150_SubreadCounts.forwardstrand.possaf.txt \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0001.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0002.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0003.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0004.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0005.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0006.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0007.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0008.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0009.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0010.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0011.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0012.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam

featureCounts -F SAF -T 18 -s 0 -O -a 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150.sorted.merged.pos.saf -o 10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150_SubreadCounts.allcounts.possaf.txt \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0001.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0002.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0003.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0004.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0005.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0006.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0007.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0008.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0009.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0010.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0011.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam \
  /Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM/10767-MB-0012.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam
```

## Filter for peaks with bidirectional transcription 
### Counts tables filtered > 0.001 & _rev !=0
```{r}

counts_strand <- read.table("/Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-FeatureCounts/10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150_SubreadCounts.forwardstrand.possaf.txt", header = TRUE)
CountColumns <- c("PeakID","Chrom","Start","End","Strand","Length","minusdT_A","minusdT_B","minusdT_C","plusdT15m_A","plusdT15m_B","plusdT15m_C","plusdT30m_A","plusdT30m_B","plusdT30m_C","plusdT1h_A","plusdT1h_B","plusdT1h_C")
colnames(counts_strand) <- CountColumns
counts_reverse <- read.table("/Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-FeatureCounts/10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150_SubreadCounts.reversestrand.possaf.txt", header = TRUE)
colnames(counts_reverse) <- CountColumns
counts_total <- read.table("/Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-FeatureCounts/10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150_SubreadCounts.allcounts.possaf.txt", header = TRUE)
colnames(counts_total) <- CountColumns
AllCounts <- counts_strand %>% left_join(counts_reverse, by = "PeakID", suffix=c("_fwd", "_rev")) %>%
  left_join(counts_total, by = "PeakID", suffix=c("_rev", "_both"))

AllCounts_filtered <- AllCounts %>%
  dplyr::filter((`minusdT_A_rev` != 0 ) & (`minusdT_A_fwd`/`minusdT_A_rev`) > 0.001 & (`minusdT_B_rev` != 0 ) & (`minusdT_B_fwd`/`minusdT_B_rev`) > 0.001 & (`minusdT_C_rev` != 0 ) & (`minusdT_C_fwd`/`minusdT_C_rev`) > 0.001 |  (`plusdT1h_A_rev` != 0 ) & (`plusdT1h_A_fwd`/`plusdT1h_A_rev`) > 0.001 & (`plusdT1h_B_rev` != 0 ) & (`plusdT1h_B_fwd`/`plusdT1h_B_rev`) > 0.001  & (`plusdT1h_C_rev` != 0 ) &  (`plusdT1h_C_fwd`/`plusdT1h_C_rev`) > 0.001 ) %>%
  glimpse()

AllCounts_filtered_bed <- AllCounts_filtered %>%
  dplyr::select("Chrom", "Start", "End", "PeakID","Length","Strand")

write_tsv(AllCounts_filtered_bed, "10767-MB-AllPeaks-PS.q0.01.fBAM.gmm.callsummits_peaks.intron.intergenic.filtered.outsideTSS750.plusminus150.greaterthan0.001.bed", col_names = FALSE)
```

## Differential analysis using NRSA normalization factors and DESeq2 
### use the counts table from both strands to perform the differential anslysis (all counts)

```{r}
countData1h <- AllCounts_filtered[,c(41:43,50:52)]

colData_1h <- data.frame(condition = factor(c("DMSO", "DMSO","DMSO", "dT1h","dT1h","dT1h"), levels = c("DMSO", "dT1h")), sample = c("minusdT_A","minusdT_B","minusdT_C", "dT1h_A","dT1h_B","dT1h_C"))

dds1h <- DESeqDataSetFromMatrix(countData = countData1h, colData = colData_1h, design = ~condition)
## these size factors are output from NRSA
sizeFactors(dds1h) <- c(1/0.721057608621947,1/0.994705249048991,1/1.288104653729,1/1.05075548011527,1/0.910653829570098,1/1.12421859250555)
#dds <- DESeq(dds)
dds1h <- estimateDispersions(dds1h)
dds1h <- nbinomWaldTest(dds1h)

##Create PCA plot of normalized data

vsd<-varianceStabilizingTransformation(dds1h)
z <- plotPCA(vsd, intgroup=c("condition"))
nudge <- position_nudge(y = 1)
z + geom_text(aes(label = vsd$sample), position = nudge, size=2)
 
normdata <- counts(dds1h,normalized=TRUE)
colnames(normdata)<-paste("norm",c("minusdT_A","minusdT_B","minusdT_C", "dT1h_A","dT1h_B","dT1h_C"),sep="_")
 
result_eRNA_1h <- results(dds1h,contrast=c("condition","dT1h","DMSO"))
result_eRNA_1h_table <- cbind(AllCounts_filtered[,c(1:5)],countData1h,normdata,result_eRNA_1h)
glimpse(result_eRNA_1h_table)
write_tsv(result_eRNA_1h_table,"10767-MB-result_eRNA_1h_table.includeintron.intergentic.excludenearTSS750.plusminus150.filtered0.001.txt")
#create a MA plot of differential changes 
plotMA(result_eRNA_1h, ylim=c(-7.5,5),main="eRNA",cex.lab=1.75, cex.axis=1.75 ,cex.main=1.75, alpha=0.05)

```
