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