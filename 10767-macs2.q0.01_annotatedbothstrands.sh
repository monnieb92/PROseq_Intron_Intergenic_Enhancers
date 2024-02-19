#!/bin/bash

ALIGNMENTS=/Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-ALIGN

BAM=/Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-MB-BAM

MACS2=/Users/Hiebertlab/anaconda3/bin/macs2 

OUTPUT=/Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-macs2.q0.01

ANNO=/Volumes/Mbomber_7_14TB/Genomes/mm39_gtfs/refGene.gtf

OUTPUT2=/Volumes/Mbomber_7_14TB/10767-MB-MEL_Mtg16_PROSeq/10767-macs2.q0.01/10767-PS-annotations

FASTA=/Volumes/Mbomber_7_14TB/Genomes/mm39_fasta/mm39.fa

# define arrays
ids=("0001" "0002" "0003" "0004" "0005" "0006" "0007" "0008" "0009" "0010" "0011" "0012")
times=("0hA" "0hB" "0hC" "15mA" "15mB" "15mC" "30mA" "30mB" "30mC" "1hA" "1hB" "1hC")

# iterate over arrays
for i in "${!ids[@]}"; do
	echo "Calling peaks for ${times[$i]}"
	$MACS2 callpeak -t ${BAM}/10767-MB-${ids[$i]}.75.adap.trail.15bp.revcomp.mm39.F4q10.sorted.bam -q 0.01 -f BAM -g mm --call-summits -n $OUTPUT/10767-MB-${times[$i]}-PS.q0.01.fBAM.gmm.callsummits
	
	echo "Annotating Peaks for ${times[$i]}"
	
	annotatePeaks.pl $OUTPUT/10767-MB-${times[$i]}-PS.q0.01.fBAM.gmm.callsummits_peaks.narrowPeak $FASTA -size 200 -strand both -gtf $ANNO > $OUTPUT2/10767-MB-${times[$i]}-PS.q0.01.fBAM.gmm.callsummits_peaks.annotated.bothstrand.txt
	
done 