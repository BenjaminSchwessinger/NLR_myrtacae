#!/bin/env bash

set -xeuo pipefail # this is "bash safe mode", which catches many errors. use it!

#This is NLR annotation analysis test on Arabidopsis genome <- the whole purpose of the script is to work on everything not just one species!!!
#Things to fix
#make it general for species
#make it general for computer set up
#let the user input the path for NLR-praser.jar files including meme file. Put everything in one location.

genome=$1
#I made this a second input file so this becomes more flexible.
blastproteins=$2

#Reuse the genome as a variable

echo ${genome}
echo ${blastproteins}

genomebase=$(basename ${genome} .fa)

echo ${genomebase}

#Chopping the genome sequence into overlapping subsequences
java -jar ~/NLR-parser/scripts/ChopSequence.jar -i ${genome} \
	-o ${genomebase}_choppedseq.fa -l 20000 -p 5000

#Searching the chopped subsequences for pre-determined NLR-associated motifs
java -jar ~/NLR-parser/scripts/NLR-Parser.jar -t 10 \
	-y ~/anaconda3/envs/NLR_Annotator/bin/meme_4.9.1/bin/mast \
	-x ~/NLR-parser/scripts/meme.xml -i ${genomebase}_choppedseq.fa \
	-c ${genomebase}_NLRparser.xml

#Generate the GFF format of NLR loci for the searched motifs
java -jar ~/NLR-parser/scripts/NLR-Annotator.jar -i ${genomebase}_NLRparser.xml \
	-g ${genomebase}_NLRparser.gff

#make a database folder for the blastdb

mkdir -p blastdb

#Generate species genome database to be used for comparison against nucleotide or protein query sequence.
makeblastdb -in ${genome} -dbtype nucl -parse_seqids \
        -out ./blastdb/${genomebase}

#Generate tblastn of protein query sequences of genes which cannot be captured using NLR-parser program such as ADR1 and ZAR1
#TBLASTN helps to compare a protein query against the all six reading frames of a nucleotide sequence database.
#Additional search fasta file on the command line is required. This should be ADR1, RPW8, ZAR1 like proteins & others.
tblastnout=${genomebase}_tblastn${blastproteins}.outfmt6
tblastn -query ${blastproteins} \
	-db ./blastdb/${genomebase} \
	-evalue 0.0001 \
        -outfmt 6 > ${tblastnout}

#Convert tblastn file into bed format. Ideally we move blast2bed.sh into the PATH.
bash blast2bed.sh ${tblastnout}

#Indexing reference sequence
samtools faidx ${genome}

#Cereate genome file. FIXED this here too to make it general.
cut -d $'\t' -f1,2 ${genome}.fai \
        > ${genomebase}.genomefile

#Generate 20kb flanking BED file for blastx file
bedtools slop -b 20000 -s -i ${tblastnout}.bed \
        -g ${genomebase}.genomefile | bedtools sort -i - | bedtools merge \
        -s -d 100 -i - > ${genomebase}_tblastn_20kbflanking.bed

#Generate 20kb flanking BED file for NLR-parser file
bedtools slop -b 20000 -s -i ${genomebase}_NLRparser.gff -g ${genomebase}.genomefile \
        | bedtools sort -i - | bedtools merge -s -d 100 -i - \
        > ${genomebase}_NLRparser_20kbflanking.bed

#Merge the two bed files (blastn.bed and NLRparser.bed into one BED file
cat ${genomebase}_tblastn_20kbflanking.bed ${genomebase}_NLRparser_20kbflanking.bed \
        | bedtools sort -i - \
        | bedtools merge -d 100 -i - \
        > ${genomebase}_all_20kbflanking.bed

#Convert the merged bed file into FASTA format
bedtools getfasta -fi ${genome} -bed ${genomebase}_all_20kbflanking.bed \
        > ${genomebase}_all_20kbflanking.fa

#Convert all the sequences in 20kb flanking fasta into uppercase
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' ${genomebase}_all_20kbflanking.fa \
	 > ${genomebase}_all_20kbflanking_upper.fa

#Gene prediction by BRAKER using extended regions around NB-ARCs by 20kb up and downstream
braker.pl --genome=${genomebase}_all_20kbflanking_upper.fa \
	--prot_seq=RefPlant_235NLR_Representatives.fa \
	--species=${genomebase} --epmode --cores=15 --softmasking --prg=ph \
	--ALIGNMENT_TOOL_PATH=~/anaconda3/envs/braker2/bin/spaln --gff3

#Sort amino acid sequences starting with methionine
cat ./braker/augustus.hints.aa | seqkit grep -s -r -i -p ^M \
        > ${genomebase}_braker_sorted_proteincoding.fa

#Search NB-ARC domain against library of Pfam
pfam_scan.pl -fasta ${genomebase}_braker_sorted_proteincoding.fa \
	-dir ~/ddatabase/PfamScan -as \
	-outfile ${genomebase}_braker_pfamscan_output.txt

#Parsing the PfamScan output to extracts all domains for each proteins and 
#removes redundant nested hits with larger e-values
# Reference: Sarris et al., 2016

perl K-parse_Pfam_domains_v3.1.pl --pfam ${genomebase}_braker_pfamscan_output.txt \
	--evalue 0.001 --output ${genomebase}_pfamscan_output.parser.verbose --verbose T

# Parsing the output of PfamScan output parser using the script "K-parse_Pfam_domains_NLR-fusions-v2.2.pl"
# The script generates 
# Summary of number of NLRs and NLR-IDs
# Summary of integrated domains with species list for each domain
# Abundace list of IDs counted once for each family
# Contingency table per ID domain for each species
# Reference: Sarris et al., 2016

mkdir -p ${genomebase}_pfam-parser

perl K-parse_Pfam_domains_NLR-fusions-v2.2.pl --indir ${genomebase}_pfam-parser \
	--evalue 0.001 -o ${genomebase}_pfam-parser -d db_descriptions.txt
