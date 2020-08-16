# NLR_myrtacae

###NLR_filter_gene_models.py

Prerequiste:

* Bedtools and pybedtools
* pandas
* Biopython
* getAnnoFasta.pl from Augustus in the PATH
* python 3.5+
* assume to be run in the folder that also contains the original PREFIX.fa genome file used for NLR annotator

usage: NLR_filter_gene_models.py [-h]
                                 anno_fn [anno_fn ...] pfam_fn [pfam_fn ...]
                                 fa_fn [fa_fn ...] pro_fn [pro_fn ...]

This program parses out the gene models on the identified NLR loci. It filters
them for NB-ARC domain containing proteins only. Works with augustus.gff3 or
braker.gtf. It needs 'gene' in the feature column. It needs 'transcript_id'
and 'gene_id' in the featur column.

positional arguments:
  anno_fn     Please provide a gtf or gff3 annotation file.
  pfam_fn     Please provide pfam output file for the corresponding proteins.
  fa_fn       Please provide NLR loci filename used for gene prediction.
  pro_fn      Please provide protein filename of the annotated gene models.

optional arguments:
  -h, --help  show this help message and exit
