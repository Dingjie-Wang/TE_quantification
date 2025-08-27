# TE_quantification
TE_quantification is a package that compares TE copy annotations from RepeatMasker with the transcript annotation file of a given sample based on genomic coordinates. Transcripts were classified into three major categories:
1. TE-alone: autonomously expressed TE transcripts, defined as transcripts overlapping TE regions by >80% with the first exon also covered by a TE (i.e., TE promoter-driven transcription).
2. TE-Gene: TE-Gene chimeric transcripts with partial overlap between TE annotations and gene regions.
3. Gene-alone: transcripts without any overlap with TE-related sequences.
TE-alone transcripts were further classified into subfamilies, families, and types (DNA transposon, LTR, LINE, SINE, and SVA) following the human TE classification system.
## Installation
1. Make sure Python 3.6 or above and conda installed on your system
2. `conda env create -f environment.yml -n TE_quantification`
3. `source activate TE_quantification & cd TE_quantification`

## Classification of TE-alone, TE-Gene and Gene-alone transcripts
```
usage: TE_quantification.py annot_select [-h] -rptmsk REPEAT_MASKER_ANNOTATION
                                         -trans_annot TRANSCRIPT_ANNOTATION -o
                                         OUTPUT -isoform_o ISOFORM_OUTPUT

optional arguments:
  -h, --help            show this help message and exit

required named arguments for TE annotation selection:
  -rptmsk REPEAT_MASKER_ANNOTATION, --repeat_masker_annotation REPEAT_MASKER_ANNOTATION
                        The path of RepeatMasker annotation file [.fa.out]
  -trans_annot TRANSCRIPT_ANNOTATION, --transcript_annotation TRANSCRIPT_ANNOTATION
                        The path of transcript annotation file [GTF]
  -o OUTPUT, --output OUTPUT
                        The path of output curated annotation [GTF]
  -isoform_o ISOFORM_OUTPUT, --isoform_output ISOFORM_OUTPUT
                        The path of output isoform TE overlapping length [TSV]
```
