# TE_quantification
## Installation
1. Make sure Python 3.6 or above and conda installed on your system
2. `conda env create -f environment.yml -n TE_quantification`
3. `source activate TE_quantification & cd TE_quantification`
## Quantify
```
usage: TE_quantification.py quantif [-h] -annot ANNOTATION -aln ALIGNMENT -o
                                    OUTPUT [--num_iterations NUM_ITERATIONS]
                                    [--only_uniq_mapping ONLY_UNIQ_MAPPING]
                                    [--infer_init_weights INFER_INIT_WEIGHTS]
                                    [--alignment_multiple_counting ALIGNMENT_MULTIPLE_COUNTING]
                                    [--mapping_threshold MAPPING_THRESHOLD]
                                    [--count_mapping COUNT_MAPPING]

optional arguments:
  -h, --help            show this help message and exit

required named arguments:
  -annot ANNOTATION, --annotation ANNOTATION
                        The path of TE annotation file [GTF]
  -aln ALIGNMENT, --alignment ALIGNMENT
                        The path of alignment file [BAM\SAM]
  -o OUTPUT, --output OUTPUT
                        The path of output directory

optional arguments:
  --num_iterations NUM_ITERATIONS
                        Number of iterations for EM algorithm [default:100]
  --only_uniq_mapping ONLY_UNIQ_MAPPING
                        Only count unique mapping reads [default:False]
  --infer_init_weights INFER_INIT_WEIGHTS
                        Infer initital weights for EM algorithm from unique
                        mapping [default:True]
  --alignment_multiple_counting ALIGNMENT_MULTIPLE_COUNTING
                        Count TE alignment multiple times [default:True]
  --mapping_threshold MAPPING_THRESHOLD
                        Mapping proportion threshold of TE/read [default:1.0]
```
## Create curated TE annotation
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
