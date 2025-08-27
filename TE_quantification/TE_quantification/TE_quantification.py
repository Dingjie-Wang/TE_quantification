import argparse
import time

from parse_annotation import parse_TE_annotation
from map_reads import map_reads
from generate_output import generate_output
from util import count_unique_multi_reads_proportion
from create_curated_annotation import create_curated_annotation

def TE_quantification(annotation_file_path,alignment_file_path,output_path,num_EM_iterations=100,infer_init_weight_from_unique_reads=True,mapping_threshold=0.5,alignment_multiple_counting=True,only_uniq_mapping=False):
    start_time = time.time()
    print('Parsing TE annotation...')
    TE_interval_tree_dict,TE_features_dict,TE_length_matrix = parse_TE_annotation(annotation_file_path)
    print('Mapping reads to TE annotation...')
    TE_features_dict,TE_TPM,TE_CPM,TE_abundance_matrix = map_reads(alignment_file_path,TE_interval_tree_dict,TE_features_dict,TE_length_matrix,num_EM_iterations,infer_init_weight_from_unique_reads,mapping_threshold,alignment_multiple_counting,only_uniq_mapping)
    generate_output(TE_features_dict,TE_TPM,TE_CPM,TE_abundance_matrix,output_path)
    print('Done in %f seconds'%(time.time() - start_time))
def parse_arguments():
    """
    Parse the arguments
    """
    parser = argparse.ArgumentParser(description="TE quantification tools",add_help=True)
    subparsers = parser.add_subparsers(help='sub-command help',dest="subparser_name")
    parser_quantif = subparsers.add_parser('quantif', help='quantify TE expression')
    parser_select_annot = subparsers.add_parser('annot_select', help='Create curated TE annotation')
    
    requiredNamed = parser_quantif.add_argument_group('required named arguments')
    requiredNamed.add_argument('-annot','--annotation', type=str, help="The path of TE annotation file [GTF]",required=True)
    requiredNamed.add_argument('-aln','--alignment', type=str, help="The path of alignment file [BAM\SAM]",required=True)
    requiredNamed.add_argument('-o','--output', type=str, help="The path of output directory",required=True)
    optional = parser_quantif.add_argument_group('optional arguments')
    optional.add_argument('--num_iterations',type=int,default=100, help="Number of iterations for EM algorithm [default:100]")
    optional.add_argument('--only_uniq_mapping',type=bool,default=False, help="Only count unique mapping reads [default:False]")
    optional.add_argument('--infer_init_weights',type=bool,default=True, help="Infer initital weights for EM algorithm from unique mapping [default:True]")
    optional.add_argument('--alignment_multiple_counting',type=bool,default=True, help="Count TE alignment multiple times [default:True]")
    optional.add_argument('--mapping_threshold',type=float,default=1.0, help="Mapping proportion threshold of TE/read [default:1.0]")
    # optional.add_argument('--count_mapping',type=bool,default=False,help="Count mapping reads proportion")
    
    
    requiredNamed_select_annot = parser_select_annot.add_argument_group('required named arguments for TE annotation selection')
    requiredNamed_select_annot.add_argument('-rptmsk','--repeat_masker_annotation', type=str, help="The path of RepeatMasker annotation file [.fa.out]",required=True)
    requiredNamed_select_annot.add_argument('-trans_annot','--transcript_annotation', type=str, help="The path of transcript annotation file [GTF]",required=True)
    requiredNamed_select_annot.add_argument('-o','--output', type=str, help="The path of output curated annotation [GTF]",required=True)
    requiredNamed_select_annot.add_argument('-isoform_o','--isoform_output', type=str, help="The path of output isoform TE overlapping length [TSV]",required=True)
    args = parser.parse_args()
    if args.subparser_name == 'annot_select':
        create_curated_annotation(args.transcript_annotation,args.repeat_masker_annotation,args.output,args.isoform_output)
    elif args.subparser_name == 'quantif':
        # if (args.count_mapping):
        #     count_unique_multi_reads_proportion(args.annotation,args.alignment,args.output)
        # else:
        TE_quantification(args.annotation,args.alignment,args.output,args.num_iterations,args.infer_init_weights,args.mapping_threshold,args.alignment_multiple_counting,args.only_uniq_mapping)
    else:
        print('Unknown command!') 
if __name__ == "__main__":
    parse_arguments()
