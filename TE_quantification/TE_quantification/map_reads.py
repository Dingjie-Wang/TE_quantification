import pysam
import numpy as np
from collections import defaultdict
from reads_assignment import assign_unique_mapping_reads,assign_multi_mapping_reads
from util import sync_reference_name

def normalize_read_counts(TE_abundance_matrix,TE_features_dict,total_num_TE_reads,total_num_genome_reads):
    TE_length_matrix = np.zeros(TE_abundance_matrix.shape) 
    for class_name in TE_features_dict:
        for family_name in TE_features_dict[class_name]:
            for subfamily_name in TE_features_dict[class_name][family_name]:
                for repeat_name in TE_features_dict[class_name][family_name][subfamily_name]:
                    TE_feature = TE_features_dict[class_name][family_name][subfamily_name][repeat_name]
                    TE_length_matrix[TE_feature['index']] = TE_feature['end'] - TE_feature['start']
#     TE_TPM = TE_FPKM = TE_abundance_matrix / TE_length_matrix
    TE_TPM = TE_CPM = TE_abundance_matrix
    try:
        TE_TPM = TE_TPM * 1e6 / sum(TE_TPM)
        TE_CPM = TE_CPM * 1e6 / total_num_genome_reads
    except:
        print(sum(TE_TPM))
        print(total_num_TE_reads)
    return TE_TPM,TE_CPM
def get_mapping_dict(alignment_file_path,TE_interval_tree_dict,mapping_strategy="TE",mapping_strategy_proportion = 0.5):
    all_reads_TE_mapping_dict = defaultdict(list)
    if alignment_file_path.endswith('.bam'):
        samfile = pysam.AlignmentFile(alignment_file_path, "rb")
    elif alignment_file_path.endswith('.sam'):
        samfile = pysam.AlignmentFile(alignment_file_path, "r")
    else:
        raise Exception('Invalid alignment file given!')
    read_names = set()
    try:
        for alignment in samfile:
            if (not alignment.is_unmapped):
                chr_name,start_pos,read_name = sync_reference_name(alignment.reference_name),alignment.reference_start,alignment.query_name
                read_names.add(read_name)
                if alignment.query_alignment_length / alignment.infer_read_length() < 0.8:
                    continue
                strand = '-' if alignment.is_reverse else '+'
                if ((chr_name not in TE_interval_tree_dict) or (strand not in TE_interval_tree_dict[chr_name])):
                    continue
                TE_index = TE_interval_tree_dict[chr_name][strand]
                blocks = alignment.get_blocks()
                TE_mappings = set()
                for block in blocks:
                    read_begin,read_end = block[0],block[1]
                    overlapped_TE = TE_index.overlap(read_begin,read_end)
                    filtered_TE = set()
                    for TE in overlapped_TE:
                        # Noted here the interval tree library I use does not include the upper bound of interval
                        TE_begin,TE_end = TE.begin,TE.end - 1
                        overlapped_length = (TE_end - TE_begin) + (read_end - read_begin) - (max(TE_end,read_end) - min(TE_begin,read_begin))
                        if ((TE_end - TE_begin) > (read_end - read_begin)):
                            mapping_strategy = "read"
                        else:
                            mapping_strategy = "TE"
                        if ((mapping_strategy == "TE" and overlapped_length >= mapping_strategy_proportion * (TE_end - TE_begin)) or (mapping_strategy == "read" and overlapped_length >= mapping_strategy_proportion * (read_end - read_begin))):
                            filtered_TE.add(TE)
                    TE_mappings = TE_mappings.union(filtered_TE)
                if (len(TE_mappings)>0):
                    all_reads_TE_mapping_dict[read_name].append(TE_mappings)
    except:
        pass
    samfile.close()         
    return all_reads_TE_mapping_dict,read_names
    
def map_reads(alignment_file_path,TE_interval_tree_dict,TE_features_dict,TE_length_matrix,num_EM_iterations,infer_init_weight_from_unique_reads,mapping_threshold,alignment_multiple_counting,only_uniq_mapping):
    all_reads_TE_mapping_dict,read_names = get_mapping_dict(alignment_file_path,TE_interval_tree_dict,"TE",mapping_threshold)
    total_num_TE_reads = len(all_reads_TE_mapping_dict)
    total_num_genome_reads = len(read_names)
    TE_unique_reads_abundance_matrix = assign_unique_mapping_reads(all_reads_TE_mapping_dict,TE_features_dict,TE_length_matrix,alignment_multiple_counting)
    if (only_uniq_mapping):
        TE_abundance_matrix = TE_unique_reads_abundance_matrix.copy()
    else:
        TE_multi_reads_abundance_matrix = assign_multi_mapping_reads(all_reads_TE_mapping_dict,TE_features_dict,TE_length_matrix,TE_unique_reads_abundance_matrix,num_EM_iterations,infer_init_weight_from_unique_reads,alignment_multiple_counting)
        TE_abundance_matrix = (TE_unique_reads_abundance_matrix + TE_multi_reads_abundance_matrix+0.5).astype(int)
    TE_TPM,TE_CPM = normalize_read_counts(TE_abundance_matrix,TE_features_dict,total_num_TE_reads,total_num_genome_reads)
    return TE_features_dict,TE_TPM,TE_CPM,TE_abundance_matrix