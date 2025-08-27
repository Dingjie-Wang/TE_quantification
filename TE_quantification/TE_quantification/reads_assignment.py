import numpy as np
from collections import defaultdict

def norm(x):
    if sum(x) == 0:
        return 0
    if (not isinstance(x,np.ndarray)):
        x = np.array(x)
    return x/sum(x)

def assign_multi_mapping_reads(all_reads_TE_mapping_dict,TE_features_dict,TE_length_matrix,TE_unique_reads_abundance_matrix,num_iterations,infer_init_weight_from_unique_reads,alignment_multiple_counting):
    num_TEs = TE_length_matrix.shape[0]
    TE_multi_reads_abundance_matrix = np.zeros(num_TEs)
    #initial estimation
    multi_reads_TE_mapping_fraction_dict = {}
    multi_reads_TE_mapping_index_dict = {}
    for read_name in all_reads_TE_mapping_dict:
        num_alignments = len(all_reads_TE_mapping_dict[read_name])
        init_weight = 1/num_alignments
        if (num_alignments > 1):
            multi_reads_TE_mapping_index_dict[read_name] = []
            for TE_mappings in all_reads_TE_mapping_dict[read_name]:
                mapped_TE_indics = []
                for TE_mapping in TE_mappings: 
                    TE_index = TE_features_dict[TE_mapping.data['class_name']][TE_mapping.data['family_name']][TE_mapping.data['subfamily_name']][TE_mapping.data['repeat_name']]['index']
                    if (alignment_multiple_counting):
                        TE_multi_reads_abundance_matrix[TE_index] += init_weight
                    else:
                        TE_multi_reads_abundance_matrix[TE_index] += init_weight/len(TE_mappings)
                    mapped_TE_indics.append(TE_index)
                multi_reads_TE_mapping_index_dict[read_name].append(mapped_TE_indics)
            # init weight
            if (infer_init_weight_from_unique_reads):
                multi_reads_TE_mapping_fraction_dict[read_name] = np.zeros((num_alignments,))
                unique_reads_count_list_per_alignment = []
                for mapped_TE_indics in multi_reads_TE_mapping_index_dict[read_name]:
                    unique_reads_count_list_per_alignment.append(TE_unique_reads_abundance_matrix[mapped_TE_indics[0]])
                num_zero_count_alignments = unique_reads_count_list_per_alignment.count(0)
                num_non_zero_count_alignments = len(unique_reads_count_list_per_alignment) - num_zero_count_alignments
                sum_non_zero_count_alignments = sum(unique_reads_count_list_per_alignment)
                for i in range(len(unique_reads_count_list_per_alignment)):
                    count = unique_reads_count_list_per_alignment[i]
                    if count == 0:
                        multi_reads_TE_mapping_fraction_dict[read_name][i] = 1 / num_zero_count_alignments
                    else:
                        multi_reads_TE_mapping_fraction_dict[read_name][i] = (1 / num_non_zero_count_alignments) * (count / sum_non_zero_count_alignments)
            else:
                multi_reads_TE_mapping_fraction_dict[read_name] = np.full((num_alignments,),init_weight)
#     TE_relative_abundance = norm(TE_multi_reads_abundance_matrix / TE_length_matrix)
    TE_relative_abundance = norm(TE_multi_reads_abundance_matrix)
    print('Multi_mapping_reads: %d'%(len(multi_reads_TE_mapping_fraction_dict)))
    if (len(multi_reads_TE_mapping_fraction_dict)==0):
        return np.zeros(num_TEs)
    #iterative estimation
    print('Start assigning multi_mapping_reads by EM algorithm')
    for i in range(num_iterations):
        for read_name in multi_reads_TE_mapping_fraction_dict:
            mapped_TE_indics = multi_reads_TE_mapping_index_dict[read_name]
            if (len(mapped_TE_indics)==0):
                continue
            # Eq 2:calculate new fractions
            multi_reads_TE_mapping_fraction_dict[read_name] = norm([sum(TE_relative_abundance[indics]) for indics in mapped_TE_indics])
            
        TE_relative_abundance = np.zeros(TE_relative_abundance.shape)
        # Eq 3:calculate relative abundance
        for read_name in multi_reads_TE_mapping_fraction_dict:
            mapped_TE_indics = multi_reads_TE_mapping_index_dict[read_name]
            for j in range(len(mapped_TE_indics)):
#                 TE_relative_abundance[mapped_TE_indics[j]] += multi_reads_TE_mapping_fraction_dict[read_name][j] / TE_length_matrix[mapped_TE_indics[j]]
                TE_relative_abundance[mapped_TE_indics[j]] += multi_reads_TE_mapping_fraction_dict[read_name][j]
        
        TE_relative_abundance = norm(TE_relative_abundance)
    TE_read_counts = (TE_relative_abundance * len(multi_reads_TE_mapping_fraction_dict) + 0.5).astype(int)

    return TE_read_counts

def assign_unique_mapping_reads(all_reads_TE_mapping_dict,TE_features_dict,TE_length_matrix,alignment_multiple_counting):
    num_TEs = TE_length_matrix.shape[0]
    TE_unique_reads_abundance_matrix = np.zeros(num_TEs)
    num_unique_reads = 0
    for read_name in all_reads_TE_mapping_dict:
        if (len(all_reads_TE_mapping_dict[read_name]) == 1):
            num_unique_reads += 1
            num_mapped_TE = len(all_reads_TE_mapping_dict[read_name][0])
            for TE_mapping in all_reads_TE_mapping_dict[read_name][0]:
                if (alignment_multiple_counting):
                    TE_unique_reads_abundance_matrix[TE_features_dict[TE_mapping.data['class_name']][TE_mapping.data['family_name']][TE_mapping.data['subfamily_name']][TE_mapping.data['repeat_name']]['index']] += 1
                else:
                    TE_unique_reads_abundance_matrix[TE_features_dict[TE_mapping.data['class_name']][TE_mapping.data['family_name']][TE_mapping.data['subfamily_name']][TE_mapping.data['repeat_name']]['index']] += 1/num_mapped_TE
#     TE_relative_abundance = TE_unique_reads_abundance_matrix / TE_length_matrix
    print('Unique_mapping_reads: %d'%(num_unique_reads))
    TE_unique_reads_abundance_matrix = (TE_unique_reads_abundance_matrix+0.5).astype(int)
    return TE_unique_reads_abundance_matrix