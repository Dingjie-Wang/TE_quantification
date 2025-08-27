from collections import defaultdict
from intervaltree import Interval, IntervalTree
import numpy as np

from util import sync_reference_name,parse_GTF_feature
# Depreciated: for repeatmasker annoation
# def parse_TE_annotation(annotation_file_path):
#     TE_interval_tree_dict = defaultdict(lambda : defaultdict(IntervalTree))
#     TE_length_dict = defaultdict(dict)
#     TE_features_dict = defaultdict(lambda : defaultdict(dict))
#     TE_names_dict = defaultdict(lambda : defaultdict(lambda : {'occurance':0}))
#     with open(annotation_file_path) as f:
#         lines = f.readlines()[3:]
#         for line in lines:
#             fields = [x for x in line.split(' ') if x !='']
#             chr_name, begin_pos, end_pos, strand, repeat_name, family_name = sync_reference_name(fields[4]), int(fields[5]), int(fields[6]), fields[8], fields[9],fields[10]
#             TE_names_dict[family_name][repeat_name]['occurance'] += 1
#             repeat_name_new = repeat_name + '.' + str(TE_names_dict[family_name][repeat_name]['occurance'])
#             TE_features_dict[family_name][repeat_name_new] = {'chr_name':chr_name,'strand':strand,'begin_pos':begin_pos,'end_pos':end_pos}
#             features = {'repeat_name':repeat_name_new,'family_name':family_name}
#             #Noted here the interval tree library I use does not include the upper bound of interval
#             #So an additional base needed after the end_pos
#             TE_interval_tree_dict[chr_name][strand].addi(begin_pos,end_pos + 1,features)
#     TE_length_list = []
#     num_TEs = 0
#     for family_name in TE_features_dict:
#         for repeat_name in TE_features_dict[family_name]:
#             TE_features_dict[family_name][repeat_name]['index'] = num_TEs
#             TE_length_list.append(TE_features_dict[family_name][repeat_name]['end_pos'] - TE_features_dict[family_name][repeat_name]['begin_pos'])
#             num_TEs += 1
#     TE_length_matrix = np.array(TE_length_list)
#     return TE_interval_tree_dict,TE_features_dict,TE_length_matrix


def parse_TE_annotation(annotation_file_path):
    TE_interval_tree_dict = defaultdict(lambda : defaultdict(IntervalTree))
    TE_length_dict = defaultdict(dict)
    TE_features_dict = defaultdict(lambda : defaultdict(lambda : defaultdict(dict)))
    with open(annotation_file_path) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            ## noted gtf uses 1-base index.Transform to 0-base here
            chr_name, start, end, strand = sync_reference_name(fields[0]), int(fields[3]) - 1, int(fields[4]), fields[6]
            feature_dict = parse_GTF_feature(fields[8])
            repeat_name, subfamily_name, family_name, class_name,consensus_length,truncated_proportion = feature_dict['transcript_id'],feature_dict['gene_id'],feature_dict['family_id'],feature_dict['class_id'],int(feature_dict['consensus_length']),float(feature_dict['truncated_proportion'])
            TE_features_dict[class_name][family_name][subfamily_name][repeat_name] = {'chr_name':chr_name,'strand':strand,'start':start,'end':end,'consensus_length':consensus_length,'truncated_proportion':truncated_proportion}
            feature_index = {'repeat_name':repeat_name,'subfamily_name':subfamily_name,'family_name':family_name,'class_name':class_name,'consensus_length':consensus_length,'truncated_proportion':truncated_proportion}
            #Noted here the interval tree library I use does not include the upper bound of interval
            #So an additional base needed after the end
            TE_interval_tree_dict[chr_name][strand].addi(start,end + 1,feature_index)
    TE_length_list = []
    num_TEs = 0
    for class_name in TE_features_dict:
        for family_name in TE_features_dict[class_name]:
            for subfamily_name in TE_features_dict[class_name][family_name]:
                for repeat_name in TE_features_dict[class_name][family_name][subfamily_name]:
                    TE_features_dict[class_name][family_name][subfamily_name][repeat_name]['index'] = num_TEs
                    TE_length_list.append(TE_features_dict[class_name][family_name][subfamily_name][repeat_name]['end'] - TE_features_dict[class_name][family_name][subfamily_name][repeat_name]['start'])
                    num_TEs += 1
    TE_length_matrix = np.array(TE_length_list)
    return TE_interval_tree_dict,TE_features_dict,TE_length_matrix