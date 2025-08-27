import re
import pysam
from collections import defaultdict
## for gtf i.e. gff ver.2 only
def create_GTF_feature(key_val_pairs):
    feature_string = ''
    for key,val in key_val_pairs:
        feature_string += '{} "{}"; '.format(key,val)
    return feature_string
def parse_GTF_feature(feature_string):
    feature_dict = {}
    key_val_pairs = [s for s in feature_string.split(";") if (s!= "" and not s.isspace()) ]
    for key_val_pair_string in key_val_pairs:
        key_val_pair = [s for s in key_val_pair_string.split(" ") if s != ""]
        key,val = key_val_pair[0:2]
        feature_dict[key] = val.replace('"','')
    return feature_dict
def sync_reference_name(ref_name):
    ref_name = ref_name.upper()
    match = re.search("(?<=CHR).*", ref_name)
    if match:
        ref_name = match.group(0)
    if ref_name == "M":
        ref_name = "MT"
    return ref_name
def count_alignments(alignment_file_path,TE_interval_tree_dict,mapping_strategy="TE",mapping_strategy_proportion = 0.5):
    all_reads_TE_mapping_dict = defaultdict(lambda: {'num_alignments':0,'mapped_to_TE':False})
    samfile = pysam.AlignmentFile(alignment_file_path, "r")
    read_names = set()
    unmapped_reads = 0
    try:
        for alignment in samfile:
            if (alignment.is_unmapped):
                unmapped_reads += 1
                continue
            chr_name,start_pos,read_name = sync_reference_name(alignment.reference_name),alignment.reference_start,alignment.query_name
            read_names.add(read_name)
            strand = '-' if alignment.is_reverse else '+'
            if ((chr_name not in TE_interval_tree_dict) or (strand not in TE_interval_tree_dict[chr_name])):
                continue
            all_reads_TE_mapping_dict[read_name]['num_alignments'] += 1
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
                    if ((mapping_strategy == "TE" and overlapped_length >= mapping_strategy_proportion * (TE_end - TE_begin)) or (mapping_strategy == "read" and overlapped_length >= mapping_strategy_proportion * (read_end - read_begin))):
                        filtered_TE.add(TE)
                TE_mappings = TE_mappings.union(filtered_TE)
            if (len(TE_mappings)>0):
                all_reads_TE_mapping_dict[read_name]['mapped_to_TE'] = True
    except:
        pass
    samfile.close()
    return all_reads_TE_mapping_dict,read_names,unmapped_reads
def count_unique_multi_reads_proportion(annotation_file_path,alignment_file_path,output_path):
    from parse_annotation import parse_TE_annotation
    TE_interval_tree_dict,TE_features_dict,TE_length_matrix = parse_TE_annotation(annotation_file_path)
    unique_read_names,multi_read_names = set(),set()
    with open(output_path,'w') as f:
#         for (strategy,proportion) in [("TE",0.5),("read",0.5),("TE",1.0),("read",1.0)]:
        for (strategy,proportion) in [("TE",0.0)]:
            all_reads_TE_mapping_dict,read_names,num_unmapped_reads = count_alignments(alignment_file_path,TE_interval_tree_dict,strategy,proportion)
            num_multi_reads = 0
            num_unique_reads = 0
            num_TE_mapped_multi_reads = 0
            num_TE_mapped_unique_reads = 0
            for read_name in all_reads_TE_mapping_dict:
                num_alignments = all_reads_TE_mapping_dict[read_name]['num_alignments']
                mapped_to_TE = all_reads_TE_mapping_dict[read_name]['mapped_to_TE']
                if (num_alignments > 1):
                    num_multi_reads += 1
                    if mapped_to_TE:
                        num_TE_mapped_multi_reads += 1
                elif (num_alignments == 1):
                    num_unique_reads += 1
                    if mapped_to_TE:
                        num_TE_mapped_unique_reads += 1
            num_TE_unmapped_reads = len(read_names) - num_TE_mapped_unique_reads - num_TE_mapped_multi_reads
            f.write("{}\t{}\t{}\t".format(num_unique_reads,num_multi_reads,num_unmapped_reads))
            f.write("{}\t{}\t{}\n".format(num_TE_mapped_unique_reads,num_TE_mapped_multi_reads,num_TE_unmapped_reads))
            