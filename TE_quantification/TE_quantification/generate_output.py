from collections import defaultdict
from pathlib import Path

def add_abundance(abundance_sum,abundance):
    abundance_sum['TPM'],abundance_sum['CPM'],abundance_sum['read_count'] = abundance_sum['TPM'] + abundance['TPM'],abundance_sum['CPM'] + abundance['CPM'],abundance_sum['read_count'] + abundance['read_count']
    return abundance_sum
def generate_output(TE_features_dict,TE_TPM,TE_CPM,TE_abundance_matrix,output_path):
    Path(output_path).mkdir(parents=True, exist_ok=True)
    TE_abundance_dict = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : {"TPM":0,"CPM":0,"read_count":0})))
    with open(output_path+'/TE_locus.out','w') as f:
        f.write('Repeat_name\tSubfamily_name\tFamily_name\tClass_name\tChr\tStrand\tStart\tEnd\tCPM\tTPM\tRead_count\tConsensus_length\tTE_Consensus_proportion\tTruncated_proportion\n')
        for class_name in TE_features_dict:
            for family_name in TE_features_dict[class_name]:
                for subfamily_name in TE_features_dict[class_name][family_name]:
                    for repeat_name in TE_features_dict[class_name][family_name][subfamily_name]:
                        TE_feature = TE_features_dict[class_name][family_name][subfamily_name][repeat_name]
                        TE_abundance_dict[class_name][family_name][subfamily_name]['TPM'] += TE_TPM[TE_feature['index']]
                        TE_abundance_dict[class_name][family_name][subfamily_name]['CPM'] += TE_CPM[TE_feature['index']]
                        TE_abundance_dict[class_name][family_name][subfamily_name]['read_count'] += TE_abundance_matrix[TE_feature['index']]
                        consensus_length = TE_feature['consensus_length']
                        truncated_proportion = TE_feature['truncated_proportion']
                        TE_Consensus_proportion = (TE_feature['end'] - TE_feature['start'])/TE_feature['consensus_length']
                        f.write('%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%f\t%f\t%d\t%d\t%f\t%f\n'%(repeat_name,subfamily_name,family_name,class_name,TE_feature['chr_name'],TE_feature['strand'],TE_feature['start'],TE_feature['end'],TE_CPM[TE_feature['index']],TE_TPM[TE_feature['index']],TE_abundance_matrix[TE_feature['index']],consensus_length,TE_Consensus_proportion,truncated_proportion))
    with open(output_path+'/TE_subfamily.out','w') as f_sub:
        with open(output_path+'/TE_family.out','w') as f_fam:
            with open(output_path+'/TE_class.out','w') as f_cla:
                f_sub.write('Subfamily_name\tFamily_name\tClass_name\tCPM\tTPM\tRead_count\n')
                f_fam.write('Family_name\tClass_name\tCPM\tTPM\tRead_count\n')
                f_cla.write('Class_name\tCPM\tTPM\tRead_count\n')
                for class_name in TE_abundance_dict:
                    abundance_class_sum = defaultdict(lambda : 0)
                    for family_name in TE_abundance_dict[class_name]:
                        abundance_family_sum = defaultdict(lambda : 0)
                        for subfamily_name in TE_abundance_dict[class_name][family_name]:
                            abundance = TE_abundance_dict[class_name][family_name][subfamily_name]
                            abundance_family_sum = add_abundance(abundance_family_sum,abundance)
                            f_sub.write('%s\t%s\t%s\t%f\t%f\t%d\n'%(subfamily_name,family_name,class_name,abundance['CPM'],abundance['TPM'],abundance['read_count']))
                        abundance_class_sum = add_abundance(abundance_class_sum,abundance_family_sum)
                        f_fam.write('%s\t%s\t%f\t%f\t%d\n'%(family_name,class_name,abundance_family_sum['CPM'],abundance_family_sum['TPM'],abundance_family_sum['read_count']))
                    f_cla.write('%s\t%f\t%f\t%d\n'%(class_name,abundance_class_sum['CPM'],abundance_class_sum['TPM'],abundance_class_sum['read_count']))