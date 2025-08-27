from intervaltree import Interval, IntervalTree
from collections import defaultdict
from util import sync_reference_name,parse_GTF_feature,create_GTF_feature
# GTF 1-BASED
def parse_isoform_annotation(isoform_annotation):
    isoform_dict = defaultdict(lambda : defaultdict(IntervalTree))
    isoform_length_dict = defaultdict(lambda: {})
    with open(isoform_annotation,'r') as f:
        for line in f.readlines():
            if not line.startswith('#'):
                fields = line.split('\t')
                if fields[2] == 'exon':
                    chrom,strand,start,end = sync_reference_name(fields[0]),fields[6],int(fields[3]),int(fields[4])
                    feature_dict = parse_GTF_feature(fields[8])
                    isoform_name,exon_number,gene = feature_dict['transcript_id'],feature_dict['exon_number'],feature_dict['gene_id']
                    isoform_length_dict[isoform_name][exon_number] = {'isoform_name':isoform_name,'gene':gene,'length':end-start+1,'exon_number':exon_number}
                    if strand == '.':
                        isoform_dict[chrom]['+'].addi(start,end + 1,(isoform_name,exon_number))
                        isoform_dict[chrom]['-'].addi(start,end + 1,(isoform_name,exon_number))
                    else:
                        isoform_dict[chrom][strand].addi(start,end + 1,(isoform_name,exon_number))
    return isoform_dict,isoform_length_dict

def intersect_interval(interval,new_interval):
    if interval.overlaps(new_interval.begin,new_interval.end):
        return Interval(max(interval.begin,new_interval.begin),min(interval.end,new_interval.end))
    else:
        return False
def calculate_overlap_TE_proportion(isoform_length_dict,overlapped_interval_dict,overlapped_TE_dict,out_isoform_path):
    f1 = open('{}_exon.tsv'.format(out_isoform_path),'w')
    overlapped_isoform_dict = defaultdict(lambda:defaultdict(lambda:[]))
    all_TE_dict = {}
    with open(out_isoform_path,'w') as f:
        f.write('isoform\tgene\tisoform_length\toverlapped_length\toverlapped_proportion\tmax_overlapped_TE\tmax_overlapped_length\tsubfamily\tfamily\tclass\n')
        for isoform_name in isoform_length_dict:
            gene = list(isoform_length_dict[isoform_name].values())[0]['gene']
            if isoform_name in overlapped_TE_dict:
                total_overlapped_length = 0
                isoform_length = 0
                max_overlapped_TE_length = 0
                max_overlapped_TE_name = ''
                max_overlapped_TE_dict = {}
                for exon_number in isoform_length_dict[isoform_name]:
                    exon_length = isoform_length_dict[isoform_name][exon_number]['length']
                    isoform_length += exon_length
                    if exon_number in overlapped_TE_dict[isoform_name]:
                        unioned = []
                        for interval,TE_name in sorted(zip(overlapped_interval_dict[isoform_name][exon_number],overlapped_TE_dict[isoform_name][exon_number]['name'])):
                            begin,end = interval.begin,interval.end
                            overlapped_TE_length = (end - 1) - begin + 1
                            overlapped_isoform_dict[TE_name][isoform_name].append((exon_number,overlapped_TE_length))
                            all_TE_dict[TE_name] = overlapped_TE_dict[isoform_name][exon_number]['dict'][TE_name]
                            if overlapped_TE_length > max_overlapped_TE_length:
                                max_overlapped_TE_length = overlapped_TE_length
                                max_overlapped_TE_name = TE_name
                                max_overlapped_TE_dict = overlapped_TE_dict[isoform_name][exon_number]['dict'][TE_name]
                            if unioned and unioned[-1].end >= begin - 1:
                                unioned[-1] = Interval(unioned[-1].begin,max(unioned[-1].end, end))
                            else:
                                unioned.append(Interval(begin, end))
                        exon_overlapped_length = 0
                        for interval in unioned:
                            exon_overlapped_length += (interval.end - 1) - interval.begin + 1
                        total_overlapped_length += exon_overlapped_length
                        f1.write('{}\t{}\t{}\t{}\n'.format(isoform_name,exon_number,exon_overlapped_length,';'.join(overlapped_TE_dict[isoform_name][exon_number]['name'])))
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(isoform_name,gene,isoform_length,total_overlapped_length,total_overlapped_length/isoform_length,max_overlapped_TE_name,max_overlapped_TE_length,max_overlapped_TE_dict['subfamily'],max_overlapped_TE_dict['family'],max_overlapped_TE_dict['class']))
            else:
                f.write('{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n'.format(isoform_name,gene))
    f1.close()
    with open('{}_TE.tsv'.format(out_isoform_path),'w') as f:
        f.write('loci\tsubfamily\tfamily\tclass\tTE_length\toverlapped_isoform\n')
        for TE_name in overlapped_isoform_dict:
            TE_dict = all_TE_dict[TE_name]
            overlapped_isoform = []
            for isoform_name in overlapped_isoform_dict[TE_name]:
                overlapped_length = sum([length for (exon,length) in overlapped_isoform_dict[TE_name][isoform_name]])
                overlapped_isoform.append((isoform_name,overlapped_length))
            overlapped_isoform = sorted(overlapped_isoform,key=lambda tup:tup[1],reverse=True)
            overlapped_isoform_str = ';'.join(['{}({:.2%})'.format(isoform_name,length/TE_dict['length']) for (isoform_name,length) in overlapped_isoform])
            f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(TE_name,TE_dict['subfamily'],TE_dict['family'],TE_dict['class'],TE_dict['length'],overlapped_isoform_str))

def create_curated_annotation(isoform_annotation,TE_annotation,out_path,out_isoform_path):
    print('Start creating curated annotation...')
    isoform_dict,isoform_length_dict = parse_isoform_annotation(isoform_annotation)
    overlapped_interval_dict = defaultdict(lambda :defaultdict(lambda:[]))
    overlapped_TE_dict = defaultdict(lambda : defaultdict(lambda:{'name':[],'dict':{}}))
    repNameIndex = defaultdict(lambda : 0)
    TE_gtf_lines = []
    # repeat masker out use 0-based indexing, covert to 1-based here
    with open(TE_annotation,'r') as f:
        for line in f.readlines()[3:]:
            if not line.startswith('#'):
                fields = [f for f in line.split(' ') if f != '']
                chrom,start,end = sync_reference_name(fields[4]),int(fields[5])+1,int(fields[6])
                TE_length = end - start + 1
                repSubfamily = fields[9]
                if '/' in fields[10]:
                    repClass,repFamily = fields[10].split('/')
                else:
                    repClass = fields[10]
                    repFamily = repSubfamily
                repeat_name = '{}_dup{}'.format(repSubfamily,repNameIndex[repSubfamily]+1)
                if fields[8] == 'C':
                    strand = '-'
                else:
                    strand = fields[8]
                if not (chrom in isoform_dict) and (strand in isoform_dict[chrom]):
                    continue
                exon_candidates = isoform_dict[chrom][strand].overlap(start,end+1)
                overlap_with_isoform = False
                isoform_candidates = defaultdict(lambda:[])
                for exon in exon_candidates:
                    isoform,exon_number = exon.data
                    isoform_candidates[isoform].append(exon)
                for isoform in isoform_candidates:
                    # intervaltree exclude right edge
                    overlapped_length = 0
                    for exon in isoform_candidates[isoform]:
                            # 1 based length: END - START + 1
                        overlapped_length += (exon.end - exon.begin + 1) + TE_length - (max(exon.end,end) - min(exon.begin,start)+1)
                    if overlapped_length > 10:
                        overlap_with_isoform = True
                        for exon in isoform_candidates[isoform]:
                            isoform_name,exon_number = exon.data
                            intersected_interval = intersect_interval(exon,Interval(start,end+1))
                            overlapped_interval_dict[isoform_name][exon_number].append(intersected_interval)
                            overlapped_TE_dict[isoform_name][exon_number]['name'].append(repeat_name)
                            overlapped_TE_dict[isoform_name][exon_number]['dict'][repeat_name] = {'length':TE_length,'subfamily':repSubfamily,'family':repFamily,'class':repClass}
                if not overlap_with_isoform:
                    continue
                repStart = abs(int(fields[11].strip('()')))
                repEnd = abs(int(fields[12]))
                repLeft = abs(int(fields[13].strip('()')))
                if strand == '+':
                    consensus_len =  repEnd + repLeft
                    truncated_proportion = (repEnd - repStart) / consensus_len
                elif strand == '-':
                    consensus_len =  repStart + repEnd
                    truncated_proportion = (repEnd - repLeft) / consensus_len
                else:
                    continue
                repNameIndex[repSubfamily] += 1
                features = [('gene_id',repSubfamily),
                            ('transcript_id',repeat_name),
                            ('family_id',repFamily),
                            ('class_id',repClass),
                            ('consensus_length',str(consensus_len)),
                           ('truncated_proportion',str(truncated_proportion))]
                feature_str = create_GTF_feature(features)
                gtf_ln = '{}\tuser_provided\texon\t{}\t{}\t.\t{}\t.\t{}\n'.format(chrom,str(start),str(end),strand,feature_str)
                TE_gtf_lines.append(gtf_ln)
    calculate_overlap_TE_proportion(isoform_length_dict,overlapped_interval_dict,overlapped_TE_dict,out_isoform_path)
#     return TE_gtf_lines,isoform_TE_overlapped_dict
    with open(out_path,'w') as f:
        for ln in TE_gtf_lines:
            f.write(ln)
    print('Done!')