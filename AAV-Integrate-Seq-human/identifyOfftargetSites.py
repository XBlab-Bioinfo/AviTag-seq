# -*- coding: utf-8 -*-
import argparse
import collections
import numpy
import os
import string
import operator
import pyfaidx
import re
import swalign
import logging
import tqdm
import time
import numpy as np
import sys

from wordcloud import WordCloud, ImageColorGenerator

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

logger = logging.getLogger('root')


# chromosomePosition defines a class to keep track of the positions.
class chromosomePosition():

    def __init__(self, reference_genome):
        self.chromosome_dict = {}
        self.chromosome_barcode_dict = {}
        self.position_summary = []
        self.index_stack = {}           # we keep track of the values by index here
        self.genome = pyfaidx.Fasta(reference_genome)

    def addPositionBarcode(self, chromosome, position, strand, barcode, primer, count):
        # Create the chromosome keyValue if it doesn't exist
        if chromosome not in self.chromosome_barcode_dict:
            self.chromosome_barcode_dict[chromosome] = {}
        # Increment the position on that chromosome if it exists, otherwise initialize it with 1
        if position not in self.chromosome_barcode_dict[chromosome]:
            self.chromosome_barcode_dict[chromosome][position] = {}
            self.chromosome_barcode_dict[chromosome][position]['+_total'] = 0
            self.chromosome_barcode_dict[chromosome][position]['+primer1_total'] = 0
            self.chromosome_barcode_dict[chromosome][position]['+primer2_total'] = 0
            self.chromosome_barcode_dict[chromosome][position]['+nomatch_total'] = 0

            self.chromosome_barcode_dict[chromosome][position]['-_total'] = 0
            self.chromosome_barcode_dict[chromosome][position]['-primer1_total'] = 0
            self.chromosome_barcode_dict[chromosome][position]['-primer2_total'] = 0
            self.chromosome_barcode_dict[chromosome][position]['-nomatch_total'] = 0

            self.chromosome_barcode_dict[chromosome][position]['+'] = collections.Counter()
            self.chromosome_barcode_dict[chromosome][position]['+primer1'] = collections.Counter()
            self.chromosome_barcode_dict[chromosome][position]['+primer2'] = collections.Counter()
            self.chromosome_barcode_dict[chromosome][position]['+nomatch'] = collections.Counter()

            self.chromosome_barcode_dict[chromosome][position]['-'] = collections.Counter()
            self.chromosome_barcode_dict[chromosome][position]['-primer1'] = collections.Counter()
            self.chromosome_barcode_dict[chromosome][position]['-primer2'] = collections.Counter()
            self.chromosome_barcode_dict[chromosome][position]['-nomatch'] = collections.Counter()
        if primer!="nomatch":
            self.chromosome_barcode_dict[chromosome][position][strand][barcode] += count
            self.chromosome_barcode_dict[chromosome][position][strand + primer][barcode] += count
            self.chromosome_barcode_dict[chromosome][position][strand + primer + '_total'] += count
            self.chromosome_barcode_dict[chromosome][position][strand + '_total'] += count

    def getSequence(self, genome, chromosome, start, end, strand="+"):
        
        if strand == "+":
            seq = self.genome[chromosome][int(start):int(end)]
        elif strand == "-":
            seq = self.genome[chromosome][int(start):int(end)].reverse.complement
        return seq

    # Generates a summary of the barcodes by position
    def SummarizeBarcodePositions(self):
        self.barcode_position_summary = [[chromosome, position,
                                          len(self.chromosome_barcode_dict[chromosome][position]['+']),
                                          len(self.chromosome_barcode_dict[chromosome][position]['-']),
                                          self.chromosome_barcode_dict[chromosome][position]['+_total'],
                                          self.chromosome_barcode_dict[chromosome][position]['-_total'],
                                          len(self.chromosome_barcode_dict[chromosome][position]['+primer1']),
                                          len(self.chromosome_barcode_dict[chromosome][position]['+primer2']),
                                          len(self.chromosome_barcode_dict[chromosome][position]['-primer1']),
                                          len(self.chromosome_barcode_dict[chromosome][position]['-primer2']),
                                          ]
                                         for chromosome in sorted(self.chromosome_barcode_dict)
                                         for position in sorted(self.chromosome_barcode_dict[chromosome])]
        return self.barcode_position_summary

    # Summarizes the chromosome, positions within a 10 bp window
    def SummarizeBarcodeIndex(self):
        last_chromosome, last_position, window_index = 0, 0, 0
        index_summary = []
        for chromosome, position, barcode_plus_count, barcode_minus_count, total_plus_count, total_minus_count, plus_primer1_count, plus_primer2_count,\
                minus_primer1_count, minus_primer2_count in self.barcode_position_summary:
            if chromosome != last_chromosome or abs(position - last_position) > 5:
                window_index += 1   # new index
                last_chromosome, last_position = chromosome, position
            if window_index not in self.index_stack:
                self.index_stack[window_index] = []
            self.index_stack[window_index].append([chromosome, int(position),
                                                   int(barcode_plus_count), int(barcode_minus_count),
                                                   int(barcode_plus_count) + int(barcode_minus_count),
                                                   int(total_plus_count), int(total_minus_count),
                                                   int(total_plus_count) + int(total_minus_count),
                                                   int(plus_primer1_count), int(plus_primer2_count),
                                                   int(minus_primer1_count), int(minus_primer2_count)
                                                   ])
        cc = 0
        for index in self.index_stack:
            cc+=1
            # print(cc)
            sorted_list = sorted(self.index_stack[index], key=operator.itemgetter(4))   # sort by barcode_count_total
            chromosome_list, position_list, \
                barcode_plus_count_list, barcode_minus_count_list, barcode_sum_list,\
                total_plus_count_list, total_minus_count_list, total_sum_list, \
                plus_primer1_list, plus_primer2_list, minus_primer1_list, minus_primer2_list\
                = zip(*sorted_list)
            
            barcode_plus = sum(barcode_plus_count_list)
            barcode_minus = sum(barcode_minus_count_list)
            total_plus = sum(total_plus_count_list)
            total_minus = sum(total_minus_count_list)
            plus_primer1 = sum(plus_primer1_list)
            plus_primer2 = sum(plus_primer2_list)
            minus_primer1 = sum(minus_primer1_list)
            minus_primer2 = sum(minus_primer2_list)
            position_std = numpy.std(position_list)
            min_position = min(position_list)
            max_position = max(position_list)
            
            barcode_sum = barcode_plus + barcode_minus
            barcode_geometric_mean = (barcode_plus * barcode_minus) ** 0.5
            total_sum = total_plus + total_minus
            total_geometric_mean = (total_plus * total_minus) ** 0.5
            primer1 = plus_primer1 + minus_primer1
            primer2 = plus_primer2 + minus_primer2
            primer_geometric_mean = (primer1 * primer1) ** 0.5
            most_frequent_chromosome = sorted_list[-1][0]
            most_frequent_position = sorted_list[-1][1]
            #BED_format_chromosome = "chr" + most_frequent_chromosome
            
            BED_format_chromosome = most_frequent_chromosome
            BED_name = BED_format_chromosome + "_" + str(most_frequent_position) + "_" + str(barcode_sum)
            # print(2222222222222222)
            # print(most_frequent_chromosome, most_frequent_position, most_frequent_position + 35)
            offtarget_sequence = self.getSequence(self.genome, most_frequent_chromosome, most_frequent_position, most_frequent_position + 35)
            # print("11111111111")
            summary_list = [str(x) for x in [index, most_frequent_chromosome, most_frequent_position, offtarget_sequence,                        # pick most frequently occurring chromosome and position
                                             BED_format_chromosome, min_position, max_position, BED_name,
                                             barcode_plus, barcode_minus, barcode_sum, barcode_geometric_mean,
                                             total_plus, total_minus, total_sum, total_geometric_mean,
                                             primer1, primer2, primer_geometric_mean, position_std]]
            # if (barcode_geometric_mean > 0 or primer_geometric_mean > 0):
            index_summary.append(summary_list)
            # index_summary.append(summary_list)
        # print(index_summary)
        return index_summary    # WindowIndex, Chromosome, Position, Plus.mi, Minus.mi,
        # BidirectionalArithmeticMean.mi, BidirectionalGeometricMean.mi,
        # Plus, Minus,
        # BidirectionalArithmeticMean, BidirectionalGeometricMean,


def alignSequences(ref_seq, query_seq):
    match = 2
    mismatch = -1
    ref_length = len(ref_seq)
    # matches_required = len(ref_seq) - 1 - 7  # allow up to 8 mismatches
    # LSA modified for SpCas9
    matches_required = len(ref_seq) - 1 - 6  # allow up to 7 mismatches
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring, gap_penalty=-100, gap_extension_penalty=-100, prefer_gap_runs=True)  # you can also choose gap penalties, etc...
    # sw = swalign.LocalAlignment(scoring, gap_penalty=-10, gap_extension_penalty=-0.5, prefer_gap_runs=True)  # you can also choose gap penalties, etc...
    forward_alignment = sw.align(ref_seq, query_seq)
    reverse_alignment = sw.align(ref_seq, reverseComplement(query_seq))
    if forward_alignment.matches >= matches_required and forward_alignment.matches > reverse_alignment.matches:
        start_pad = forward_alignment.r_pos
        start = forward_alignment.q_pos - start_pad
        end_pad = ref_length - forward_alignment.r_end
        end = forward_alignment.q_end + end_pad
        strand = "+"
        return [forward_alignment.query[start:end], ref_length - forward_alignment.matches - 1, end - start, strand, start, end]
    elif reverse_alignment.matches >= matches_required and reverse_alignment.matches > forward_alignment.matches:
        start_pad = reverse_alignment.r_pos
        start = reverse_alignment.q_pos - start_pad
        end_pad = ref_length - reverse_alignment.r_end
        end = reverse_alignment.q_end + end_pad
        strand = "-"
        return [reverse_alignment.query[start:end], ref_length - reverse_alignment.matches - 1, end - start, strand, start, end]
    else:
        return ["", "", "", "", "", ""]


"""
annotation is in the format:

"""
def analyze(sam_filename, reference_genome, outfile, annotations,out_cons,outimg):
    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    logger.info("Processing SAM file %s", sam_filename)
    file = open(sam_filename, 'rU')
    __, filename_tail = os.path.split(sam_filename)
    chromosome_position = chromosomePosition(reference_genome)
    mm = 0
    file_name = sam_filename[:-4].replace("aligned/","")
    ft = open("{}id.txt".format(file_name),'w')
    print(sam_filename)
    for line in file:
        fields = line.split('\t')
        if len(fields) >= 16:
            # These are strings--need to be cast as ints for comparisons.
            full_read_name, sam_flag, chromosome, position, mapq, cigar, name_of_mate, position_of_mate, template_length, read_sequence, read_quality = fields[:11]
            NM,MD,MC,AS = fields[11].split(":")[-1],fields[12].split(":")[-1],fields[13].split(":")[-1],fields[14].split(":")[-1]
            
            num,num_s,S_lis,M_lis = "",0,[],[]
            for s in cigar:
                if s == "S":
                    num_s = int(num)
                    S_lis.append(int(num))
                    num= ""
                if s == "M":
                    M_lis.append(int(num))
                    num = ""
                if not s.isalpha():
                    num+=s
            mm+=1
            
            if int(mapq) >= 60 and int(sam_flag) > 128 and "S" not in MC and int(NM)==0 and MD.isdigit() and 115<=int(MD)<=118 and -1000< int(template_length) < 0:
            # if int(mapq) >= 60:
            # if int(mapq) >= 60 and (int(sam_flag) ==163 or int(sam_flag) == 147) and (len(M_lis) == 1 and 50 <=sum(M_lis)):
            # if int(mapq) >= 60 and (int(sam_flag) ==163 or int(sam_flag) == 147):
             
                barcode, count = parseReadName(full_read_name)
                primer = assignPrimerstoReads(read_sequence, sam_flag)
                if int(template_length) < 0:  # Reverse read
                    read_position = int(position_of_mate) + abs(int(template_length)) - 1
                    strand = "-"
                    if primer == "primer1" or primer == "primer2":
                        ft.write("id:"+full_read_name+"\tfinal_pos:"+str(read_position)+"\tchromosome:"+chromosome+"\n")
                    chromosome_position.addPositionBarcode(chromosome, read_position, strand, barcode, primer, count)
                elif int(template_length) > 0:  # Forward read
                    read_position = int(position)
                    strand = "+"
                    if primer == "primer1" or primer == "primer2":
                        ft.write("id:"+full_read_name+"\tfinal_pos:"+str(read_position)+"\tchromosome:"+chromosome+"\n")
                    chromosome_position.addPositionBarcode(chromosome, read_position, strand, barcode, primer, count)
    ft.close()
    # Generate barcode position summary
    stacked_summary = chromosome_position.SummarizeBarcodePositions()
    with open(outfile, 'w') as f:
        # Write header
        f.write('\t'.join(['#BED Chromosome', 'BED Min.Position',
                           'BED Max.Position', 'Gene_id','Gene_Symbol','BED Name', 'Filename', 'WindowIndex', 'Chromosome', 'Position', 'Sequence', '+.mi', '-.mi', 'bi.sum.mi', 'bi.geometric_mean.mi', '+.total',
                           '-.total', 'total.sum', 'total.geometric_mean', 'primer1.mi', 'primer2.mi', 'primer.geometric_mean',
                           'position.stdev', 'Off-Target Sequence', 'Mismatches', 'Length', 'BED off-target Chromosome', 'BED off-target start', 'BED off-target end', 'BED off-target name', 'BED Score', 'Strand', 'Cells', 'Targetsite', 'Target Sequence']) + '\n')

        # Output summary of each window
        summary = chromosome_position.SummarizeBarcodeIndex()
        target_sequence = annotations["Sequence"]
        annotation = [annotations['Description'],
                      annotations['Targetsite'],
                      annotations['Sequence']]
        allinfo,allbed = [],{}
        for row in tqdm.tqdm(summary):
            window_sequence = row[3]

            if target_sequence:
                sequence, mismatches, length, strand, target_start_relative, target_end_relative = alignSequences(target_sequence, window_sequence)
                BED_chromosome = row[4]
                BED_name = row[7]
                BED_score = 1
                if strand == "+":
                    target_start_absolute = target_start_relative + int(row[2]) - 25
                    target_end_absolute = target_end_relative + int(row[2]) - 25
                elif strand == "-":
                    target_start_absolute = int(row[2]) + 25 - target_end_relative
                    target_end_absolute = int(row[2]) + 25 - target_start_relative
                else:
                    BED_chromosome, target_start_absolute, target_end_absolute, BED_score, BED_name = [""] * 5
                f.write('\t'.join(row[4:8] + [filename_tail] + row[0:4] + row[8:] +
                                  [str(x) for x in sequence, mismatches, length, BED_chromosome, target_start_absolute,
                                   target_end_absolute, BED_name, BED_score, strand] + [str(x) for x in annotation] + ['\n']))
            else:
            
                gene_id = ""
                '''make anno.bed'''
                line = "{}\t{}\t{}\n".format(row[4],row[5],row[6])
                allbed["{}+{}".format(row[4],row[5])]=None
                with open("anno.bed",'a') as bed:
                    bed.write(line)
                    bed.close()
                allinfo.append([str(x) for x in row[4:7]+[row[7]] + [filename_tail] + row[0:3] +[""+row[3].lower()]+ row[8:] + [""] * 9 + annotation] + ['\n'])
        bedtools = "bedtools intersect -a gencode.v29.primary_assembly.annotation.gtf -b anno.bed -u > result.txt"
        os.system(bedtools)
        result = [] 
        record_symbol = {}
        record_symbol["Nomatch"]="None"
        with open("result.txt",'r') as rr:
            old_geneid = ""
            con = rr.readline()
            chrr,start,end = con.split("\t")[0],con.split("\t")[3],con.split("\t")[4]
            
            old_geneid = con.split("\t")[8].split(" ")[1].strip(";\"")

            while con:
                if old_geneid not in record_symbol:
                    frot = ""
                    for ge in con.split("\t")[8].split(" "):
                        if frot =="gene_name":
                            record_symbol["{}".format(old_geneid)]=ge.strip('";')
                            break
                        frot = ge
                con_lis = con.split("\t")
                if old_geneid!=con_lis[8].split(" ")[1].strip(";\""):
                    result.append([chrr,start,end,old_geneid])
                    old_geneid = con_lis[8].split(" ")[1].strip(";\"")
                    chrr,start,end =con_lis[0], con_lis[3],con_lis[4]
                else:
                    if int(con_lis[3])<int(start):
                        start = con_lis[3]
                    if int(con_lis[4])>int(end):
                        end = con_lis[4]
                con = rr.readline()
            rr.close()
        # print(record_symbol)
        for bed in tqdm.tqdm(allbed):
            chr_,pos = bed.split("+")[0],bed.split("+")[1]
            flag_bed = False
            for res in result:
                if chr_ == res[0] and int(res[1])<=int(pos)<=int(res[2]):
                    allbed["{}".format(bed)]=res[3]
                    flag_bed = True
                    break
            if flag_bed == False:
                allbed["{}".format(bed)]="Nomatch"
        for info in tqdm.tqdm(allinfo):
            geneid = allbed["{}+{}".format(info[0],info[1])]
            gene_symbol = record_symbol["{}".format(geneid)].strip('";\n')
            if geneid != "Nomatch":
                f.write('\t'.join(info[:3]+[geneid]+[gene_symbol]+info[3:]))
        f.close
    # next_collapse_allinfo = []
    os.remove("anno.bed")
    os.remove("result.txt")
    with open(out_cons, 'w') as wr:
        wr.write(
            '\t'.join(['Chromosome', 'Position', 'Gene_id',  'Gene_Symbol', 'Sites_Num','consolidated', 'Sequence']) + '\n')
        with open(outfile, 'r') as co:
            head = co.readline()  # 取第一行
            if head.split("\t")[4] != "Gene_Symbol":
                print("the file error.")
                sys.exit()
            main_line_ = co.readline()  # 取第二行
            main_line, collapse, sites_num = main_line_.split("\t"), {}, {}  # 分开
            front_pos, front_geneid, front_symbol, front_chr, front_seq, front_pri = main_line[2], main_line[3], main_line[
                4], main_line[8], main_line[10], int(main_line[13])
            collapse["{}".format(front_symbol)] = front_pri
            if front_pri == 0:
                sites_num["{}".format(front_symbol)] = 0
            else:
                sites_num["{}".format(front_symbol)] = 1
            main_line_ = co.readline()
            while main_line_:
                main_lis = main_line_.split("\t")
                if main_lis[3] == "Nomatch":
                    main_line_ = co.readline()
                    continue
                if main_lis[4] not in collapse:
                    if collapse["{}".format(front_symbol)] != 0:
                        wr.write('\t'.join(
                            [front_chr, front_pos, front_geneid, front_symbol, str(sites_num["{}".format(front_symbol)]),
                            str(collapse["{}".format(front_symbol)]), front_seq]) + '\n')
                    # next_collapse_allinfo.append([front_pos,front_geneid,front_symbol,collapse["{}".format(front_symbol)],front_chr,front_seq])
                    front_pos, front_geneid, front_symbol, front_chr, front_seq, front_pri = main_lis[2], main_lis[3], \
                                                                                            main_lis[4], main_lis[8], \
                                                                                            main_lis[10], int(main_lis[13])
                    collapse["{}".format(main_lis[4])] = front_pri
                    sites_num["{}".format(main_lis[4])] = 1
                else:
                    collapse["{}".format(main_lis[4])] += int(main_lis[13])
                    if int(main_lis[13]) != 0:
                        sites_num["{}".format(main_lis[4])] += 1
                main_line_ = co.readline()
            if collapse["{}".format(front_symbol)] != 0:
                wr.write('\t'.join(
                    [front_chr, front_pos, front_geneid, front_symbol, str(sites_num["{}".format(front_symbol)]),
                    str(collapse["{}".format(front_symbol)]), front_seq]) + '\n')
            co.close()
        wr.close()
        draw_cloud(collapse,outimg)
def draw_cloud(dic,out):
    image = Image.open('bg.png')  # 作为背景轮廓图
    graph = np.array(image)
    # 参数分别是指定字体、背景颜色、最大的词的大小、使用给定图作为背景形状
    wc = WordCloud( background_color='#FAFAD2', max_words=200, mask=graph)#font_path='simkai.ttf'
    
    wc.generate_from_frequencies(dic)  # 根据给定词频生成词云
    image_color = ImageColorGenerator(graph)
    #plt.imshow(wc)
    #plt.axis("off")  # 不显示坐标轴
    #plt.show()
    wc.to_file(out)  # 图片命名




def assignPrimerstoReads(read_sequence, sam_flag):
    # Get 20-nucleotide sequence from beginning or end of sequence depending on orientation
    if int(sam_flag) & 16:
        readstart = reverseComplement(read_sequence[-20:])
    else:
        readstart = read_sequence[:20]
    # if readstart == "TTGAGTTGTCATATGTTAAT":
    #     return "primer1"
    # elif readstart == "ACATATGACAACTCAATTAA":
    #     return "primer2"
    # XB 2020-10-20: Change the above 4 lines to
    # if readstart == "AGTTGTCATATGTTAATAAC":
        # return "primer1"
    # elif readstart == "ACATATGACAACTCAATTAA":
        # return "primer2"
    if readstart == "GAGTTGGCCACTCCCTCTCTG":
        return "primer1"
    elif readstart == "TCTCTGCGCGCTCGCTCG":
        return "primer2"
    else:
        return "nomatch"


def loadFileIntoArray(filename):
    with open(filename, 'rU') as f:
        keys = f.readline().rstrip('\r\n').split('\t')[1:]
        data = collections.defaultdict(dict)
        for line in f:
            filename, rest = processLine(line)
            line_to_dict = dict(zip(keys, rest))
            data[filename] = line_to_dict
    return data


def parseReadName(read_name):
    m = re.search(r'([ACGTN]{7}_[ACGTN]{6}_[ACGTN]{6})_([0-9]*)', read_name)#Should be improved to randomly match 1-9 characters
    if m:
        molecular_index, count = m.group(1), m.group(2)
        return molecular_index, int(count)
    else:
        # print read_name
        return None, None


def processLine(line):
    fields = line.rstrip('\r\n').split('\t')
    filename = fields[0]
    rest = fields[1:]
    return filename, rest


def reverseComplement(sequence):
    transtab = string.maketrans("ACGTacgt", "TGCATGCA")
    return sequence.translate(transtab)[::-1]


def main():
    # This sets up the command line components of the program.
    parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
    parser.add_argument('--ref', help='Reference Genome Fasta', required=True)
    parser.add_argument('--samfile', help='SAM file', nargs='*')
    parser.add_argument('--outfile', help='File to output identified sites to.', required=True)
    parser.add_argument('--demo')
    parser.add_argument('--target', default='')

    args = parser.parse_args()

    annotations = {'Description': 'test description', 'Targetsite': 'dummy targetsite', 'Sequence': args.target}
    analyze(args.samfile[0], args.ref, args.outfile, annotations)


if __name__ == "__main__":
    # Run main program
    main()
