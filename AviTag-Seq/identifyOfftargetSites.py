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

        self.chromosome_barcode_dict[chromosome][position][strand][barcode] += count
        self.chromosome_barcode_dict[chromosome][position][strand + primer][barcode] += count
        self.chromosome_barcode_dict[chromosome][position][strand + primer + '_total'] += count
        self.chromosome_barcode_dict[chromosome][position][strand + '_total'] += count

    def getSequence(self, genome, chromosome, start, end, strand="+"):
        if strand == "+":
            seq = self.genome[chromosome][int(start):int(end)]
        elif strand == "-":
            seq = self.genome[chromosome][int(start):int(end)].reverse.complement
        
        # print("==============={}".format(seq))
        return str(seq).replace("\r","")

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
            if chromosome != last_chromosome or abs(position - last_position) > 10:
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
        for index in self.index_stack:
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
            primer_geometric_mean = (primer1 * primer2) ** 0.5
            most_frequent_chromosome = sorted_list[-1][0]
            most_frequent_position = sorted_list[-1][1]
            BED_format_chromosome = "chr" + most_frequent_chromosome
            BED_name = BED_format_chromosome + "_" + str(most_frequent_position) + "_" + str(barcode_sum)
            offtarget_sequence = self.getSequence(self.genome, most_frequent_chromosome, most_frequent_position - 25, most_frequent_position + 25)

            summary_list = [str(x) for x in [index, most_frequent_chromosome, most_frequent_position, offtarget_sequence,                        # pick most frequently occurring chromosome and position
                                             BED_format_chromosome, min_position, max_position, BED_name,
                                             barcode_plus, barcode_minus, barcode_sum, barcode_geometric_mean,
                                             total_plus, total_minus, total_sum, total_geometric_mean,
                                             primer1, primer2, primer_geometric_mean, position_std]]

            if (barcode_geometric_mean > 0 or primer_geometric_mean > 0):
                index_summary.append(summary_list)
        return index_summary    # WindowIndex, Chromosome, Position, Plus.mi, Minus.mi,

import time
def alignSequences(ref_seq, query_seq):
    # print(len(ref_seq),len(query_seq))
    # time.sleep(60)
    match = 2
    mismatch = -1
    ref_length = len(ref_seq)
    # matches_required = len(ref_seq) - 1 - 9  # allow up to 8 mismatches
    # LSA modified for SpCas9
    # matches_required = len(ref_seq) - 1 - 6  # allow up to 7 mismatches
    matches_required = len(ref_seq) - 1 - 5  # allow up to 6 mismatches
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring, gap_penalty=-100, gap_extension_penalty=-100, prefer_gap_runs=True)  # you can also choose gap penalties, etc...
    forward_alignment = sw.align(ref_seq, query_seq)
    reverse_alignment = sw.align(ref_seq, reverseComplement(query_seq))
    if forward_alignment.matches >= matches_required and forward_alignment.matches > reverse_alignment.matches:
        start_pad = forward_alignment.r_pos
        start = forward_alignment.q_pos - start_pad
        end_pad = ref_length - forward_alignment.r_end
        end = forward_alignment.q_end + end_pad
        strand = "+"
        if end>len(forward_alignment.query):
            return ["", ref_length - forward_alignment.matches - 1, end - start, strand, start, end]
        return [forward_alignment.query[start:end], ref_length - forward_alignment.matches - 1, end - start, strand, start, end]
    elif reverse_alignment.matches >= matches_required and reverse_alignment.matches > forward_alignment.matches:
        start_pad = reverse_alignment.r_pos
        start = reverse_alignment.q_pos - start_pad
        end_pad = ref_length - reverse_alignment.r_end
        end = reverse_alignment.q_end + end_pad
        strand = "-"
        if end>len(reverse_alignment.query):
            return ["", ref_length - reverse_alignment.matches - 1, end - start, strand, start, end]
        return [reverse_alignment.query[start:end], ref_length - reverse_alignment.matches - 1, end - start, strand, start, end]
    else:
        return ["", "", "", "", "", ""]

def analyze(sam_filename, reference_genome, outfile, annotations):
    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    logger.info("Processing SAM file %s", sam_filename)
    file = open(sam_filename, 'rU')
    __, filename_tail = os.path.split(sam_filename)
    chromosome_position = chromosomePosition(reference_genome)

    # total_reads and use_reads
    total_reads = 0
    use_reads = 0
    primer11,primer12,primer21,primer22,nomatchh = 0,0,0,0,0
    ttargget = 0


    for line in file:
        fields = line.split('\t')
        if len(fields) >= 10:
            total_reads +=1
            # These are strings--need to be cast as ints for comparisons.
            full_read_name, sam_flag, chromosome, position, mapq, cigar, name_of_mate, position_of_mate, template_length, read_sequence, read_quality = fields[:11]
            if int(mapq) >= 50 and int(sam_flag) & 128 and not int(sam_flag) & 2048:
                if 'GCATTTCAGGGCAAGGTTTAA' in read_sequence:
                    ttargget +=1
                # Second read in pair
                barcode, count = parseReadName(full_read_name)
                primer,stas_primer = assignPrimerstoReads(read_sequence, sam_flag)
                use_reads += count

                if stas_primer == 'primer11':
                    primer11 += 1
                elif stas_primer == 'primer12':
                    primer12 += 1
                elif stas_primer == 'primer21':
                    primer21 += 1
                elif stas_primer == 'primer22':
                    primer22 += 1
                elif stas_primer =='nomatch':
                    nomatchh += 1

                if int(template_length) < 0:  # Reverse read
                    read_position = int(position_of_mate) + abs(int(template_length)) - 1
                    strand = "-"
                    chromosome_position.addPositionBarcode(chromosome, read_position, strand, barcode, primer, count)
                elif int(template_length) > 0:  # Forward read
                    read_position = int(position)
                    strand = "+"
                    chromosome_position.addPositionBarcode(chromosome, read_position, strand, barcode, primer, count)

    # Generate barcode position summary
    stacked_summary = chromosome_position.SummarizeBarcodePositions()

    with open(outfile, 'w') as f:
        f.write('\t'.join(['#BED Chromosome', 'BED Min.Position',
                           'BED Max.Position', 'BED Name', 'Filename', 'WindowIndex', 'Chromosome', 'Position', 'Sequence', '+.mi', '-.mi', 'bi.sum.mi', 'bi.geometric_mean.mi', '+.total',
                           '-.total', 'total.sum', 'total.geometric_mean', 'primer1.mi', 'primer2.mi', 'primer.geometric_mean',
                           'position.stdev', 'Off-Target Sequence', 'Mismatches', 'Length', 'BED off-target Chromosome', 'BED off-target start', 'BED off-target end', 'BED off-target name', 'BED Score', 'Strand', 'Cells', 'Targetsite', 'Target Sequence']) + '\n')

        # Output summary of each window
        summary = chromosome_position.SummarizeBarcodeIndex()
        target_sequence = annotations["Sequence"]
        annotation = [annotations['Description'],
                      annotations['Targetsite'],
                      annotations['Sequence']]
        for row in summary:
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
                                  [str(x) for x in (sequence, mismatches, length, BED_chromosome, target_start_absolute,
                                   target_end_absolute, BED_name, BED_score, strand)] + [str(x) for x in annotation] + ['\n']))
            else:
                f.write('\t'.join([str(x) for x in row[4:8] + [filename_tail] + row[0:4] + row[8:] + [""] * 9 + annotation] + ['\n']))

    print("primer11:{},primer12:{},primer21:{},primer22:{},nomatch:{}".format(primer11,primer12,primer21,primer22,nomatchh))
    print("ttagget:{}".format(ttargget))
    return total_reads//2,use_reads//2,primer11,primer12,primer21,primer22,nomatchh

def assignPrimerstoReads(read_sequence, sam_flag):
    if int(sam_flag) & 16:
        # readstart = reverseComplement(read_sequence[-20:])
        readstart = reverseComplement(read_sequence)
    else:
        # readstart = read_sequence[:20]
        readstart = read_sequence

    # old primer11
    if "AGTTGTCATATGTTAATAAC" in readstart:
        return "primer1", "primer11"
    # new(2024-1-26) primer12
    elif "AGGAACCCCTAGTGATGGAG" in readstart:
        return "primer1", "primer12"
    
    # old primer21
    elif "ACATATGACAACTCAATTAA" in readstart:
        return "primer2", "primer21"
    # new(2024-1-26) primer22
    elif "GAGTTGGCCACTCCCTCT" in readstart:
        return "primer2", "primer22"
    else:
        return "nomatch", "nomatch"


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
    m = re.search(r'([ACGTN]{8}_[ACGTN]{6}_[ACGTN]{6})_([0-9]*)', read_name)
    if m:
        molecular_index, count = m.group(1), m.group(2)
        return molecular_index, int(count)
    else:
        return None, None


def processLine(line):
    fields = line.rstrip('\r\n').split('\t')
    filename = fields[0]
    rest = fields[1:]
    return filename, rest


def reverseComplement(sequence):
    transtab = str.maketrans("ACGTacgt", "TGCATGCA")
    return sequence.translate(transtab)[::-1]
