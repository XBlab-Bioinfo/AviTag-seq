import svgwrite
import sys
import os

boxWidth = 10
box_size = 15
v_spacing = 3

colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3', 'V': '#B3B3B3'}


def parseSitesFile(infile):
    offtargets = []
    chrs = []
    with open(infile, 'r') as f:
        f.readline()
        seq_dic = {}
        for line in f:
            line_items = line.split('\t')
            offtarget_sequence = line_items[21]
            offtarget_reads = line_items[11]
            ref_seq_ = line_items[32]
            dif_chr = line_items[0][3:]
            if offtarget_sequence != '':
                key = '{}_{}'.format(dif_chr.strip(), offtarget_sequence.strip())
                if key not in seq_dic:
                    seq_dic[key] = int(offtarget_reads.strip())
                else:
                    seq_dic[key] += int(offtarget_reads.strip())
        for key, value in seq_dic.items():
            offtargets.append({'seq': key.split('_')[1], 'reads': value})
            chrs.append({'seq': key.split('_')[0] + '-' + key.split('_')[1], 'reads': value})

    # 自动检测分割依据字符（'V' 或 'N'）
    split_char = detect_split_char(ref_seq_)

    offtargets = sorted(offtargets, key=lambda x: x['reads'], reverse=True)
    chrs = sorted(chrs, key=lambda x: x['reads'], reverse=True)
    chrs = onTargetsTop_chr(chrs, ref_seq_, split_char)  # 传递 split_char
    chrs = [item['seq'] for item in chrs]
    return offtargets, ref_seq_, chrs


def detect_split_char(ref_seq):
    """自动检测分割依据字符，优先使用 'N'，如果没有则使用 'V'"""
    if 'N' in ref_seq:
        return 'N'
    elif 'V' in ref_seq:
        return 'V'
    else:
        raise ValueError("Reference sequence does not contain 'N' or 'V' for splitting")

def onTargetsTop_chr(chrs, ref_seq, split_char):
    tempSeq_front = ref_seq.split(split_char)[0]
    tempSeq_end = ref_seq.split(split_char)[1]
    chr_off, chr_on = [], []
    for chr in chrs:
        seq = chr["seq"].split('-')[1]
        if seq[:len(tempSeq_front)] == tempSeq_front and seq[len(tempSeq_front)+1:] == tempSeq_end:
            chr_on.append(chr)
        else:
            chr_off.append(chr)
    return chr_on + chr_off

def onTargetsTop(offtargets, ref_seq, split_char):
    tempSeq_front = ref_seq.split(split_char)[0]
    tempSeq_end = ref_seq.split(split_char)[1]
    offtarget_off, offtarget_on = [], []
    for target in offtargets:
        if target["seq"][:len(tempSeq_front)] == tempSeq_front and target["seq"][len(tempSeq_front)+1:] == tempSeq_end:
            offtarget_on.append(target)
        else:
            offtarget_off.append(target)
    return offtarget_on + offtarget_off

def visualizeOfftargets(infile, outfile, title=None):
    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    offtargets, ref_seq, chrs = parseSitesFile(infile)

    # 自动检测分割字符
    split_char = detect_split_char(ref_seq)

    offtargets = onTargetsTop(offtargets, ref_seq, split_char)

    dwg = svgwrite.Drawing(outfile + '.svg', profile='full')

    if title is not None:
        x_offset = 20
        y_offset = 50
        dwg.add(dwg.text(title, insert=(x_offset, 30), style="font-size:20px; font-family:Arial"))
    else:
        x_offset = 20
        y_offset = 20

    tick_locations = [1, len(ref_seq)]
    tick_locations += range(len(ref_seq) + 1)[::10][1:]
    for x in tick_locations:
        dwg.add(dwg.text(str(x), insert=(x_offset + (x - 1) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Arial"))

    for i, c in enumerate(ref_seq):
        y = y_offset
        x = x_offset + i * box_size
        dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
        dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Arial"))

    dwg.add(dwg.text('Reads', insert=(x_offset + box_size * len(ref_seq) + 16, y_offset + box_size - 3), style="font-size:15px; font-family:Arial"))

    y_offset += 10
    for j, seq in enumerate(offtargets):
        y = y_offset + j * box_size
        for i, c in enumerate(seq['seq']):
            x = x_offset + i * box_size
            if c == ref_seq[i]:
                dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Arial"))
            else:
                dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Arial"))

        reads_text = dwg.text(str(seq['reads']), insert=(box_size * (len(ref_seq) + 1) + 20, y_offset + box_size * (j + 2) - 2), fill='black', style="font-size:15px; font-family:Arial")
        dwg.add(reads_text)

    dwg.save()

    off_on_target_reads = sum(int(item['reads']) for item in offtargets)
    f = open(outfile + '.txt', 'w')
    f.write('chr-seq')
    f.write('\t')
    f.write(ref_seq)
    f.write('\t')
    f.write('reads')
    f.write('\n')
    count = 0
    for value in offtargets:
        f.write(chrs[count])
        f.write('\t')
        f.write(str(value['reads']))
        f.write('\n')
        count += 1
    f.close()

    return off_on_target_reads

def main():
    if len(sys.argv) >= 3:
        if len(sys.argv) == 4:
            title = sys.argv[3]
        else:
            title = None
        visualizeOfftargets(sys.argv[1], sys.argv[2], title=title)
    else:
        print('Usage: python visualization.py INFILE OUTFILE [TITLE]')

if __name__ == '__main__':
    main()
