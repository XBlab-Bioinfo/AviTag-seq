# -*- coding:utf-8 -*-
#分流，合并PCR产物
import os
import gzip
import time
import yaml
import datetime
import subprocess

import datetime
import glob
import os.path


def consolidate(file_umi):
    '''
    把分流出来的数据（也就是umi标记的）合并PCR产物。
    里面利用Q值保留了同PCR产物中Q值最高的数据。
    '''
    pcr_con_read = 0
    if not os.path.exists("./consolidated"):
        os.mkdir("./consolidated")
    outf = file_umi.replace("umitagged.", "consolidated.").replace("umitagged", "consolidated")
    # outf = os.path.join(out_dir, of)
    outfile = open(outf, 'a')
    with open(file_umi, 'r') as f:
        s1 = f.readline()
        s2 = f.readline()
        s3 = f.readline()
        s4 = f.readline()
        front_umi_id, count = "", 0
        # all_q = sum([(ord(x) - 33) for x in s4.replace("\n", "")])

        while s1:
            cur_umi = s1.split(" ")[-1].replace("\n", "")
            if front_umi_id == "":
                front_seq, front_q, all_q = s2, s4, sum([(ord(x) - 33) for x in s4.replace("\n", "")])
                front_umi_id = s1.split(" ")[-1].replace("\n", "")
                count += 1
                s1 = f.readline()
                s2 = f.readline()
                s3 = f.readline()
                s4 = f.readline()
                continue
            if cur_umi != front_umi_id:
                pcr_con_read+=1
                outfile.write("@" + front_umi_id + "_{}\n".format(count) + front_seq + s3 + front_q)
                front_umi_id = cur_umi
                front_seq, front_q, all_q, count = s2, s4, sum([(ord(x) - 33) for x in s4.replace("\n", "")]), 1
            else:
                cur_q = sum([(ord(x) - 33) for x in s4.replace("\n", "")])
                if cur_q > all_q:
                    front_seq, front_q, all_q = s2, s4, cur_q
                count += 1
            s1 = f.readline()
            s2 = f.readline()
            s3 = f.readline()
            s4 = f.readline()
        f.close()
    outfile.close()
    return pcr_con_read
def dep(file1, file2, config,out_dirumi):
    '''
    利用adaptor、itr、barcode分流。
    '''
    outfiles_r1 = {}
    outfiles_r2 = {}
    # outfiles_i1 = {}
    # outfiles_i2 = {}
    samples_name = {}
    for sample, value in config["samples"].items():
        if sample != "control":
            samples_name["{}".format(config["samples"]["{}".format(sample)]["barcode1"])] = sample
            outfiles_r1[sample] = open(os.path.join(out_dirumi, '%s.r1.tempumitagged.fastq' % sample), 'a')
            outfiles_r2[sample] = open(os.path.join(out_dirumi, '%s.r2.tempumitagged.fastq' % sample), 'a')
            # outfiles_i1[sample] = open(os.path.join(out_dir, '%s.i1.umitagged.fastq' % sample), 'a')
            # outfiles_i2[sample] = open(os.path.join(out_dir, '%s.i2.fastq' % sample), 'a')

    # read file31   取 barcode  【0：8】     umi：【21：29】   R1：【49：151】
    # 先不考虑NP的
    adaptor = "GCTCCCTCGCCT"
    itr = "CTGCGCGCTCGCTCGC"  
    all_reads = 0
    umi_reads=0
    for file31,file32 in zip(file1,file2):
        f31 = gzip.open(file31, 'rb')
        f32 = gzip.open(file32, 'rb')
        s1_1 = f31.readline().decode('utf-8')
        s1_2 = f31.readline().decode('utf-8')
        s1_3 = f31.readline().decode('utf-8')
        s1_4 = f31.readline().decode('utf-8')

        s2_1 = f32.readline().decode('utf-8')
        s2_2 = f32.readline().decode('utf-8')
        # .replace("\n", "").replace("A", "B").replace("G", "D").replace("T",
        #                                                                                                     "A").replace(
        #     "C", "G").replace("B", "T").replace("D", "C")[::-1] + "\n"
        s2_3 = f32.readline().decode('utf-8')
        s2_4 = f32.readline().decode('utf-8')
        # .replace("\n", "")[::-1] + "\n"
        while s1_1:
            all_reads+=1
            if adaptor in s1_2[34:49] and s1_2[:8] in samples_name.keys() and itr in s2_2[0:40]:
                umi_reads+=1
                outfiles_r1['{}'.format(samples_name["{}".format(s1_2[:8])])].write(
                    s1_1.replace("\n", " ") + s1_2[26:34] + "_" + s1_2[46:52] + "_" + s2_2[:6] + "\n" + s1_2[
                                                                                                        46:] + s1_3 + s1_4[
                                                                                                                      46:])
                outfiles_r2['{}'.format(samples_name["{}".format(s1_2[:8])])].write(
                    s2_1.replace("\n", " ") + s1_2[26:34] + "_" + s1_2[46:52] + "_" + s2_2[:6] + "\n" + s2_2 + s2_3 + s2_4)
            elif adaptor not in s1_2[34:49] and adaptor in s2_2[34:49] and s2_2[:8] in samples_name.keys() and itr in s1_2[0:40]:
                umi_reads+=1
                outfiles_r1['{}'.format(samples_name["{}".format(s2_2[:8])])].write(
                    s2_1.replace("\n", " ") + s2_2[26:34] + "_" + s2_2[46:52] + "_" + s1_2[:6] + "\n" + s2_2[
                                                                                                        46:] + s2_3 + s2_4[
                                                                                                                      46:])
                outfiles_r2['{}'.format(samples_name["{}".format(s2_2[:8])])].write(
                    s1_1.replace("\n", " ") + s2_2[26:34] + "_" + s2_2[46:52] + "_" + s1_2[:6] + "\n" + s1_2 + s1_3 + s1_4)
                # outfiles_i1['{}'.format(samples_name["{}".format(s1_2[:8])])].write(
                #     s1_1 + s1_2[:8] + "\n" + s1_3 + s1_4[:8] + "\n")
                # outfiles_i2['{}'.format(samples_name["{}".format(s1_2[:8])])].write(
                #     s1_1 + s1_2[22:30] + "\n" + s1_3 + s1_4[22:30] + "\n")
            s1_1 = f31.readline().decode('utf-8')
            s1_2 = f31.readline().decode('utf-8')
            s1_3 = f31.readline().decode('utf-8')
            s1_4 = f31.readline().decode('utf-8')

            s2_1 = f32.readline().decode('utf-8')
            s2_2 = f32.readline().decode('utf-8')
            # .replace("\n", "").replace("A", "B").replace("G", "D").replace("T",
            #                                                                                                     "A").replace(
            #     "C", "G").replace("B", "T").replace("D", "C")[::-1] + "\n"
            s2_3 = f32.readline().decode('utf-8')
            s2_4 = f32.readline().decode('utf-8')
            # .replace("\n", "")[::-1] + "\n"
        f31.close()
        f32.close()
    for sample, value in config["samples"].items():
        if sample != "control":
            samples_name["{}".format(config["samples"]["{}".format(sample)]["barcode1"])] = sample
            outfiles_r1[sample].close()
            outfiles_r2[sample].close()
            # outfiles_i1[sample].close()
            # outfiles_i2[sample].close()

    for sample, value in config["samples"].items():#按照首字母顺序排序
        if sample != "control":
            r1_umitagged_unsorted_file = os.path.join(out_dirumi, '%s.r1.tempumitagged.fastq' % sample)
            r2_umitagged_unsorted_file = os.path.join(out_dirumi, '%s.r2.tempumitagged.fastq' % sample)
            read1_out = os.path.join(out_dirumi, '%s.r1.umitagged.fastq' % sample)
            read2_out = os.path.join(out_dirumi, '%s.r2.umitagged.fastq' % sample)
            cmd = 'cat ' + r1_umitagged_unsorted_file + ' | paste - - - - | sort -k3,3 -k1,1 | tr "\t" "\n" >' + read1_out
            subprocess.check_call(cmd, shell=True, env=os.environ.copy())
            cmd = 'cat ' + r2_umitagged_unsorted_file + ' | paste - - - - | sort -k3,3 -k1,1 | tr "\t" "\n" >' + read2_out
            subprocess.check_call(cmd, shell=True, env=os.environ.copy())
            os.remove(r1_umitagged_unsorted_file)
            os.remove(r2_umitagged_unsorted_file)

    return all_reads,umi_reads
def main(config,data1,data2):
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    out_dirumi = "./umitagged"
    if not os.path.exists(out_dirumi):
        os.mkdir(out_dirumi)
    
    file_1 = data1
    file_2 = data2

    all_reads_,all_umi_reads=dep(file_1, file_2, config,out_dirumi)
    fileall = glob.glob("./umitagged/*.fastq")
    all_pcr_reads = 0
    for file in fileall:
        all_pcr_reads += consolidate(file)
    with open("reads.txt","w") as ff:#统计每个阶段的reads
        ff.write("all reads: {}   (1)\numi marked: {}   ({})\npcr combined: {}   ({})\n".format(all_reads_,all_umi_reads,round(all_umi_reads/all_reads_,6),all_pcr_reads,round(all_pcr_reads/all_reads_,6)))
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
