# -*- coding: utf-8 -*-
"""

guideseq.py
===========
serves as the wrapper for all guideseq pipeline

"""

from itertools import count
import os
import sys
import yaml
import argparse
import traceback
import log
import tagged

logger = log.createCustomLogger('root')

from alignReads import alignReads
import identifyOfftargetSites
import statisticalTopSites

DEFAULT_DEMULTIPLEX_MIN_READS = 10000
MAX_MISMATCHES = 6

CONSOLIDATE_MIN_QUAL = 15
CONSOLIDATE_MIN_FREQ = 0.9


class GuideSeq:

    def __init__(self):
        pass

    def parseManifest(self, manifest_path):
        logger.info('Loading manifest...')

        self.manifest_data = yaml.safe_load(open(manifest_path, 'r'))

        # with open(manifest_path, 'r') as f:
        #     manifest_data = yaml.load(f)

        try:
            # Validate manifest data
            # validation.validateManifest(manifest_data)

            self.BWA_path = self.manifest_data['bwa']
            self.bedtools = self.manifest_data['bedtools']
            self.reference_genome = self.manifest_data['reference_genome']
            self.output_folder = self.manifest_data['output_folder']
            # self.undemultiplexed = self.manifest_data['undemultiplexed']
            self.samples = self.manifest_data['samples']
            # self.ori_data = self.manifest_data['ori_data']

            self.data1 = []
            self.data2 = []
            for filename in os.listdir(self.output_folder):
                if filename.endswith('fq.gz'):
                    if '_1.' in filename or '_1_' in filename:
                        self.data1.append(filename)
                    elif '_2.' in filename or '_2_' in filename:
                        self.data2.append(filename)
            print("data1:",self.data1)
            print("data2:",self.data2)

        except Exception as e:
            logger.error(
                'Incorrect or malformed manifest file. Please ensure your manifest contains all required fields.')
            sys.exit()

        # Allow the user to specify min reads for demultiplex if they want
        if 'demultiplex_min_reads' in self.manifest_data:
            self.demultiplex_min_reads = self.manifest_data['demultiplex_min_reads']
        else:
            self.demultiplex_min_reads = DEFAULT_DEMULTIPLEX_MIN_READS

        # Make sure the user has specified a control barcode
        if 'control' not in self.samples.keys():
            raise AssertionError('Your manifest must have a control sample specified.')

        # Make sure the user has both a sample and a control
        if len(self.samples) < 2:
            raise AssertionError('Your manifest must have at least one control and one treatment sample.')

        logger.info('Successfully loaded manifest.')

    def dataTagged(self):
        logger.info('Tagged reads...')
        try:
            tagged.main(self.manifest_data,self.data1,self.data2)
            logger.info('Finished tagging reads.')

        except Exception as e:
            logger.error('Error aligning')
            logger.error(traceback.format_exc())
            quit()
  
    def consolidate(self, min_freq=CONSOLIDATE_MIN_FREQ, min_qual=CONSOLIDATE_MIN_QUAL):
        logger.info('Consolidating reads...')

        try:
            self.consolidated = {}

            for sample in self.samples:
                if sample != "control":
                
                    self.consolidated[sample] = {}
                    self.consolidated[sample]['read1'] = os.path.join(self.output_folder, 'consolidated',
                                                                    sample + '.r1.consolidated.fastq')
                    self.consolidated[sample]['read2'] = os.path.join(self.output_folder, 'consolidated',
                                                                    sample + '.r2.consolidated.fastq')

                    consolidate.consolidate(self.umitagged[sample]['read1'], self.consolidated[sample]['read1'], min_qual,
                                            min_freq)
                    consolidate.consolidate(self.umitagged[sample]['read2'], self.consolidated[sample]['read2'], min_qual,
                                        min_freq)

            logger.info('Successfully consolidated reads.')
        except Exception as e:
            logger.error('Error umitagging')
            logger.error(traceback.format_exc())
            quit()

    def alignReads(self):
        logger.info('Aligning reads...')

        try:
            self.aligned = {}
            for sample in self.samples:
                if sample != "control":
                    sample_alignment_path = os.path.join(self.output_folder, 'aligned', sample + '.sam')
                    alignReads(self.BWA_path,
                             self.reference_genome,
                             self.consolidated[sample]['read1'],
                             self.consolidated[sample]['read2'],
                             sample_alignment_path)
                    self.aligned[sample] = sample_alignment_path
                    logger.info('Finished aligning reads to genome.')

        except Exception as e:
            logger.error('Error aligning')
            logger.error(traceback.format_exc())
            quit()

    def identifyOfftargetSites(self):
        logger.info('Identifying offtarget sites...')

        try:
            self.identified = {}
            self.iden_consolidated = {}
            # Identify offtarget sites for each sample
            for sample in self.samples:
                if sample != "control":
                    # Prepare sample annotations
                    sample_data = self.samples[sample]
                    annotations = {}
                    annotations['Description'] = sample_data['description']
                    annotations['Targetsite'] = sample

                    if sample is 'control':
                        annotations['Sequence'] = ''
                    else:
                        annotations['Sequence'] = sample_data['target']

                    samfile = self.aligned[sample]

                    self.identified[sample] = os.path.join(self.output_folder, 'identified',
                                                        sample + '_identifiedOfftargets.txt')
                    self.iden_consolidated[sample] = os.path.join(self.output_folder, 'identified',
                                                        sample + '_consolidated.txt')
                    outimg = os.path.join(self.output_folder, 'identified',
                                                        sample + '_wordcloud.png')

                    identifyOfftargetSites.analyze(samfile, self.reference_genome, self.identified[sample], annotations,self.iden_consolidated[sample],outimg)

            logger.info('Finished identifying offtarget sites.')

        except Exception as e:
            logger.error('Error identifying offtarget sites.')
            logger.error(traceback.format_exc())
            quit()


    def checkTheSameSites(self):
        '''
        假设有8和BC
        先从每个consolidated获取到前15位置的整合基因。每个BC做一个excel，15*n行9列，第一列是（基因+位点），第二列是这个BC的，其余七列是其他的。
        前提是位点要一样，从identifiedofftargets.txt里取位点（取成字典，key是id+pos，value是num）
        如果同一个基因有两个位点，就都要考虑。

        同时记录每个BC的整合中0.5%的位点，拿出来单独记录。
        :return:
        '''

        for sample in self.samples:
            if sample != "control":
                # identified.txt: self.identified[sample]
                # consolidated: self.iden_consolidated[sample]
                statisticalTopSites.topSites(self.samples,sample,self.identified[sample],self.iden_consolidated[sample],self.output_folder)
        logger.info('Finished statistical top sites.')


def parse_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(description='Individual Step Commands',
                                       help='Use this to run individual steps of the pipeline',
                                       dest='command')

    all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
    all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)

    return parser.parse_args()


def main():
    args = parse_args()

    if args.command == 'all':
        g = GuideSeq()
        g.parseManifest(args.manifest)#加载yaml参数
        
        g.dataTagged()
        g.consolidate()#保留这个函数只是为了保留原来文件路径，没实际作用
        g.alignReads()#基因组比对，生成sam文件
        
        g.identifyOfftargetSites()#满病毒整合位点识别，也就是LTR整合位置分析，从sam文件中计算出整合位置
        g.checkTheSameSites()
    
    else:
        logger.error('Program parameter error.')


if __name__ == '__main__':
    main()
