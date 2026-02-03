# -*- coding: utf-8 -*-
# 9.23更新去掉 most_frequent_position - 25, most_frequent_position + 25误差外的数据
# 9.23更新ontarget放上边
# 2023-10-25添加自动获取fq.gz功能
# 2023-10-25添加findoff子命令
# 2024-1-26添加统计每个BC所有reads,有用reads,on or off target reads数
# 2024-1-26添加了新的引物
# 2024-4-15修改添加统计每个BC所有reads的bug
# 2024-4-30允许错错配从7变6
import os
import sys
import yaml
import argparse
import traceback
import subprocess


import log

logger = log.createCustomLogger('root')

from alignReads import alignReads
from visualization import visualizeOfftargets
import identifyOfftargetSites
import tagged


class OliTagSeq:

    def __init__(self):
        pass

    def parseManifest(self, manifest_path):
        logger.info('Loading manifest...')
        self.manifest_data = yaml.safe_load(open(manifest_path, 'r'))
        try:
            self.BWA_path = self.manifest_data['bwa']
            self.bedtools = self.manifest_data['bedtools']
            self.reference_genome = self.manifest_data['reference_genome']
            self.output_folder = self.manifest_data['output_folder']
            self.samples = self.manifest_data['samples']
            
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

    def alignReads(self, flag=True):
        logger.info('Aligning reads...')

        try:
            self.aligned = {}
            for sample in self.samples:
                sample_alignment_path = os.path.join(self.output_folder, 'aligned', sample + '.sam')

                consolidated1 = os.path.join(self.output_folder, 'consolidated',
                                             sample + '.r1.consolidated.fastq')
                consolidated2 = os.path.join(self.output_folder, 'consolidated',
                                             sample + '.r2.consolidated.fastq')
                if flag:
                    alignReads(self.BWA_path,
                            self.reference_genome,
                            consolidated1,
                            consolidated2,
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
            # 用于统计每个BC里的reads
            self.statics = []
            for sample in self.samples:
                sample_data = self.samples[sample]
                annotations = {}
                annotations['Description'] = sample_data['description']
                annotations['Targetsite'] = sample

                annotations['Sequence'] = sample_data['target']

                samfile = self.aligned[sample]

                self.identified[sample] = os.path.join(self.output_folder, 'identified',
                                                       sample + '_identifiedOfftargets.txt')

                self.total_reads,self.use_reads,self.primer11,self.primer12,self.primer21,self.primer22,self.nomatch =  identifyOfftargetSites.analyze(samfile, self.reference_genome, self.identified[sample], annotations)
                
                self.statics.append({sample:[self.total_reads,self.use_reads,self.primer11,self.primer12,self.primer21,self.primer22]})

            logger.info('Finished identifying offtarget sites.')

        except Exception as e:
            logger.error('Error identifying offtarget sites.')
            logger.error(traceback.format_exc())
            quit()


    def visualize(self):
        logger.info('Visualizing off-target sites')
        try:
            self.visua_reads = []
            for sample in self.samples:
                if sample != 'control':
                    infile = self.identified[sample]
                    outfile = os.path.join(self.output_folder, 'visualization', sample + '_offtargets')

                    # dict_with_on1 = next(item for item in self.statics if sample in item)
                    # total_reads,use_reads,primer11,primer12,primer21,primer22 = dict_with_on1[sample]
                    # visualizeOfftargets(total_reads,use_reads,primer11,primer12,primer21,primer22,infile, outfile, title=sample)
                    self.off_on_target_reads = visualizeOfftargets(infile, outfile, title=sample)
                    self.visua_reads.append({sample:self.off_on_target_reads})

            # 将统计信息统一写入新文件
            f = open('reads_type.txt','w')
            f.write('BC'+'\t' + 'total_reads' + '\t' + 'use_reads' + '\t' + 'primer11' + '\t' + 'primer12' + '\t' + 'primer21' + '\t' + 'primer22' + '\t'+'off_on_target_reads')
            f.write('\n')
            for sample in self.samples:
                if sample != 'control':
                    dict_with_on1 = next(item for item in self.statics if sample in item)
                    total_reads,use_reads,primer11,primer12,primer21,primer22 = dict_with_on1[sample]
                    dict_with_on2 = next(item for item in self.visua_reads if sample in item)
                    off_on_target_reads = dict_with_on2[sample]
                    f.write(sample+'\t'+str(total_reads)+'\t'+str(use_reads)+'\t'+str(primer11)+'\t'+str(primer12)+'\t'+str(primer21)+'\t'+str(primer22)+'\t'+str(off_on_target_reads))
                    f.write('\n')

            f.close

            logger.info('Finished visualizing off-target sites')

        except Exception as e:
            logger.error('Error visualizing off-target sites.')
            logger.error(traceback.format_exc())


def parse_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(description='Individual Step Commands',
                                       help='Use this to run individual steps of the pipeline',
                                       dest='command')

    all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
    all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    all_parser.add_argument('--identifyAndFilter', action='store_true', default=False)

    part_parser = subparsers.add_parser('findoff', help='Run all steps except 数据分流')
    part_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)

    return parser.parse_args()


def main():
    args = parse_args()

    if args.command == 'all':
        g = OliTagSeq()
        g.parseManifest(args.manifest)
        g.dataTagged()
        g.alignReads()
        g.identifyOfftargetSites()
        g.visualize()

    elif args.command == 'findoff':
        g = OliTagSeq()
        g.parseManifest(args.manifest)
        g.alignReads(flag=False)
        g.identifyOfftargetSites()
        g.visualize()
    else:
        logger.error('Program parameter error.')


if __name__ == '__main__':
    main()
