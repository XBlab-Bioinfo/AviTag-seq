# -*- coding:utf-8 -*-
# @Time    :2023/8/18 17:28
# @Author  :ZZK
# @ File   :statisticalTopSites.py
# Description:
import os
from xml.sax import ContentHandler
import xlwt


def topSites(samples, currentSample, file_identified, file_idenConsolidated, output_folder):
    '''

        假设有8和BC
        先从每个consolidated获取到前15位置的整合基因。每个BC做一个excel，15*n行9列，第一列是（基因+位点），第二列是这个BC的，其余七列是其他的。
        前提是位点要一样，从identifiedofftargets.txt里取位点（取成字典，key是id+pos，value是num）
        如果同一个基因有两个位点，就都要考虑。

        同时记录每个BC的整合中0.5%的位点，拿出来单独记录。
    :param file_identified: 所有位点txt
    :param file_idenConsolidated: 合并相同基因的位点txt
    :return:
    '''
    # sortSamples = sorted([sp.split("-")[1] for sp in samples])
    print(file_idenConsolidated)
    workbook = xlwt.Workbook(encoding='utf-8')
    sheet = workbook.add_sheet('Sheet1')
    head_n = 2
    for _, data_ in enumerate(samples):  # write table head
        if data_ != "control":
            if data_ == currentSample:
                sheet.write(0, 1, data_.split("-")[1])
            else:
                sheet.write(0, head_n, data_.split("-")[1])
                head_n += 1
    topFifteenGene = []  # key:gene name   value:[siteNum,total abundance]
    # paddingGeneSites = []  # gene name+
    tableCurrentSample = []  # key:gene name-position  value:siteNum
    tableCurrSiteNum = []
    cur_conAll, contentKeyName = open(file_idenConsolidated, 'r').readlines()[1:], []
    cur_identifidedAll = open(file_identified, 'r').readlines()[1:]

    contentKeySitesNum, contentKeyabundance = [], []
    # contIdenKey = {}  #key:gene name-position  value:siteNum

    for content in cur_conAll:
        contentKeyName.append(content.split("\t")[3])
        contentKeySitesNum.append(content.split("\t")[4])
        contentKeyabundance.append(int(content.split("\t")[5]))

    totalAbundance = sum(contentKeyabundance)  # 总丰度
    # print("11111",contentKeyabundance)
    # print(max(contentKeyabundance),contentKeyabundance.index(max(contentKeyabundance)))
    for _ in range(15):
        max_index = contentKeyabundance.index(max(contentKeyabundance))
        # print(contentKeyName[max_index])
        topFifteenGene.append(contentKeyName[max_index])
        contentKeyName.pop(max_index)
        contentKeySitesNum.pop(max_index)
        contentKeyabundance.pop(max_index)
    # print(topFifteenGene)
    cont1,cont2,cont3 = [],[],[]
    for content_ in cur_identifidedAll:  # 从top15的基因中取具体的基因+位点
        cont1.append(content_.split("\t")[4])
        cont2.append(content_.split("\t")[9])
        cont3.append(content_.split("\t")[19])

    for gene in topFifteenGene:
        while True:
            try:
                gene_index = cont1.index(gene)
            except:
                break
            tableCurrentSample.append(cont1[gene_index]+"-"+cont2[gene_index])
            tableCurrSiteNum.append(int(cont3[gene_index]))
            cont1.pop(gene_index)
            cont2.pop(gene_index)
            cont3.pop(gene_index)
    # print(tableCurrentSample.keys())

    i = 1
    print()
    for key, value in zip(tableCurrentSample,tableCurrSiteNum):  # 一次写入一行
        totalGeneSite = float(value)
        row_data = [key,float(value)]
                    # '{}/{}={}%'.format(value, totalAbundance, 100 * round(float(value) / float(totalAbundance), 4))]
        for sample in samples:
            if sample != "control" and sample != currentSample:
                sampleFileIdentified = os.path.join(output_folder, 'identified',
                                                    sample + '_identifiedOfftargets.txt')  # 获取：总丰度、name、pos、sitesNum

                sampleFileIdenConsolidated = os.path.join(output_folder, 'identified',
                                                          sample + '_consolidated.txt')  # 为了获取总丰度
                thiscontentKeyabundance = []
                for cont in open(sampleFileIdenConsolidated, 'r').readlines()[1:]:  # 取这个sample的总丰度
                    thiscontentKeyabundance.append(int(cont.split("\t")[5]))
                thisSampleAllAbundance = sum(thiscontentKeyabundance)

                thisidentifidedAll = open(sampleFileIdentified, 'r').readlines()[1:]
                flag_None = False
                for content_ in thisidentifidedAll:  #
                    if key == content_.split("\t")[4] + "-" + content_.split("\t")[9]:
                        totalGeneSite += float(content_.split("\t")[19])
                        flag_None = True
                        row_data.append(float(content_.split("\t")[19]))
                        # ('{}/{}={}%'.format(int(content_.split("\t")[19]), thisSampleAllAbundance,
                        #                                100 * round(float(content_.split("\t")[19]) / float(
                        #                                    thisSampleAllAbundance), 4)))
                if not flag_None:
                    row_data.append("0.00%")
        print(row_data)
        for col, data in enumerate(row_data):
            try:
                data.isalpha()
                sheet.write(i, col, data)
            except:
                sheet.write(i, col, '{}/{}={}%'.format(data, totalGeneSite, 100 * round(data / totalGeneSite, 4)))
        i += 1
    workbook.save(os.path.join(output_folder, 'identified', currentSample + '_topFifteen.xls'))
    return True
