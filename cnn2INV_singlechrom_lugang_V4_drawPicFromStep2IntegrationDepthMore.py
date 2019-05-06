#!/usr/bin/python
# -*- coding: UTF-8 -*-
# cd /usr/python3/bin
# python3.5 /mnt/hde/gao/wuzhongjia/deepLearning/code/cnn2INV_singlechrom_lugang_V4_drawPicFromStep2Integration.py

import numpy as np
import os
import pysam
import re
from PIL import Image, ImageDraw

def grabFeatures_ori(bamFileName, chr_id, candidates, resultDataRoot):

    fetchFeatureList = []
    depth = os.popen('samtools depth ' + bamFileName +' -r ' + chr_id + ':' + str(candidates[0][0]) + '-' + str(candidates[0][0]) + ' | less').readlines()  
    try:
        depthInt = int(depth[0].split('\t')[-1].replace('\n', ''))
    except:
        depthInt = 0
    if depthInt < 110:
        lPointFetchFeatures = fetchReads(bamFileName, chr_id, int(candidates[0][0])-500, int(candidates[0][0])+800, resultDataRoot)
        rPointFetchFeatures = fetchReads(bamFileName, chr_id, int(candidates[0][1])-500, int(candidates[0][1])+800, resultDataRoot)
        fetchFeatureList.append([lPointFetchFeatures, rPointFetchFeatures])
    i = 1
    for i in range(len(candidates)-1):
        candidateOneScore = int(candidates[i][2].split(':')[-1])
        if (candidateOneScore > 100) or (int(candidates[i+1][0]) - int(candidates[i][0]) < 1000) or (int(candidates[i][0]) - int(candidates[i-1][0]) < 1000): # 修改 2019年1月16日10:44:17
            depth = os.popen('samtools depth ' + bamFileName +' -r ' + chr_id + ':' + str(candidates[i][0]) + '-' + str(candidates[i][0]) + ' | less').readlines()  
            try:
                depthInt = int(depth[0].split('\t')[-1].replace('\n', ''))
            except:
                depthInt = 0
            if depthInt > 100:
                continue
        lPointFetchFeatures = fetchReads(bamFileName, chr_id, int(candidates[i][0])-500, int(candidates[i][0])+800, resultDataRoot)
        rPointFetchFeatures = fetchReads(bamFileName, chr_id, int(candidates[i][1])-500, int(candidates[i][1])+800, resultDataRoot)
        fetchFeatureList.append([lPointFetchFeatures, rPointFetchFeatures])
    depth = os.popen('samtools depth ' + bamFileName +' -r ' + chr_id + ':' + str(candidates[-1][0]) + '-' + str(candidates[-1][0]) + ' | less').readlines()  
    try:
        depthInt = int(depth[0].split('\t')[-1].replace('\n', ''))
    except:
        depthInt = 0
    if depthInt < 100:
        lPointFetchFeatures = fetchReads(bamFileName, chr_id, int(candidates[-1][0])-500, int(candidates[-1][0])+800, resultDataRoot)
        rPointFetchFeatures = fetchReads(bamFileName, chr_id, int(candidates[-1][1])-500, int(candidates[-1][1])+800, resultDataRoot)
        fetchFeatureList.append([lPointFetchFeatures, rPointFetchFeatures])
    return fetchFeatureList

def mkdir(path, sampleChromList):

    for sample in sampleChromList:
        resPath = path + sample + '/'
        command = 'mkdir ' + resPath
        os.system(command)
        command = 'mkdir ' + resPath + sample
        os.system(command)
        p1 = resPath + sample + '/0/'
        command = 'mkdir ' + p1
        os.system(command)
        p2 = resPath + sample + '/1/'
        command = 'mkdir ' + p2
        os.system(command)

        resPath = path + sample + '/' + 'except' + sample + '/'
        command = 'mkdir ' + resPath
        os.system(command)
        p1 = resPath + '/0/'
        command = 'mkdir ' + p1
        os.system(command)
        p2 = resPath + '/1/'
        command = 'mkdir ' + p2
        os.system(command)

        resPath = path + sample + '/' + 'lrPic' + sample + '/'
        command = 'mkdir ' + resPath
        os.system(command)
        p1 = resPath + '/l_1/'
        command = 'mkdir ' + p1
        os.system(command)
        p2 = resPath + '/r_1/'
        command = 'mkdir ' + p2
        os.system(command)

        resPath = path + sample + '/' + 'fetch' + sample + '/'
        command = 'mkdir ' + resPath
        os.system(command)


def getSampleVCFPoints(rootPath, vcfFile, sampleChromList):

    # wholeVCFFile = '/home/wuzhongjia/vcf/INV.allsamples.ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.getChromPoint'
    # print(rootPath, vcfFile, sampleChromList)
    wholeVCFPoints = np.loadtxt(vcfFile, str)
    for aPoint in wholeVCFPoints:
        if aPoint[0] in sampleChromList:
            items = '"' + str(aPoint[1]) + ' '  + str(aPoint[2]) + ' '  + str(aPoint[3].split(';END=')[1].split(';EUR_AF')[0]) + '"' 
            sample = aPoint[0]
            sampleVcfFile = rootPath + sample + '/' + sample + '_vcf'
            command = 'echo ' + items + ' >> ' + sampleVcfFile
            # print(command)
            os.system(command)

def fetchReads(bamFileName, chr_id, lPoint, rPoint, resultDataRoot):

    resFileName = resultDataRoot + bamFileName.split('/')[-1] + '_' + str(lPoint) + '_' + str(rPoint) + '_fetchFeatures_V4'
    bamfile = pysam.AlignmentFile(bamFileName, "rb")
    readFeatures = []
    for read in bamfile.fetch(chr_id, lPoint, rPoint):
        readFeatures.append([read.cigarstring, read.is_proper_pair, read.mapping_quality, read.is_reverse == read.mate_is_reverse, read.reference_start, read.query_name, read.is_reverse, read.next_reference_start, read.template_length, read.is_paired, read.query_sequence])
    readFeatures = np.array(readFeatures)
    np.savetxt(resFileName, readFeatures, newline='\n',fmt='%s')
    bamfile.close()
    return resFileName

def getMaxDepth(candidate, bamFileName, chr_id):

    lPoint = candidate[0]
    rPoint = candidate[1]
    depthInt = 0
    for i in range(2):
        depthInt = 0
        tempI = int(candidate[i])-500
        for j in range(40):
            depth = os.popen('samtools depth ' + bamFileName +' -r ' + chr_id + ':' + str(tempI) + '-' + str(tempI) + ' | less').readlines()
            try:
                depthIntTemp = int(depth[0].split('\t')[-1].replace('\n', ''))
            except:
                depthIntTemp = 0
            if(depthInt < depthIntTemp):
                depthInt = depthIntTemp
            if(depthInt > 100):
                return depthInt
            tempI += 30
    return depthInt


def grabFeatures(bamFileName, chr_id, candidates, resultDataRoot):

    fetchFeatureList = []
    depthInt = getMaxDepth(candidates[0], bamFileName, chr_id)
    if depthInt < 100:
        lPointFetchFeatures = fetchReads(bamFileName, chr_id, int(candidates[0][0])-500, int(candidates[0][0])+800, resultDataRoot)
        rPointFetchFeatures = fetchReads(bamFileName, chr_id, int(candidates[0][1])-500, int(candidates[0][1])+800, resultDataRoot)
        fetchFeatureList.append([lPointFetchFeatures, rPointFetchFeatures])
    i = 1
    for i in range(len(candidates)-1):
        if((int(candidates[i][1]) - int(candidates[i][0])) > 1000000):
            continue
        candidateOneScore = int(candidates[i][2].split(':')[-1])
        # if((candidateOneScore > 100) or (int(candidates[i+1][0]) - int(candidates[i][0]) < 1000) or (int(candidates[i][0]) - int(candidates[i-1][0]) < 1000)): # 修改 2019年1月16日10:44:17
        #     depthInt = getMaxDepth(candidates[i], bamFileName, chr_id)
        #     if depthInt > 100:
        #         continue
        depthInt = getMaxDepth(candidates[i], bamFileName, chr_id)
        if depthInt > 100:
            continue
        lPointFetchFeatures = fetchReads(bamFileName, chr_id, int(candidates[i][0])-500, int(candidates[i][0])+800, resultDataRoot)
        rPointFetchFeatures = fetchReads(bamFileName, chr_id, int(candidates[i][1])-500, int(candidates[i][1])+800, resultDataRoot)
        fetchFeatureList.append([lPointFetchFeatures, rPointFetchFeatures])

    depthInt = getMaxDepth(candidates[-1], bamFileName, chr_id)
    if depthInt < 100:
        lPointFetchFeatures = fetchReads(bamFileName, chr_id, int(candidates[-1][0])-500, int(candidates[-1][0])+800, resultDataRoot)
        rPointFetchFeatures = fetchReads(bamFileName, chr_id, int(candidates[-1][1])-500, int(candidates[-1][1])+800, resultDataRoot)
        fetchFeatureList.append([lPointFetchFeatures, rPointFetchFeatures])
    return fetchFeatureList

def drawFeatureReads(fetchReadsFileName):

    if os.path.getsize(fetchReadsFileName) > 999999: # 文件过大不做图 5999999
        return
    # print(fetchReadsFileName)
    readFeatures = np.loadtxt(fetchReadsFileName, dtype=bytes).astype(str)
    if len(readFeatures)==0:
        return
    allStartPos = int(readFeatures[0][4])
    lastReadId = -1
    while(readFeatures[lastReadId][0]=='None'):
        lastReadId -= 1
    cigarIntList = re.findall(r"[0-9]+", readFeatures[lastReadId][0])
    lastReadLen = 0
    for partLen in cigarIntList:
        lastReadLen += int(partLen)
    try:
        allEndPos = int(readFeatures[-1][4]) + lastReadLen
    except:
        return # fetchRead 只有一条

    width = allEndPos - allStartPos
    height = -1
    col = row = red = green = blue = 100
    endPos = [-1]*100
    image = Image.new ("RGB", (width,1000),(255,255,255))
    draw = ImageDraw.Draw(image)
    for read in readFeatures:

        readStartPos = int(read[4]) - allStartPos
        cigarString = read[0]
        cigarTypeList = re.findall(r"[A-Z]+", cigarString)
        cigarIntList = re.findall(r"[0-9]+", cigarString)

        if (read[3] == 'False') and (read[1] == 'True') and (len(cigarTypeList) == 1) and (int(read[2]) >= 60): # 这一步的改善可以直接在fetch candidates reads 那一步就改善
            continue # 需要注释掉这两条语句重新执行一次

        readLen = 0
        for partLen in cigarIntList:
            readLen += int(partLen)

        for i in range(len(endPos)):
            if readStartPos > endPos[i]:
                endPos[i] = readStartPos + readLen
                break
        yIdex = i
        height = max(height, yIdex)
        r = 0
        g = 0
        b = 0     
        if len(cigarIntList) > 0:
            if read[3] == 'True': # 同向
                r += 120
            else:
                r += 50
            if read[1] == 'False': # ISPE异常
                r += 120
            else:
                r += 50
            if int(read[2]) < 60: # 质量异常
                g += 50
            if int(read[2]) < 30: # 质量异常
                g += 50
            if int(read[2]) < 20: # 质量异常
                g += 50
            for i in range(len(cigarTypeList)):
                
                if cigarTypeList[i] == 'S': # soft-clip
                    b = 200
                elif cigarTypeList[i] == 'M':
                    b = 50
                else:
                    b = 100
                draw.rectangle((readStartPos, yIdex*10, readStartPos + int(cigarIntList[i]), yIdex*10+10), (r, g, b), (r, g, b))
                readStartPos += int(cigarIntList[i])
        else: # cigar == 'None' (single read)
            g = 255
            r = b = 0
            draw.rectangle((readStartPos, yIdex*10, readStartPos + 100, yIdex*10+10), (r, g, b), (r, g, b))
            endPos[i] = readStartPos + 100
    image = image.crop((0, 0, width, height*10+10))
    return image

def mergeLRPointImage(lPointPNG, rPointPNG, resFileName):

    width = max(lPointPNG.size[0], rPointPNG.size[0])
    height = lPointPNG.size[1] + rPointPNG.size[1]
    mergePNG = Image.new('RGB', (width, height), (255,255,255))
    box = (0, 0, lPointPNG.size[0], lPointPNG.size[1]) 
    region = lPointPNG.crop(box)  
    mergePNG.paste(region, box)
    box = (0, 0, rPointPNG.size[0], rPointPNG.size[1])
    region = rPointPNG.crop(box)
    box = (0, lPointPNG.size[1], rPointPNG.size[0], lPointPNG.size[1]+rPointPNG.size[1])
    mergePNG.paste(region, box) 
    mergePNG.save(resFileName, "PNG")

def drawSamplesPic(sampleChromList, rootPath, bamFilePath, candidatePath, lrPicPath, candidatesThreScore): # 提炼出 lrPicPath 作为一个参数

    sampleCount = 0
    chr_idCount = 0
    for sampleName in sampleChromList:
        print(sampleName, 'start')
        sampleCount += 1
        for chr_id in sampleChromList.get(sampleName):
            chr_idCount += 1
            vcfSamplePoints = np.loadtxt(rootPath + sampleName + '/' + sampleName + '_vcf', int)
            vcfPoints = []
            for vcf in vcfSamplePoints:
                if chr_id == vcf[0]:
                    vcfPoints.append(vcf)
            chr_id = str(chr_id)
            bamFileName = bamFilePath + sampleName + '/' + sampleName + '.' + chr_id + '.high.bam'
            invCandidatesFileName = candidatePath + sampleName + '.' + chr_id + '.high.bam_candidates_V4_treScore' + str(candidatesThreScore)
            candidates = np.loadtxt(invCandidatesFileName, str).reshape(-1,3)
            fetchRoot = rootPath + sampleName + '/' + 'fetch' + sampleName + '/'
            fetchFeatureList = grabFeatures(bamFileName, chr_id, candidates, fetchRoot)
            print(sampleName, chr_id, '提取特征完成')
            lrPicRoot = rootPath + sampleName + '/' + 'lrPic' + sampleName + '/'
            for fetchFeature in fetchFeatureList:
                lPointPNG = drawFeatureReads(fetchFeature[0])
                rPointPNG = drawFeatureReads(fetchFeature[1])
                if lPointPNG == None or rPointPNG == None:
                    continue
                bothBool = 0
                for vcf in vcfPoints:
                    if chr_id == str(vcf[0]):
                        if (int(fetchFeature[0].split('.bam_')[1].split('_')[0]) < vcf[1]) and (int(fetchFeature[0].split('.bam_')[1].split('_')[1]) > vcf[1]):
                            bothBool = 2
                            lPointPNGName = lrPicRoot + 'l_1/' + sampleName + '.' + chr_id + '.high.bam_' + fetchFeature[0].split('.bam_')[1].split('_')[0] + '-' + fetchFeature[0].split('.bam_')[1].split('_')[1] + '_' + chr_id + '_lPoint_filterNormal_V4__1.png'
                            lPointPNG.save(lPointPNGName)
                            os.system('cp ' + lPointPNGName + ' ' + lrPicPath + 'l_1' + lPointPNGName.split('l_1')[1])
                            rPointPNGName = lrPicRoot + 'r_1/' + sampleName + '.' + chr_id + '.high.bam_' + fetchFeature[0].split('.bam_')[1].split('_')[0] + '-' + fetchFeature[0].split('.bam_')[1].split('_')[1] + '_' + chr_id + '_rPoint_filterNormal_V4__1.png' # 这个名字取的有问题，取成左点的名字了
                            rPointPNG.save(rPointPNGName)
                            os.system('cp ' + rPointPNGName + ' ' + lrPicPath + 'r_1' + rPointPNGName.split('r_1')[1])
                            bothPointPNGName = lrPicRoot.split('lrPic')[0] + sampleName + '/1/' + sampleName + '.' + chr_id + '.high.bam_' + fetchFeature[0].split('.bam_')[1].split('_')[0] + '-' + fetchFeature[0].split('.bam_')[1].split('_')[1] + '_' + fetchFeature[1].split('.bam_')[1].split('_')[0] +'-' +  fetchFeature[1].split('.bam_')[1].split('_')[1] + '_' + chr_id + '_filterNormal_V4__1.png'
                            mergeLRPointImage(lPointPNG, rPointPNG, bothPointPNGName)
                            break
                        if (int(fetchFeature[1].split('.bam_')[1].split('_')[0]) < vcf[2]) and (int(fetchFeature[1].split('.bam_')[1].split('_')[1]) > vcf[2]):
                            bothBool = 3
                            lPointPNGName = lrPicRoot + 'l_1/' + sampleName + '.' + chr_id + '.high.bam_' + fetchFeature[0].split('.bam_')[1].split('_')[0] + '-' + fetchFeature[0].split('.bam_')[1].split('_')[1] + '_' + chr_id + '_lPoint_filterNormal_V4__1.png'
                            lPointPNG.save(lPointPNGName)
                            os.system('cp ' + lPointPNGName + ' ' + lrPicPath + 'l_1' + lPointPNGName.split('l_1')[1])
                            rPointPNGName = lrPicRoot + 'r_1/' + sampleName + '.' + chr_id + '.high.bam_' + fetchFeature[0].split('.bam_')[1].split('_')[0] + '-' + fetchFeature[0].split('.bam_')[1].split('_')[1] + '_' + chr_id + '_rPoint_filterNormal_V4__1.png'
                            rPointPNG.save(rPointPNGName)
                            os.system('cp ' + rPointPNGName + ' ' + lrPicPath + 'r_1' + rPointPNGName.split('r_1')[1])
                            bothPointPNGName = lrPicRoot.split('lrPic')[0] + sampleName + '/1/' + sampleName + '.' + chr_id + '.high.bam_' + fetchFeature[0].split('.bam_')[1].split('_')[0] + '-' + fetchFeature[0].split('.bam_')[1].split('_')[1] + '_' + fetchFeature[1].split('.bam_')[1].split('_')[0] +'-' +  fetchFeature[1].split('.bam_')[1].split('_')[1] + '_' + chr_id + '_filterNormal_V4__1.png'
                            mergeLRPointImage(lPointPNG, rPointPNG, bothPointPNGName)
                            break

                if bothBool == 0:
                    bothPointPNGName = lrPicRoot.split('lrPic')[0] + sampleName + '/0/' + sampleName + '.' + chr_id + '.high.bam_' + fetchFeature[0].split('.bam_')[1].split('_')[0] + '-' + fetchFeature[0].split('.bam_')[1].split('_')[1] + '_' + fetchFeature[1].split('.bam_')[1].split('_')[0] +'-' +  fetchFeature[1].split('.bam_')[1].split('_')[1] + '_' + chr_id + '_filterNormal_V4__0.png'
                    mergeLRPointImage(lPointPNG, rPointPNG, bothPointPNGName)

                else:
                    bothBool = 0

            print(sampleCount, sampleName, chr_id, len(vcfPoints), 'end')
    print('染色体总数：', chr_idCount)

if __name__ == '__main__':

    print('cnn2INV_singlechrom_lugang_V4_drawPicFromStep2Integration.py start')

    # sampleChromList = {'NA18525':[1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 14, 16, 21],\
    # 'NA19017':[1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 14, 15, 16, 18, 21], 'NA19238':[1, 2, 3, 4, 6, 7, 9, 10, 11, 12, 14, 16, 18, 21, 22],\
    # 'NA19239':[1, 2, 3, 4, 6, 7, 9, 10, 11, 12, 14, 16, 17, 18, 21], 'NA19625':[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20, 21],\
    # 'NA19648':[1, 2, 3, 5, 7, 9, 10, 11, 12, 16, 18, 21], 'NA20502':[1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 14, 16, 18, 21, 22],\
    # 'NA20845':[1, 2, 3, 5, 6, 7, 9, 10, 11, 12, 14, 16, 18, 20, 21],\
    # 'HG00096':[1, 2, 4, 5, 7, 9, 10, 12, 16, 18, 21], 'HG00268':[1, 2, 4, 5, 7, 8, 9, 10, 12, 14, 16, 17, 20, 21],\
    # 'HG00419':[1, 2, 5, 7, 8, 9, 10, 12, 14, 16, 18, 20, 21], 'HG00759':[1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 14, 16, 18, 21],\
    # 'HG01051':[1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 17, 18, 21], 'HG01112':[1, 2, 3, 5, 7, 10, 11, 12, 14, 16, 17, 21],\
    # 'HG01500':[1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 14, 16, 20, 21], 'HG01565':[1, 2, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 21],\
    # 'HG01583':[1, 2, 5, 7, 9, 10, 11, 12, 14, 16, 21], 'HG02568':[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 17, 18, 20, 21],\
    # 'HG03642':[1, 2, 4, 5, 7, 9, 10, 12, 14, 16, 18, 20, 21], 'HG03742':[1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 21]} # 285
    
    sampleChromList = {'HG01112':[16], 'HG01583':[16]}

    candidatesThreScore = 5
    rootPath = '/mnt/xiaoju/cnnResult/'
    bamFilePath = '/mnt/xiaohei/highBam/'
    candidatePath = '/mnt/xiaolan/wzj/realBamResult/step2_Integration/step2_V4_threScore' + str(candidatesThreScore) + '_merge/' #'/mnt/hdb/wzj/realBamResult/20samplesCandites/'
    lrPicPath = rootPath + 'lrPic/'

    # os.system('mkdir ' + lrPicPath)
    # os.system('mkdir ' + lrPicPath + 'l_1')
    # os.system('mkdir ' + lrPicPath + 'r_1')
    # os.system('mkdir ' + lrPicPath + 'rangeMergeLR')
    # mkdir(rootPath, sampleChromList)
    # print('mkdir end')

    # vcfFile = '/mnt/xiaolan/wzj/realBamResult/INV.allsamples.ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.getChromPoint'
    # getSampleVCFPoints(rootPath, vcfFile, sampleChromList) # sampleVcfFile = rootPath + sample + '/' + sample + '_vcf'
    # print('getSampleVCFPoints end')

    drawSamplesPic(sampleChromList, rootPath, bamFilePath, candidatePath, lrPicPath, candidatesThreScore)
    # # print('generate PNG end')

    print('cnn2INV_singlechrom_lugang_V4_drawPicFromStep2Integration.py end')
