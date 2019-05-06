#!/usr/bin/python
# -*- coding: UTF-8 -*-
# cd /usr/python3/bin
# python3.5 /mnt/hde/gao/wuzhongjia/deepLearning/code/cnn2INV_singlechrom_lugang_V4_step2_Integration.py
import numpy as np
import pysam
from PIL import Image, ImageDraw
import re
import time
import multiprocessing
import os

# get all inversion candidates' points from same orientation candidates' reads
def getInvCandidates(step1Dir, step2Dir, sameOrienFileName, candidatesThreScore):

    oriPoints = np.loadtxt(step1Dir + sameOrienFileName, str)
    if len(oriPoints)==0:
        return
    resList = []

    candidateLPoint = int(oriPoints[0][0][1:-1])
    candidateRPoint = int(oriPoints[0][1][0:-1])

    oriPoints = oriPoints[1:]
    score = 0
    for point in oriPoints:
        lPoint = int(point[0][1:-1])
        rPoint = int(point[1][0:-1])
        if lPoint>=candidateLPoint and lPoint<=candidateRPoint:
            candidateRPoint = rPoint
            score += 1
        elif score>candidatesThreScore:
            resList.append([candidateLPoint, candidateRPoint, score+1])
            score = 0
        else:
            candidateLPoint = lPoint
            candidateRPoint = rPoint
            score = 0
            
    if score>candidatesThreScore:
        resList.append([candidateLPoint, candidateRPoint, score+1])
        score = 0

    invCanditesFileName = step2Dir + sameOrienFileName.split('soft_sameOrien')[0]+'candidates_V4_treScore' + str(candidatesThreScore) # candidatePath + sampleName + '.' + chr_id + '.high.bam_candidates_V3' HG01500.1.high.bam_candidates_V3
    resList = np.array(resList)
    np.savetxt(invCanditesFileName, resList, newline='\n',fmt='%s')
    print('inversion candites:', invCanditesFileName)
    return invCanditesFileName

# 把ispe对于比上一条的read的ispe大于2000的candidate拆分成两个candidate(为了断点更精确)
def getInvCanNear(step1Dir, step2Dir, sameOrienFileName, candidatesThreScore):

    oriPoints = np.loadtxt(step1Dir + sameOrienFileName, str)
    if len(oriPoints)==0:
        return
    resList = []

    candidateLPoint = int(oriPoints[0][0][1:-1])
    candidateRPoint = int(oriPoints[0][1][0:-1])
    candidateISPE = int(oriPoints[0][5][0:-1])

    oriPoints = oriPoints[1:]
    score = 0
    for point in oriPoints:
        lPoint = int(point[0][1:-1])
        rPoint = int(point[1][0:-1])
        ispePoint = int(point[5][0:-1])
        absIspe = abs(candidateISPE - ispePoint)
        if(absIspe>2000):
            resList.append([candidateLPoint, candidateRPoint, score+1])
            score = 0
            candidateLPoint = lPoint
            candidateISPE = ispePoint
            continue
        if lPoint>=candidateLPoint and lPoint<=candidateRPoint:
            candidateRPoint = rPoint
            score += 1
        elif score>candidatesThreScore:
            resList.append([candidateLPoint, candidateRPoint, score+1])
            score = 0
            candidateISPE = ispePoint
        else:
            candidateLPoint = lPoint
            candidateRPoint = rPoint
            score = 0
            candidateISPE = ispePoint
            
    if score>candidatesThreScore:
        resList.append([candidateLPoint, candidateRPoint, score+1])
        score = 0

    invCanditesFileName = step2Dir + sameOrienFileName.split('soft_sameOrien')[0]+'candidates_V4_treScore' + str(candidatesThreScore) # candidatePath + sampleName + '.' + chr_id + '.high.bam_candidates_V3' HG01500.1.high.bam_candidates_V3
    resList = np.array(resList)
    np.savetxt(invCanditesFileName, resList, newline='\n',fmt='%s')
    print('inversion candites:', invCanditesFileName)
    return invCanditesFileName

def miFilReviseOriens(oriSameOrienDir, reviseOrienDir):
    for oriSameOrienFileNamePart in os.listdir(oriSameOrienDir):
        print(oriSameOrienFileNamePart)
        oriSameOrienFileName = oriSameOrienDir + oriSameOrienFileNamePart
        reviseOrienFileName = reviseOrienDir + oriSameOrienFileNamePart
        oriSameOrienList = np.loadtxt(oriSameOrienFileName, str)
        reviseOrienList = []
        for line in oriSameOrienList:
            ispeMi = int(line[-1][0:-1])
            if(ispeMi<6):
                reviseOrienList.append(line)
        reviseOrienList = np.array(reviseOrienList)
        np.savetxt(reviseOrienFileName, reviseOrienList, newline='\n',fmt='%s')
        reviseOrienList = []

def getMaxMi(oriKeySameOrienList):
    keyDict = {}
    for line in oriKeySameOrienList:
        tempKey = int(line[-1].replace(']', ''))
        for key in range(tempKey-1, tempKey + 2):
            if key in keyDict:
                keyDict[key] = keyDict[key] + 1
            else:
                keyDict[key] = 1
        keyDict[tempKey] = keyDict[tempKey] + 1
    return max(keyDict, key=keyDict.get)

def miReviseOriens(oriSameOrienDir, oriCandidateDir, reviseOrienDir, candidatesThreScore):
    print(oriSameOrienDir, oriCandidateDir, reviseOrienDir)
    for oriSameOrienFileNamePart in os.listdir(oriSameOrienDir):
        print(oriSameOrienFileNamePart)
        oriSameOrienFileName = oriSameOrienDir + oriSameOrienFileNamePart
        oriCandidateFileName = oriCandidateDir + oriSameOrienFileNamePart.split('_')[0]+'_candidates_V4_treScore' + str(candidatesThreScore)
        # reviseOrienFileName = reviseOrienDir + oriSameOrienFileNamePart.split('_')[0]+'_soft_sameOrien_1_V4_revise'
        reviseOrienFileName = reviseOrienDir + oriSameOrienFileNamePart
        oriSameOrienList = np.loadtxt(oriSameOrienFileName, str)
        oriCandidateList = np.loadtxt(oriCandidateFileName, int)
        reviseOrienList = []
        curOrienLoc = 0
        curCandiLoc = 0
        oriSameOrienListSize = len(oriSameOrienList)
        oriCandidateListSize = len(oriCandidateList)
        while((curOrienLoc < oriSameOrienListSize) and (curCandiLoc < oriCandidateListSize)):
            if(oriSameOrienList[curOrienLoc][0][1:-1] == str(oriCandidateList[curCandiLoc][0])):
                # print(oriSameOrienFileNamePart, curOrienLoc, oriSameOrienListSize, oriCandidateList[curCandiLoc])
                oriKeySameOrienList = oriSameOrienList[curOrienLoc:curOrienLoc+oriCandidateList[curCandiLoc][2]]
                maxMi = getMaxMi(oriKeySameOrienList)
                # print(maxMi)
                curOrienLoc += oriCandidateList[curCandiLoc][2]
                curCandiLoc += 1
                for line in oriKeySameOrienList:
                    tempKey = int(line[-1].replace(']', ''))
                    if((tempKey >= (maxMi-1)) and (tempKey <= (maxMi+1))):
                        reviseOrienList.append(line)
            else:
                curOrienLoc += 1
        
        # print(reviseOrienList)
        reviseOrienList = np.array(reviseOrienList)
        np.savetxt(reviseOrienFileName, reviseOrienList, newline='\n',fmt='%s')
        reviseOrienList = []

def getVariance(keySameOrienList):

    oriKeySameOrienNP = np.array(keySameOrienList)
    sum1=oriKeySameOrienNP.sum()
    narray2=oriKeySameOrienNP*oriKeySameOrienNP
    sum2=narray2.sum()
    mean=sum1/len(oriKeySameOrienNP)
    variance=sum2/len(oriKeySameOrienNP)-mean**2
    # print('oriKeySameOrienNP', oriKeySameOrienNP)
    return variance

def getMaxVariance(oriKeySameOrienList):

    keyDict = {}
    for line in oriKeySameOrienList:
        tempKey = int(line[-1].replace(']', ''))
        for key in range(tempKey-1, tempKey + 2):
            if key in keyDict:
                keyDict[key].append(int(line[5].replace(',', '')))
            else:
                keyDict[key] = []
                keyDict[key].append(int(line[5].replace(',', '')))

    oriVarianceDict = {}
    varianceDict = {}
    for key in keyDict:
        value = getVariance(keyDict[key])
        oriVarianceDict[key] = value
        if(value != 0):
            varianceDict[key] = value

    if(len(varianceDict) == 0): # 当方差为0时，说明每组数据的长度都是1.或每组数据的组内元素完全一样
        minKey = 9999999
        for key in oriVarianceDict:
            if(key < minKey):
                minKey = key
        minKeyValue = oriVarianceDict[minKey]
    else:
        minKey = min(varianceDict, key=varianceDict.get) # 也有可能不一定是方差小的是最凝聚的，会不会因为数量多导致方差变大？
        minKeyValue = varianceDict[minKey]
    minCount = 0
    for key in varianceDict:
        if(varianceDict[key] == minKeyValue):
            minCount += 1

    # print('minKey', minKey, 'minKeyValue:', minKeyValue)
    if(minCount == 1):
        # print(minKey)
        return minKey
    else:
        # print(minKey+1)
        return (minKey+1)

def varianceReviseOriens(oriSameOrienDir, oriCandidateDir, reviseOrienDir, candidatesThreScore):
    print(oriSameOrienDir, oriCandidateDir, reviseOrienDir)

    for oriSameOrienFileNamePart in os.listdir(oriSameOrienDir):
        print(oriSameOrienFileNamePart)
        oriSameOrienFileName = oriSameOrienDir + oriSameOrienFileNamePart
        oriCandidateFileName = oriCandidateDir + oriSameOrienFileNamePart.split('_')[0]+'_candidates_V4_treScore' + str(candidatesThreScore)
        # reviseOrienFileName = reviseOrienDir + oriSameOrienFileNamePart.split('_')[0]+'_soft_sameOrien_1_V4_variance_revise'
        reviseOrienFileName = reviseOrienDir + oriSameOrienFileNamePart
        oriSameOrienList = np.loadtxt(oriSameOrienFileName, str)
        oriCandidateList = np.loadtxt(oriCandidateFileName, int)
        reviseOrienList = []
        curOrienLoc = 0
        curCandiLoc = 0
        oriSameOrienListSize = len(oriSameOrienList)
        oriCandidateListSize = len(oriCandidateList)
        while((curOrienLoc < oriSameOrienListSize) and (curCandiLoc < oriCandidateListSize)):
            if(oriSameOrienList[curOrienLoc][0][1:-1] == str(oriCandidateList[curCandiLoc][0])):
                oriKeySameOrienList = oriSameOrienList[curOrienLoc:curOrienLoc+oriCandidateList[curCandiLoc][2]]
                minVarianceMi = getMaxVariance(oriKeySameOrienList)
                curOrienLoc += oriCandidateList[curCandiLoc][2]
                curCandiLoc += 1
                for line in oriKeySameOrienList:
                    tempKey = int(line[-1].replace(']', ''))
                    if((tempKey >= (minVarianceMi-1)) and (tempKey <= (minVarianceMi+1))):
                        reviseOrienList.append(line)
            else:
                curOrienLoc += 1
        
        reviseOrienList = np.array(reviseOrienList)
        np.savetxt(reviseOrienFileName, reviseOrienList, newline='\n',fmt='%s')
        print(reviseOrienFileName, 'finished')
        reviseOrienList = []

def getReadDistance(oriKeySameOrienList):
    count = 0
    posDict = {}
    mposDict = {}
    for read in oriKeySameOrienList:
        count += 1
        distancePos = 0
        distanceMpos = 0
        for read1 in oriKeySameOrienList:
            if((read[0] != read1[0]) or (read[1] != read1[1])):
                distancePos += abs(int(read[0].replace('[', '').replace(',', '')) - int(read1[0].replace('[', '').replace(',', '')))
                distanceMpos += abs(int(read[1].replace(',', '')) - int(read1[1].replace(',', '')))
                # print(read1[0], read1[1])
        # print(count, '****', read, distancePos, distanceMpos)
        if(read[0] in posDict):
            if(distancePos < posDict[read[0]]):
                posDict[read[0]] = distancePos
        else:
            posDict[read[0]] = distancePos
        if(read[1] in mposDict):
            if(distanceMpos < mposDict[read[1]]):
                mposDict[read[1]] = distanceMpos
        else:
            mposDict[read[1]] = distanceMpos
    minPos = min(posDict, key=posDict.get)[1:-1]
    minMPos = min(mposDict, key=mposDict.get)[:-1]
    # print('min############################', minPos, minMPos)
    return minPos, minMPos


def getDistance(oriSameOrienDir, oriCandidateDir, reviseCandidateDir, candidatesThreScore):
    print(oriSameOrienDir, oriCandidateDir, reviseCandidateDir)
    for oriSameOrienFileNamePart in os.listdir(oriSameOrienDir):
        print(oriSameOrienFileNamePart)
        oriSameOrienFileName = oriSameOrienDir + oriSameOrienFileNamePart
        oriCandidateFileName = oriCandidateDir + oriSameOrienFileNamePart.split('_')[0]+'_candidates_V4_treScore' + str(candidatesThreScore)
        # reviseOrienFileName = reviseCandidateDir + oriSameOrienFileNamePart.split('_')[0]+'_soft_sameOrien_1_V4_revise'
        reviseOrienFileName = reviseCandidateDir + oriSameOrienFileNamePart.split('_')[0]+'_candidates_V4_treScore' + str(candidatesThreScore)
        oriSameOrienList = np.loadtxt(oriSameOrienFileName, str)
        oriCandidateList = np.loadtxt(oriCandidateFileName, int)
        reviseOrienList = []
        curOrienLoc = 0
        curCandiLoc = 0
        oriSameOrienListSize = len(oriSameOrienList)
        oriCandidateListSize = len(oriCandidateList)
        while((curOrienLoc < oriSameOrienListSize) and (curCandiLoc < oriCandidateListSize)):
            # 注意：以下这句(oriCandidateList[curCandiLoc][2] < 500)是过滤得分高于500的情况（如HG00096.1.[46385800 46419189 15432]，这个需要另外慎重考虑会不会miss一些
            if((oriCandidateList[curCandiLoc][2] < 500) and (oriSameOrienList[curOrienLoc][0][1:-1] == str(oriCandidateList[curCandiLoc][0]))):
                oriKeySameOrienList = oriSameOrienList[curOrienLoc:curOrienLoc+oriCandidateList[curCandiLoc][2]] # 由candidate获得了一组read
                minPos, minMPos = getReadDistance(oriKeySameOrienList)
                # print(minPos, minMPos, oriCandidateList[curCandiLoc][2])
                reviseOrienList.append([minPos, minMPos, oriCandidateList[curCandiLoc][2]])
                curOrienLoc += oriCandidateList[curCandiLoc][2]
                curCandiLoc += 1
            else:
                curOrienLoc += 1
        
        reviseOrienList = np.array(reviseOrienList)
        np.savetxt(reviseOrienFileName, reviseOrienList, newline='\n',fmt='%s')
        reviseOrienList = []

def getReadScore(oriKeySameOrienList, ispeTime, meanISPE):

    clearKeySameOrienList = []
    posDict = {}
    mposDict = {}
    posDistanceDict = {}
    mposDistanceDict = {}
    for read in oriKeySameOrienList:
        clearKeySameOrienList.append([int(read[0].replace('[', '').replace(',', '')), int(read[1].replace(',', ''))])
    limitScore = len(clearKeySameOrienList)
    loc = 0
    for read in clearKeySameOrienList:
        # print(read)
        keyPosLeft = read[0] - ispeTime*meanISPE
        keyPosRight = read[0] + ispeTime*meanISPE
        keyPosScore = 0

        keyMposLeft = read[1] - ispeTime*meanISPE
        keyMposRight = read[1] + ispeTime*meanISPE
        keyMposScore = 0

        lrCount = 0
        posDistanceDict[read[0]] = 0
        mposDistanceDict[read[1]] = 0
        for i in range(1, limitScore):
            if(lrCount > 1):
                break
            nowLocLeft = loc - i
            if((nowLocLeft > 0) and (clearKeySameOrienList[nowLocLeft][0]>keyPosLeft) and (clearKeySameOrienList[nowLocLeft][0]<keyPosRight)):
                keyPosScore += 1
                posDistanceDict[read[0]] += abs(read[0]-clearKeySameOrienList[nowLocLeft][0])
            else:
                lrCount += 1
            if((nowLocLeft > 0) and (clearKeySameOrienList[nowLocLeft][1]>keyMposLeft) and (clearKeySameOrienList[nowLocLeft][1]<keyMposRight)):
                keyMposScore += 1
                mposDistanceDict[read[1]] += abs(read[1]-clearKeySameOrienList[nowLocLeft][1])
            else:
                lrCount += 1

        lrCount = 0
        for i in range(1, limitScore):
            if(lrCount > 1):
                break
            nowLocRight = loc + i
            if((nowLocRight < limitScore) and (clearKeySameOrienList[nowLocRight][0]>keyPosLeft) and (clearKeySameOrienList[nowLocRight][0]<keyPosRight)):
                keyPosScore += 1
                posDistanceDict[read[0]] += abs(read[0]-clearKeySameOrienList[nowLocLeft][0])
            else:
                lrCount += 1
            if((nowLocRight < limitScore) and (clearKeySameOrienList[nowLocRight][1]>keyMposLeft) and (clearKeySameOrienList[nowLocRight][1]<keyMposRight)):
                keyMposScore += 1
                mposDistanceDict[read[1]] += abs(read[1]-clearKeySameOrienList[nowLocLeft][1])
            else:
                lrCount += 1

        posDict[read[0]] = keyPosScore
        mposDict[read[1]] = keyMposScore
        loc += 1

    # print('posDict:', posDict)
    # print('mposDict:', mposDict)
    # print('posDistanceDict', posDistanceDict)
    # print('mposDistanceDict', mposDistanceDict)
    maxPos = max(posDict, key=posDict.get)
    maxPosValue = posDict[maxPos]
    if(maxPosValue != 0):
        minDistance = 999999999
        for pos in posDict:
            if(posDict[pos] == maxPosValue):
                if(posDistanceDict[pos] < minDistance):
                    minDistance = posDistanceDict[pos]
                    maxPos = pos

        # 以下两个for是为了把更精确可以在上面得到的maxPos往上推，找到在距离这个[maxPos - ispeTime*meanISPE, maxPos]间离maxPos最远的点，但不用找此时的maxPos应该也符合误差区间
        for i in range(limitScore):
            if(clearKeySameOrienList[i][0] == maxPos):
                posLoc = i
                break
        for i in range(posLoc):
            curPosLoc = posLoc-1-i
            if(curPosLoc<0): 
                break
            if((clearKeySameOrienList[curPosLoc][0] <= maxPos) and (clearKeySameOrienList[curPosLoc][0] >= (maxPos-ispeTime*meanISPE))):# (maxPos-ispeTime*meanISPE)考虑此处ispeTime换一个参数名或取固定的3？
                maxPos = clearKeySameOrienList[curPosLoc][0]

    else:    # 当（posDict+mposDict+posDistanceDict+mposDistanceDict）的value都为0时，应该取原有的candidate
        maxPos = clearKeySameOrienList[0][0]

    maxMPos = max(mposDict, key=mposDict.get)
    maxMPosValue = mposDict[maxMPos]
    if(maxMPosValue != 0):
        minDistance = 999999999
        for mpos in mposDict:
            if(mposDict[mpos] == maxMPos):
                if(mposDistanceDict[mpos] < minDistance):
                    minDistance = mposDistanceDict[mpos]
                    maxMPos = mpos
        # 以下两个for是为了更精确可以在上面得到的maxMPos往下推，找到在距离这个[maxMPos, maxMPos + ispeTime*meanISPE]间离maxMPos最远的点，但不用找此时的maxMPos应该也符合误差区间
        for i in range(limitScore):
            mposLoc = limitScore - 1 - i
            if(clearKeySameOrienList[mposLoc][1] == maxMPos):
                mposLoc = mposLoc
                break
        for i in range(1, limitScore):
            curMPosLoc = mposLoc + i
            if(curMPosLoc >= limitScore):
                break
            if((clearKeySameOrienList[curMPosLoc][1] >= maxMPos) and (clearKeySameOrienList[curMPosLoc][1] <= (maxMPos+ispeTime*meanISPE))):
                maxMPos = clearKeySameOrienList[curMPosLoc][1]

    else:    # 当（posDict+mposDict+posDistanceDict+mposDistanceDict）的value都为0时，应该取取原有的candidate
        maxMPos = clearKeySameOrienList[-1][1]

    if(maxPos>maxMPos): # 当筛选后的右断点大于左断点时，取原有的candidate
        maxPos = clearKeySameOrienList[0][0]
        maxMPos = clearKeySameOrienList[-1][1]

    # print('-----------------------------maxPos, maxMPos:', maxPos, posDict[maxPos], posDistanceDict[maxPos], ';;;', maxMPos, mposDict[maxMPos], mposDistanceDict[maxMPos])
    return maxPos, maxMPos


def getScore(oriSameOrienDir, oriCandidateDir, ispeTime, reviseCandidateDir, candidatesThreScore):
    print(oriSameOrienDir, oriCandidateDir, reviseCandidateDir)
    for oriSameOrienFileNamePart in os.listdir(oriSameOrienDir):
        print(oriSameOrienFileNamePart)
        meanISPE = int(np.loadtxt('/mnt/xiaolan/wzj/realBamResult/pindel/config/config_'+oriSameOrienFileNamePart.split('_')[0], dtype=bytes).astype(str)[1]) # config_HG03642.12.high.bam
        oriSameOrienFileName = oriSameOrienDir + oriSameOrienFileNamePart
        oriCandidateFileName = oriCandidateDir + oriSameOrienFileNamePart.split('_')[0]+'_candidates_V4_treScore' + str(candidatesThreScore)
        # reviseOrienFileName = reviseCandidateDir + oriSameOrienFileNamePart.split('_')[0]+'_soft_sameOrien_1_V4_revise'
        reviseOrienFileName = reviseCandidateDir + oriSameOrienFileNamePart.split('_')[0]+'_candidates_V4_treScore' + str(candidatesThreScore)
        oriSameOrienList = np.loadtxt(oriSameOrienFileName, str)
        oriCandidateList = np.loadtxt(oriCandidateFileName, int)
        reviseOrienList = []
        curOrienLoc = 0
        curCandiLoc = 0
        oriSameOrienListSize = len(oriSameOrienList)
        oriCandidateListSize = len(oriCandidateList)
        while((curOrienLoc < oriSameOrienListSize) and (curCandiLoc < oriCandidateListSize)):
            # print('candidate: ', oriCandidateList[curCandiLoc])
            # 注意：以下这句(oriCandidateList[curCandiLoc][2] < 500)是过滤得分高于500的情况（如HG00096.1.[46385800 46419189 15432]，这个需要另外慎重考虑会不会miss一些
            if((oriCandidateList[curCandiLoc][2] < 500) and (oriSameOrienList[curOrienLoc][0][1:-1] == str(oriCandidateList[curCandiLoc][0]))):
                oriKeySameOrienList = oriSameOrienList[curOrienLoc:curOrienLoc+oriCandidateList[curCandiLoc][2]] # 由candidate获得了一组read
                minPos, minMPos = getReadScore(oriKeySameOrienList, ispeTime, meanISPE)
                # print(minPos, minMPos, oriCandidateList[curCandiLoc][2])
                reviseOrienList.append([minPos, minMPos, oriCandidateList[curCandiLoc][2]])
                curOrienLoc += oriCandidateList[curCandiLoc][2]
                curCandiLoc += 1
            else:
                curOrienLoc += 1
        
        reviseOrienList = np.array(reviseOrienList)
        np.savetxt(reviseOrienFileName, reviseOrienList, newline='\n',fmt='%s')
        reviseOrienList = []

if __name__ == '__main__':

    # 用不同的策略生成candidate start
    print('cnn2INV_singlechrom_lugang_V4_step2_Integration.py start')
    commenDir = '/mnt/xiaolan/wzj/realBamResult/'
    step1CommenDir = commenDir + 'step1_Integration/' # 同向read的文件夹
    step2CommenDir = commenDir + 'step2_Integration/' # candidate的文件夹
    candidatesThreScore = 5 # 1 2 5

    # # get 未过滤 candidates start 1  由同向read获取candidate 
    # step1Dir = step1CommenDir + 'step1out/'
    # step2Dir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '/'
    # for sameOrienFileName in os.listdir(step1Dir):
    #     print(step1Dir + sameOrienFileName)
    #     invCandidatesFileName = getInvCandidates(step1Dir, step2Dir, sameOrienFileName, candidatesThreScore) # ccc 1
    # # get 未过滤 candidates end

    # # 过滤ispe位数大于5的 start
    # oriSameOrienDir = step1CommenDir + 'step1out/'
    # reviseOrienDir = step1CommenDir + 'step1out_miFil_threScore' + str(candidatesThreScore) + '_revise/' # step1out_miFil_revise
    # miFilReviseOriens(oriSameOrienDir, reviseOrienDir)

    # # get candidates start 2  由同向read获取candidate
    # step1Dir = step1CommenDir + 'step1out_miFil_threScore' + str(candidatesThreScore) + '_revise/' # step1out step1out_mi_revise step1out_variance_revise 
    # step2Dir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_miFil_revise/' # step2_V4_threScore1 step2_V4_threScore1_mi_revise step2_V4_threScore1_variance_revise
    # for sameOrienFileName in os.listdir(step1Dir):
    #     print(step1Dir + sameOrienFileName)
    #     invCandidatesFileName = getInvCandidates(step1Dir, step2Dir, sameOrienFileName, candidatesThreScore) # ccc 2
    # # get candidates end
    # # 过滤ispe位数大于5的 end

    # # # get candidates start   由同向read获取candidate
    # # 把与上一条ispe的距离大于2000的拆分开 start
    # step1Dir = step1CommenDir + 'step1out_miFil_threScore' + str(candidatesThreScore) + '_revise/' # step1out step1out_mi_revise step1out_variance_revise 
    # step2Dir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_miFilNear_revise/' # step2_V4_threScore1 step2_V4_threScore1_mi_revise step2_V4_threScore1_variance_revise
    # for sameOrienFileName in os.listdir(step1Dir):
    #     print(step1Dir + sameOrienFileName)
    #     invCandidatesFileName = getInvCanNear(step1Dir, step2Dir, sameOrienFileName, candidatesThreScore) # ccc 3
    # # 把与上一条ispe的距离大于2000的拆分开 end
    # # # get candidates end

    # # 幂次过滤修正candidates start
    # # mi revise orientation start 由ISPE幂次修正同向的read,(输入是ori_同向read与ori_candidate)
    # oriSameOrienDir = step1CommenDir + 'step1out/'
    # oriCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '/'
    # reviseOrienDir = step1CommenDir + 'step1out_mi_threScore'+ str(candidatesThreScore) +'_revise/'
    # miReviseOriens(oriSameOrienDir, oriCandidateDir, reviseOrienDir, candidatesThreScore)
    # # mi revise orientation end

    # # get candidates start  由同向read获取candidate
    # step1Dir = step1CommenDir + 'step1out_mi_threScore'+ str(candidatesThreScore) +'_revise/' # step1out step1out_mi_revise step1out_variance_revise 
    # step2Dir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_mi_revise/' # step2_V4_threScore1 step2_V4_threScore1_mi_revise step2_V4_threScore1_variance_revise
    # for sameOrienFileName in os.listdir(step1Dir):
    #     print(step1Dir + sameOrienFileName)
    #     invCandidatesFileName = getInvCandidates(step1Dir, step2Dir, sameOrienFileName, candidatesThreScore) # ccc 4
    # # get candidates end
    # # 幂次过滤修正candidates end

    # # 方差过滤修正candidates start
    # # variance revise orientation start 由方差修正同向read
    # oriSameOrienDir = step1CommenDir + 'step1out/'
    # oriCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '/'
    # reviseOrienDir = step1CommenDir + 'step1out_variance_threScore'+ str(candidatesThreScore) +'_revise/'
    # varianceReviseOriens(oriSameOrienDir, oriCandidateDir, reviseOrienDir, candidatesThreScore)
    # # variance revise orientation end
    
    # # get candidates start  由同向read获取candidate
    # step1Dir = step1CommenDir + 'step1out_variance_threScore'+ str(candidatesThreScore) +'_revise/' # step1out step1out_mi_revise step1out_variance_revise 
    # step2Dir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_variance_revise/' # step2_V4_threScore1 step2_V4_threScore1_mi_revise step2_V4_threScore1_variance_revise
    # for sameOrienFileName in os.listdir(step1Dir):
    #     print(step1Dir + sameOrienFileName)
    #     invCandidatesFileName = getInvCandidates(step1Dir, step2Dir, sameOrienFileName, candidatesThreScore) # ccc 5
    # # get candidates end
    # # 方差过滤修正candidates end


    # innerDistance start
    # # innerDistance修正原始candidates start
    # oriSameOrienDir = step1CommenDir + 'step1out/'
    # oriCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '/'
    # reviseCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_innerDistance_revise/'
    # getDistance(oriSameOrienDir, oriCandidateDir, reviseCandidateDir, candidatesThreScore) # ccc 6
    # # innerDistance修正原始candidates end

    # # innerDistance修正用mi修正过的candidates start
    # oriSameOrienDir = step1CommenDir + 'step1out_mi_threScore' + str(candidatesThreScore) + '_revise/' 
    # oriCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_mi_revise/'
    # reviseCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_mi_revise_innerDistance_revise/' # 为什么有一些文件大小为0？？？
    # getDistance(oriSameOrienDir, oriCandidateDir, reviseCandidateDir, candidatesThreScore) # ccc 7
    # # innerDistance修正用mi修正过的candidates end

    # # innerDistance修正用方差修正过的candidates start
    # oriSameOrienDir = step1CommenDir + 'step1out_variance_threScore' + str(candidatesThreScore) + '_revise/'
    # oriCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_variance_revise/' 
    # reviseCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_variance_revise_innerDistance_revise/' # 为什么有一些文件大小为0？？？
    # getDistance(oriSameOrienDir, oriCandidateDir, reviseCandidateDir, candidatesThreScore) # ccc 8
    # # innerDistance修正用方差修正过的candidates end
    # # innerDistance end

    # # getScore start 
    # ispeTime = 20

    # # getScore修正原始candidates start
    # oriSameOrienDir = step1CommenDir + 'step1out/'
    # oriCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '/'
    # reviseCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_innerScore_revise/'
    # getScore(oriSameOrienDir, oriCandidateDir, ispeTime, reviseCandidateDir, candidatesThreScore) # ccc 9
    # # getScore修正原始candidates end

    # # getScore修正用mi修正过的candidates start
    # oriSameOrienDir = step1CommenDir + 'step1out_mi_threScore' + str(candidatesThreScore) + '_revise/'
    # oriCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_mi_revise/'
    # reviseCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_mi_revise_innerScore_revise/'
    # getScore(oriSameOrienDir, oriCandidateDir, ispeTime, reviseCandidateDir, candidatesThreScore) # ccc 10
    # # getScore修正用mi修正过的candidates end

    # # getScore修正用方差修正过的candidates start
    # oriSameOrienDir = step1CommenDir + 'step1out_variance_threScore' + str(candidatesThreScore) + '_revise/'
    # oriCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_variance_revise/'
    # reviseCandidateDir = step2CommenDir + 'step2_V4_threScore' + str(candidatesThreScore) + '_variance_revise_innerScore_revise/'
    # getScore(oriSameOrienDir, oriCandidateDir, ispeTime, reviseCandidateDir, candidatesThreScore) # ccc 11
    # # getScore修正用方差修正过的candidates end
    # # getScore end

    # 用不同的策略生成candidate end
    print('cnn2INV_singlechrom_lugang_V4_step2_Integration.py end') #/mnt/xiaolan/wzj/realBamResult/step2_Integration/
