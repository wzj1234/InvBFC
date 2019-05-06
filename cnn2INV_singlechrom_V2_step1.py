# -*- coding: UTF-8 -*-
# cd /usr/local/lib/python3.6/dist-packages
import numpy as np
import pysam
from PIL import Image, ImageDraw
import re
import time
import multiprocessing
import os

def getMinInvLenReadInform(bamFileName, chr_id, lPoint, rPoint, resFileName): 

    samfile = pysam.AlignmentFile(bamFileName, "rb")
    reads = []
    for read in samfile.fetch(chr_id, lPoint, rPoint):
        
        cigarstring = read.cigarstring
        if cigarstring==None:
            continue
        template_length = read.template_length 
        if template_length<0: # 同向的read只记录了lPoint，所以也解决了多线程断点read不会因被切分漏读的问题
            continue
        mate_is_reverse = read.mate_is_reverse
        query_name = read.query_name
        is_reverse = read.is_reverse
        isSameOrien = (mate_is_reverse == is_reverse)
        reference_start = read.reference_start
        mpos = read.mpos

        if isSameOrien==True and mpos>reference_start:
            try: 
                mate = samfile.mate(read)
            except:
                continue
            reads.append([reference_start, mpos, isSameOrien, is_reverse, query_name, template_length, cigarstring, read.query_length])  # len(read.query_alignment_qualities), read.is_paired, mate_is_reverse, 
            if len(reads) >= 100:
                # if os.path.getsize(file_dir + file) > 500000:
                #     print(file_dir + file)
                f = open(resFileName,"a")
                for read in reads:
                    f.write(str(read)+'\n')
                f.close()
                reads = []

    samfile.close() 

    f = open(resFileName,"a")
    for read in reads:
        f.write(str(read)+'\n')
    f.close()

def multiProcese(bamFileName, chr_id, lPoint, rPoint, resultDataRoot):

    anInterval = (rPoint - lPoint)/10
    resGetMinInvLenReadInform = []
    commonPath = resultDataRoot + bamFileName.split('/')[-1].split('.bam')[0]+'_'
    for i in range(1, 11):
        resGetMinInvLenReadInform.append(commonPath+str(lPoint+(i-1)*anInterval)+'-'+str(lPoint+i*anInterval)+'_soft_sameOrien_' + chr_id + '_V2.txt')

    try:

        p1 = multiprocessing.Process(target=getMinInvLenReadInform,args=(bamFileName, chr_id, lPoint, lPoint + anInterval, resGetMinInvLenReadInform[0]))
        p2 = multiprocessing.Process(target=getMinInvLenReadInform,args=(bamFileName, chr_id, lPoint + anInterval, lPoint + 2*anInterval, resGetMinInvLenReadInform[1]))
        p3 = multiprocessing.Process(target=getMinInvLenReadInform,args=(bamFileName, chr_id, lPoint + 2*anInterval, lPoint + 3*anInterval, resGetMinInvLenReadInform[2]))
        p4 = multiprocessing.Process(target=getMinInvLenReadInform,args=(bamFileName, chr_id, lPoint + 3*anInterval, lPoint + 4*anInterval, resGetMinInvLenReadInform[3]))
        p5 = multiprocessing.Process(target=getMinInvLenReadInform,args=(bamFileName, chr_id, lPoint + 4*anInterval, lPoint + 5*anInterval, resGetMinInvLenReadInform[4]))
        p6 = multiprocessing.Process(target=getMinInvLenReadInform,args=(bamFileName, chr_id, lPoint + 5*anInterval, lPoint + 6*anInterval, resGetMinInvLenReadInform[5]))
        p7 = multiprocessing.Process(target=getMinInvLenReadInform,args=(bamFileName, chr_id, lPoint + 6*anInterval, lPoint + 7*anInterval, resGetMinInvLenReadInform[6]))
        p8 = multiprocessing.Process(target=getMinInvLenReadInform,args=(bamFileName, chr_id, lPoint + 7*anInterval, lPoint + 8*anInterval, resGetMinInvLenReadInform[7]))
        p9 = multiprocessing.Process(target=getMinInvLenReadInform,args=(bamFileName, chr_id, lPoint + 8*anInterval, lPoint + 9*anInterval, resGetMinInvLenReadInform[8]))
        p10 = multiprocessing.Process(target=getMinInvLenReadInform,args=(bamFileName, chr_id,lPoint + 9*anInterval, rPoint, resGetMinInvLenReadInform[9]))
        
        p1.start()
        p2.start()
        p3.start()
        p4.start()
        p5.start()
        p6.start()
        p7.start()
        p8.start()
        p9.start()
        p10.start()
         
        p1.join()
        p2.join()
        p3.join()
        p4.join()
        p5.join()
        p6.join()
        p7.join()
        p8.join()
        p9.join()
        p10.join()
         
    except:
        print("Error: unable to start thread")

    sameOrienFileName = commonPath + str(lPoint) + '-' + str(rPoint) + '_soft_sameOrien_' + chr_id + '_V2'
    for aRes in resGetMinInvLenReadInform: # merge all same orientation candidates' reads
        os.system('cat ' + aRes + ' >> ' + sameOrienFileName)
        os.system('rm ' + aRes)
    print('same orientation reads: ', sameOrienFileName)
    return sameOrienFileName


if __name__ == '__main__':

    print('cnn2INV_singlechrom_V2_step1.py start')
    start = time.time()
    # os.system('')

# NA12878 ['1', '2', '5', '7', '8', '9', '10', '12', '16', '17', '21']
# NA18525 ['1', '2', '3', '4', '5', '7', '9', '10', '11', '12', '14', '16', '21']
# NA19017 ['1', '2', '3', '4', '5', '6', '7', '9', '10', '11', '12', '14', '15', '16', '18', '21']
# NA19238 ['1', '2', '3', '4', '6', '7', '9', '10', '11', '12', '14', '16', '18', '21', '22']
# NA19239 ['1', '2', '3', '4', '6', '7', '9', '10', '11', '12', '14', '16', '17', '18', '21']
# NA19625 ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '14', '16', '20', '21']
# NA19648 ['1', '2', '3', '5', '7', '9', '10', '11', '12', '16', '18', '21']
# NA20502 ['1', '2', '3', '5', '7', '8', '9', '10', '11', '12', '14', '16', '18', '21', '22']
# NA20845 ['1', '2', '3', '5', '6', '7', '9', '10', '11', '12', '14', '16', '18', '20', '21']
# cd /mnt/4T/gao/1000_genome_bam/
# cailei  NA12878  NA18525  NA19017  NA19238  NA19239  NA19240  NA19625  NA19648  NA20502  NA20845  touch  vcf  view_chr1_chr2.sh
# cd /mnt/hdb/Xunlei/
# NA12878  NA18525  NA19017  NA19238  NA19240  NA19625  NA20502  NA20845
# 1    LN:249250621
# 2    LN:243199373
# 3    LN:198022430
# 4    LN:191154276
# 5    LN:180915260
# 6    LN:171115067
# 7    LN:159138663
# 8    LN:146364022
# 9    LN:141213431
# 10   LN:135534747
# 11   LN:135006516
# 12   LN:133851895
# 13   LN:115169878
# 14   LN:107349540
# 15   LN:102531392
# 16   LN:90354753
# 17   LN:81195210
# 18   LN:78077248
# 19   LN:59128983
# 20   LN:63025520
# 21   LN:48129895
# 22   LN:51304566

    resultDataRoot = "/home/wuzhongjia/deepLearning/data/candidates_orien_png_V6/"

    # filter same orientation start
    bamFileName = "/mnt/new4T/highBam/NA18525/NA18525.1.high.bam" # /mnt/4T/BAM/high/NA19648.11.high.bam /mnt/hdb/NA20845.11/NA20845.11.high.bam /mnt/4T/BAM/high/NA18525.3.high.bam
    chr_id = "1"
    lPoint = 0
    rPoint = 250000000
    print(bamFileName, chr_id, rPoint)
    sameOrienFileName = multiProcese(bamFileName, chr_id, lPoint, rPoint, resultDataRoot)
    # filter same orientation end

    end = time.time()
    print('cnn2INV_singlechrom_V2_step1.py running time:', str(round(end-start,3))+'s')
    print('cnn2INV_singlechrom_V2_step1.py end')



