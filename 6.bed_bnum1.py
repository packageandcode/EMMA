########!/~/miniconda3/bin/python
# -*-coding:utf-8 -*-

#######

#file format:
#sample_name  chr start end ...
#colnames are splited by tab
#file have been sorted by " sort -t $'\t' -k 1,1  -k  2,2 "

#output file format:
#chr start end num...
#colnames are splited by tab

#bed_bnum.py input_file.bed step output_file.bed

import time
import sys
import numpy as np

def Checktime():
    return time.asctime( time.localtime(time.time()) )

def Readfile(rootdir):
    file = []
    with open(rootdir, "r") as lines:
        for line in lines:
            line = line.rstrip()
            arr = line.split("\t")
            ins = list(arr)
            file.append(ins)
    return file


#file_ampl=Readfile("E:/cooperation/cfDNA/tissue/CNV/cg.regions/num.merge/wgs.ampl.events.res.nochrY.bed")
#file_dele=Readfile("E:/cooperation/cfDNA/tissue/CNV/cg.regions/num.merge/wgs.dele.events.res.nochrY.bed")
file_input_path = sys.argv[1]
bed_list = Readfile(file_input_path)
step=sys.argv[2]
asctime = Checktime()
print(str(asctime)+" : "+str(file_input_path)+" have been read, which had "+str(len(bed_list))+" rows. " +
      " step: "+str(step))


#######
def bed_merge_num(bed_list,step):
    # merge same cnv in adjacent regions in same sample
    bed_list=bed_list
    step=int(step)
    bed_hash={}
    keys=""
    for i in range(len(bed_list)):
        keys = bed_list[i][0] + "," + str(bed_list[i][1])  # set keys as "sam,chr"
        if (keys in bed_hash.keys()):
            bed_hash[keys] = bed_hash[keys]+list([bed_list[i][2],bed_list[i][3]]) # if keys have been exsited, append current value to Previous values
        else:
            bed_hash[keys] = list([bed_list[i][2],bed_list[i][3]]) #{sam,chr:start,end,start,end,start,end.....}

    uniq_list=[]
    uniq_hash={}
    k = list(bed_hash.keys())
    for i in range(len(k)):
        uniq_list=[int(int(x) / step) for x in bed_hash[k[i]]]
        for j in range(len(uniq_list)-2,-1,-1): #-1使倒序排列
            if (uniq_list[j] == uniq_list[j - 1]):
                uniq_list.remove(uniq_list[j])
                uniq_list.remove(uniq_list[j - 1])   # 01sam chr1 [0,1,1,2,8,10] →  01sam chr1 [0,2,8,10]
        keys=k[i].split(",")[1]
        if (keys in uniq_hash.keys()):
            uniq_hash[keys] = uniq_hash[keys]+uniq_list # if keys have been exsited, append current value to Previous values
        else:
            uniq_hash[keys] = uniq_list #{chr:start,end,start,end,start,end.....}

    #sum num of sample which is in region
    res_list=[]
    kq = list(uniq_hash.keys())
    bed_chr=[]
    for i in range(len(kq)):
            k_h=uniq_hash[kq[i]]
            start_list = k_h[::2]  #从0开始步长为2取奇数位
            end_list = k_h[1::2] #从1开始步长为2取奇数位
            bed_idx=[]
            bed_idx=[0 for x in range(max(end_list)+1)]
            for j in range(len(start_list)):
                start=start_list[j]
                end=end_list[j]
                for site in range(start, end+1):
                    bed_idx[site]=bed_idx[site]+1
            num_list=[]  #Delete adjacent bins which have same number, leaving start and end
            num_list.append(list([0,bed_idx[0]]))
            for idx in range(1,len(bed_idx)):
                num = bed_idx[idx]
                num_1 = bed_idx[idx - 1]
                if num != num_1:
                    num_list.append(list([idx-1, num_1]))
                    num_list.append(list([idx, num]))
                elif idx == (len(bed_idx) - 1):
                    num_list.append(list([idx, num]))
            for ni in range(1,len(num_list)):     #merge the same number of bins
                if (num_list[ni][1]==num_list[ni-1][1])&(num_list[ni][1]!=0):  #chr,start,end,num
                    res_list.append(list([kq[i],num_list[ni-1][0]*step,num_list[ni][0]*step,num_list[ni][1]]))
    return res_list


res_list=bed_merge_num(bed_list,step)

#time2=time.asctime( time.localtime(time.time()) )
asctime = Checktime()
print(str(asctime)+" : "+" Counting the occurrences of multiple samples within the same genomic region has been completed.")

out_path = sys.argv[3]
#out_path = 'E:/cooperation/cfDNA/tissue/CNV/cg.regions/num.merge/ampl.chr.regions.num.res.txt'
result=res_list
f = open(out_path, 'w')
for i in range(len(result)):
    for j in range(len(result[i])):
        f.write(str(result[i][j]))
        if j < len(result[i]) - 1:
            f.write('\t')
    f.write('\n')
f.close()

#time3=time.asctime( time.localtime(time.time()) )
asctime = Checktime()
print(str(asctime)+" : "+str(out_path)+" have been output, which had "+str(len(res_list))+" rows.")
