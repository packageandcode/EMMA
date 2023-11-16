import pandas as pd 
import numpy as np 
from pybedtools import BedTool
from optparse import OptionParser
import time

schema_region = {"start":np.int64,"end":np.int64}
column_region = ["chr","start","end"]
schema_CpGcontext = {"start":np.int64,"end":np.int64,"status":np.str_}
column_CpGcontext = ["readId","r1-r2","chr","loci","status"]
schema_input = {"start":np.int64,"end":np.int64}
column_input = ["chrom","start","end","readId","status1","status2"]
schema_reference = {"start":np.int64,"end":np.int64}

def qlda(x,p1):
    x.loc[:, ['readId']] = range(x.shape[0])
    return_v = {
        "chr":x.iloc[0].chrom,
        "start":x.start.min(),
        "end":x.end_context.max(),
        "readId":x.name,
        "dQB":x.QB_p.sum()-2*np.log(1-p1),
        "dQM":x.QM_p.sum()-2*np.log(p1),
    }
    return pd.Series(return_v)


if __name__=="__main__":

    parser = OptionParser()
    parser.add_option("-i","--input",dest="input_file",default="error",type = "string")
    parser.add_option("-o","--output",dest="output_file",default="error",type = "string")
    parser.add_option("--reference",dest="reference_file",type = "string")
    parser.add_option("--region",dest="region_file",type = "string")
    (options, args) = parser.parse_args()

    star = time.time()
    reference = pd.read_csv(options.reference_file,sep="\t",dtype=schema_reference)
    elapsed = (time.time() - star)
    print("R Time used:",elapsed)
    star = time.time()
    CpGcontext = pd.read_csv(options.input_file,sep="\t",header=None,names=column_input,dtype=schema_input)
    elapsed = (time.time() - star)
    print("C Time used:",elapsed)
    star = time.time()    
    CpGcontext = CpGcontext.drop_duplicates(["chrom","start","readId"])
    elapsed = (time.time() - star)
    print("dedup Time used:",elapsed)
    star = time.time()
    joined = pd.merge(left=CpGcontext,right=reference,on=["chrom","start"],suffixes=('_context','_refercence'))
    elapsed = (time.time() - star)
    print("merged Time used:",elapsed)
    star = time.time()
    joined['joined_status3'] = joined.status1.map({"+":1,"-":0})
    joined['QB_p'] = ((joined['joined_status3'] - joined['Benign_mean'])/(joined['Benign_sd']+0.01))**2 + ((joined['Benign_sd']+0.01)**2).apply(np.log)
    joined['QM_p'] = ((joined['joined_status3'] - joined['Malignant_mean'])/(joined['Malignant_sd']+0.01))**2 + ((joined['Malignant_sd']+0.01)**2).apply(np.log)
    print(joined.head(1))
    print(joined.shape[0])
    print(joined.shape[1])
    elapsed = (time.time() - star)
    print("joined Time used:",elapsed)
    star = time.time()
    
    joined = joined.loc[:, ['chrom','start','end_context','readId','QB_p','QM_p']]
    dqda_all_reads = joined.groupby("readId").apply(qlda,p1=0.2)
    elapsed = (time.time() - star)
    print("all Time used:",elapsed)
    star = time.time()

    dqda_m_reads = dqda_all_reads[dqda_all_reads['dQM']<dqda_all_reads['dQB']]
    elapsed = (time.time() - star)
    print("m Time used:",elapsed)
    star = time.time()

    dqda_areads_bed = BedTool.from_dataframe(dqda_all_reads)
    dqda_mreads_bed = BedTool.from_dataframe(dqda_m_reads)
    elapsed = (time.time() - star)
    print("dqda Time used:",elapsed)
    star = time.time()
    region = BedTool(options.region_file)
    elapsed = (time.time() - star)
    print("region Time used:",elapsed) 
    star = time.time()
    mdf = region.coverage(dqda_mreads_bed).to_dataframe()
    adf = region.coverage(dqda_areads_bed).to_dataframe()
    pd.DataFrame({'mcount':mdf['name'],'acount':adf['name']}).to_csv(options.output_file,index=False,sep="\t")
    elapsed = (time.time() - star)
    print("out Time used:",elapsed)
    




