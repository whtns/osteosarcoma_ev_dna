#!/usr/bin/python

import sys
import glob
import pandas as pd
import shutil

bam_dir = sys.argv[1]

maverix_bams = glob.glob(bam_dir+"*bam*")

maverix_metadata = pd.read_csv("~/os_exosome_pipeline/data/Phase1_OS_controls_CHL_HLOH_sent_SBI.csv")

maverix_metadata["maverix_ids"] = maverix_metadata["maverix_ids"].str.replace("Sample_", "")

id_dict = dict(zip(maverix_metadata["maverix_ids"], maverix_metadata["study_id"]))

def multipleReplace(text, wordDict):
    for key in wordDict:
        text = text.replace(key+".", wordDict[key]+".")
    return text

maverix_bams_clean = [multipleReplace(i, id_dict) for i in maverix_bams]

file_name_dict = dict(zip(maverix_bams, maverix_bams_clean))

for key in file_name_dict:
	if (not key == file_name_dict[key]):
		shutil.copyfile(key, file_name_dict[key])
