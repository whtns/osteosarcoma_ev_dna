#!/usr/bin/env python

import os
import re

pattern = r'\.bam.*\.bai.*'
pattern2 = r'\.sorted.*\.bam.*'

replacement = '.bam.bai'
replacement2 = '.sorted.bam'

file_path = '/media/thor/storage/os_exosome_pipeline/chla_normal_vs_os/'
for file in os.listdir(file_path):
	if not file.endswith('.bai'):
		print file
		f = re.sub(pattern2, replacement2 , file)
		os.rename(os.path.join(file_path, file), os.path.join(file_path, f))


