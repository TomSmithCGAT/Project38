'''
Purpose:
rename fastqs 

To do:
Generalise!
Assert to avoid error

This script will rename fastq files based on a sample_id_map.csv file 
in the format below.
The script assumes the original fastq filenames end with 
_1_sequence.txt.gz and _2_sequence.txt.gz for the paired ends.


File    patient_id      sample_type
WTCHG_37419_04  56JCV   Adeno
WTCHG_37419_02  56JCV   Control
WTCHG_37419_05  FFHUZ   Control
WTCHG_37419_06  FFHUZ   Adeno

Note: if some files have been renamed already, script will through an error!
-re-work to remove this with an assert statement

'''
import os
os.chdir("/ifs/projects/proj038/readqc")
sample_map=open("sample_id_map.csv")
sample_map.next()
for line in sample_map:
    old=line.split()[0]
    new="-".join(line.split()[1:])
    os.rename(old+"_1_sequence.txt.gz",new+"-1.fastq.1.gz")
    os.rename(old+"_2_sequence.txt.gz",new+"-1.fastq.2.gz")                
   




