
# This script takes all Fastq files from a folder and does the following:
# 1) breaks the sequence into fragments (user defined)
# 2) determines what fragments pass the qscore threshold (user defined)
# 3) keeps sequences that have at least one fragment passing threshold, currently...
#... the sequences also need to be in a row: If fragment one and three are good but...
#... fragment two is bad then we only keep fragment one.

# run 'pip list' on your terminal to verify these packages are installed
# if not please install or create an enviroment

# import packages
from Bio import SeqIO
from Bio import Seq
import numpy as np
import glob
import os
import json
import textwrap
import time
import subprocess


#----------------------------------------------------------------
#----------------------------------------------------------------
##### Users should modulate these parametes to run this script

length_of_sequence  = 20 # how many bp do you want to use 
length_of_fragments = 10 # fragment size, make length_of_sequence divisible by the fragment size
threshold_qscore = 20    # Is the mean of the fragments qscore > threshold_qscore

fastq_folder = "/Users/raul/Documents/GitHub/Barcode_cut_length/test_data/" # "/path/to/FastqFolder/"
output_folder = "/Users/raul/Documents/GitHub/Barcode_cut_length/test_data/" # this folder should exist  "/path/to/output_folder/"


what_read = "R1" # R1 or R2 what read have the sequence

file_name = ["Test1_gDNA_L001_R1_001.fastq","Test2_gDNA_L001_R1_001.fastq"] # leave this empty [] if you want to run all the fastqs in your folder


#### END OF USER DEFINED PARAMETERS ##### 









#----------------------------------------------------------------
#----------------------- functions -----------------------
# to kee track of itme 
class Timer(object):
	def __init__(self, name=None):
		self.name = name
	def __enter__(self):
		self.tstart = time.time()
		if self.name:
			print('		%s' % self.name,)

	def __exit__(self, type, value, traceback):
		print('			    Elapsed: %.4f sec ' % (time.time() - self.tstart))
#----------------------------------------------------------

print("")
print(' Running...')
print("")


all_R = glob.glob(fastq_folder  + "*" + what_read + "*.fastq")



#unzip all files
print(" Unzip...")
gunzipCommand = ['gunzip', '-r', fastq_folder ]
subprocess.call(gunzipCommand)
print(" ")


if file_name != []:
	hold_R = []
	dir_name = os.path.dirname(all_R[0])
	for files_R in file_name:
		hold_R.append(dir_name+'/'+files_R)
	all_R = hold_R

with open(output_folder+"summary.txt","w") as write_file:
	write_file.write("Summary of Fastq_sequence_fragment_and_clean.py" + "\n")
	write_file.write("" + "\n")

	for n,path in enumerate(all_R):

		with Timer("fastq: " + path):
			

			new_records = []
			length_of_new_seq = []
			for record in SeqIO.parse(path, "fastq"):
				sequence = str(record.seq)
				letter_annotations = record.letter_annotations

				# You first need to empty the existing letter annotations
				record.letter_annotations = {}


				seq_i =  str(record.seq)[0:length_of_sequence]
				fragment_seqs = textwrap.wrap(seq_i, length_of_fragments)
				

				qscore_i = json.loads(str(letter_annotations["phred_quality"]));
				qscore_i = np.array(qscore_i);
				qscore_i = qscore_i[0:length_of_sequence]
				fragment_qscore = np.split(qscore_i, length_of_sequence/length_of_fragments)


				new_sequence = ""
				new_qscore = []
				cnt = 0
				for i,qs in enumerate(fragment_qscore):
				 	if np.mean(qs) > threshold_qscore:
				 		cnt = 1
				 		new_sequence = new_sequence + fragment_seqs[i]
				 		new_qscore.append(qs)
				 	else:
				 		break
				 	
				
				if cnt > 0:
					new_qscore = np.concatenate(new_qscore)

					record.seq = Seq.Seq(new_sequence)
					record.letter_annotations = {'phred_quality': new_qscore}

					new_records.append(record)

					length_of_new_seq.append(len(new_sequence))
				cnt = 0
				

			print("              Sequence [bp] =   mean:%.2f   median:%.2f   std:%.2f"%(np.mean(length_of_new_seq),np.median(length_of_new_seq),np.std(length_of_new_seq)))
			print("              Qscore =   mean:%.2f   median:%.2f   std:%.2f"%(np.mean(new_qscore),np.median(new_qscore),np.std(new_qscore))+ "\n")
			write_file.write("fastq: " + path + "\n")
			write_file.write("       Sequence [bp] =   mean:%.2f   median:%.2f   std:%.2f"%(np.mean(length_of_new_seq),np.median(length_of_new_seq),np.std(length_of_new_seq))+ "\n")
			write_file.write("              Qscore =   mean:%.2f   median:%.2f   std:%.2f"%(np.mean(new_qscore),np.median(new_qscore),np.std(new_qscore))+ "\n")
			write_file.write("" + "\n")

		with open(output_folder+"modified_"+path.split('/')[-1], 'w') as output_handle:
		    SeqIO.write(new_records, output_handle, "fastq")




print("")
print(' Done :D')
print("")
