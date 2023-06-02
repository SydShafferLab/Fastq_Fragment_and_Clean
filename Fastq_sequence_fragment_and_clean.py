
"""
This script takes all Fastq files from a folder and does the following:

1. Breaks the sequence into fragments (user defined)
2. Determines what fragments pass the qscore threshold (user defined)
3. Keeps sequences that have at least one fragment passing threshold, currently...
... the sequences also need to be in a row: If fragment one and three are good but...
... fragment two is bad then we only keep fragment one.

run 'pip list' on your terminal to verify these packages are installed
if not please install or create an enviroment


You can also subset your fastq in the terminal with:
cat file.fastq | head -n 100 > subset.fastq

or

tail -n 52 file.fastq > subset.fastq


 Inputs

input_folder -> The path to the input fastq file.
output_folder -> The path to the output fastq file.
what_read -> R1 or R2, what file has the sequence of interest.
file_name -> The names of the files (Ex. file1.fastq.gz file2.fastq.gz), use none if you want to run all the fastqs in your folder.
length_of_sequence -> The total bp of the sequences to use.
length_of_fragments -> The size of the fragments to create.
qscore_threshold -> The minimum qscore for a fragment to be kept.



"""

# import packages
from Bio import SeqIO, bgzf, Seq
import gzip
import numpy as np
import glob
import os
import json
import textwrap
import time
import subprocess




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


def fragment_and_filter_reads(input_folder, output_folder, what_read, file_name, length_of_sequence, length_of_fragments, qscore_threshold):
	print("")
	print(' Running...')
	print("")

	# Get all the Fastq files in the folder.
	all_R = glob.glob(input_folder  + "*" + what_read + "*.fastq.gz")

	# If file_name is not empty, then only process the specified files.
	if file_name[0] != 'none':

		# Get the list of files to process.
		hold_R = []
		dir_name = os.path.dirname(all_R[0])

		for files_R in file_name:
			hold_R.append(dir_name+'/'+files_R)
		# Set all_R to the list of files to process.
		all_R = hold_R

		# Get file names
	else:
		file_name= []
		for filepath in all_R:
			basename = os.path.basename(filepath)
			file_name.append(basename)

	# Create a summary file.
	with open(output_folder+"summary.txt","w") as write_file:
		write_file.write("Summary of Fastq_sequence_fragment_and_clean.py" + "\n")
		write_file.write("" + "\n")

		for n,path in enumerate(all_R):

			with Timer("fastq: " + file_name[n]):
				
				# Create a list of new records.
				new_records = []
				length_of_new_seq = []
				all_qscores = []

				# Iterate over the records.
				with gzip.open(path, "rt") as f:
					for record in SeqIO.parse(f, "fastq"):

						# Get the sequence and letter annotations.
						sequence = str(record.seq)
						letter_annotations = record.letter_annotations

						# You first need to empty the existing letter annotations
						record.letter_annotations = {}

						# Break the sequence into fragments.
						seq_i =  str(record.seq)[0:length_of_sequence]
						fragment_seqs = textwrap.wrap(seq_i, length_of_fragments)

						# Get the qscores for each fragment.
						qscore_i = json.loads(str(letter_annotations["phred_quality"]));
						qscore_i = np.array(qscore_i);
						qscore_i = qscore_i[0:length_of_sequence]
						fragment_qscore = np.split(qscore_i, length_of_sequence/length_of_fragments)

		                # Create a new sequence and qscore.
						new_sequence = ""
						new_qscore = []
						cnt = 0
						for i,qs in enumerate(fragment_qscore):

						 	if np.mean(qs) > qscore_threshold:
						 		cnt = 1
						 		new_sequence = new_sequence + fragment_seqs[i]
						 		new_qscore.append(qs)
						 	else:
						 		break

			            # If there is at least one fragment that passes the threshold, then keep the record.
						if cnt > 0:
							new_qscore = np.concatenate(new_qscore)

							# Add q-score to the record
							all_qscores.append(new_qscore)

							# Set the sequence and qscore for the new record.
							record.seq = Seq.Seq(new_sequence)
							record.letter_annotations = {'phred_quality': new_qscore}

							# Add the record to the list of new records.
							new_records.append(record)
							length_of_new_seq.append(len(new_sequence))

						cnt = 0

					all_qscores = np.concatenate(all_qscores)
					if len(all_qscores) > 0:
						print("              Sequence [bp] =   mean:%.2f   median:%.2f   std:%.2f"%(np.mean(length_of_new_seq),np.median(length_of_new_seq),np.std(length_of_new_seq)))
						print("              Qscore =   mean:%.2f   median:%.2f   std:%.2f"%(np.mean(all_qscores),np.median(all_qscores),np.std(all_qscores))+ "\n")
						write_file.write("fastq: " + file_name[n]  + "\n")
						write_file.write("       Sequence [bp] =   mean:%.2f   median:%.2f   std:%.2f"%(np.mean(length_of_new_seq),np.median(length_of_new_seq),np.std(length_of_new_seq))+ "\n")
						write_file.write("              Qscore =   mean:%.2f   median:%.2f   std:%.2f"%(np.mean(all_qscores),np.median(all_qscores),np.std(all_qscores))+ "\n")
						write_file.write("" + "\n")
					else:
						print("              No sequences passed the threshold"+"\n")
						write_file.write("              No sequences passed the threshold"+"\n")
						write_file.write("" + "\n")

	        
	        # Write the new records to a zipped file.
			with bgzf.BgzfWriter(output_folder+"modified_"+path.split('/')[-1], "wb") as output_handle:
				SeqIO.write(sequences=new_records, handle=output_handle, format="fastq")


if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser()
	parser.add_argument('--input_folder', type=str, help="The path to the input fastq file.",required=True)
	parser.add_argument('--output_folder', type=str, help="The path to the output fastq file.",required=True)
	parser.add_argument('--what_read', type=str, help="R1 or R2, what file has the sequence of interest.",required=True)
	parser.add_argument('--file_name', type=str,nargs="*", help="The names of the files (Ex. file1.fastq.gz file2.fastq.gz), use none if you want to run all the fastqs in your folder.",required=True)
	parser.add_argument('--length_of_sequence', type=int, default=20, help="The total bp of the sequences to use.")
	parser.add_argument('--length_of_fragments', type=int, default=100, help="The size of the fragments to create.")
	parser.add_argument('--qscore_threshold', type=int, default=20, help="The minimum qscore for a fragment to be kept.")

	args = parser.parse_args()

	fragment_and_filter_reads(args.input_folder, args.output_folder, args.what_read, args.file_name, args.length_of_sequence, args.length_of_fragments, args.qscore_threshold)



