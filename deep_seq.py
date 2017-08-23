#=coding:utf-8
import sys
import os
import numpy as np 
import pandas as pd 
import Bio
import Bio.SeqIO
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
import math 
import re
from pandas import DataFrame
from Bio.SeqUtils import GC


def statistics_reads(design_filename,out_filename):
	design_frame = pd.read_csv(design_filename,delimiter = "\t", header =0)
	oligo_table = pd.read_csv('/home/comp/rmlab/maxim/sumithra_crispr/TF/human/development/QC/jun20_data/oligo_design_all',delimiter = "\t",header = 0)

	x = design_frame['oligo_name'].value_counts()
	x_df = pd.DataFrame(x).reset_index()
	x_df.columns = ['oligo_name','read_count']
	print x_df['read_count'].sum(),x_df['read_count'].mean(),x_df['read_count'].median(),x_df['read_count'].max(),x_df['read_count'].min()

	x_new = pd.merge(x_df,oligo_table, on = 'oligo_name', how = 'inner')

	x_new['Gene'] = x_new['oligo_name'].str.split('_').str[0]
	x_new['transcript_id'] = x_new['oligo_name'].str.split("_").str[1]+"_"+x_new['oligo_name'].str.split("_").str[2]

	x_new['GC_content'] = x_new['oligo_sequence'].apply(GC)
	x_new.index = x_new['Gene']
	x_new = x_new.sort()
	del x_new['Gene']

	x_new.to_csv(out_filename,sep = '\t')

def library_characterization(library_file):

	# readin the library design file
	library_frame = pd.read_csv(library_file,delimiter = '\t',header =0)
	
	# make columns for Gene and transcript ids
	library_frame['Gene'] = library_frame['oligo_name'].str.split("_").str[0]
	library_frame['transcript_id'] = library_frame['oligo_name'].str.split("_").str[1]+"_"+library_frame['oligo_name'].str.split("_").str[2]

	#throw out controls
	only_gRNA = library_frame[~library_frame['oligo_name'].str.contains('_cntrl')]

	x = only_gRNA['transcript_id'].value_counts()
	x_df = pd.DataFrame(x).reset_index()
	x_df.columns = ['transcript_id','gRNA_count']
	x_df.index = x_df['transcript_id']
	del x_df['transcript_id']
        x_df.to_csv("gRNApertranscript_library_histogram",sep = "\t")

	x = only_gRNA['Gene'].value_counts()
	x_df = pd.DataFrame(x).reset_index()
	x_df.columns = ['Gene','gRNA_count']
	x_df.index = x_df['Gene']
	del x_df['Gene']
        x_df.to_csv("gRNApergene_library_histogram",sep = "\t")

def chosen_oligo_set(inputfile,outputfile):
	reads_summary_df = pd.read_csv(inputfile,delimiter = "\t",header = 0)
	reads_good = reads_summary_df[reads_summary_df['dist_median_measure'] <= 10]
	reads_good.to_csv(outputfile,sep = "\t",index = False)

def hamming(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def quality_control(input_file,qc_file):
	downstream_read2_start = "aacttgctaTTTCtagctc"
	qc_sequences = []

#QC reads for exact match (upto 3 mm) with downstream_read2_start
	for record in SeqIO.parse(open(input_file,"rU"),"fastq"):
      		s2 = str(record.seq[55:74])
        	if hamming(downstream_read2_start.upper(),s2) <= 3 :
			print downstream_read2_start,s2
                	qc_sequences.append(record)


	output_handle = open(qc_file,"w")
	SeqIO.write(qc_sequences,output_handle,"fastq")
	output_handle.close()

def main(): 
	
	library_characterization('/home/comp/rmlab/maxim/sumithra_crispr/TF/human/development/QC/jun20_data/oligo_design_all')

	inp_file = sys.argv[1]
	quality_control(inp_file,"Reads2_qc.fastq")

	#map back reads to library
	os.system('bowtie --suppress 6,7 -v 3 -5 80 -3 150 /home/comp/rmlab/maxim/sumithra_crispr/TF/human/development/QC/human_tf_bowtie_indexes/human_tf /home/comp/rmlab/maxim/sumithra_crispr/TF/human/development/QC/jun20_data/Reads2_qc.fastq  /home/comp/rmlab/maxim/sumithra_crispr/TF/human/development/QC/jun20_data/Reads2.bowtie') 

	#save bowtie file as dataframe and create coloumn ids and rewrite into same bowtie file
	reads2_df = pd.read_csv('/home/comp/rmlab/maxim/sumithra_crispr/TF/human/development/QC/jun20_data/Reads2.bowtie',delimiter = "\t",names = ['read_id','ref_strand','oligo_name','offset','sequence','mismatches'],header = None) 
	reads2_df.index = reads2_df['read_id']
	del reads2_df['read_id']
	reads2_df.to_csv("/home/comp/rmlab/maxim/sumithra_crispr/TF/human/development/QC/jun20_data/Reads2.bowtie",sep = '\t')

	#from new bowtie file, cut the oligo_names which aligned
	os.system('cut -f 3 /home/comp/rmlab/maxim/sumithra_crispr/TF/human/development/QC/jun20_data/Reads2.bowtie > /home/comp/rmlab/maxim/sumithra_crispr/TF/human/development/QC/jun20_data/Reads2.txt')

	#now count the number of reads per oligo, normalize and calculate distance median measure
	statistics_reads('/home/comp/rmlab/maxim/sumithra_crispr/TF/human/development/QC/jun20_data/Reads2.txt','/home/comp/rmlab/maxim/sumithra_crispr/TF/human/development/QC/jun20_data/reads_summary.txt')
	
if __name__ == '__main__': 
	main() 

######end of main code####

