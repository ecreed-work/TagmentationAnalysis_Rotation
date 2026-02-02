import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import NGSanalysis_module as ngs
import subprocess
import tempfile
import os
import re
import numpy as np


#configuration:
#samtools
samtools_bin = '/Users/ekello73/src/samtools-1.22.1/bin/bin/'
#bbmap
bbmap_bin = '/Users/ekello73/src/bbmap/'

st.title("Kellogg Lab Tagmentation Analysis Toolkit")


# initialize
fastq_filepath = 0
referencegenome_filepath = 0
primer = 0
first_ten_in_donor = 0
sgRNA_sequence = 0

#example inputs:
# fastq_filepath = "/Users/ekello73/projects/NGS/PrashantTagmentation/PD1_merged.fastq"
# primer = "CTGAAAAACAACCACCACGACATTAATTTGCGAATAACGACACTAAATTGCGAAAAGCGACATTTAATTTGCGAATGTACA" # This is actually the primer -> end of LE
# referencegenome_filepath = "/Users/ekello73/projects/NGS/genomes/cJP003_assembly.fasta"
# first_ten_in_donor = "TATATAATGG"
# sgRNA_sequence = "AGCCCAAAAAAACGGGTATGGAGA"

fastq_filepath = st.file_uploader("fastq file")
primer = st.text_input("primer mapping to donor:")
referencegenome_filepath = st.file_uploader("reference genome fasta")
first_ten_in_donor = st.text_input("first ten basepairs of donor")
sgRNA_sequence = st.text_input("sgRNA sequence")

# Filepath for Output files from bbmap -- change to the absolute filepath where you want them exported to, end a dir with "/"
#will change to a default path and not make this an option
#path_to_output_dir = st.text_input("")
path_to_output_dir= "/Users/ekello73/projects/NGS/TagmentationAnalysis/test/"

# Adds this to the beginning of each output file for simplicity
#sample_abbreviation = st.text_input("output_prefix")
sample_abbreviation = "test"

submit = st.button("Run Analysis")

if submit:
	if not fastq_filepath or referencegenome_filepath == 0:
		st.warning("Please complete all fields before running the analysis.")
	else:
		st.success(f"Running analysis...")

        # Useful variables
		genome_sequence = referencegenome_filepath.read().decode("utf-8")
		genome_sequence = re.sub(r'^>.*\n?','',genome_sequence,flags=re.MULTILINE)
		genome_sequence = genome_sequence.replace('\n','').replace('\r','')
		ccdb_gene = "ATGCAGTTTAAGGTTTACACCTATAAAAGAGAGAGCCGTTATCGTCTGTTTGTGGATGTACAGAGTGATATTATTGACACGCCCGGGCGACGGATGGTGATCCCCCTGGCCAGTGCACGTCTGCTGTCAGATAAAGTCTCCCGTGAACTTTACCCGGTGGTGCATATCGGGGATGAAAGCTGGCGCATGATGACCACCGATATGGCCAGTGTGCCGGTCTCCGTTATCGGGGAAGAAGTGGCTGATCTCAGCCACCGCGAAAATGACATCAAAAACGCCATTAACCTGATGTTTTGGGGAATATAA"
		ccdb_start_site = genome_sequence.find(ccdb_gene)
        
		print("ccdb_start_site: "+str(ccdb_start_site)+" genome size: "+str(len(genome_sequence)))
		#from the new jupyter notebook... is this needed??
		# Length of the CAST's PAM Sequence
		PAM_len = 3
		# Is this a naturally occuring CAST system?
		natural_system = True # Set to false if Novel CAST -- integration site relative to the PAM is different than in natural systems
		sgRNA_start_site, PAM_start_site, PAM_end_site, is_rc = ngs.sgRNA_finder(sgRNA_sequence,genome_sequence, PAM_len, natural_system)
		print("sgRNA start site is: " + str(sgRNA_start_site))

		highlight = [
			(sgRNA_start_site, sgRNA_start_site + len(sgRNA_sequence), "red"),
			(PAM_start_site, PAM_end_site, "black"),
		]

    
		#tmpdir = tempfile.TemporaryDirectory()
		tmpdir = 'test_dir_webapp'
		fastq_tempfile = ngs.save_to_temp(fastq_filepath,tmpdir)
		genome_tempfile = ngs.save_to_temp(referencegenome_filepath,tmpdir)

		[forward_insertions, reverse_insertions] = ngs.process_reads(fastq_tempfile,
                                                                     primer,
                                                                     first_ten_in_donor,
                                                                     genome_tempfile,
                                                                     bbmap_bin,
                                                                     samtools_bin,
                                                                     path_to_output_dir)


		st.download_button(
			label="Download forward insertion reads",
			data="\n".join(map(str, forward_insertions)),
			file_name="forward_insertions.txt",
			mime="text",
		)

		st.download_button(
			label="Download reverse insertion reads",
			data="\n".join(map(str, reverse_insertions)),
			file_name="reverse_insertions.txt",
			mime="text",
		)
		# Plot Integration profiles in genome and in CCDB gene
		fig1 = ngs.plot_insertion_histogram(forward_insertions,reverse_insertions,region_start=0,region_end=len(genome_sequence),title="Genome",highlight_regions=highlight)
		st.pyplot(fig1.gcf())
		fig2 = ngs.plot_insertion_histogram(forward_insertions,reverse_insertions,region_start=ccdb_start_site-50,region_end=ccdb_start_site+len(ccdb_gene),title="CCDB Gene",highlight_regions=highlight)
		st.pyplot(fig2.gcf())
        	# Plot Integration Profiles
		modified_forward_insertions = ngs.integration_finder(forward_insertions,ccdb_gene,ccdb_start_site,PAM_start_site,PAM_end_site,natural_system,is_rc)
		modified_reverse_insertions = ngs.integration_finder(reverse_insertions,ccdb_gene,ccdb_start_site,PAM_start_site,PAM_end_site,natural_system,is_rc)
		combined_integrations = modified_forward_insertions + modified_reverse_insertions
		fig3 = ngs.plot_integration_profile(combined_integrations)
		st.pyplot(fig3.gcf())

		# --- Calculate integration stats ---
		fwd_total = len(forward_insertions)
		rev_total = len(reverse_insertions)
		fwd_on_target = len(modified_forward_insertions)
		rev_on_target = len(modified_reverse_insertions)
		total = fwd_total + rev_total
		on_target_total = fwd_on_target + rev_on_target

		fwd_pct = 100 * fwd_on_target / fwd_total if fwd_total else 0
		rev_pct = 100 * rev_on_target / rev_total if rev_total else 0
		total_pct = 100 * on_target_total / total if total else 0

		all_lengths = modified_forward_insertions + modified_reverse_insertions
		mean_len = np.mean(all_lengths) if all_lengths else 0
		median_len = np.median(all_lengths) if all_lengths else 0
		std_len = np.std(all_lengths) if all_lengths else 0

		# Percent Integration in ccdb
		report_text = f"""Integration Analysis Report
==================================
Total insertions:
  Forward strand: {fwd_total}
  Reverse strand: {rev_total}
  Combined total: {total}

On-target insertions (within CCDB region):
  Forward strand: {fwd_on_target} ({fwd_pct:.2f}%)
  Reverse strand: {rev_on_target} ({rev_pct:.2f}%)
  Combined on-target: {on_target_total} ({total_pct:.2f}%)

Integration length statistics (bp):
  Mean length:   {mean_len:.2f}
  Median length: {median_len:.2f}
  Std deviation: {std_len:.2f}
"""
		st.write(report_text)
