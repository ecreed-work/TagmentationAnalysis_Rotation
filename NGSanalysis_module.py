from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
from pathlib import Path
import matplotlib.pyplot as plt
import subprocess
import tempfile
import os

def fasta_to_string(fasta_path):
    return str(next(SeqIO.parse(fasta_path, "fasta")).seq)

def trim_reads_by_primer(fastq_filepath, primer):
    total_reads = []
    
    trimmed_reads = []

    primer = primer.upper()

    primer_len = len(primer)
    for rec in SeqIO.parse(fastq_filepath, "fastq"):
        seq = str(rec.seq)
        total_reads.append(seq)
        primer_idx = seq.find(primer)
        if primer_idx != -1:
            trimmed_seq = seq[primer_idx + primer_len:]
            trimmed_qual = rec.letter_annotations["phred_quality"][primer_idx + primer_len:]
            trimmed_reads.append((rec.id, trimmed_seq, trimmed_qual))
    
    if len(trimmed_reads) > 0:
        print(f"Total reads before trimming: {len(total_reads)}")
        return trimmed_reads

    else:
        print("Could not trim reads -- double check your primer sequence and offset.")


def write_trimmed_fastq(trimmed_reads, output_path):
    directory = os.path.dirname(output_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(output_path, "w") as f:
        for read_id, seq, qual in trimmed_reads:
            qual_str = ''.join([chr(q + 33) for q in qual])
            f.write(f"@{read_id}\n{seq}\n+\n{qual_str}\n")

def extract_insertion_sites(bam_file, strand="f"):
    result = subprocess.run(["samtools", "view", bam_file],
                            capture_output=True, text=True, check=True)
    positions = []
    for line in result.stdout.strip().split("\n"):
        if not line or line.startswith('@'):
            continue
        fields = line.split("\t")      # 20250905 one base off FIX
        pos1 = int(fields[3])          # SAM POS is 1-based
        start0 = pos1 - 1              # convert to 0-based for fix
        read_len = len(fields[9])
        if strand == "f":
            positions.append(start0)
        else:
            positions.append(start0 + read_len - 1)
    return positions
	
def sgRNA_finder(sgRNA, genome, PAM_len, natural_system):
    sgRNA = sgRNA.upper()
    genome = genome.upper()
    
    if natural_system:
        if genome.find(sgRNA) != -1:
            return genome.find(sgRNA), genome.find(sgRNA), genome.find(sgRNA) - PAM_len, False
        elif genome.find(str(Seq(sgRNA).reverse_complement())) != -1:
            return genome.find(str(Seq(sgRNA).reverse_complement())), genome.find(str(Seq(sgRNA).reverse_complement()))  + len(sgRNA), genome.find(str(Seq(sgRNA).reverse_complement())) + len(sgRNA) + PAM_len, True
        else:
            print("sgRNA Sequence Not Found in Genome")
    
    else:
        if genome.find(sgRNA) != -1:
            return genome.find(sgRNA), genome.find(sgRNA) + len(sgRNA), genome.find(sgRNA) + len(sgRNA) + PAM_len, False
        elif genome.find(str(Seq(sgRNA).reverse_complement())) != -1:
            return genome.find(str(Seq(sgRNA).reverse_complement())), genome.find(str(Seq(sgRNA).reverse_complement())), genome.find(str(Seq(sgRNA).reverse_complement())) - PAM_len, True
        else:
            print("sgRNA Sequence Not Found in Genome")


def save_to_temp(uploadedFile,tmpdir):
#	with tempfile.NamedTemporaryFile(dir=tmpdir.name,delete=False, suffix=".fasta") as tmp_file:
#	directory = os.path.dirname(tmpdir)
	if not os.path.exists(tmpdir):
		os.makedirs(tmpdir)
	with tempfile.NamedTemporaryFile(dir=tmpdir,delete=False, suffix=".fasta") as tmp_file:
		tmp_file.write(uploadedFile.getvalue())
		temp_path = tmp_file.name  # Full path to the saved file
		print("temporary file: " + temp_path)
	return temp_path

def integration_finder(insertion_idx_list, ccdb_gene, ccdb_start_site, PAM_start_site, PAM_end_site, natural_system, is_rc):
    dist_list = []

    if natural_system:
        if is_rc:
            def dist_calc(idx): return (PAM_start_site - idx)

        else:
            def dist_calc(idx): return (idx - PAM_start_site)
    
    else:
        if is_rc:
            def dist_calc(idx): return (PAM_end_site - idx)

        else:
            def dist_calc(idx): return (idx - PAM_end_site)

    for idx in insertion_idx_list:
        if ccdb_start_site + len(ccdb_gene) >= idx >= ccdb_start_site:
            modified_idx = dist_calc(idx)
            dist_list.append(modified_idx)

    return dist_list

#the main function
def process_reads(fastq,primer,first_ten_in_donor,refgenome,bbmap_bin,samtools_bin,path_to_output_dir):
	# 1. Trim reads by primer
	trimmed_reads = trim_reads_by_primer(fastq, primer)

	# 2. Remove donor reads
	donor_removed_reads = [read for read in trimmed_reads if not read[1].startswith(first_ten_in_donor)]
	print(f"Total reads after primer trimming: {len(trimmed_reads)}")
	print(f"Total reads after donor removal: {len(donor_removed_reads)}")

	# 3. Write output FASTQ with original qualities preserved
	trimmed_fastq = Path(f"{path_to_output_dir}/trimmed.fastq")
	write_trimmed_fastq(donor_removed_reads, trimmed_fastq)

	# 4. BBMap indexing and alignment
	subprocess.run([f"{bbmap_bin}/bbmap.sh", f"ref={refgenome}"], check=True)
	mapped_bam = trimmed_fastq.with_suffix(".bam")
	subprocess.run([
	    f"{bbmap_bin}/bbmap.sh", f"in={trimmed_fastq}", f"outm={mapped_bam}",
	    "minid=0.9", "ambig=random"
	    ], check=True)

	# 5. Separate forward/reverse BAM
	forward_bam = mapped_bam.with_name(mapped_bam.stem + "_f.bam")
	reverse_bam = mapped_bam.with_name(mapped_bam.stem + "_r.bam")
	subprocess.run(f"{samtools_bin}/samtools view -F 0x10 -h {mapped_bam} | {samtools_bin}/samtools view -bS - > {forward_bam}", shell=True, check=True)
	subprocess.run(f"{samtools_bin}/samtools view -f 0x10 -h {mapped_bam} | {samtools_bin}/samtools view -bS - > {reverse_bam}", shell=True, check=True)

	# 6. Proper insertion site extraction
	forward_insertions = extract_insertion_sites(forward_bam, samtools_bin, "f")
	reverse_insertions = extract_insertion_sites(reverse_bam, samtools_bin, "r")

	return [forward_insertions,reverse_insertions]


def plot_insertion_histogram(
    forward_positions, reverse_positions,
    region_start=0, region_end=None,
    bins=100, title="Insertion Sites",
    highlight_regions=None,
    filename=None
):

    if region_end is None:
        region_end = max(forward_positions + reverse_positions)

    plt.figure(figsize=(10, 4))
    ax = plt.gca()

    # Plot histograms
    ax.hist(forward_positions, bins=bins, range=(region_start, region_end),
            alpha=0.5, label='Top Strand', color='black')
    ax.hist(reverse_positions, bins=bins, range=(region_start, region_end),
            alpha=0.5, label='Bottom Strand', color='grey')

    # Plot horizontal red bars
    if highlight_regions:
        ymin, _ = ax.get_ylim()
        y = ymin - 0.01 * (ax.get_ylim()[1] - ymin)
        for region in highlight_regions:
            if len(region) == 3:
                x_start, x_end, color = region
            else:
                x_start, x_end = region
                color = 'red'
            ax.hlines(y=2, xmin=x_start, xmax=x_end, colors=color, linewidth=5)

    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("Insertion Count")
    ax.set_title(title)
    ax.set_xlim(region_start, region_end)
    ax.legend()
    plt.tight_layout()

    if filename:
        plt.savefig(filename, dpi=300)
        plt.close()
    else:
        plt.show()
    return plt

def plot_integration_profile(data, title="Integration Distance From PAM", xlabel="Distance From PAM (bp)", ylabel="Count", 
                   xlabel_fontsize=12, ylabel_fontsize=12, title_fontsize=14, 
                   hist_color='grey', edge_color='black', linewidth=5, dpi=300, xlim=(20,100), filename=None):
    
    plt.clf()
    plt.hist(data, edgecolor=edge_color, color=hist_color, linewidth=linewidth)
    plt.title(title, fontsize=title_fontsize)
    plt.xlabel(xlabel, fontsize=xlabel_fontsize)
    plt.ylabel(ylabel, fontsize=ylabel_fontsize)
    
    if xlim:
        plt.xlim(xlim)
    
    plt.gcf().set_dpi(dpi)
    
    plt.tick_params(axis='both', labelsize=10)

    if filename:
        plt.savefig(filename, dpi=dpi)
        plt.close()
    else:
        plt.show()
    return plt
