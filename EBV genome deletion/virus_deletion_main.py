#! /data/home/hanbw/software/anaconda2/bin/python2.7
import os
import sys
import optparse
import time
import threading

LIB_PATH = './'
sys.path.append(LIB_PATH)
import virus_deletion_subs

# PATH to software
Trim_galore_PATH = "~/bin/trim_galore"
BWA_PATH = "~/bin/bwa-0.7.12/bwa"
Samtools_PATH = "~/bin/samtools-1.2/samtools"
SortSam_PATH = "~/bin/picard-tools-1.119/SortSam.jar"
Picard_PATH = "~/bin/picard/1.127.11/dist/picard.jar"
AdapterRemoval_PATH = "~/bin/AdaptorRemoval/2.1.7/bin/AdapterRemoval"
HISAT_PATH = "~/bin/hisat2"
bedtools_intersectBed_path = "~/bin/bedtools/2.23.0/bin/intersectBed"
bedtools_path = "~/bin/bedtools/2.23.0/bin/bedtools"

# PATH to annotation
hg19_PATH = "~/database/hg19/ucsc_hg19"
hg19_EBV_PATH = "~/database/hg19_EBV/hg19_ebv_all"
EBV_NCBI = "~/database/hg19_EBV/ebv_NC_007605.1"
EBV_PATH = "~/hg19_EBV/EBV_genome_common"
repeatmasker_annotation_file_path = "~/database/hg19/hg19.repeatmasker.anno"

# get arguments
optParser = optparse.OptionParser(

	usage="%prog  -g [host_genome] -v [virus_genome] -s [virus_standard_genome] -c [combined_genome] -1 <fastq_file_1> -2 <fastq_file_2>",

	description=
	"This script analyse virus integration sites in both virus and host genome",

	epilog=
	"Written by Han Bowei (hanbw@sysucc.org), SYSUCC, Dec 2016 \n"
)

optParser.add_option("-1", "--fastq1", type="string", dest="fastq1", default="",
					 help="input fastq file 1")

optParser.add_option("-2", "--fastq2", type="string", dest="fastq2", default="",
					 help="input fastq file 2")

optParser.add_option("-g", "--host_genome", type="string", dest="host_genome", default=hg19_PATH,
					 help="pathway to host genome (default: hg19)")

optParser.add_option("-v", "--virus_genome", type="string", dest="virus_genome", default=EBV_PATH,
					 help="pathway to virus genome (default: EBV)")

optParser.add_option("-s", "--virus_standard_genome", type="string", dest="virus_standard_genome", default=EBV_NCBI,
					 help="pathway to standard virus genome (default: EBV_NC_007605.1)")

optParser.add_option("-c", "--combined_genome", type="string", dest="combined_genome", default=hg19_EBV_PATH,
					 help="pathway to combined genome (default: hg19 + EBV)")

optParser.add_option("-p", "--project_prefix", type="string", dest="prefix", default="",
					 help="sample prefix")

optParser.add_option("-m", "--merge", type="choice", dest="merge_reads",
					 choices = ( "True", "T", "False", "F" ), default="False",
					 help="merge paired reads into a long read (True/False, default: False)")

optParser.add_option("-t", "--thread", type="string", dest="thread", default="8",
					 help="thread to use (default: 8)")

optParser.add_option("-r", "--run_merging_result", type="choice", dest="run_merge",
					 choices = ("True", "T", "False", "F"), default="True",
					 help="Run steps to merge result (default: True)")

optParser.add_option("-i", "--insertion_length", type="int", dest="insertion_length",
					 default=450, help="Max insertion length of sequencing library (default: 450)")

# optParser.add_option("-t", "--type", type="choice", dest="vcf_type",
# 					 choices = ( "RNA", "DNA" ), default="RNA",
# 					 help="VCF file generate from RNA-seq or DNA-seq (default: RNA)")


if len(sys.argv) == 1:
	optParser.print_help()
	sys.exit(1)

(opts, args) = optParser.parse_args()

if opts.fastq1 == "" or opts.fastq2 == "":
	if opts.prefix != "" and os.path.exists(opts.prefix + "step01.sig"):
		opts.fastq1 = opts.prefix + "_R1.fastq.gz"
		opts.fastq2 = opts.prefix + "_R2.fastq.gz"
	else:
		sys.stderr.write(sys.argv[0] + ": Error: Please provide fastq file or file prefix.\n")
		sys.stderr.write("  Call with '-h' to get usage information.\n")
		sys.exit(1)

if opts.prefix == "":
	prefix = opts.fastq1.split("/")[-1].split("_R1")[0]
else:
	prefix = opts.prefix

# Step 1. Data cleaning
current_time = time.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write("Step 1. Data cleaning " + current_time +"\n")
if os.path.exists(prefix+"_Vdl_step01.sig"):
	sys.stdout.write("Step 1 result files exist, pass.\n")
	fastq_clean_1 = opts.fastq1.split(".")[0] + "_val_1.fq"
	fastq_clean_2 = opts.fastq2.split(".")[0] + "_val_2.fq"
else:
	fastq_clean_1, fastq_clean_2 = virus_deletion_subs.data_clean(Trim_galore_PATH, opts.fastq1, opts.fastq2)
	print opts.fastq1 + " and " + opts.fastq2
	with open(prefix+"_Vdl_step01.sig", "w") as sig_file:
		current_time = time.strftime('%Y-%m-%d %H:%M:%S')
		sig_file.write(current_time)

# Step 2. Mapping to hg19 and EBV
current_time = time.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write("Step 2. Mapping to EBV " + current_time +"\n")
if os.path.exists(prefix + "_Vdl_step02.sig"):
	sys.stdout.write("Step 2 result files exist, pass. \n")
	fastq_EBV_1 = fastq_clean_1.split(".")[0] + "_read1.EBV.fq"
	fastq_EBV_2 = fastq_clean_2.split(".")[0] + "_read2.EBV.fq"
else:
	outsam_hg19_1 = prefix + "_read1.hg19.sam"
	outsam_hg19_2 = prefix + "_read2.hg19.sam"
	postfix = ".hg19.fq"
	virus_deletion_subs.mapping_single_end(BWA_PATH, prefix, hg19_PATH, fastq_clean_1, outsam_hg19_1, opts.thread)
	virus_deletion_subs.mapping_single_end(BWA_PATH, prefix, hg19_PATH, fastq_clean_2, outsam_hg19_2, opts.thread)
	fastq_hg19_1, fastq_hg19_2 = virus_deletion_subs.removing_mapped_intersection(outsam_hg19_1, outsam_hg19_2, fastq_clean_1, fastq_clean_2, postfix)

	outsam_EBV_1 = prefix + "_read1.EBV.sam"
	outsam_EBV_2 = prefix + "_read2.EBV.sam"
	postfix = ".EBV.fq"
	virus_deletion_subs.mapping_single_end(BWA_PATH, prefix, EBV_PATH, fastq_hg19_1, outsam_EBV_1, opts.thread)
	virus_deletion_subs.mapping_single_end(BWA_PATH, prefix, EBV_PATH, fastq_hg19_2, outsam_EBV_2, opts.thread)
	fastq_EBV_1, fastq_EBV_2 = virus_deletion_subs.removing_mapped_intersection(outsam_EBV_1, outsam_EBV_2, fastq_hg19_1, fastq_hg19_2, postfix)
	os.system("rm -f %s %s %s %s" % (outsam_EBV_1, outsam_EBV_2, outsam_hg19_1, outsam_hg19_2))
	with open(prefix+"_Vdl_step02.sig", "w") as sig_file:
		current_time = time.strftime('%Y-%m-%d %H:%M:%S')
		sig_file.write(current_time)

# Step 3. Collapsing reads and remapping
current_time = time.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write("Step 3. Collapsing reads " + current_time +"\n")
if os.path.exists(prefix + "_Vdl_step03.sig"):
	sys.stdout.write("Step 3 result files exist, pass. \n")
	outsam_merge_reads = prefix + "_merge.sam"
	outsam_single_read = prefix + "_single.sam"
	if opts.merge_reads == "True" or "T":
		merge_reads, single_reads_1, single_reads_2 = (prefix + ".collapsed", prefix + ".pair1.truncated", prefix + ".pair2.truncated")
		outsam_merge_reads, outsam_single_read = (prefix + "_merge.sam", prefix + "_single.sam")
	else:
		merge_reads, single_reads_1, single_reads_2 = ("", fastq_EBV_1, fastq_EBV_2)
		outsam_merge_reads, outsam_single_read = ("", prefix + "single.sam")

else:
	virus_deletion_subs.fastq_tab(fastq_EBV_1)
	virus_deletion_subs.fastq_tab(fastq_EBV_2)
	try:
		virus_deletion_subs.fastq_tab(fastq_hg19_1)
		virus_deletion_subs.fastq_tab(fastq_hg19_2)
	except:
		pass

	if opts.merge_reads == "True" or "T":  # reads go through insert fragments (e.g. MISEQ for targeting sequencing)
		merge_reads, single_reads_1, single_reads_2 = virus_deletion_subs.collapsing_reads(AdapterRemoval_PATH, prefix, fastq_EBV_1, fastq_EBV_2, opts.thread)
		outsam_merge_reads = prefix + "_merge.sam"
		outsam_single_read = prefix + "_single.sam"

		virus_deletion_subs.hisat_mapping_single_end(HISAT_PATH, EBV_PATH, merge_reads, outsam_merge_reads, opts.thread)
		virus_deletion_subs.hisat_mapping_paired_ends(HISAT_PATH, EBV_PATH, single_reads_1, single_reads_2, outsam_single_read, opts.thread)

	else:  # reads not go through insert fragments (e.g. HISEQ for WGS)
		outsam_merge_reads, merge_reads = ("", "")
		outsam_single_read = prefix + "single.sam"
		single_reads_1,single_reads_2  = fastq_EBV_1, fastq_EBV_2
		virus_deletion_subs.hisat_mapping_paired_ends(HISAT_PATH, EBV_PATH, fastq_EBV_1, fastq_EBV_2, outsam_single_read, opts.thread)

	with open(prefix + "_Vdl_step03.sig", "w") as sig_file:
		current_time = time.strftime('%Y-%m-%d %H:%M:%S')
		sig_file.write(current_time)

# Step 4. Picking gap reads
current_time = time.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write("Step 4. Picking gap reads " + current_time +"\n")
merge_clipping_reads_fq = merge_reads + ".clipping.fq"
single_clipping_reads_fq_1 = prefix + ".clipping_read1.fq"
single_clipping_reads_fq_2 = prefix + ".clipping_read2.fq"
if os.path.exists(prefix + "_Vdl_step04.sig"):
	sys.stdout.write("Step 4 result files exist, pass. \n")
else:
	if opts.merge_reads == "True" or "T":
		merge_clipping_reads = virus_deletion_subs.analyze_merged_deletion_reads(outsam_merge_reads)
		virus_deletion_subs.delete_reads_ori(merge_reads, merge_clipping_reads_fq, merge_clipping_reads)
		#del merge_clipping_reads

		single_clipping_reads = virus_deletion_subs.analyze_single_deletion_reads(outsam_single_read)
		virus_deletion_subs.delete_reads_ori(single_reads_1, single_clipping_reads_fq_1, single_clipping_reads)
		virus_deletion_subs.delete_reads_ori(single_reads_2, single_clipping_reads_fq_2, single_clipping_reads)
		#del single_clipping_reads

		# os.system("rm -f %s %s " % (outsam_merge_reads, outsam_single_read))

	else:
		single_clipping_reads = virus_deletion_subs.analyze_single_deletion_reads(outsam_single_read)
		virus_deletion_subs.delete_reads_ori(single_reads_1, single_clipping_reads_fq_1, single_clipping_reads)
		virus_deletion_subs.delete_reads_ori(single_reads_2, single_clipping_reads_fq_2, single_clipping_reads)
		#del single_clipping_reads
		# os.system("rm -f %s" % (outsam_single_read))

#	os.system("rm -f %s %s %s " % (merge_reads, single_reads_1, single_reads_2))

	with open(prefix + "_Vdl_step04.sig", "w") as sig_file:
		current_time = time.strftime('%Y-%m-%d %H:%M:%S')
		sig_file.write(current_time)

# Step 5. Remapping gap reads
current_time = time.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write("Step 5. Remapping gap reads " + current_time +"\n")
all_pos_list = prefix + "_all_del_pos.txt"
if os.path.exists(prefix + "_Vdl_step05.sig"):
	sys.stdout.write("Step 5 result files exist, pass. \n")
else:
	outsam_clipping_merge_reads = prefix + ".merge_clipping.sam"
	outsam_clipping_single_reads = prefix + ".single_clipping.sam"
	library_max_size = 2000  # sequencing library insertion size
	if opts.merge_reads == "True" or "T":
		virus_deletion_subs.hisat_mapping_single_end(HISAT_PATH, EBV_NCBI, merge_clipping_reads_fq, outsam_clipping_merge_reads, opts.thread)
		virus_deletion_subs.hisat_mapping_paired_ends(HISAT_PATH, EBV_NCBI, single_clipping_reads_fq_1, single_clipping_reads_fq_2, outsam_clipping_single_reads, opts.thread)
		merge_pos_list = virus_deletion_subs.output_merged_deletion_reads(outsam_clipping_merge_reads)
		single_pos_list = virus_deletion_subs.output_single_deletion_reads(outsam_clipping_single_reads, library_max_size)
		os.system("cat %s %s > %s" % (single_pos_list, merge_pos_list, all_pos_list))
		os.system("rm -f %s %s %s %s " % (outsam_clipping_merge_reads, outsam_clipping_single_reads, single_pos_list, merge_pos_list))

	else:
		virus_deletion_subs.hisat_mapping_paired_ends(HISAT_PATH, EBV_NCBI, single_clipping_reads_fq_1, single_clipping_reads_fq_2, outsam_clipping_single_reads, opts.thread)
		single_pos_list = virus_deletion_subs.output_single_deletion_reads(outsam_clipping_single_reads, library_max_size)
		os.system("mv %s %s" % (single_pos_list, all_pos_list))
		os.system("rm -f %s " % (outsam_clipping_single_reads))

	with open(prefix + "_Vdl_step05.sig", "w") as sig_file:
		current_time = time.strftime('%Y-%m-%d %H:%M:%S')
		sig_file.write(current_time)

# Step 6. Repeatmasker and GC filter
current_time = time.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write("Step 6. Repeatmasker and GC filter " + current_time +"\n")
if os.path.exists(prefix + "_Vdl_step06.sig"):
	sys.stdout.write("Step 6 result files exist, pass. \n")
	all_pos_list_filter = prefix + "_all_del_pos_filter.txt"
	pos_dict = virus_deletion_subs.read_pos_list(all_pos_list_filter)
else:
	# read position list
	pos_dict = virus_deletion_subs.read_pos_list(all_pos_list)

	# bedtools filtering repeatmasker region
	bedtools_tmp_name = prefix + "_repeatmasker_tmp.bed"
	pos_dict = virus_deletion_subs.remove_repeat_masker_bedtools(pos_dict, bedtools_tmp_name, bedtools_intersectBed_path, repeatmasker_annotation_file_path)

	# bedtools filtering GC percentage
	bedtools_tmp_name = prefix + "_GC_tmp.bed"
	pos_dict = virus_deletion_subs.gc_adjust(pos_dict, bedtools_tmp_name, bedtools_path, EBV_NCBI)

	# output result
	all_pos_list_filter = prefix + "_all_del_pos_filter.txt"
	with open(all_pos_list_filter, "w") as pos_list_filter:
		pos_list_filter.write("read_ID\tread_start\tgap_start\tgap_end\tread_end\tdel_size\tother_info\n")
		for bam_id, pos_list in pos_dict.iteritems():
			pos_list_filter.write("\t".join(pos_list) + "\n")

	with open(prefix + "_Vdl_step06.sig", "w") as sig_file:
		current_time = time.strftime('%Y-%m-%d %H:%M:%S')
		sig_file.write(current_time)

# Step 7. Group positions into windows and write positions with read counts
if opts.run_merge == "True":
	current_time = time.strftime('%Y-%m-%d %H:%M:%S')
	sys.stdout.write("Step 7. Grouping positions into windows and writing positions with read counts " + current_time +"\n")
	if os.path.exists(prefix + "_Vdl_step07.sig"):
		sys.stdout.write("Step 7 result files exist, pass. \n")
	else:
		output_file_name = prefix + "_del_distribution.txt"
		virus_deletion_subs.group_pos(pos_dict, output_file_name)

		with open(prefix + "_Vdl_step07.sig", "w") as sig_file:
			current_time = time.strftime('%Y-%m-%d %H:%M:%S')
			sig_file.write(current_time)


