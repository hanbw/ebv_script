#! /data/home/hanbw/software/anaconda2/bin/python2.7
import os
import sys
import optparse
import time
import threading

LIB_PATH = '/data/home/hanbw/pipeline/Virus_Integration/bin'
sys.path.append(LIB_PATH)
import virus_integration_subs

# PATH to software
Trim_galore_PATH = "~/bin/trim_galore_v0.4.0/trim_galore"
BWA_PATH = "~/bin/bwa-0.7.12/bwa"
Samtools_PATH = "~/bin/samtools-1.2/samtools"
SortSam_PATH = "~/bin/picard-tools-1.119/SortSam.jar"
Picard_PATH = "~/bin/picard/1.127.11/dist/picard.jar"
AdapterRemoval_PATH = "~/bin/AdaptorRemoval/2.1.7/bin/AdapterRemoval"
bedtools_intersectBed_path = "~/bin/bedtools/2.23.0/bin/intersectBed"
bedtools_path = "~/bin/bedtools/2.23.0/bin/bedtools"

# PATH to annotation
hg19_PATH = "~/database/hg19/ucsc_hg19"
hg19_EBV_PATH = "~/reference/combined_ucsc-hg19_EBV"
EBV_NCBI = "~/database/hg19_EBV/ebv_NC_007605.1"
EBV_PATH = "~/database/hg19_EBV/EBV_genome_common"
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
					 help="pathway to combined genome (default: hg19 + EBV)")

optParser.add_option("-m", "--merge", type="choice", dest="merge_reads",
					 choices = ( "True", "T", "False", "F" ), default="False",
					 help="merge paired reads into a long read (True/False, default: False)")

optParser.add_option("-t", "--thread", type="string", dest="thread", default="8",
					 help="thread to use (default: 8)")

optParser.add_option("-r", "--run_merging_result", type="choice", dest="run_merge",
					 choices = ("True", "T", "False", "F"), default="True",
					 help="Run steps to merge result (default: True)")

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
if os.path.exists(prefix+"_Vig_step01.sig"):
	sys.stdout.write("Step 1 result files exist, pass.\n")
	fastq_clean_1 = opts.fastq1.split(".")[0] + "_val_1.fq"
	fastq_clean_2 = opts.fastq2.split(".")[0] + "_val_2.fq"

else:
	fastq_clean_1, fastq_clean_2 = virus_integration_subs.data_clean(Trim_galore_PATH, opts.fastq1, opts.fastq2)
	print opts.fastq1 + " and " + opts.fastq2
	with open(prefix+"_Vig_step01.sig", "w") as sig_file:
		current_time = time.strftime('%Y-%m-%d %H:%M:%S')
		sig_file.write(current_time)

# Step 2. Mapping to hg19
current_time = time.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write("Step 2. Mapping to hg19 " + current_time + "\n")
if os.path.exists(prefix+"_Vig_step02.sig"):
	sys.stdout.write("Step 2 result files exist, pass. \n")
	fastq_hg19_1 = fastq_clean_1.split(".")[0] + "_read1.hg19.fq"
	fastq_hg19_2 = fastq_clean_2.split(".")[0] + "_read2.hg19.fq"
else:
	outsam_hg19_1 = prefix + "_read1.hg19.sam"
	outsam_hg19_2 = prefix + "_read2.hg19.sam"
	postfix = ".hg19.fq"
	virus_integration_subs.mapping_single_end(BWA_PATH, prefix, hg19_PATH, fastq_clean_1, outsam_hg19_1, opts.thread)
	virus_integration_subs.mapping_single_end(BWA_PATH, prefix, hg19_PATH, fastq_clean_2, outsam_hg19_2, opts.thread)
	fastq_hg19_1, fastq_hg19_2 = virus_integration_subs.removing_mapped_intersection(outsam_hg19_1, outsam_hg19_2, fastq_clean_1, fastq_clean_2, postfix)
	os.system("rm -f %s %s %s %s" % (fastq_clean_1, fastq_clean_2, outsam_hg19_1, outsam_hg19_2))
	with open(prefix+"_Vig_step02.sig", "w") as sig_file:
		current_time = time.strftime('%Y-%m-%d %H:%M:%S')
		sig_file.write(current_time)

# Step 3. Mapping to EBV
current_time = time.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write("Step 3. Mapping to EBV " + current_time +"\n")
if os.path.exists(prefix + "_Vig_step03.sig"):
	sys.stdout.write("Step 3 result files exist, pass. \n")
	fastq_EBV_1 = fastq_hg19_1.split(".")[0] + "_read1.EBV.fq"
	fastq_EBV_2 = fastq_hg19_2.split(".")[0] + "_read2.EBV.fq"
else:
	outsam_EBV_1 = prefix + "_read1.EBV.sam"
	outsam_EBV_2 = prefix + "_read2.EBV.sam"
	postfix = ".EBV.fq"
	virus_integration_subs.mapping_single_end(BWA_PATH, prefix, EBV_PATH, fastq_hg19_1, outsam_EBV_1, opts.thread)
	virus_integration_subs.mapping_single_end(BWA_PATH, prefix, EBV_PATH, fastq_hg19_2, outsam_EBV_2, opts.thread)
	fastq_EBV_1, fastq_EBV_2 = virus_integration_subs.removing_mapped_intersection(outsam_EBV_1, outsam_EBV_2, fastq_hg19_1, fastq_hg19_2, postfix)
	os.system("rm -f %s %s %s %s" % (fastq_hg19_1, fastq_hg19_2, outsam_EBV_1, outsam_EBV_2))
	with open(prefix+"_Vig_step03.sig", "w") as sig_file:
		current_time = time.strftime('%Y-%m-%d %H:%M:%S')
		sig_file.write(current_time)

# Step 4. Collapsing reads and remapping
current_time = time.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write("Step 4. Collapsing reads " + current_time +"\n")
if os.path.exists(prefix + "_Vig_step04.sig"):
	sys.stdout.write("Step 4 result files exist, pass. \n")
	outsam_merge_reads_ebv = prefix + "_merge_ebv.sam"
	outsam_merge_reads_hg19 = prefix + "_merge_hg19.sam"
	outsam_single_reads_ebv = prefix + "_single_ebv.sam"
	outsam_single_reads_hg19 = prefix + "_single_hg19.sam"
	if opts.merge_reads == "True" or "T":
		merge_reads, single_reads_1, single_reads_2 = (prefix + ".collapsed", prefix + ".pair1.truncated", prefix + ".pair2.truncated")
		outsam_merge_reads, outsam_single_read = (prefix + "_merge.sam", prefix + "_single.sam")
	else:
		merge_reads, single_reads_1, single_reads_2 = ("", fastq_EBV_1, fastq_EBV_2)
		outsam_merge_reads, outsam_single_read = ("", prefix + "single.sam")
else:
	if opts.merge_reads == "True" or "T":  # reads go through insert fragments (e.g. MISEQ for targeting sequencing)
		merge_reads, single_reads_1, single_reads_2 = virus_integration_subs.collapsing_reads(AdapterRemoval_PATH, prefix, fastq_EBV_1, fastq_EBV_2, opts.thread)
		outsam_merge_reads_ebv = prefix + "_merge_ebv.sam"
		outsam_merge_reads_hg19 = prefix + "_merge_hg19.sam"
		outsam_single_reads_ebv = prefix + "_single_ebv.sam"
		outsam_single_reads_hg19 = prefix + "_single_hg19.sam"

		virus_integration_subs.mapping_single_end(BWA_PATH, prefix, EBV_PATH, merge_reads, outsam_merge_reads_ebv, opts.thread)
		virus_integration_subs.mapping_single_end(BWA_PATH, prefix, hg19_PATH, merge_reads, outsam_merge_reads_hg19, opts.thread)
		virus_integration_subs.mapping_paired_ends(BWA_PATH, prefix, EBV_PATH, single_reads_1, single_reads_2, outsam_single_reads_ebv, opts.thread)
		virus_integration_subs.mapping_paired_ends(BWA_PATH, prefix, hg19_PATH, single_reads_1, single_reads_2, outsam_single_reads_hg19, opts.thread)

		os.system("rm -f %s %s" % (fastq_EBV_1, fastq_EBV_2))

	else:  # reads not go through insert fragments (e.g. HISEQ for WGS)
		outsam_single_reads_ebv = prefix + "_single_ebv.sam"
		outsam_single_reads_hg19 = prefix + "_single_hg19.sam"
		outsam_merge_reads_hg19, outsam_merge_reads_ebv, merge_reads = ("","", "")
		single_reads_1, single_reads_2 = fastq_EBV_1, fastq_EBV_2
		virus_integration_subs.mapping_paired_ends(BWA_PATH, prefix, EBV_PATH, fastq_EBV_1, fastq_EBV_2, outsam_single_reads_ebv, opts.thread)
		virus_integration_subs.mapping_paired_ends(BWA_PATH, prefix, hg19_PATH, fastq_EBV_1, fastq_EBV_2, outsam_single_reads_hg19, opts.thread)

	with open(prefix + "_Vig_step04.sig", "w") as sig_file:
		current_time = time.strftime('%Y-%m-%d %H:%M:%S')
		sig_file.write(current_time)

# Step 5. Picking clipping reads
current_time = time.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write("Step 5. Picking clipping reads " + current_time + "\n")
merge_clipping_reads_fq = merge_reads + ".clipping.fq"
single_clipping_reads_fq_1 = prefix + ".clipping_read1.fq"
single_clipping_reads_fq_2 = prefix + ".clipping_read2.fq"
if os.path.exists(prefix + "_Vig_step05.sig"):
	sys.stdout.write("Step 5 result files exist, pass. \n")
else:
	if opts.merge_reads == "True" or "T":
		merge_clipping_reads = virus_integration_subs.analyze_integration_reads(outsam_merge_reads_hg19, outsam_merge_reads_ebv)
		virus_integration_subs.delete_reads_ori(merge_reads, merge_clipping_reads_fq, merge_clipping_reads)
		del merge_clipping_reads

		single_clipping_reads = virus_integration_subs.analyze_integration_reads(outsam_single_reads_hg19, outsam_single_reads_ebv)
		virus_integration_subs.delete_reads_ori(single_reads_1, single_clipping_reads_fq_1, single_clipping_reads)
		virus_integration_subs.delete_reads_ori(single_reads_2, single_clipping_reads_fq_2, single_clipping_reads)
		del single_clipping_reads
		os.system("rm -f %s %s %s %s" % (outsam_merge_reads_hg19, outsam_merge_reads_ebv, outsam_single_reads_hg19, outsam_single_reads_ebv))

	else:
		single_clipping_reads = virus_integration_subs.analyze_integration_reads(outsam_single_reads_hg19, outsam_single_reads_ebv)
		virus_integration_subs.delete_reads_ori(single_reads_1, single_clipping_reads_fq_1, single_clipping_reads)
		virus_integration_subs.delete_reads_ori(single_reads_2, single_clipping_reads_fq_2, single_clipping_reads)
		del single_clipping_reads
		os.system("rm -f %s %s" % (outsam_single_reads_hg19, outsam_single_reads_ebv))

	os.system("rm -f %s %s %s " % (merge_reads, single_reads_1, single_reads_2))

	with open(prefix + "_Vig_step05.sig", "w") as sig_file:
		current_time = time.strftime('%Y-%m-%d %H:%M:%S')
		sig_file.write(current_time)

# Step 6. Remapping clipping reads
current_time = time.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write("Step 6. Remapping clipping reads " + current_time +"\n")
all_pos_list = prefix + "_all_int_pos.txt"
if os.path.exists(prefix + "_Vig_step06.sig"):
	sys.stdout.write("Step 6 result files exist, pass. \n")
else:
	outsam_clipping_merge_reads_hg19 = prefix + ".merge_clipping_hg19.sam"
	outsam_clipping_merge_reads_ebv = prefix + ".merge_clipping_ebv.sam"
	outsam_clipping_single_reads_r1_hg19 = prefix + "_r1.single_clipping_r1_hg19.sam"
	outsam_clipping_single_reads_r1_ebv = prefix + "_r1.single_clipping_r1_ebv.sam"
	outsam_clipping_single_reads_r2_hg19 = prefix + "_r2.single_clipping_r2_hg19.sam"
	outsam_clipping_single_reads_r2_ebv = prefix + "_r2.single_clipping_r2_ebv.sam"

	if opts.merge_reads == "True" or "T":
		virus_integration_subs.mapping_single_end(BWA_PATH, prefix, hg19_PATH, merge_clipping_reads_fq, outsam_clipping_merge_reads_hg19, opts.thread)
		virus_integration_subs.mapping_single_end(BWA_PATH, prefix, EBV_NCBI, merge_clipping_reads_fq, outsam_clipping_merge_reads_ebv, opts.thread)

		merge_pos_list = virus_integration_subs.output_merged_integration_reads(outsam_clipping_merge_reads_hg19, outsam_clipping_merge_reads_ebv)

	virus_integration_subs.mapping_single_end(BWA_PATH, prefix, hg19_PATH, single_clipping_reads_fq_1, outsam_clipping_single_reads_r1_hg19, opts.thread)
	virus_integration_subs.mapping_single_end(BWA_PATH, prefix, EBV_NCBI, single_clipping_reads_fq_1, outsam_clipping_single_reads_r1_ebv, opts.thread)
	virus_integration_subs.mapping_single_end(BWA_PATH, prefix, hg19_PATH, single_clipping_reads_fq_2, outsam_clipping_single_reads_r2_hg19, opts.thread)
	virus_integration_subs.mapping_single_end(BWA_PATH, prefix, EBV_NCBI, single_clipping_reads_fq_2, outsam_clipping_single_reads_r2_ebv, opts.thread)

	single_pos_list, ebv_r1_matched, ebv_r2_matched, human_r1_matched, human_r2_matched = virus_integration_subs.output_paired_non_clipping_integration_reads(outsam_clipping_single_reads_r1_hg19, outsam_clipping_single_reads_r2_hg19, outsam_clipping_single_reads_r1_ebv, outsam_clipping_single_reads_r2_ebv)

	single_r1_pos_list = virus_integration_subs.output_merged_integration_reads(outsam_clipping_single_reads_r1_hg19, outsam_clipping_single_reads_r1_ebv)
	single_r2_pos_list = virus_integration_subs.output_merged_integration_reads(outsam_clipping_single_reads_r2_hg19, outsam_clipping_single_reads_r2_ebv)

	library_max_size = 2000  # sequencing library insertion size
	virus_integration_subs.analyse_paired_one_clipping_integration_reads(single_r1_pos_list, ebv_r2_matched, human_r2_matched, library_max_size)
	virus_integration_subs.analyse_paired_one_clipping_integration_reads(single_r2_pos_list, ebv_r1_matched, human_r1_matched, library_max_size)
	single_r1_pos_list_p = single_r1_pos_list.split(".")[0] + "_integration_single_pos.txt"
	single_r2_pos_list_p = single_r2_pos_list.split(".")[0] + "_integration_single_pos.txt"

	virus_integration_subs.analyse_paired_two_clipping_integration_reads(single_r1_pos_list, single_r2_pos_list, merge_pos_list)

	os.system("cat %s %s %s > %s" % (merge_pos_list, single_r1_pos_list_p, single_r2_pos_list_p, all_pos_list))
	os.system("rm -f %s %s %s %s %s %s" % (outsam_clipping_merge_reads_hg19, outsam_clipping_merge_reads_ebv, outsam_clipping_single_reads_r1_hg19, outsam_clipping_single_reads_r1_ebv, outsam_clipping_single_reads_r2_hg19, outsam_clipping_single_reads_r2_ebv))
	os.system("rm -f %s %s %s" % (merge_pos_list, single_r1_pos_list, single_r2_pos_list))

	with open(prefix + "_Vig_step06.sig", "w") as sig_file:
		current_time = time.strftime('%Y-%m-%d %H:%M:%S')
		sig_file.write(current_time)

# Step 7. Repeatmasker and GC filter
current_time = time.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write("Step 7. Repeatmasker and GC filter " + current_time + "\n")
all_pos_list_filter = prefix + "_all_int_pos_filter.txt"
if os.path.exists(prefix + "_Vig_step07.sig"):
	sys.stdout.write("Step 7 result files exist, pass. \n")
	pos_dict = virus_integration_subs.read_pos_list(all_pos_list_filter)
else:
	# read position list
	pos_dict = virus_integration_subs.read_pos_list(all_pos_list)

	# bedtools filtering repeatmasker region
	bedtools_tmp_name = prefix + "_repeatmasker_tmp.bed"
	pos_dict = virus_integration_subs.remove_repeat_masker_bedtools(pos_dict, bedtools_tmp_name, bedtools_intersectBed_path, repeatmasker_annotation_file_path)

	# bedtools filtering GC percentage
	bedtools_tmp_name = prefix + "_GC_tmp.bed"
	pos_dict = virus_integration_subs.gc_adjust(pos_dict, bedtools_tmp_name, bedtools_path, hg19_EBV_PATH)

	# output result
	with open(all_pos_list_filter, "w") as pos_list_filter:
		pos_list_filter.write("read_ID\tebv_start\tebv_gap\tebv_strand\thg19_gap\thg19_end\thg19_strand\thg_19_chr\toverlap\tother_info\n")
		for bam_id, pos_list in pos_dict.iteritems():
			pos_list_filter.write("\t".join(pos_list) + "\n")

	with open(prefix + "_Vig_step07.sig", "w") as sig_file:
		current_time = time.strftime('%Y-%m-%d %H:%M:%S')
		sig_file.write(current_time)

# Step 8. Group positions into windows and write positions with read counts
if opts.run_merge == "True":
	current_time = time.strftime('%Y-%m-%d %H:%M:%S')
	sys.stdout.write("Step 8. Grouping positions into windows and writing positions with read counts " + current_time + "\n")
	if os.path.exists(prefix + "_Vig_step08.sig"):
		sys.stdout.write("Step 8 result files exist, pass. \n")
	else:
		output_file_name_EBV = prefix + "_int_distribution_EBV.txt"
		output_file_name_hg19 = prefix + "_int_distribution_hg19.txt"

		virus_integration_subs.group_pos(pos_dict, output_file_name_EBV, output_file_name_hg19)
		with open(prefix + "_Vig_step08.sig", "w") as sig_file:
			current_time = time.strftime('%Y-%m-%d %H:%M:%S')
			sig_file.write(current_time)

