#! /usr/bin/env python
import os
import sys
import threading
import Queue
import HTSeq
import itertools

queue_1 = Queue.Queue(maxsize=3)

def data_clean(Trim_galore_PATH, fastq1, fastq2):
	adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
	adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
	min_length = "60"
	os.system("%s --phred33 -a %s -a2 %s --length %s -q 20 --paired --dont_gzip --stringency 3 %s %s"
			  % (Trim_galore_PATH, adapter1, adapter2, min_length, fastq1, fastq2))
	fastq1_out = fastq1.split(".")[0] + "_val_1.fq"
	fastq2_out = fastq2.split(".")[0] + "_val_2.fq"
	return fastq1_out, fastq2_out

def mapping_single_end(BWA_PATH, prefix, genome_PATH, fastq, outsam, thread):
	os.system("%s mem -t %s -M -R \"@RG\\tID:%s\\tLB:%s\\tPL:ILLUMINA\\tSM:%s\" %s %s > %s"
			  % (BWA_PATH, thread, prefix, prefix, prefix, genome_PATH, fastq, outsam))


def hisat_mapping_single_end(HISAT_PATH, genome_PATH, fastq, outsam, thread):
	os.system("%s -p %s -x %s -U %s -S %s --no-unal" % (HISAT_PATH, thread, genome_PATH, fastq, outsam))


def hisat_mapping_paired_ends(HISAT_PATH, genome_PATH, fastq1, fastq2, outsam, thread):
	os.system("%s -p %s -x %s -1 %s -2 %s -S %s --no-unal" % (HISAT_PATH, thread, genome_PATH, fastq1, fastq2, outsam))


def analyse_reads(sam_file):
	with open(sam_file, "r") as sam_lines:
		mapped_reads = set()
		for sam_line in sam_lines:
			if sam_line.startswith("@"):
				pass
			elif ("H" in sam_line.split("\t")[5]) or ("S" in sam_line.split("\t")[5]) or ("*" in sam_line.split("\t")[5]) or ("D" in sam_line.split("\t")[5]):
				pass
			else:
				# mapped_reads.add("@" + sam_line.split("\t")[0])
				mapped_reads.add(sam_line.split("\t")[0])
	queue_1.put((mapped_reads, sam_file))
	# return mapped_reads


def get_all_reads(fastq_in):
	with open(fastq_in, "r") as fastq_in_lines:
		seq_header = fastq_in_lines.readline().split(":")[0]

	with open(fastq_in, "r") as fastq_in_lines:
		all_reads = set()
		for fastq_line in fastq_in_lines:
			if fastq_line.startswith(seq_header):
				fastq_line = fastq_line.replace("@", "")
				all_reads.add(fastq_line.split(" ")[0].split("\t")[0])
	return all_reads


def delete_reads_ori(fastq_in, fastq_out, reads_to_left):
	with open(fastq_in, "r") as fastq_in_lines:
		seq_header = fastq_in_lines.readline().split(":")[0]
		# print seq_header
	with open(fastq_in, "r") as fastq_in_lines:
		with open(fastq_out, "w") as fastq_out_file:
			for fastq_line in fastq_in_lines:
				if fastq_line.startswith(seq_header):
					if fastq_line.split(" ")[0] in reads_to_left:
						fastq_out_file.write(fastq_line)
						fastq_out_file.write(fastq_in_lines.next())
						fastq_out_file.write(fastq_in_lines.next())
						fastq_out_file.write(fastq_in_lines.next())


def fastq_tab(fastq_in):
	with open(fastq_in, "r") as fastq_in_lines:
		seq_header = fastq_in_lines.readline().split(":")[0]
	fastq_out = fastq_in + ".tmp"

	with open(fastq_in, "r") as fastq_in_lines:
		with open(fastq_out, "w") as fastq_out_line:
			for fastq_line in fastq_in_lines:
				if fastq_line.startswith(seq_header):
					fastq_line = fastq_line.replace("\t", " ")
				fastq_out_line.write(fastq_line)

	os.system("rm -f %s" % fastq_in)
	os.system("mv %s %s" % (fastq_out, fastq_in))


def delete_reads(fastq_in, fastq_out, reads_to_left):
	reads_to_left_list = fastq_in + ".list"
	with open(reads_to_left_list, "w") as reads_to_left_list_file:
		reads_to_left_list_file.write("\n".join(reads_to_left))
	os.system("seqtk subseq %s %s > %s" %(fastq_in, reads_to_left_list, fastq_out))
	fastq_tab(fastq_out)
	os.system("rm -f %s" % reads_to_left_list)


def removing_mapped_intersection(sam1, sam2, fastq1, fastq2, postfix):
	fastq_rm_1 = fastq1.split(".")[0] + "_read1" + postfix
	fastq_rm_2 = fastq2.split(".")[0] + "_read2" + postfix

	queue_result = list()
	analyse_fastq1 = threading.Thread(target=analyse_reads, args=(sam1,))
	analyse_fastq2 = threading.Thread(target=analyse_reads, args=(sam2,))
	analyse_fastq1.start()
	analyse_fastq2.start()
	analyse_fastq1.join()
	analyse_fastq2.join()
	while not queue_1.empty():
		queue_result.append(queue_1.get())
	for item in queue_result:
		if item[1] == sam1:
			read_ana_1 = item[0]
		elif item[1] == sam2:
			read_ana_2 = item[0]

	# intersection_reads = list(read_ana_1.intersection(read_ana_2))
	intersection_reads = read_ana_1.intersection(read_ana_2)
	read_all = get_all_reads(fastq1)
	read_left = list(read_all.difference(intersection_reads))
	del(read_ana_1, read_ana_2, intersection_reads, read_all)

	delete_fastq1 = threading.Thread(target=delete_reads, args=(fastq1, fastq_rm_1, read_left))
	delete_fastq2 = threading.Thread(target=delete_reads, args=(fastq2, fastq_rm_2, read_left))
	delete_fastq1.start()
	delete_fastq2.start()
	delete_fastq1.join()
	delete_fastq2.join()

	return fastq_rm_1, fastq_rm_2


def collapsing_reads(AdapterRemoval_PATH, prefix, fastq1, fastq2, threads):
	merge_reads = prefix + ".collapsed"
	single_reads_1 = prefix + ".pair1.truncated"
	single_reads_2 = prefix + ".pair2.truncated"
	os.system("%s --file1 %s --file2 %s --collapse --basename %s --threads %s --adapter1 TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG --adapter2 GTCTCGTGGGCTCGGAGATGTGTATAAGAGAGAG --qualitymax 50"
			  % (AdapterRemoval_PATH, fastq1, fastq2, prefix, threads))
	return merge_reads, single_reads_1, single_reads_2


def analyze_merged_deletion_reads(input_sam_file):
	input_sam = HTSeq.SAM_Reader(input_sam_file)
	Clipping_reads = set()
	Mapping_reads = set()

	for sam_line in input_sam:
		if sam_line.aligned:
			clipping = 0
			for cigar_line in sam_line.cigar:
				if cigar_line.type == "N":
					Clipping_reads.add("@" + sam_line.get_sam_line().split("\t")[0])
					# insert_size = cigar_line.size
					clipping += 1
				elif cigar_line.type == "D" and cigar_line.size > 2:
					Clipping_reads.add("@" + sam_line.get_sam_line().split("\t")[0])
					clipping += 1
			if clipping == 0:
				Mapping_reads.add("@" + sam_line.get_sam_line().split("\t")[0])

	kept_reads = Clipping_reads.difference(Mapping_reads)
	return kept_reads


def analyze_single_deletion_reads(input_sam_file):
	input_sam = HTSeq.SAM_Reader(input_sam_file)
	input_sam = HTSeq.pair_SAM_alignments(input_sam)
	Clipping_reads_1 = set()
	Mapping_reads_1 = set()
	Clipping_reads_2 = set()
	Mapping_reads_2 = set()
	for sam_line in input_sam:
		if (sam_line[0] is not None and sam_line[0].aligned) and (sam_line[1] is not None and sam_line[1].aligned):
			clipping_1 = 0
			for cigar_line_1 in sam_line[0].cigar:
				if cigar_line_1.type == "N":
					Clipping_reads_1.add("@" + sam_line[0].get_sam_line().split("\t")[0])
					clipping_1 += 1
				elif cigar_line_1.type == "D" and cigar_line_1.size > 2:
					Clipping_reads_1.add("@" + sam_line[0].get_sam_line().split("\t")[0])
					clipping_1 += 1
			if clipping_1 == 0:
				Mapping_reads_1.add("@" + sam_line[0].get_sam_line().split("\t")[0])

			clipping_2 = 0
			for cigar_line_2 in sam_line[1].cigar:
				if cigar_line_2.type == "N":
					Clipping_reads_2.add("@" + sam_line[1].get_sam_line().split("\t")[0])
					clipping_2 += 1
				elif cigar_line_2.type == "D" and cigar_line_2.size > 2:
					Clipping_reads_2.add("@" + sam_line[1].get_sam_line().split("\t")[0])
					clipping_2 += 1
			if clipping_2 == 0:
				Mapping_reads_2.add("@" + sam_line[1].get_sam_line().split("\t")[0])

	Clipping_reads = Clipping_reads_1.union(Clipping_reads_2)
	Mapping_reads = Mapping_reads_1.intersection(Mapping_reads_2)
	kept_reads = Clipping_reads.difference(Mapping_reads)
	return kept_reads


def cigar_analyse(sam_line):
	if sam_line.aligned:
		clipping = 0
		read_start = -1
		mapped_size = 0
		for cigar_id in range(0, len(sam_line.cigar)):
			if sam_line.cigar[cigar_id].type == "M":
				mapped_size += sam_line.cigar[cigar_id].size
				if clipping == 0:
					if read_start == -1:
						read_start = sam_line.cigar[cigar_id].ref_iv.start
						read_start_clip = sam_line.cigar[cigar_id].ref_iv.end
					else:
						read_start_clip = sam_line.cigar[cigar_id].ref_iv.end
				elif clipping >= 1:
					if sam_line.cigar[cigar_id - 1].type == "N" or (
							sam_line.cigar[cigar_id - 1].type == "D" and sam_line.cigar[cigar_id - 1].size > 3):
						read_end_clip = sam_line.cigar[cigar_id].ref_iv.start
						read_end = sam_line.cigar[cigar_id].ref_iv.end
					elif (sam_line.cigar[cigar_id - 1].type == "I" or sam_line.cigar[cigar_id - 1].type == "D") and (sam_line.cigar[cigar_id - 2].type == "N" or (
							sam_line.cigar[cigar_id - 2].type == "D" and sam_line.cigar[cigar_id - 2].size > 3)):
						read_end_clip = sam_line.cigar[cigar_id].ref_iv.start
						read_end = sam_line.cigar[cigar_id].ref_iv.end
					else:
						read_end = sam_line.cigar[cigar_id].ref_iv.end

			elif sam_line.cigar[cigar_id].type == "N" or (
					sam_line.cigar[cigar_id].type == "D" and sam_line.cigar[cigar_id].size > 3):
				if clipping == 0:
					insert_size = sam_line.cigar[cigar_id].size
				else:
					insert_size += sam_line.cigar[cigar_id].size
				clipping += 1

		if clipping > 0:
			return (clipping, read_start, read_start_clip, read_end_clip, read_end, insert_size, mapped_size)
		else:
			return (clipping, read_start, 0, 0, read_start_clip, 0, mapped_size)
	else:
		return (-1, 0, 0, 0, 0, 0, 0)


def output_merged_deletion_reads(input_sam_file):
	output_file = input_sam_file.split(":")[0] + "merged_pos.txt"
	input_sam = HTSeq.SAM_Reader(input_sam_file)
	with open(output_file, "w") as output_list:
		# output_list.write("read_ID\tread_start\tgap_start\tgap_end\tread_end\tdel_size\tother_info\n")
		for sam_line in input_sam:

			(clipping, read_start, read_start_clip, read_end_clip, read_end, insert_size, mapped_size) = cigar_analyse(sam_line)

			if clipping > 0:
				if (read_start_clip - read_start) >= 30 and (read_end - read_end_clip) >= 30:
					output_list.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(
						sam_line.get_sam_line().split("\t")[0], str(read_start), str(read_start_clip), str(read_end_clip), str(read_end), str(insert_size), "merged"))
				elif clipping > 1 and (mapped_size + read_start + read_end_clip - read_end - read_start_clip) >= 30:
					output_list.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(
						sam_line.get_sam_line().split("\t")[0], str(read_start), str(read_start_clip), str(read_end_clip), str(read_end), str(insert_size), "merged,multi_clip"))

	return output_file


def output_single_deletion_reads(input_sam_file, library_max_size):
	output_file = input_sam_file.split(":")[0] + "single_pos.txt"
	input_sam = HTSeq.SAM_Reader(input_sam_file)
	input_sam = HTSeq.pair_SAM_alignments(input_sam)
	with open(output_file, "w") as output_list:
		output_list.write("read_ID\tread_start\tgap_start\tgap_end\tread_end\tdel_size\tother_info\n")
		for sam_line in input_sam:
			if (sam_line[0] is not None and sam_line[0].aligned) and (sam_line[1] is not None and sam_line[1].aligned):
				(clipping_1, read_start_1, read_start_clip_1, read_end_clip_1, read_end_1, insert_size_1, mapped_size_1) = cigar_analyse(sam_line[0])
				(clipping_2, read_start_2, read_start_clip_2, read_end_clip_2, read_end_2, insert_size_2, mapped_size_2) = cigar_analyse(sam_line[1])
				read_start = min(read_start_1, read_start_2)
				read_end = max(read_end_1, read_end_2)

				if read_end - read_start - insert_size_1 - insert_size_2 < library_max_size:
					if clipping_1 + clipping_2 > 0:  # at least one clip read
						if clipping_1 * clipping_2 == 0:  # only one clip read
							if clipping_1 > 0:  # read1 clipped
								insert_size = insert_size_1
								read_start_clip = read_start_clip_1
								read_end_clip = read_end_clip_1
								mapped_size = mapped_size_1
							else:  # read2 clipped
								insert_size = insert_size_2
								read_start_clip = read_start_clip_2
								read_end_clip = read_end_clip_2
								mapped_size = mapped_size_2

							if (read_start_clip - read_start) >= 30 and (read_end - read_end_clip) >= 30:
								output_list.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
									sam_line[0].get_sam_line().split("\t")[0], str(read_start), str(read_start_clip),
									str(read_end_clip), str(read_end), str(insert_size), "paired"))
							elif clipping_1 + clipping_2 > 1 and (mapped_size + read_start + read_end_clip - read_end - read_start_clip) >= 30:
								output_list.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
									sam_line[0].get_sam_line().split("\t")[0], str(read_start), str(read_start_clip),
									str(read_end_clip), str(read_end), str(insert_size), "paired,multi"))

						else:  # both are clip reads
							insert_size = insert_size_1 + insert_size_2
							read_start_clip = min(read_start_clip_1, read_start_clip_2)
							read_end_clip = max(read_end_clip_1, read_end_clip_2)
							mapped_size = mapped_size_1 + mapped_size_2

							if (read_start_clip - read_start) >= 30 and (read_end - read_end_clip) >= 30:
								output_list.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
									sam_line[0].get_sam_line().split("\t")[0], str(read_start), str(read_start_clip),
									str(read_end_clip), str(read_end), str(insert_size), "paired_clip"))
							elif (clipping_1 > 1 or clipping_2 > 1) and (mapped_size + read_start + read_end_clip - read_end - read_start_clip) >= 30:
								output_list.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
									sam_line[0].get_sam_line().split("\t")[0], str(read_start), str(read_start_clip),
									str(read_end_clip), str(read_end), str(insert_size), "paired_clip,multi"))

	return output_file


def read_pos_list(all_pos_list):
	pos_dict = {}

	with open(all_pos_list, "r") as input_pos_list:
		try:
			next(input_pos_list)
			repeat_pos_set = set()
			pos_dict_tmp = {}
			for pos_line in input_pos_list:
				pos_info_list = pos_line.rstrip().split("\t")
				bam_id = pos_info_list[0]
				if bam_id != "read_ID":
					if pos_dict.get(bam_id):
						repeat_pos_set.add(bam_id)
					else:
						pos_dict_tmp[bam_id] = pos_info_list  # [bam_id, read_start, gap_start, gap_end, read_end, del_size, other_info]

			for bam_id, pos_info_list in pos_dict_tmp.iteritems():
				if bam_id not in repeat_pos_set:
					pos_dict[bam_id] = pos_info_list

		except:
			sys.stderr.write("No available result for further analyze, please check post list file.\n")
			sys.exit(1)

	return pos_dict


def remove_repeat_masker_bedtools(pos_dict, bedtools_tmp_name, bedtools_intersectBed_path, repeatmasker_annotation_file_path):
	repeat_pos = list()
	bedtools_tmp_name2 = bedtools_tmp_name + "2"

	with open(bedtools_tmp_name, "w") as bedtools_tmp:
		for bam_id, pos_info_list in pos_dict.iteritems():
			read_start = pos_info_list[1]
			gap_start = pos_info_list[2]
			bedtools_tmp.write("%s\t%s\t%s\t%s\n" % ("EBV", read_start, gap_start, bam_id))

			gap_end = pos_info_list[3]
			read_end = pos_info_list[4]
			bedtools_tmp.write("%s\t%s\t%s\t%s\n" % ("EBV", gap_end, read_end, bam_id))

	os.system("%s -a %s -b %s -wa -f 0.50 > %s" % (bedtools_intersectBed_path, bedtools_tmp_name, repeatmasker_annotation_file_path, bedtools_tmp_name2))

	with open(bedtools_tmp_name2, "r") as bedtools_tmp2:
		for bedtools_line in bedtools_tmp2:
			bam_id_info = bedtools_line.rstrip().split("\t")[3]
			if pos_dict.get(bam_id_info):
				pos_dict.pop(bam_id_info)

	os.system("rm -f %s %s" %(bedtools_tmp_name2, bedtools_tmp_name))

	return pos_dict


def gc_calculate(sequence):
	gc_count = 0
	sequence = sequence.rstrip()
	read_length = len(sequence)
	for i in sequence:
		if i in ['G', 'C', 'g', 'c']:
			gc_count += 1
		elif i == "N":
			read_length -= 1
	gc_content = float(gc_count) / read_length
	return gc_content


def gc_adjust(pos_dict, bedtools_tmp_name, bedtools_path, EBV_NCBI):
	bedtools_tmp_name2 = bedtools_tmp_name + ".fa"
	gc_content_rate = 0.8
	gc_content_rate_2 = 1 - gc_content_rate

	with open(bedtools_tmp_name, "w") as bedtools_tmp:
		for bam_id, pos_info_list in pos_dict.iteritems():
			read_start = pos_info_list[1]
			gap_start = pos_info_list[2]
			bedtools_tmp.write("%s\t%s\t%s\t%s\n" % ("EBV", read_start, gap_start, bam_id))

			gap_end = pos_info_list[3]
			read_end = pos_info_list[4]
			bedtools_tmp.write("%s\t%s\t%s\t%s\n" % ("EBV", gap_end, read_end, bam_id))

	os.system("%s getfasta -fi %s.fa -bed %s -fo %s -name -fullHeader" %
			  (bedtools_path, EBV_NCBI, bedtools_tmp_name, bedtools_tmp_name2))

	with open(bedtools_tmp_name2, "r") as read_seq_fa:
		for read_seq in read_seq_fa:
			if read_seq.startswith(">"):
				bam_id = read_seq.rstrip().replace(">", "")
			else:
				gc_content_seq = gc_calculate(read_seq)
				if gc_content_seq > gc_content_rate or gc_content_seq < gc_content_rate_2:
					if pos_dict.get(bam_id):
						pos_dict.pop(bam_id)

	os.system("rm -f %s %s" % (bedtools_tmp_name2, bedtools_tmp_name))

	return pos_dict


def group_pos(pos_dict, output_file_name):
	pos_EBV = {}
	for bam_id, pos_info_list in pos_dict.iteritems():
		gap_start = pos_info_list[2]
		gap_end = pos_info_list[3]

		pos_num_1 = int(gap_start) / 1000  # 1,000nt for a single window (EBV)
		pos_num_2 = int(gap_end) / 1000  # 1,000nt for a single window (EBV)

		del_list = range(pos_num_1, pos_num_2 + 1)

		for del_pos in del_list:
			del_pos = '0' * (3 - len(str(del_pos))) + str(del_pos)
			if not pos_EBV.get(del_pos):
				pos_EBV[del_pos] = []
			gap_id = bam_id + "(" + gap_start + "-" + gap_end + ")"
			pos_EBV[del_pos].append(gap_id)

	with open(output_file_name, "w") as output_file:
		output_file.write("EBV_position\tdeletion_counts\tdeletion_reads\n")
		for del_pos in sorted(pos_EBV.keys()):
			output_file.write(del_pos + "\t" + str(len(pos_EBV[del_pos])) + "\t" + ",".join(pos_EBV[del_pos]) + "\n")


def merge_left_and_right(left_fastq, right_fastq):
	out_fastq = left_fastq.split("_left.fq")[0] + "_left_right.fq"
	with open(out_fastq, "w") as out_seq:
		with open(left_fastq, "r") as left_seq:
			seq_header = left_seq.readline().split(":")[0]
		with open(left_fastq, "r") as left_seq:
			with open(right_fastq, "r") as right_seq:
				for left_line, right_line in zip(left_seq, right_seq):
					if left_line.startswith(seq_header):
						out_seq.write(left_line)
					elif left_line == "+\n":
						out_seq.write(left_line)
					else:
						out_seq.write(left_line.rstrip() + right_line)


def fastq_cut_left_and_right(fastq_in):
	left_fastq = fastq_in.split(".fq")[0] + "_left.fq"
	right_fastq = fastq_in.split(".fq")[0] + "_right.fq"
	os.system("seqtk trimfq -B 50 %s > %s" % (fastq_in, left_fastq))
	os.system("seqtk trimfq -E 50 %s > %s" % (fastq_in, right_fastq))
	merge_left_and_right(left_fastq, right_fastq)


def filter_gap_length(pos_dict, gap_length):
	pos_dict_out = {}
	for bam_id, pos_info_list in pos_dict.iteritems():
		if int(pos_info_list[5]) > gap_length:
			pos_dict_out[bam_id] = pos_info_list
	return pos_dict_out






