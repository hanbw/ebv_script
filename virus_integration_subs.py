#! /usr/bin/env python
import os
import sys
import threading
import Queue
import HTSeq

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


def mapping_paired_ends(BWA_PATH, prefix, genome_PATH, fastq1, fastq2, outsam, thread):
	os.system("%s mem -t %s -M -R \"@RG\\tID:%s\\tLB:%s\\tPL:ILLUMINA\\tSM:%s\" %s %s %s > %s"
			  % (BWA_PATH, thread, prefix, prefix, prefix, genome_PATH, fastq1, fastq2, outsam))


def analyse_reads(sam_file):
	with open(sam_file, "r") as sam_lines:
		mapped_reads = set()
		for sam_line in sam_lines:
			if sam_line.startswith("@"):
				pass
			elif ("H" in sam_line.split("\t")[5]) or ("S" in sam_line.split("\t")[5]) or ("*" in sam_line.split("\t")[5]):
				pass
			else:
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
	read_left = read_all.difference(intersection_reads)
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


def analyze_integration_reads(input_sam_file_1, input_sam_file_2):
	input_sam_1 = HTSeq.SAM_Reader(input_sam_file_1)
	input_sam_2 = HTSeq.SAM_Reader(input_sam_file_2)

	sam1_reads = set()
	sam2_reads = set()

	for sam_line_1 in input_sam_1:
		if sam_line_1.aligned:
			sam1_reads.add("@" + sam_line_1.get_sam_line().split("\t")[0])

	for sam_line_2 in input_sam_2:
		if sam_line_2.aligned:
			sam2_reads.add("@" + sam_line_2.get_sam_line().split("\t")[0])

	kept_reads = sam1_reads.intersection(sam2_reads)
	return kept_reads


def analyse_read_interval(query_read):
	read_length = len(query_read.get_sam_line().split("\t")[9])
	for i in query_read.cigar:
		part = 0
		if i.type == "M":
			if part == 0:
				chrom = i.ref_iv.chrom
				ref_start = i.ref_iv.start
				ref_end = i.ref_iv.end
				part += 1
			elif part > 0:
				ref_end = i.ref_iv.end
				part += 1
			else:
				print "Err in analyse_read_interval"
			if i.ref_iv.strand == "+":
				if part == 1:
					strand = "+"
					query_start = i.query_from
					query_end = i.query_to
				elif part > 1:
					query_end = i.query_to

			elif i.ref_iv.strand == "-":
				if part == 1:
					strand = "-"
					query_end = read_length - i.query_from
					query_start = read_length - i.query_to
				elif part > 1:
					query_start = read_length - i.query_to
			else:
				print "Err in analyse_read_interval (strand)"
	return [query_start, query_end, chrom, ref_start, ref_end, strand, read_length, ""]


def analyse_sam_interval(input_sam):
	pos_dict = {}
	pos_dict_second = {}
	for line in input_sam:
		read_id = line.get_sam_line().split("\t")[0]
		second = 0
		if line.aligned and not line.not_primary_alignment:
			pos_dict[read_id] = analyse_read_interval(line)
		elif line.aligned and line.not_primary_alignment == "True":
			second += 1
			second_align = analyse_read_interval(line)
			if second_align[0] >= pos_dict[read_id][0] and second_align[1] <= pos_dict[read_id][1]:  # second alignment is covered by first alignment
				pos_dict[read_id][7] = "multi_cov_%s_%s_%s_%s_%s_%s" % (
					pos_dict[read_id][0], pos_dict[read_id][1], pos_dict[read_id][2],
					pos_dict[read_id][3], pos_dict[read_id][4], pos_dict[read_id][5])

			else:
				if second == 1:
					pos_dict_second[read_id] = []
				pos_dict_second[read_id].append(second_align)
				pos_dict[read_id][7] = "multi_uncov"

				# gap model
				if pos_dict[read_id][2] == second_align[2] and pos_dict[read_id][5] == second_align[5]:  # same chr and same strand
					read_length = len(line.get_sam_line().split("\t")[9])
					if min(pos_dict[read_id][0], second_align[0]) <= 1 and max(pos_dict[read_id][1], second_align[1]) >= (read_length - 1):  # cover start and end of read (1bp gap)
						if pos_dict[read_id][1] - pos_dict[read_id][0] + second_align[1] - second_align[0] < read_length:  # with a gap inside read
							pos_gap = ["gap", min(pos_dict[read_id][0], second_align[0]),
													 min(pos_dict[read_id][1], second_align[1]),
													 max(pos_dict[read_id][0], second_align[0]),
													 max(pos_dict[read_id][1], second_align[1]),
													 second_align[2], min(pos_dict[read_id][3], second_align[3]),
													 min(pos_dict[read_id][4], second_align[4]),
													 max(pos_dict[read_id][3], second_align[3]),
													 max(pos_dict[read_id][4], second_align[4]), second_align[5]]
							pos_dict_second[read_id].append(pos_gap)

	return pos_dict, pos_dict_second


def identify_clipping(ebv_pos, human_pos):
	judgement = 0
	gap_length = human_pos[6]
	if ebv_pos[1] - ebv_pos[0] < ebv_pos[6] - 30 and  human_pos[1] - human_pos[0] < human_pos[6] - 30:  # each species must remain at least 30 bp unmapped
		if min(ebv_pos[0], human_pos[0]) <= 1 and max(ebv_pos[1], human_pos[1]) >= human_pos[6] - 1:  # cover start and end of read (1bp gap)
			gap_length = max(ebv_pos[0], human_pos[0]) - min(ebv_pos[1], human_pos[1])
			if gap_length < 10:  # gap of two map result < 10 bp
				judgement = 1
	return judgement, gap_length


def identify_fully_matched(genome_pos):
	judgement = 0
	if genome_pos[0] < 2 and genome_pos[1] > genome_pos[6] - 2:
		judgement = 1
	return judgement


def get_single_read_output(ebv_pos, human_pos):
	ebv_strand = ebv_pos[5]
	human_strand = human_pos[5]

	if ebv_pos[0] < 2:
		if ebv_strand == "+":
			ebv_start = ebv_pos[3]
			ebv_gap = ebv_pos[4]
		elif ebv_strand == "-":
			ebv_start = ebv_pos[4]
			ebv_gap = ebv_pos[3]
		if human_strand == "+":
			human_gap = human_pos[3]
			human_end = human_pos[4]
		elif human_strand == "-":
			human_gap = human_pos[4]
			human_end = human_pos[3]
	else:
		if ebv_strand == "+":
			ebv_start = ebv_pos[4]
			ebv_gap = ebv_pos[3]
		elif ebv_strand == "-":
			ebv_start = ebv_pos[3]
			ebv_gap = ebv_pos[4]
		if human_strand == "+":
			human_end = human_pos[3]
			human_gap = human_pos[4]
		elif human_strand == "-":
			human_end = human_pos[4]
			human_gap = human_pos[3]

	human_chr = human_pos[2]

	return [ebv_start, ebv_gap, ebv_strand, human_gap, human_end, human_strand, human_chr]


def output_merged_integration_reads(outsam_clipping_merge_reads_hg19, outsam_clipping_merge_reads_ebv):
	output_file = outsam_clipping_merge_reads_hg19.split(":")[0] + "_integration_merged_pos.txt"
	with open(output_file, "w") as output_list:
		output_list.write("read_ID\tebv_start\tebv_gap\tebv_strand\thg19_gap\thg19_end\thg19_strand\thg_19_chr\toverlap\tother_info\n")

		input_sam_hg19 = HTSeq.SAM_Reader(outsam_clipping_merge_reads_hg19)
		input_sam_ebv = HTSeq.SAM_Reader(outsam_clipping_merge_reads_ebv)

		hg19_pos_dict, hg19_pos_dict_second = analyse_sam_interval(input_sam_hg19)
		ebv_pos_dict, ebv_pos_dict_second = analyse_sam_interval(input_sam_ebv)

		for sam_id,sam_line in hg19_pos_dict.iteritems():
			if ebv_pos_dict.get(sam_id):
				judgement, gap_length = identify_clipping(ebv_pos_dict[sam_id], sam_line)
				if judgement > 0:
					read_output = get_single_read_output(ebv_pos_dict[sam_id], sam_line)  # in (query_start, query_end, chrom, ref_start, ref_end, strand, read_length, "")
					read_output = map(str, read_output)
					output_list.write(sam_id + "\t" + "\t".join(read_output) + "\t" + str(gap_length) + "\tmerge\n")

				else:
					# min_gap = sam_line[6]
					min_gap = 10
					min_read_hg19 = ""
					min_read_ebv = ""
					gap_model = 0
					if hg19_pos_dict_second.get(sam_id):
						for sam_line_2nd in hg19_pos_dict_second[sam_id]:
							if sam_line_2nd[0] != "gap":
								judgement, gap_length = identify_clipping(ebv_pos_dict[sam_id], sam_line_2nd)
								if judgement > 0:
									if gap_length < min_gap:
										min_read_hg19 = sam_line_2nd
										min_read_ebv = ebv_pos_dict[sam_id]
							else:
								gap_model += 1

					if ebv_pos_dict_second.get(sam_id):
						for sam_line_2nd in ebv_pos_dict_second[sam_id]:
							if sam_line_2nd[0] != "gap":
								judgement, gap_length = identify_clipping(sam_line_2nd, hg19_pos_dict[sam_id])
								if judgement > 0:
									if gap_length < min_gap:
										min_read_hg19 = hg19_pos_dict[sam_id]
										min_read_ebv = sam_line_2nd
								else:
									gap_model += 1

					if hg19_pos_dict_second.get(sam_id) and ebv_pos_dict_second.get(sam_id):
						for sam_line_2nd_hg19 in hg19_pos_dict_second[sam_id]:
							if sam_line_2nd_hg19[0] != "gap":
								for sam_line_2nd_ebv in ebv_pos_dict_second[sam_id]:
									if sam_line_2nd_ebv[0] != "gap":
										judgement, gap_length = identify_clipping(sam_line_2nd_ebv, sam_line_2nd_hg19)
										if judgement > 0:
											if gap_length < min_gap:
												min_read_hg19 = sam_line_2nd_hg19
												min_read_ebv = sam_line_2nd_ebv

					# if gap_length < 10:
					if min_read_hg19 != "":
						read_output = get_single_read_output(min_read_ebv, min_read_hg19)  # query_start, query_end, chrom, ref_start, ref_end, strand, read_length, ""
						read_output = map(str, read_output)
						output_list.write(sam_id + "\t" + "\t".join(read_output) + str(gap_length) + "\tmerge_secondary\n")

					elif gap_model > 0:
						print sam_id + "gap_model\n"

	return output_file


def output_paired_non_clipping_integration_reads(outsam_clipping_single_reads_r1_hg19, outsam_clipping_single_reads_r2_hg19, outsam_clipping_single_reads_r1_ebv, outsam_clipping_single_reads_r2_ebv):
	output_file = outsam_clipping_single_reads_r1_hg19.split(".")[0] + "_integration_single_pos.txt"

	with open(output_file, "w") as output_list:
		output_list.write("read_ID\tebv_start\tebv_gap\tebv_strand\thg19_gap\thg19_end\thg19_strand\thg_19_chr\toverlap\tother_info\n")

		input_sam_r1_hg19 = HTSeq.SAM_Reader(outsam_clipping_single_reads_r1_hg19)
		input_sam_r1_ebv = HTSeq.SAM_Reader(outsam_clipping_single_reads_r1_ebv)
		input_sam_r2_hg19 = HTSeq.SAM_Reader(outsam_clipping_single_reads_r2_hg19)
		input_sam_r2_ebv = HTSeq.SAM_Reader(outsam_clipping_single_reads_r2_ebv)

		hg19_pos_dict_r1, hg19_pos_dict_second_r1 = analyse_sam_interval(input_sam_r1_hg19)
		hg19_pos_dict_r2, hg19_pos_dict_second_r2 = analyse_sam_interval(input_sam_r2_hg19)
		ebv_pos_dict_r1, ebv_pos_dict_second_r1 = analyse_sam_interval(input_sam_r1_ebv)
		ebv_pos_dict_r2, ebv_pos_dict_second_r2 = analyse_sam_interval(input_sam_r2_ebv)

		# identify (r1-human + r2-ebv) or (r1-ebv + r2-human) pairs

		ebv_r1_matched = {}
		for sam_id, sam_line in ebv_pos_dict_r1.iteritems():
			if identify_fully_matched(sam_line) > 0:  # fully matched
				ebv_r1_matched[sam_id] = sam_line
				if hg19_pos_dict_r2.get(sam_id):  # paired-end mapped
					if identify_fully_matched(hg19_pos_dict_r2[sam_id]) > 0 and ((not ebv_pos_dict_r2.get(sam_id)) or (ebv_pos_dict_r2.get(sam_id) and identify_fully_matched(ebv_pos_dict_r2[sam_id]) == 0)):
						sam_line = map(str, sam_line)
						hg19_pos_dict_r2_list = map(str, hg19_pos_dict_r2[sam_id])
						if sam_line[5] == "+" and hg19_pos_dict_r2_list[5] == "-":
							output_list.write(sam_id + "\t" + sam_line[3] + "\t" + sam_line[4] + "\t" + sam_line[5] + "\t" +
								hg19_pos_dict_r2_list[3] + "\t" + hg19_pos_dict_r2_list[4] + "\t" + hg19_pos_dict_r2_list[5] + "\t" + hg19_pos_dict_r2_list[2] + "\tNA\tpaired\n")
						elif sam_line[5] == "-" and hg19_pos_dict_r2_list[5] == "+":
							output_list.write(sam_id + "\t" + sam_line[4] + "\t" + sam_line[3] + "\t" + sam_line[5] + "\t" +
								hg19_pos_dict_r2_list[4] + "\t" + hg19_pos_dict_r2_list[3] + "\t" + hg19_pos_dict_r2_list[5] + "\t" + hg19_pos_dict_r2_list[2] + "\tNA\tpaired\n")
						elif sam_line[5] == "+" and hg19_pos_dict_r2_list[5] == "+":
							output_list.write(sam_id + "\t" + sam_line[3] + "\t" + sam_line[4] + "\t" + sam_line[5] + "\t" +
								hg19_pos_dict_r2_list[4] + "\t" + hg19_pos_dict_r2_list[3] + "\t" + hg19_pos_dict_r2_list[5] + "\t" + hg19_pos_dict_r2_list[2] + "\tNA\tpaired\n")
						elif sam_line[5] == "-" and hg19_pos_dict_r2_list[5] == "-":
							output_list.write(sam_id + "\t" + sam_line[4] + "\t" + sam_line[3] + "\t" + sam_line[5] + "\t" +
								hg19_pos_dict_r2_list[3] + "\t" + hg19_pos_dict_r2_list[4] + "\t" + hg19_pos_dict_r2_list[5] + "\t" + hg19_pos_dict_r2_list[2] + "\tNA\tpaired\n")

		ebv_r2_matched = {}
		for sam_id, sam_line in ebv_pos_dict_r2.iteritems():
			if identify_fully_matched(sam_line) > 0:
				ebv_r2_matched[sam_id] = sam_line
				if hg19_pos_dict_r1.get(sam_id):  # paired-end mapped
					if identify_fully_matched(hg19_pos_dict_r1[sam_id]) > 0 and ((not ebv_pos_dict_r1.get(sam_id)) or (ebv_pos_dict_r1.get(sam_id) and identify_fully_matched(ebv_pos_dict_r1[sam_id]) == 0)):
						sam_line = map(str, sam_line)
						hg19_pos_dict_r1_list = map(str, hg19_pos_dict_r1[sam_id])
						if sam_line[5] == "+" and hg19_pos_dict_r1_list[5] == "-":
							output_list.write(sam_id + "\t" + sam_line[3] + "\t" + sam_line[4] + "\t" + sam_line[5] + "\t" +
								hg19_pos_dict_r1_list[3] + "\t" + hg19_pos_dict_r1_list[4] + "\t" + hg19_pos_dict_r1_list[5] + "\t" + hg19_pos_dict_r1_list[2] + "\tNA\tpaired\n")
						elif sam_line[5] == "-" and hg19_pos_dict_r1[sam_id][5] == "+":
							output_list.write(sam_id + "\t" + sam_line[4] + "\t" + sam_line[3] + "\t" + sam_line[5] + "\t" +
								hg19_pos_dict_r1_list[4] + "\t" + hg19_pos_dict_r1_list[3] + "\t" + hg19_pos_dict_r1_list[5] + "\t" + hg19_pos_dict_r1_list[2] + "\tNA\tpaired\n")
						elif sam_line[5] == "+" and hg19_pos_dict_r1[sam_id][5] == "+":
							output_list.write(sam_id + "\t" + sam_line[3] + "\t" + sam_line[4] + "\t" + sam_line[5] + "\t" +
								hg19_pos_dict_r1_list[4] + "\t" + hg19_pos_dict_r1_list[3] + "\t" + hg19_pos_dict_r1_list[5] + "\t" + hg19_pos_dict_r1_list[2] + "\tNA\tpaired\n")
						elif sam_line[5] == "-" and hg19_pos_dict_r1[sam_id][5] == "-":
							output_list.write(sam_id + "\t" + sam_line[4] + "\t" + sam_line[3] + "\t" + sam_line[5] + "\t" +
								hg19_pos_dict_r1_list[3] + "\t" + hg19_pos_dict_r1_list[4] + "\t" + hg19_pos_dict_r1_list[5] + "\t" + hg19_pos_dict_r1_list[2] + "\tNA\tpaired\n")

		human_r1_matched = {}
		for sam_id, sam_line in hg19_pos_dict_r1.iteritems():
			if identify_fully_matched(sam_line) > 0:
				human_r1_matched[sam_id] = sam_line

		human_r2_matched = {}
		for sam_id, sam_line in hg19_pos_dict_r2.iteritems():
			if identify_fully_matched(sam_line) > 0:
				human_r2_matched[sam_id] = sam_line

		return output_file, ebv_r1_matched, ebv_r2_matched, human_r1_matched, human_r2_matched


def analyse_paired_one_clipping_integration_reads(single_pos_list, ebv_matched, human_matched, library_max_size):
	output_file = single_pos_list.split(".")[0] + "_integration_single_pos.txt"
	with open(output_file, "a") as output_list:
		with open(single_pos_list, "r") as input_list:
			next(input_list)
			for read_line in input_list:  # sam_id0, ebv_start1, ebv_gap2, ebv_strand3, human_gap4, human_end5, human_strand6, human_chr7, gap_length8, other_info9
				read_line_list = read_line.split("\t")
				read_id = read_line_list[0]
				ebv_start = int(read_line_list[1])
				ebv_gap = int(read_line_list[2])
				human_gap = int(read_line_list[4])
				human_end = int(read_line_list[5])

				# if ebv_start < ebv_gap:  # read start with ebv, mate pair read is human
				if human_matched.get(read_id):  # human_matched[read_id] = [query_start, query_end, chrom, ref_start, ref_end, strand, read_length, ""]
					human_line = human_matched[read_id]
					if read_line.split("\t")[7] == human_line[2] and read_line_list[6] != human_line[5]:  # same chromosome and different strand

						if human_end > human_gap and human_line[4] - human_gap < library_max_size and human_line[3] - human_gap >= 0:

							# "read_ID\tebv_start\tebv_gap\tebv_strand\thg19_gap\thg19_end\thg19_strand\thg_19_chr\toverlap\tother_info\n"
							human_line = map(str, human_line)
							read_line_list = map(str, read_line_list)
							output_list.write(read_id + "\t" + str(ebv_start) + "\t" + str(ebv_gap) + "\t" + read_line_list[3] + "\t" + str(human_gap) + "\t" +
											  human_line[4] + "\t" + read_line_list[6] + "\t" + read_line_list[7] +  "\t" + read_line_list[8] + "\tclip_paired\n")
						elif human_end < human_gap and human_gap - human_line[3] < library_max_size and human_gap - human_line[4] >= 0:

							human_line = map(str, human_line)
							read_line_list = map(str, read_line_list)
							output_list.write(read_id + "\t" + str(ebv_start) + "\t" + str(ebv_gap) + "\t" + read_line_list[3] + "\t" + str(human_gap) + "\t" +
											  human_line[3] + "\t" + read_line_list[6] + "\t" + read_line_list[7] +  "\t" + read_line_list[8] + "\tclip_paired\n")

				# elif ebv_start > ebv_gap:  # read start with human, mate pair read is ebv
				if ebv_matched.get(read_id):
					ebv_line = ebv_matched[read_id]
					if read_line_list[3] != ebv_line[5]:  # same chromosome and different strand

						if ebv_start > ebv_gap and ebv_line[4] - ebv_gap < library_max_size and ebv_line[3] - ebv_gap >= 0:
							ebv_line = map(str, ebv_line)
							read_line_list = map(str, read_line_list)
							# "read_ID\tebv_start\tebv_gap\tebv_strand\thg19_gap\thg19_end\thg19_strand\thg_19_chr\toverlap\tother_info\n"
							output_list.write(read_id + "\t" + ebv_line[4] + "\t" + str(ebv_gap) + "\t" + read_line_list[3] + "\t" + str(human_gap) + "\t" +
											  read_line_list[5] + "\t" + read_line_list[6] + "\t" + read_line_list[7] + "\t" + read_line_list[8] + "\tclip_paired\n")
						elif ebv_start < ebv_gap and ebv_gap - ebv_line[3] < library_max_size and ebv_gap - ebv_line[4] >= 0:
							ebv_line = map(str, ebv_line)
							read_line_list = map(str, read_line_list)
							output_list.write(read_id + "\t" + ebv_line[3] + "\t" + str(ebv_gap) + "\t" + read_line_list[3] + "\t" + str(human_gap) + "\t" +
											  read_line_list[5] + "\t" + read_line_list[6] + "\t" + read_line_list[7] + "\t" + read_line_list[8] + "\tclip_paired\n")



def analyse_paired_two_clipping_integration_reads(single_pos_list_r1, single_pos_list_r2, output_file):
#	output_file = single_pos_list_r1.split(".")[0] + "_integration_single_pos.txt"
	with open(output_file, "a") as output_list:
		with open(single_pos_list_r1, "r") as pos_list_1:
			next(pos_list_1)
			for pos_line in pos_list_1:
				read_line_list_1 = pos_line.split("\t")
				read_id = read_line_list_1[0]
				with open(single_pos_list_r2, "r") as pos_list_2:
					next(pos_list_2)
					for pos_line_2 in pos_list_2:
						read_line_list_2 = pos_line_2.split("\t")
						if read_id == read_line_list_2[0]:
							read_line_list_1[-1] = "two_clip\n"
							read_line_list_2[-1] = "two_clip\n"
							output_list.write("\t".join(read_line_list_1))
							output_list.write("\t".join(read_line_list_2))


def read_pos_list(all_pos_list):
	pos_dict = {}

	with open(all_pos_list, "r") as input_pos_list:
		try:
			repeat_pos_set = set()
			pos_dict_tmp = {}
			for pos_line in input_pos_list:
				pos_info_list = pos_line.rstrip().split("\t")
				bam_id = pos_info_list[0]
				if bam_id != "read_ID":
					if pos_dict.get(bam_id):
						repeat_pos_set.add(bam_id)
					else:
						pos_dict_tmp[bam_id] = pos_info_list  # "read_ID0\tebv_start1\tebv_gap2\tebv_strand3\thg19_gap4\thg19_end5\thg19_strand6\thg_19_chr7\toverlap8\tother_info9\n"

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
			ebv_start = pos_info_list[1]
			ebv_gap = pos_info_list[2]
			ebv_min = str(min(int(ebv_start), int(ebv_gap)))
			ebv_max = str(max(int(ebv_start), int(ebv_gap)))
			bedtools_tmp.write("%s\t%s\t%s\t%s\n" % ("EBV", ebv_min, ebv_max, bam_id))

			human_chr = pos_info_list[7].replace("chr", "")
			human_gap = pos_info_list[4]
			human_end = pos_info_list[5]
			human_min = str(min(int(human_gap), int(human_end)))
			human_max = str(max(int(human_gap), int(human_end)))
			bedtools_tmp.write("%s\t%s\t%s\t%s\n" % (human_chr, human_min, human_max, bam_id))

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


def gc_adjust(pos_dict, bedtools_tmp_name, bedtools_path, hg19_EBV_PATH):
	bedtools_tmp_name2 = bedtools_tmp_name + ".fa"
	gc_content_rate = 0.8
	gc_content_rate_2 = 1 - gc_content_rate

	with open(bedtools_tmp_name, "w") as bedtools_tmp:
		for bam_id, pos_info_list in pos_dict.iteritems():
			ebv_start = pos_info_list[1]
			ebv_gap = pos_info_list[2]
			ebv_min = str(min(int(ebv_start), int(ebv_gap)))
			ebv_max = str(max(int(ebv_start), int(ebv_gap)))
			bedtools_tmp.write("%s\t%s\t%s\t%s\n" % ("chrEBV", ebv_min, ebv_max, bam_id))

			human_chr = pos_info_list[7]
			human_gap = pos_info_list[4]
			human_end = pos_info_list[5]
			human_min = str(min(int(human_gap), int(human_end)))
			human_max = str(max(int(human_gap), int(human_end)))
			bedtools_tmp.write("%s\t%s\t%s\t%s\n" % (human_chr, human_min, human_max, bam_id))

	os.system("%s getfasta -fi %s.fa -bed %s -fo %s -name -fullHeader" % (bedtools_path, hg19_EBV_PATH, bedtools_tmp_name, bedtools_tmp_name2))

	with open(bedtools_tmp_name2, "r") as read_seq_fa:
		for read_seq in read_seq_fa:
			if read_seq.startswith(">"):
				bam_id = read_seq.rstrip().replace(">", "")
			else:
				gc_content_seq = gc_calculate(read_seq)
				if gc_content_seq > gc_content_rate or gc_content_seq < gc_content_rate_2:
					if pos_dict.get(bam_id):
						pos_dict.pop(bam_id)

	os.system("rm -f %s %s" % (bedtools_tmp_name, bedtools_tmp_name2))

	return pos_dict


def group_pos(pos_dict, output_file_name_EBV, output_file_name_hg19):
	pos_EBV = {}
	pos_human = {}

	for bam_id, pos_info_list in pos_dict.iteritems():
		ebv_gap = pos_info_list[2]
		human_gap = pos_info_list[4]
		human_chr = pos_info_list[7]

		pos_num_ebv = int(ebv_gap) / 1000  # 1,000nt for a single window (EBV)
		pos_num_ebv = '0' * (3 - len(str(pos_num_ebv))) + str(pos_num_ebv)
		if not pos_EBV.get(pos_num_ebv):
			pos_EBV[pos_num_ebv] = []
		ebv_id = bam_id + "(" + ebv_gap + ")"
		pos_EBV[pos_num_ebv].append(ebv_id)

		pos_num_human = int(human_gap) / 100000  # 100,000nt for a single window (human)
		pos_num_human = '0' * (5 - len(str(pos_num_human))) + str(pos_num_human)
		pos_num_human = human_chr + ":" + pos_num_human
		if not pos_human.get(pos_num_human):
			pos_human[pos_num_human] = []
		human_id = bam_id + "(" + human_gap + ")"
		pos_human[pos_num_human].append(human_id)

	with open(output_file_name_EBV, "w") as output_file:
		output_file.write("EBV_position\tintegration_counts\tintegration_reads\n")
		for del_pos in sorted(pos_EBV.keys()):
			output_file.write(
				del_pos + "\t" + str(len(pos_EBV[del_pos])) + "\t" + ",".join(pos_EBV[del_pos]) + "\n")

	with open(output_file_name_hg19, "w") as output_file:
		output_file.write("hg19_position\tintegration_counts\tintegration_reads\n")
		for del_pos in sorted(pos_human.keys()):
			output_file.write(
				del_pos + "\t" + str(len(pos_human[del_pos])) + "\t" + ",".join(pos_human[del_pos]) + "\n")


