import os
import sys

input_stretcher_file = sys.argv[1]
clw_out_file = os.path.basename(input_stretcher_file) + ".clw"

clw_header = "CLUSTAL - ~/software/EMBOSS/bin/stretcher -asequence NC_007605.fa -bsequence *.fa -outfile *.out -aformat pair"

with open(clw_out_file, "w") as vcf_out:
	vcf_out.write(clw_header + "\n\n")

	with open(input_stretcher_file) as input_stretcher_pair:
		input_stretcher_line = input_stretcher_pair.readlines()
		input_stretcher_length = len(input_stretcher_line)
		info_end = 0
		REF_NAME = "EBV"
		for i in range(input_stretcher_length):
			if input_stretcher_line[i].startswith("# 1:"):
				REF_NAME = input_stretcher_line[i].rstrip().split()[2]
				ref_len = len(REF_NAME)
			elif input_stretcher_line[i].startswith("# 2:"):
				ALT_NAME = input_stretcher_line[i].rstrip().split()[2]
				alt_len = len(ALT_NAME)
			elif input_stretcher_line[i].startswith(REF_NAME):
				input_stretcher_info = input_stretcher_line[i].rstrip().split()
				REF_seq = input_stretcher_info[2]
				ALT_seq = input_stretcher_line[i+2].rstrip().split()[2]

				if info_end == 0:
					info_start = input_stretcher_line[i].index(REF_seq)
					info_end = info_start + len(REF_seq)
				Pair_info = input_stretcher_line[i+1][info_start:info_end]
				Pair_info = Pair_info.replace("|", "*").replace(".", " ")

				# Pair_info_list = list(Pair_info)
				# for j in range(len(REF_seq)):
				# 	if REF_seq[j] == "-" and ALT_seq[j] == "N":
				# 		Pair_info_list[j] = "*"
				# Pair_info = "".join(Pair_info_list)

				vcf_out.write(REF_NAME + " " * (alt_len - ref_len + 6) + REF_seq + "\n")
				vcf_out.write(ALT_NAME + " " * 6 + ALT_seq + "\n")
				vcf_out.write(" " * (alt_len + 6) + Pair_info + "\n\n")







