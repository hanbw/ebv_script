import os
import sys

# USAGE: python get_vcf_from_alignment_stretcher.py <input_fasta_file>
input_fasta_file = sys.argv[1]
# output_vcf = os.path.basename(input_alignment_file).split(".")[0] + ".vcf"

perl_script = "./fa2vcf_2seq_V3"
stretcher_bin = "~/bin/stretcher"
stretcher2clw_bin = "./stretcher2clw.py"
python_bin = "~/bin/python"

ref_EBV_fasta = "/data/home/hanbw/project/NKT_EBV/9-DATABASE/NC_007605.fa"

with open(input_fasta_file, "r") as input_fasta_data:
	fasta_tag_list = []
	for fasta_lines in input_fasta_data:
		if fasta_lines.startswith(">"):
			fasta_tag = fasta_lines.rstrip().split(" ")[0].replace(">", "")
			fasta_name = fasta_tag + ".ebv.fa"
			fasta_tag_list.append(fasta_name)
			os.system("touch %s" % (fasta_name))
			with open(fasta_name, "a") as output_fa:
				output_fa.write(fasta_lines)
		else:
			with open(fasta_name, "a") as output_fa:
				fasta_lines = fasta_lines.replace("N", "")
				#print fasta_lines
				output_fa.write(fasta_lines)

for fasta_file in fasta_tag_list:
	stretcher_file = fasta_file.replace("fa", "stretcher")
	#os.system("%s -in %s -out %s.clw -clw -clwstrict" % (muscle_bin, fasta_file, fasta_file))
	os.system("%s -asequence %s -bsequence %s -outfile %s -aformat pair" % (stretcher_bin, ref_EBV_fasta, fasta_file, stretcher_file))
	os.system("%s %s %s" % (python_bin, stretcher2clw_bin, stretcher_file))
	os.system("cat %s.clw | perl %s > %s.vcf" % (stretcher_file, perl_script, fasta_file.replace(".ebv.fa", "")))
