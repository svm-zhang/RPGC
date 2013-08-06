# Author: Simo Zhang
# Email: simozhan@indiana.edu
# Usage: python post_mapping_pipeline.py -h #

# Notes: 1. v1.00 ready to use #
#		 2. complete -s and -RG options
#		 3. add in fansier logger
#		 4. hope to add in break-and-continue function so that no need to re-run entire pipeline everytime something come up at later steps

import os, sys
import re
import argparse
from datetime import datetime
import shlex, subprocess
from subprocess import Popen, PIPE

# call Picard CleanSam #
def Clean_SAM_Files(raw_samfile, picard_prog) :
	cleaned_samfile = '.'.join(re.split('\.', raw_samfile)[0:3]) + ".cleaned.sam"
	picard_clean_sam_cmd = "java -Xmx2g -jar %s INPUT=%s OUTPUT=%s VALIDATION_STRINGENCY=STRICT" %(os.path.join(picard_prog, "CleanSam.jar"), raw_samfile, cleaned_samfile)
	sys.stdout.write("%s\n" %(picard_clean_sam_cmd))
	proc = Popen(shlex.split(picard_clean_sam_cmd), stderr=PIPE, stdout=PIPE)
	proc_stdout, proc_stderr = proc.communicate()
	if proc.returncode != 0 or re.search("(Error|error)", proc_stderr) or re.search("(Error|error)", proc_stdout) :
		sys.stderr.write(TimeStamper() + " [Program Checker]	Error: Picard CleanSam terminates unexpectedly\n")
		sys.stderr.write(TimeStamper() + " [Program Checker]	Error:%s\t%s\n" %(proc_stdout, proc_stderr))
		sys.exit(1)
	else :
		if proc_stderr != "" :
			sys.stderr.write("%s\n" %(proc_stderr))
	return cleaned_samfile

# call Picard SamFormatConverter to convert sam to bam #
def SAM_to_BAM(cleaned_samfile, picard_prog) :
	cleaned_bamfile = '.'.join(re.split('\.', cleaned_samfile)[0:3]) + ".cleaned.bam"
	picard_samtobam_cmd = "java -Xmx2g -jar %s INPUT=%s OUTPUT=%s VALIDATION_STRINGENCY=LENIENT" %(os.path.join(picard_prog, "SamFormatConverter.jar"), cleaned_samfile, cleaned_bamfile)
	sys.stdout.write("%s\n" %(picard_samtobam_cmd))
	proc = Popen(shlex.split(picard_samtobam_cmd), stderr=PIPE, stdout=PIPE)
	proc_stdout, proc_stderr = proc.communicate()
	if proc.returncode != 0 or re.search("(Error|error)", proc_stderr) or re.search("(Error|error)", proc_stdout) :
		sys.stderr.write(TimeStamper() + " [Program Checker]	Error: Picard SamFormatConverter terminates unexpectedly\n")
		sys.exit(1)
	proc = Popen(["rm", "%s" %(cleaned_samfile)], stderr=PIPE, stdout=PIPE)
	proc_stdout, proc_stderr = proc.communicate()
	if proc.returncode != 0 or re.search("(Error|error)", proc_stderr) or re.search("(Error|error)", proc_stdout) :
		sys.stderr.write(TimeStamper() + " [Program Checker]	Error: Failure to remove %s\n" %(cleaned_samfile))
		sys.exit(1)
	else :
		if proc_stderr != "" :
			sys.stderr.write("%s\n" %(proc_stderr))
	return cleaned_bamfile

# call Picard MergeSamFiles to merge bam files #
def Merge_Bams(list_bamfiles, outfile, picard_prog) :
	input_samfiles = ""
	for each_file in list_bamfiles :
		if input_samfiles == "" :
			input_samfiles += "INPUT=%s" %(each_file)
		else :
			input_samfiles += " INPUT=%s" %(each_file)
	picard_mergesam_cmd = "java -Xmx2g -jar %s %s OUTPUT=%s SORT_ORDER=coordinate USE_THREADING=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true" %(os.path.join(picard_prog, "MergeSamFiles.jar"), input_samfiles, outfile)
	sys.stdout.write("%s\n" %(picard_mergesam_cmd))
	proc = Popen(shlex.split(picard_mergesam_cmd), stderr=PIPE, stdout=PIPE)
	proc_stdout, proc_stderr = proc.communicate()
	if proc.returncode != 0 or re.search("(Error|error)", proc_stderr) or re.search("(Error|error)", proc_stdout) :
		sys.stderr.write(TimeStamper() + " [Program Checker]	Error: Picard MergeSamFiles terminates unexpectedly\n")
		sys.exit(1)
	else :
		if proc_stderr != "" :
			sys.stderr.write("%s\n" %(proc_stderr))

# call Picard AddOrReplaceReadGroups to put in @RG tag #
def AddRGs(bamfile, rg_def, outfile, picard_prog) :
	picard_addRG_cmd = "java -Xmx2g -jar %s INPUT=%s OUTPUT=%s %s SORT_ORDER=coordinate CREATE_INDEX=true" %(os.path.join(picard_prog, "AddOrReplaceReadGroups.jar"), bamfile, outfile, re.sub('\:', '=', ' '.join(re.split(r'\\t', rg_def)[1:])))
	sys.stdout.write("%s\n" %(picard_addRG_cmd))
	proc = Popen(shlex.split(picard_addRG_cmd), stderr=PIPE, stdout=PIPE)
	proc_stdout, proc_stderr = proc.communicate()
	if proc.returncode != 0 or re.search("(Error|error)", proc_stderr) or re.search("(Error|error)", proc_stdout) :
		sys.stderr.write(TimeStamper() + " [Program Checker]	Error: Picard AddOrReplaceReadGroups terminates unexpectedly\n")
		sys.exit(1)
	else :
		if proc_stderr != "" :
			sys.stderr.write("%s\n" %(proc_stderr))

# get rid of the duplicates #
def Mark_Duplicates(cleaned_bamfile, picard_prog) :
	sorted_bamfile = '.'.join(re.split('\.', cleaned_bamfile)[0:3]) + ".cleaned.sort.bam"
	picard_sort_cmd = "java -Xmx2g -jar %s I=%s O=%s SO=coordinate VALIDATION_STRINGENCY=LENIENT" %(os.path.join(picard_prog, "SortSam.jar"), cleaned_bamfile, sorted_bamfile)
	print picard_sort_cmd
	sys.stdout.write("%s [SortSam.jar]: %s\n" %(TimeStamper(), picard_sort_cmd))
	proc = Popen(shlex.split(picard_sort_cmd), stderr=PIPE, stdout=PIPE)
	proc_stdout, proc_stderr = proc.communicate()
	if proc.returncode != 0 or re.search("(Error|error)", proc_stderr) or re.search("(Error|error)", proc_stdout) :
		sys.stderr.write(TimeStamper() + " [Program Checker]	Error: Picard SortSam terminates unexpectedly\n")
		sys.exit(1)
	else :
		if proc_stderr != "" :
			sys.stderr.write("%s\n" %(proc_stderr))
	dup_log_file = '.'.join(re.split('\.', cleaned_bamfile)[0:3]) + ".dup.report"
	rmdup_bamfile = '.'.join(re.split('\.', sorted_bamfile)[0:3]) + ".cleaned.sorted.rmdup.bam"
	picard_markdup_cmd = "java -Xmx2g -jar %s REMOVE_DUPLICATES=true I=%s O=%s VALIDATION_STRINGENCY=LENIENT M=%s" %(os.path.join(picard_prog, "MarkDuplicates.jar"), sorted_bamfile, rmdup_bamfile, dup_log_file)
	sys.stdout.write("%s [MarkDuplicates.jar]: %s\n" %(TimeStamper(), picard_markdup_cmd))
	proc = Popen(shlex.split(picard_markdup_cmd), stderr=PIPE, stdout=PIPE)
	proc_stdout, proc_stderr = proc.communicate()
	if proc.returncode != 0 or re.search("(Error|error)", proc_stderr) or re.search("(Error|error)", proc_stdout) :
		sys.stderr.write(TimeStamper() + " [Program Checker]	Error: Picard MarkDuplicates terminates unexpectedly\n")
		sys.exit(1)
	else :
		if proc_stderr != "" :
			sys.stderr.write("%s\n" %(proc_stderr))
	return rmdup_bamfile

# check whether or not file and directory exist #
def Path_Checker(path, type) :
	if type == 'f' :
		tmp_path = re.split(' ', path)
		for i in range(len(tmp_path)) :
			if not os.path.exists(tmp_path[i]) :
				sys.stderr.write(TimeStamper() + " [Path Checker]	Error: cannnot reach the file: %s\n" %(tmp_path[i]))
				sys.exit(1)
			else :
				sys.stdout.write(TimeStamper() + " [Path Checker]	%s ... [checked]\n" %(tmp_path[i]))
	elif type == 'd' :
		if not os.path.exists(path) :
			os.makedirs(path)
			sys.stdout.write(TimeStamper() + " [Path Checker]	Creating directory: %s\n" %(path))

# check RG definition #
def RG_Def_Checker(rg_def) :
	tmp_def = re.split(r'\\t', rg_def)
	num_component = 0
	for i in range(1, len(tmp_def)) :
		if len(re.split(':', tmp_def[i])) == 2 :
			tmp_component = re.split(':', tmp_def[i])
			if tmp_component[0] == "" or not tmp_component[0] in ["ID", "SM", "LB", "PL", "PU"] :
				sys.stderr.write(TimeStamper() + " [@RG Definition Checker]	Error: missing @RG key %s\n") %(tmp_def[i])
				sys.stderr.write(TimeStamper() + " [@RG Definition Checker]	Example: @RG\tID:cs1\tSM:cs1\tLB:FRAG\tPL:Illumina\tPU:RUN\n")
				sys.exit(1)
			else :
				if tmp_component[1] == "" :
					sys.stderr.write(TimeStamper() + " [@RG Definition Checker]	Error: missing value for @RG_%s %s\n") %(tmp_component[0], tmp_def[i])
					sys.exit(1)
				else :
					num_component += 1
		else :
			sys.stderr.write(TimeStamper() + " [@RG Definition Checker]	Error: incorrect @RG format %s\n") %(tmp_def[i])
			sys.stderr.write(TimeStamper() + " [@RG Definition Checker]	Example: @RG\tID:cs1\tSM:cs1\tLB:FRAG\tPL:Illumina\tPU:RUN\n")
			sys.exit(1)
	if num_component != 5 :
		sys.stderr.write(TimeStamper() + " [@RG Definition Checker]	Error: missing @RG keys in the definition %s\n" %(rg_def))
		sys.exit(1)

# run Picard pipeline #
def Run_PICARD_Pipe(config_dict, merge_off, outprefix, picard_prog) :
	bamfiles_per_RGID = []
	bamfiles_all = []
	for rg_def, samfiles in config_dict.iteritems() :
		rg_id = ""
		tmp_def = re.split(r'\\t', rg_def)
		for i in range(1, len(tmp_def)) :
			tmp_component = re.split(':', tmp_def[i])
			if tmp_component[0] == "ID" :
				rgid = tmp_component[1]
		sys.stdout.write("\n*******************************\n%s [RGID]: %s\n" %(TimeStamper(), rgid))
		for each_raw_samfile in re.split('\s+', samfiles) :
			sys.stdout.write("%s [INPUT SAM]: %s\n" %(TimeStamper(), each_raw_samfile))
			sys.stdout.write("%s [CleanSAM.jar]: " %(TimeStamper()))
			each_cleaned_samfile = Clean_SAM_Files(each_raw_samfile, picard_prog)
			sys.stdout.write("%s [INPUT SAM]: %s\n" %(TimeStamper(), each_cleaned_samfile))
			sys.stdout.write("%s [SamFormatConverter.jar]: " %(TimeStamper()))
			each_cleaned_bamfile = SAM_to_BAM(each_cleaned_samfile, picard_prog)
			sys.stdout.write("%s [INPUT SAM]: %s\n" %(TimeStamper(), each_cleaned_bamfile))
			each_rmdup_bamfile = Mark_Duplicates(each_cleaned_bamfile, picard_prog)
			bamfiles_per_RGID.append(each_rmdup_bamfile)
		processed_bamfile = os.path.join(each_rmdup_bamfile.strip(re.split('\/', each_rmdup_bamfile)[-1]), rgid+".processed.bam")
		if len(bamfiles_per_RGID) > 1 :
			sample_bamfile = os.path.join(each_rmdup_bamfile.strip(re.split('\/', each_rmdup_bamfile)[-1]), rgid+".merged.bam")
			timestamper = TimeStamper()
			for bamfile in bamfiles_per_RGID :
				sys.stdout.write("%s [INPUT SAM]: %s\n" %(timestamper, bamfile))
			sys.stdout.write("%s [MergeSamFiles.jar]: " %(TimeStamper()))
			Merge_Bams(bamfiles_per_RGID, sample_bamfile, picard_prog)
			sys.stdout.write("%s [INPUT SAM]: %s\n" %(TimeStamper(), sample_bamfile))
			sys.stdout.write("%s [AddOrReplaceReadGroups.jar]: " %(TimeStamper()))
			AddRGs(sample_bamfile, rg_def, processed_bamfile, picard_prog)
			bamfiles_all.append(processed_bamfile)
		else :
			sys.stdout.write("%s [INPUT SAM]: %s\n" %(TimeStamper(), each_rmdup_bamfile))
			sys.stdout.write("%s [AddOrReplaceReadGroups.jar]: " %(TimeStamper()))
			AddRGs(each_rmdup_bamfile, rg_def, processed_bamfile, picard_prog)
			bamfiles_all.append(processed_bamfile)
		sys.stdout.write("%s [Finish]\n" %(TimeStamper()))
		sys.stdout.write("*********************************\n")
		bamfiles_per_RGID = []
	if len(bamfiles_all) > 1 :
		if not merge_off :
			sys.stdout.write("\n%s [Merging all BAM files with different @RG]\n" %(TimeStamper()))
			timestamper = TimeStamper()
			for bamfile in bamfiles_all :
				sys.stdout.write("%s [INPUT SAM]: %s\n" %(timestamper, bamfile))
			sys.stdout.write("%s [MergeSamFiles.jar]: " %(TimeStamper()))
			Merge_Bams(bamfiles_all, outprefix+".all.processed.bam", picard_prog)

# parse config file and read in each to-be-processed SAM file #
def Config_Parser(config_file) :
	line_id, rg_id, rg_lb, rg_pl, rg_pu, rg_sm = "", "", "", "", "", ""			# line id is defined by users, and should be consistent with the each sam filename #
	num_diff_RGs = 0				# denote how many RG tags defined in the file #
	raw_samfiles = []
	config_dict = {}
	fCon = open(config_file, 'r')
	for line in fCon :
		if not line.startswith('#') :
			tmp_line = re.split('\t', line.strip())
			rg_def = tmp_line[0]
			RG_Def_Checker(rg_def)
			if len(tmp_line) > 2 :
				raw_samfiles = ' '.join(tmp_line[1:])
			else :
				raw_samfiles = tmp_line[1]
			Path_Checker(raw_samfiles, 'f')
			if not config_dict.has_key(rg_def) :
				config_dict[rg_def] = raw_samfiles
			else :
				sys.stdout.write(TimeStamper() + " [Config File Checker]	Warning: repeated @RG definition %s for the following SAM files" %(rg_def))
				sys.stdout.write(TimeStamper() + " [Config File CHecker]	Warning: existing %s" %(' '.join(config_dict[rg_def])))
				sys.stdout.write(TimeStamper() + " [Config File CHecker]	Warning: new %s" %(' '.join(raw_samfiles)))
				config_dict[rg_def] += raw_samfiles
	fCon.close()
	return config_dict

def Get_SAMs_from_CMD(samfiles) :
	'''
	fix this function
	'''
	contig_dict = {}

	return contig_dict

# generate dash-separated timestamp string
def TimeStamper() :
	return datetime.now().strftime('%Y-%m-%d-%H:%M:%S')

if __name__ == "__main__" :

	# parse the options and arguments #
	parser = argparse.ArgumentParser(description="Calling PICARD to post-process mapping results in sam files. Java 1.6+ by default is assumed to be installed")
	parser.add_argument("-config", metavar="FILE", dest="config_file", help="a config file providing paths to sam files and @RG tag information. Note: either -config or -s can be specified at a time. In cases where they are specified together, only sam files in the config file will be processed")
	parser.add_argument("-s", metavar="FILE", nargs='*', dest="samfiles", help="a list of sam files separated by space. Note: either -config or -s can be specified at a time. [NOT READY]")
	parser.add_argument("-RG", metavar="STR", nargs='*', dest="rgs", help="@RG definitions, e.g. @RG\tID:foo\tSM:foo\tLB:FRAG\tPL:Illumina\tPU:RUN. Note: if multiple sam files provided, then the same number of @RG definitions provided. If only one @RG defition provided, then the @RG will be added in each of the sam files. [NOT READY]")
	parser.add_argument("-picard", metavar="DIR", dest="picard_prog", help="the full path to the directory where PICARD is installed. Note: if the path to PICARD is not under your system path, please specify")
	parser.add_argument("-mergeOFF", action='store_true', dest="merge_off", help="toggle off the mergeSamFiles function")
	parser.add_argument("-p", metavar="FILE", dest="out_prefix",  help="the prefix of the output file. Note: if -mergeOFF is provided, no output will be created no matter what you provide after -p")

	args = parser.parse_args()

	out_prefix, outdir = "", ""
	if not args.merge_off :
		if not args.out_prefix :
			sys.stderr.write("[Argument Checker]	Error: no prefix for output files provided\n")
			sys.exit(1)
		else :
			if re.search('\/', args.out_prefix) :
				outdir = args.out_prefix.strip(re.split('\/', args.out_prefix)[-1])
				Path_Checker(outdir, 'd')
				out_prefix = args.out_prefix
			else :
				outdir = os.getcwd()
				out_prefix = os.path.join(outdir, args.out_prefix)

	contig_dict = {}
	if args.config_file and (args.samfiles or not args.samfiles) :
		if args.samfiles :
			sys.stdout.write("[Argument Checker]	WARNING: all SAM file provided after -s will be ignored")
			args.samfiles = []
		Path_Checker(args.config_file, 'f')
		config_dict = Config_Parser(args.config_file)
		Run_PICARD_Pipe(config_dict, args.merge_off, out_prefix, args.picard_prog)
	else :
		if not args.rgs :
			sys.stderr.write("[Argument Checker]	Error: no @RG definition found for your input SAM files\n")
			sys.exit(1)
		else :
			if len(args.rgs) > 1 :
				if len(args.rgs) != len(args.samfiles, args.rgs) :
					sys.stderr.write("[Argument Checker]	Error: number of @RG definition is more than the number of SAM files provided\n")
					sys.exit(1)
				else :
					config_dict = Get_SAMs_from_CMD(args.samfiles)
			else :
				if len(args.rgs) < len(args.samfiles) and len(args.rgs) == 1 :
					sys.stdout.write("[Argument Checker]	Warning: one @RG definition for all %d SAM files\n" %(len(args.samfiles)))
					RG_Def_Checker(args.rgs[0])
					contig_dict[args.rgs[0]] = args.samfiles
