import os, sys
import re
import argparse

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
				raw_samfiles = tmp_line[1:]
			else :
				raw_samfiles = [tmp_line[1]]
			print raw_samfiles
			Path_Checker(raw_samfiles, 'f')
			if not config_dict.has_key(rg_def) :
				config_dict[rg_def] = raw_samfiles
			else :
				print "[@RG Definition Checker]	Warning: repeated @RG definition %s for the following SAM files" %(rg_def)
				print "[@RG Definition CHecker]	Warning: existing %s" %(' '.join(config_dict[rg_def]))
				print "[@RG Definition CHecker]	Warning: new %s" %(' '.join(raw_samfiles))
				config_dict[rg_def] += raw_samfiles
	fCon.close()
	print config_dict
	return config_dict

def clean_sam_files(raw_samfile, picard_prog) :
	cleaned_samfile = '.'.join(re.split('\.', raw_samfile)[0:3]) + ".cleaned.sam"
	picard_clean_sam_cmd = "java -Xmx4g -jar %s/CleanSam.jar INPUT=%s OUTPUT=%s VALIDATION_STRINGENCY=STRICT" %(picard_prog, raw_samfile, cleaned_samfile)
	print picard_clean_sam_cmd
	os.system(picard_clean_sam_cmd)
	return cleaned_samfile

def sam_to_bam(cleaned_samfile, picard_prog) :
	cleaned_bamfile = '.'.join(re.split('\.', cleaned_samfile)[0:3]) + ".cleaned.bam"
	picard_samtobam_cmd = "java -Xmx4g -jar %s/SamFormatConverter.jar INPUT=%s OUTPUT=%s VALIDATION_STRINGENCY=LENIENT" %(picard_prog, cleaned_samfile, cleaned_bamfile)
	print picard_samtobam_cmd
	os.system(picard_samtobam_cmd)
	os.system("rm %s" %(cleaned_samfile))
	return cleaned_bamfile

def merge_bams(list_bamfiles, outfile, picard_prog) :
	input_samfiles = ""
	for each_file in list_bamfiles :
		if input_samfiles == "" :
			input_samfiles += "INPUT=%s" %(each_file)
		else :
			input_samfiles += " INPUT=%s" %(each_file)
	picard_mergesam_cmd = "java -Xmx4g -jar %s/MergeSamFiles.jar %s OUTPUT=%s SORT_ORDER=coordinate USE_THREADING=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true" %(picard_prog, input_samfiles, outfile)
	print picard_mergesam_cmd
	os.system(picard_mergesam_cmd)

def addRG_tag(bamfile, rg_def, outfile, picard_prog) :
	picard_addRG_cmd = "java -Xmx4g -jar %s/AddOrReplaceReadGroups.jar INPUT=%s OUTPUT=%s %s SORT_ORDER=coordinate CREATE_INDEX=true" %(picard_prog, bamfile, outfile, re.sub('\:', '=', ' '.join(re.split(r'\\t', rg_def)[1:])))
	print picard_addRG_cmd
	os.system(picard_addRG_cmd)

# get rid of the duplicates #
def mark_duplicates(cleaned_bamfile, picard_prog) :
	sorted_bamfile = '.'.join(re.split('\.', cleaned_bamfile)[0:3]) + ".cleaned.sort.bam"
	picard_sort_cmd = "java -Xmx4g -jar %s/SortSam.jar I=%s O=%s SO=coordinate VALIDATION_STRINGENCY=LENIENT" %(picard_prog, cleaned_bamfile, sorted_bamfile)
	print picard_sort_cmd
	os.system(picard_sort_cmd)
	dup_log_file = '.'.join(re.split('\.', cleaned_bamfile)[0:3]) + ".dup.report"
	dedup_bamfile = '.'.join(re.split('\.', sorted_bamfile)[0:3]) + ".cleaned.sorted.dedup.bam"
	picard_markdup_cmd = "java -Xmx4g -jar %s/MarkDuplicates.jar REMOVE_DUPLICATES=true I=%s O=%s VALIDATION_STRINGENCY=LENIENT M=%s" %(picard_prog, sorted_bamfile, dedup_bamfile, dup_log_file)
	print picard_markdup_cmd
	os.system(picard_markdup_cmd)
	return dedup_bamfile

def Path_Checker(path, type) :
	if type == 'f' :
		for i in range(len(path)) :
			print path[i]
			if not os.path.exists(path[i]) :
				sys.stderr.write("[Path_Checker]	Error: cannnot reach the file: %s" %(path[i]))
				sys.exit(1)
			else :
				print "[Path_Checker]	%s ... [checked]" %(path[i])
	elif type == 'd' :
		print "checking path to directory"

def RG_Def_Checker(rg_def) :
	tmp_def = re.split(r'\\t', rg_def)
	num_component = 0
	for i in range(1, len(tmp_def)) :
		if len(re.split(':', tmp_def[i])) == 2 :
			tmp_component = re.split(':', tmp_def[i])
			if tmp_component[0] == "" or not tmp_component[0] in ["ID", "SM", "LB", "PL", "PU"] :
				sys.stderr.write("Error: missing @RG key %s\n") %(tmp_def[i])
				sys.stderr.write("Example: @RG\tID:cs1\tSM:cs1\tLB:FRAG\tPL:Illumina\tPU:RUN\n")
				sys.exit(1)
			else :
				if tmp_component[1] == "" :
					sys.stderr.write("Error: missing value for @RG_%s %s\n") %(tmp_component[0], tmp_def[i])
					sys.exit(1)
				else :
					num_component += 1
		else :
			sys.stderr.write("Error: incorrect @RG format %s\n") %(tmp_def[i])
			sys.stderr.write("Example: @RG\tID:cs1\tSM:cs1\tLB:FRAG\tPL:Illumina\tPU:RUN\n")
			sys.exit(1)
	if num_component != 5 :
		sys.stderr.write("Error: missing @RG keys in the definition %s\n" %(rg_def))
		sys.exit(1)

def Run_PICARD_Pipe(config_dict, merge_or_not, outprefix, picard_prog) :
	outdir = outprefix.strip(re.split('\/', outprefix)[-1])
	bamfiles_per_sample = []
	bamfiles_multi_samples = []
	for rg_def, samfiles in config_dict.iteritems() :
		rg_id = ""
		tmp_def = re.split(r'\\t', rg_def)
		for i in range(1, len(tmp_def)) :
			tmp_component = re.split(':', tmp_def[i])
			if tmp_component[0] == "ID" :
				rgid = tmp_component[1]
		for each_raw_samfile in samfiles :
			each_cleaned_samfile = clean_sam_files(each_raw_samfile, picard_prog)
			each_cleaned_bamfile = sam_to_bam(each_cleaned_samfile, picard_prog)
			each_dedup_bamfile = mark_duplicates(each_cleaned_bamfile, picard_prog)
			bamfiles_per_sample.append(each_dedup_bamfile)
		processed_bamfile = os.path.join(each_dedup_bamfile.strip(re.split('\/', each_dedup_bamfile)[-1]), rgid+".processed.bam")
		if len(bamfiles_per_sample) > 1 :
			sample_bamfile = os.path.join(each_dedup_bamfile.strip(re.split('\/', each_dedup_bamfile)[-1]), rgid+".merged.bam")
			merge_bams(bamfiles_per_sample, sample_bamfile, picard_prog)
			addRG_tag(sample_bamfile, rg_def, processed_bamfile, picard_prog)
			bamfiles_multi_samples.append(processed_bamfile)
		else :
			addRG_tag(each_dedup_bamfile, rg_def, processed_bamfile, picard_prog)
			bamfiles_multi_samples.append(processed_bamfile)
		bamfiles_per_sample = []
	if len(bamfiles_multi_samples) > 1 :
		if not merge_or_not :
			merged_multi_sample_bamfile = merge_bams(bamfiles_multi_samples, outprefix+".all.processed.bam", picard_prog)

if __name__ == "__main__" :

	# parse the options and arguments #
	parser = argparse.ArgumentParser(description="calling PICARD to post-process mapping results in sam files. Java by default is assumed installed")
	parser.add_argument("-config", metavar="FILE", dest="config_file", help="a config file providing paths to sam files and @RG tag information. Note: either -config or -s can be specified at a time. In cases where they are specified together, only sam files in the config file will be processed")
	parser.add_argument("-s", metavar="FILE", nargs='*', dest="samfiles", help="a list of sam files separated by space. Note: either -config or -s can be specified at a time.")
	parser.add_argument("-RG", metavar="STR", nargs='*', dest="rgs", help="@RG definitions, e.g. @RG\tID:foo\tSM:foo\tLB:FRAG\tPL:Illumina\tPU:RUN. Note: if multiple sam files provided, then the same number of @RG definitions provided. If only one @RG defition provided, then the @RG will be added in each of the sam files")
	parser.add_argument("-picard", metavar="DIR", dest="picard_prog", help="the full path to the directory where PICARD is installed. Note: if the path to PICARD is not under your system path, please specify")
	parser.add_argument("-mergeOFF", action='store_true', dest="merge_or_not", help="toggle off the mergeSamFiles function")
	parser.add_argument("-p", metavar="FILE", dest="out_prefix",  help="the prefix of the output file. Note: if -mergeOFF is provided, no output will be created no matter what you provide after -p")

	args = parser.parse_args()

	if not args.merge_or_not :
		if not args.out_prefix :
			sys.stderr.write("Error: no prefix for output files provided\n")
			sys.exit(1)
		else :
			outdir = args.out_prefix.strip(re.split('\/', args.out_prefix)[-1])
			if not os.path.exists(outdir) :
				os.makedirs(outdir)

	contig_dict = {}
	if args.config_file and (args.samfiles or not args.samfiles) :
		if args.samfiles :
			print "WARNING: all sam file provided after -s will be ignored"
			args.samfiles = []
		Path_Checker([args.config_file], 'f')
		config_dict = Config_Parser(args.config_file)
		Run_PICARD_Pipe(config_dict, args.merge_or_not, args.out_prefix, args.picard_prog)
	else :
		if not args.rgs :
			sys.stderr.write("[@RG Definition Checker]	Error: no @RG definition found for your input SAM files\n")
			sys.exit(1)
		else :
			if len(args.rgs) > 1 :
				if len(args.rgs) != len(args.samfiles, args.rgs) :
					sys.stderr.write("[@RG Definition Checker]	Error: number of @RG definition is more than the number of SAM files provided\n")
					sys.exit(1)
				else :
					Get_SAMs_from_CMD(args.samfiles)
			else :
				if len(args.rgs) < len(args.samfiles) and len(args.rgs) == 1 :
					sys.stdout.write("[@RG Definition Checker]	Warning: one @RG definition for all %d SAM files\n" %(len(args.samfiles)))
					RG_Def_Checker(args.rgs[0])
					contig_dict[args.rgs[0]] = args.samfiles
		print args.samfiles
