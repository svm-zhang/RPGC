import re
import yaml
from sys import exit
from os import mkdir, remove
from os.path import exists, basename, dirname, join
from shlex import split
from argparse import ArgumentParser
from subprocess import Popen, PIPE

from logger import logger_init
from io import check_file_existence, make_dir_if_not_exist, time_stamper

def indexing_ref(index_cmd, db, ref_in_fasta, logger) :
	"""
	index the give reference
	"""
	make_dir_if_not_exist(dirname(db))
	print index_subprogram
	index_cmd += "%s -p %s %s" %(index_subprogram, db, ref_in_fasta)
	run(index_cmd, logger)

def merge_command_setting(default_customized_combo_dict, customized_combo_dict) :
	"""
	merge each subprogram's default setting with user customized setting
	"""
	if isinstance(default_customized_combo_dict, dict) and isinstance(customized_combo_dict, dict) :
		for subprogram, settings in default_customized_combo_dict.iteritems() :
			if subprogram not in customized_combo_dict :
				customized_combo_dict[subprogram] = settings
			else :
				customized_combo_dict[subprogram] = merge_command_setting(settings, customized_combo_dict[subprogram])
	return customized_combo_dict

def parse_config(config_file, logger, debug) :
	"""
	Parse config file
	"""
	try:
		params = yaml.load(file(config_file, 'r'))
	except yaml.YAMLError, exc:
		if hasattr(exc, "problem_mark") :
			mark = exc.problem_mark
			logger.error("[CONFIG] config syntax error at position: (%s,%s)" %(mark.line, mark.column+1))
			exit()

	# get mapping program #
	mapper_prog = params["PROGRAM"]

	# get command line default setting #
	default_combo_dict = {}
	for customized_combo, default_setting in params["DEFAULT_SETTING"].iteritems() :
		default_combo_dict[customized_combo] = default_setting

	# get output info #
	out_root_dir = params["OUTPUT_FOLDER"]
	dated_output = int(params["DATED_OUTPUT"])			# this creates timestamped subfolders under the out_root_dir if set True #
	structure = params["STRUCTURE"]
	if dated_output :
		out_root_dir = join(out_root_dir, time_stamper())
	make_dir_if_not_exist(out_root_dir)

	# get indexing info #
	ref_db, ref_fasta, skip, index_subprogram = "", "", -1, ""
	for k, v in params["INDEXING"].iteritems() :
		if k == "DB" :
			ref_db = v
		elif k == "FASTA" :
			ref_fasta = v
		elif k == "SKIP" :
			skip = v
		elif k == "SUBPROGRAM" :
			index_subprogram = mapper_prog + " %s" %(v)
	if not skip :
		check_file_existence(logger, debug, ref_fasta)
		indexing_ref(index_subprogram, ref_db, ref_fasta, logger)

	# get data and customized the command line setting #
	for each_group_id, each_group_info in params["EXECUTION_GROUPS"].iteritems() :
		tmp_group_outdir = out_root_dir
		if structure == "by_group" :
			tmp_group_outdir = join(out_root_dir, each_group_id)

		customized_combo_dict = {}
		customized_cmd_dict = {}
		for customized_combo, customized_combo_setting in each_group_info["CUSTOMIZED_SETTING"].iteritems() :
			if customized_combo in default_combo_dict :
				if customized_combo_setting is not None :
					customized_combo_dict[customized_combo] = merge_command_setting(default_combo_dict[customized_combo], customized_combo_setting)
				else :
					customized_combo_dict[customized_combo] = default_combo_dict[customized_combo]
				customized_cmd_dict[customized_combo] = make_cmd(customized_combo_dict[customized_combo], mapper_prog)
		for indv, data in each_group_info["DATA"].iteritems() :
			tmp_indv_outdir = join(tmp_group_outdir, indv)
			make_dir_if_not_exist(tmp_indv_outdir)
			se_fastqs, pe_fastqs = [], []
			for data_type, fastq_file in data.iteritems() :
				check_file_existence(logger, debug, fastq_file)
				if data_type == "SE" :
					se_fastqs.append(fastq_file)
				else :
					pe_fastqs.append(fastq_file)

			for combo, combo_setting in customized_cmd_dict.iteritems() :
				if combo == "BWA_SHORT" :
					kickoff_se_bwa_aln(ref_db, indv, se_fastqs, tmp_indv_outdir, combo_setting["aln"], combo_setting["samse"], logger)
					kickoff_pe_bwa_aln(ref_db, indv, pe_fastqs, tmp_indv_outdir, combo_setting["aln"], combo_setting["sampe"], logger)
				if combo == "BWA_MEM" :
					kickoff_bwa_mem(ref_db, indv, se_fastqs, tmp_indv_outdir, combo_setting["mem"], logger, "se")
					kickoff_bwa_mem(ref_db, indv, pe_fastqs, tmp_indv_outdir, combo_setting["mem"], logger, "pe")

def make_cmd(customized_combo_dict, mapper_prog) :
	"""
	make command line string for each called subprogram
	"""
	customized_cmd_dict = {}
	for subprogram, subprogram_setting in customized_combo_dict.iteritems() :
		if subprogram_setting is not None :
			for argument, value in subprogram_setting.iteritems() :
				if not customized_cmd_dict.has_key(subprogram) :
					if value is not None :
						customized_cmd_dict[subprogram] = mapper_prog + " %s %s %s" %(subprogram, argument, value)
					else :
						customized_cmd_dict[subprogram] = mapper_prog + " %s %s" %(subprogram, argument)
				else :
					if value is not None :
						customized_cmd_dict[subprogram] += " %s %s" %(argument, value)
					else :
						customized_cmd_dict[subprogram] += " %s" %(argument)
		else :
			if not customized_cmd_dict.has_key(subprogram) :
				customized_cmd_dict[subprogram] = mapper_prog + " %s" %(subprogram)
	#customized_cmd_dict = tmp_cmd_dict
	return customized_cmd_dict

def kickoff_bwa_mem(ref_db, indv, fastqs, outdir, mem_cmd, logger, data_type) :
	if data_type == "se" :
		for fastq in fastqs :
			sam_file = join(outdir, "%s_se.bwa_mem.sam" %(indv))
			mem_cmd += " %s %s" %(ref_db, fastq)
			logger.debug(mem_cmd)
			run(mem_cmd, logger, sam_file)
	elif data_type == "pe" :
		tmp_fastqs = ""
		for i in range(len(fastqs)) :
			tmp_fastqs += fastqs[i] + " "
			if i % 2 == 1 :
				mem_cmd += " %s %s" %(ref_db, tmp_fastqs)
				logger.debug(mem_cmd)
				sam_file = join(outdir, "%s_pe.bwa_mem.sam" %(indv))
				run(mem_cmd, logger, sam_file)
				tmp_fastqs = ""

def kickoff_pe_bwa_aln(ref_db, indv, fastqs, outdir, aln_cmd, sampe_cmd, logger) :
	sai_files = []
	for i in range(len(fastqs)) :
		sai_file = join(outdir, basename(fastqs[i]).strip(".fastq").strip(".fq")+".sai")
		tmp_aln_cmd = aln_cmd + " %s %s" %(ref_db, fastqs[i])
		logger.debug(tmp_aln_cmd)
		run(tmp_aln_cmd, logger, sai_file)
		sai_files.append(sai_file)
		if i % 2 == 1 :
			sam_file = join(outdir, "%s_pe.bwa_aln.sam" %(indv))
			sampe_cmd += " %s %s %s" %(ref_db, ' '.join(sai_files), ' '.join(fastqs))
			logger.debug(sampe_cmd)
			run(sampe_cmd, logger, sam_file)
			sai_files = []

def kickoff_se_bwa_aln(ref_db, indv, fastqs, outdir, aln_cmd, samse_cmd, logger) :
	for fastq in fastqs :
		sai_file = join(outdir, basename(fastq).strip(".fastq").strip(".fq")+".sai")
		aln_cmd += " %s %s" %(ref_db, fastq)
		logger.debug(aln_cmd)
		run(aln_cmd, logger, sai_file)
		sam_file = join(outdir, "%s_se.bwa_aln.sam" %(indv))
		samse_cmd += " %s %s %s" %(ref_db, sai_file, fastq)
		logger.debug(samse_cmd)
		run(samse_cmd, logger, sam_file)

def run(cmd, logger, outfile = None) :
	if outfile is not None :
		fH = file(outfile, 'w')
		proc = Popen(split(cmd), stdout=fH, stderr=PIPE)
		proc_stderr = proc.communicate()[1]
	else :
		proc = Popen(split(cmd), stdout=PIPE, stderr=PIPE)
		proc_stdout, proc_stderr = proc.communicate()
	if proc.returncode != 0 :
		logger.error(proc_stderr)
		exit()
	else :
		if proc_stderr != "" :
			logger.info(proc_stderr)

if __name__ == "__main__" :
	parser = ArgumentParser(description="Call BWA/NovoAlign to map reads to a given reference")
	parser.add_argument("-config", metavar="FILE", dest="config_file", help="a configuration file with paths to reference, fastq files.")
	parser.add_argument("-mapper_log", metavar="FILE", dest="mapper_log_file", help="log file generated by any mapper programs")
	parser.add_argument("--debug", action="store_true", dest="debug", help="for debug purpose, more information is spit out")

	args = parser.parse_args()

	# initiate logger #
	logger = logger_init(args.mapper_log_file)

	check_file_existence(logger, args.debug, args.config_file)

	# parse the config file #
	parse_config(args.config_file, logger, args.debug)
