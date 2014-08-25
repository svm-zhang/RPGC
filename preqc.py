import sys
import os
import re
import subprocess
import shlex
import multiprocessing
import argparse
import shutil

def clean_fastq(job_queue, db, nthread, tmp_dir):

	tmp_fasta_dir, tmp_decontam_dir, tmp_rmadapter_dir = tmp_dir
	while True:
		try:
			lib, fq1, fq2, blastskip, cutadaptskip, a1, a2 = job_queue.get()
			basename1 = os.path.basename(fq1)[:os.path.basename(fq1).index('.')]
			basename2 = os.path.basename(fq2)[:os.path.basename(fq2).index('.')]
			tmp_fasta_1 = os.path.join(tmp_fasta_dir, "%s.fasta" %(basename1))
			tmp_fasta_2 = os.path.join(tmp_fasta_dir, "%s.fasta" %(basename2))

			'''
				check the existence of fasta file
				if exists, validate the file
				if not, convert the one needs to be converted
			'''
			exist_1 = check_file_if_exists(tmp_fasta_1)
			exist_2 = check_file_if_exists(tmp_fasta_2)
			if exist_1 and exist_2:
				validate_1 = validate_fasta(fq1, tmp_fasta_1)
				validate_2 = validate_fasta(fq2, tmp_fasta_2)
#				if validate_1 and validate_2:
#					pass
				if not validate_1 and not validate_2:
					fastq2fasta({fq1:tmp_fasta_1, fq2:tmp_fasta_2})
				else:
					if not validate_1:
						fastq2fasta({fq1:tmp_fasta_1})
					if not validate_2:
						fastq2fasta({fq2:tmp_fasta_2})

			elif not exist_1 and not exist_2:
				fastq2fasta({fq1:tmp_fasta_1, fq2:tmp_fasta_2})
			elif not exist_1:
					fastq2fasta({fq1:tmp_fasta_1})
			elif not exist_2:
					fastq2fasta({fq2:tmp_fasta_2})

			tmp_blast_1 = os.path.join(tmp_decontam_dir, "%s.decontam.blast" %(basename1))
			tmp_blast_2 = os.path.join(tmp_decontam_dir, "%s.decontam.blast" %(basename2))
			if not blastskip:
				fail_blast = run_blast(db, {tmp_fasta_1:tmp_blast_1, tmp_fasta_2:tmp_blast_2}, nthread)
			else:
				sys.stdout.write("skip blast for %s\n" %(lib))
				if check_file_if_exists(tmp_blast_1) and check_file_if_exists(tmp_blast_2):
					fail_blast = 0
				else:
					fail_blast = 1
			if fail_blast:
				sys.stderr.write("Error: %s failed on blast for lib %s\n" %(multiprocessing.current_process().name, lib))
			else:
				decontam_fastq_1 = os.path.join(tmp_decontam_dir, "%s.decontam.fastq" %(basename1))
				decontam_fastq_2 = os.path.join(tmp_decontam_dir, "%s.decontam.fastq" %(basename2))
				rmadapter_fastq_1 = os.path.join(tmp_rmadapter_dir, "%s.decontam.rmadapter.fastq" %(basename1))
				rmadapter_fastq_2 = os.path.join(tmp_rmadapter_dir, "%s.decontam.rmadapter.fastq" %(basename2))
				info_1 = os.path.join(tmp_rmadapter_dir, "%s.cutadapt.info" %(basename1))
				info_2 = os.path.join(tmp_rmadapter_dir, "%s.cutadapt.info" %(basename2))
				sum_1 = os.path.join(tmp_rmadapter_dir, "%s.cutadapt.summary" %(basename1))
				sum_2 = os.path.join(tmp_rmadapter_dir, "%s.cutadapt.summary" %(basename2))
				if not cutadaptskip:
					remove_contamination(tmp_blast_1, fq1, decontam_fastq_1)
					remove_contamination(tmp_blast_2, fq2, decontam_fastq_2)
					remove_adapters(a1, info_1, sum_1, decontam_fastq_1, rmadapter_fastq_1)
					remove_adapters(a2, info_2, sum_2, decontam_fastq_2, rmadapter_fastq_2)
				else:
					sys.stdout.write("skip cutadapt for %s\n" %(lib))
				final_fastq_1 = os.path.join(os.path.join(tmp_rmadapter_dir, os.pardir), "%s.cleaned.fastq" %(basename1))
				final_fastq_2 = os.path.join(os.path.join(tmp_rmadapter_dir, os.pardir), "%s.cleaned.fastq" %(basename2))
				final_fastq_se = os.path.join(os.path.join(tmp_rmadapter_dir, os.pardir), "%s_se.cleaned.fastq" %(lib))
				fix_pairing_brutal(rmadapter_fastq_1, rmadapter_fastq_2, final_fastq_1, final_fastq_2, final_fastq_se)
		finally:
			job_queue.task_done()

def fastq2fasta(fq2fa_dict):
	for fq, fa in sorted(fq2fa_dict.iteritems()):
		sys.stdout.write("%s\tconverting %s to %s ...\n" %(multiprocessing.current_process().name, os.path.basename(fq), os.path.basename(fa)))
		fFA = open(fa, 'w')
		i = 0
		with open(fq, 'r') as fFQ:
			for line in fFQ:
				if i % 4 == 0:
					fFA.write(">%s" %(line[1:]))
					i += 1
				elif i % 4 == 1:
					fFA.write("%s" %(line))
					i += 1
				elif i % 4 == 2:
					i += 1
				else:
					i = 0
				fFA.flush()
		fFA.close()

def validate_fasta(fq, fa):
	cmd = "wc -l %s" %(fq)
	p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
	nreads_in_fq = int(p.communicate()[0].split(' ')[0])/4
	cmd = "grep -c \">\" %s" %(fa)
	p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
	nreads_in_fa = int(p.communicate()[0])
	# the two number must mach up if the conversion finishes successfully
	if nreads_in_fq != nreads_in_fa:
		return 0
	else:
		return 1

def run_blast(db, fa2blast_dict, nthread):
	returncode = 0
	for fa, out_blast in sorted(fa2blast_dict.iteritems()):
		run_blast_cmd = "blastn -db %s -query %s -outfmt 6 -num_threads %d -out %s " %(db, fa, nthread, out_blast)
		sys.stdout.write("%s\n" %(run_blast_cmd))
		p = subprocess.Popen(shlex.split(run_blast_cmd))
		p.wait()
		if p.returncode != 0:
			returncode = 1
	return returncode

def remove_adapters(adapter, info_file, sum_file, tmp_decontam_fastq, tmp_rmadapter_fastq):
	cutadapt_cmd = " cutadapt -b %s -O %d -m %d --info-file %s -o %s %s " %(adapter, len(adapter)/2, len(adapter)+1, info_file, tmp_rmadapter_fastq, tmp_decontam_fastq)
	sys.stdout.write(multiprocessing.current_process().name + "\t" + cutadapt_cmd + "\n")
	p = subprocess.Popen(shlex.split(cutadapt_cmd), stdout=open(sum_file, 'w'))
	p.wait()
	return p.returncode

def remove_contamination(blast_file, fastq_file, out_decontam_fastq):
	blast_hits = get_blast_hits(blast_file)
	sys.stdout.write(multiprocessing.current_process().name+"\t"+blast_file+"\t"+"%d\n" %(len(blast_hits)))
	sys.stdout.flush()
	fOUT = open(out_decontam_fastq, 'w')
	i = 0
	tmp = 0
	with open(fastq_file, 'r') as fFASTQ:
		for line in fFASTQ:
			if i % 4 == 0:
				if len(line.split(' ')) == 1:
					reads_id = line.strip('@')
				else:
					reads_id = line.split(' ')[0].strip('@')
				rm = False
				if reads_id in blast_hits:
					tmp += 1
					rm = True
					del blast_hits[reads_id]
				else:
					fOUT.write(line)
				i += 1
			elif i % 4 != 0:
				if not rm:
					fOUT.write(line)
				if i == 3:
					i = 0
				else:
					i += 1
			fOUT.flush()
	fOUT.close()
	sys.stdout.write(multiprocessing.current_process().name+"\t"+fastq_file+"\t"+"%s\t%s\n" %(len(blast_hits), tmp))
	sys.stdout.flush()

def get_blast_hits(blast_file):
	cmd_1 = "awk \'{print $1}\' %s" %(blast_file)
	cmd_2 = "uniq -"

	p1 = subprocess.Popen(shlex.split(cmd_1), stdout = subprocess.PIPE)
	p2 = subprocess.Popen(shlex.split(cmd_2), stdin = p1.stdout, stdout = subprocess.PIPE)
	hits = {}
	for hit in p2.communicate()[0].splitlines():
		hits[hit] = 1

	return hits			# return a dictionary of reads IDs who have blast hits (these reads probably come from contamination DNA)

def fix_pairing_brutal(fq1, fq2, oq1, oq2, ose):
	''' fix reads pairing '''
	sys.stdout.write("%s\tfix pairing for %s %s\n" %(multiprocessing.current_process().name, os.path.basename(fq1), os.path.basename(fq2)))
	fq1_dict = {}
	str_except_id = ""
	with open(fq1, 'r') as fQ1:
		for i, line in enumerate(fQ1):
			if i % 4 == 0:
				tmp_line = re.split("[ |\/]", line)
				id = tmp_line[0]
				str_except_id = tmp_line[1]
			else:
				str_except_id += line
				if i % 4 == 3:
					fq1_dict[id] = str_except_id
					str_except_id = ""

	fOQ1 = open(oq1, 'w')
	fOQ2 = open(oq2, 'w')
	fORPHAN = open(ose, 'w')
	str_except_id = ""
	num_paired, num_orphan = 0, 0
	with open(fq2, 'r') as fQ2:
		for i, line in enumerate(fQ2):
			if i % 4 == 0:
				tmp_line = re.split("[ |\/]", line)
				id = tmp_line[0]
				str_except_id = tmp_line[1]
			else:
				str_except_id += line
				if i % 4 == 3:
					if id in fq1_dict:
						fOQ1.write(id + " " + fq1_dict[id])
						fOQ2.write(id + " " + str_except_id)
						num_paired += 1
						fq1_dict.pop(id)
					else:
						fORPHAN.write(id + " " + str_except_id)
						num_orphan += 1
					str_except_id = ""
			fOQ1.flush()
			fOQ2.flush()
			fORPHAN.flush()

	for k, v in fq1_dict.iteritems():
		fORPHAN.write(k + " " + v)
		fORPHAN.flush()
		num_orphan += 1

	fOQ1.close()
	fOQ2.close()
	fORPHAN.close()

	sys.stdout.write("%s\tPAIRED: %d\tORPHAN: %d\n" %(multiprocessing.current_process().name, num_paired, num_orphan))

def add_jobs(configs, job_queue):
	for fCONFIG in configs:
		blastskip = -1
		cutadaptskip = -1
		for line in fCONFIG:
			if not line.startswith('#') or not line.strip():		# skip comments
				if line.startswith('['):
					lib = line.strip().lstrip('[').rstrip(']')
				elif line.startswith("dir"):
					dir = line.split('=')[1].strip()
				elif line.startswith('fastq'):
					fq1 = os.path.join(dir, line.split('=')[1].strip().split(',')[0])
					fq2 = os.path.join(dir, line.split('=')[1].strip().split(',')[1])
					if not check_file_if_exists(fq1):
						sys.stderr.write("[Error] Cannot find %s\n" %(fq1))
						sys.exit()
					if not check_file_if_exists(fq2):
						sys.stderr.write("[Error] Cannot find %s\n" %(fq2))
						sys.exit()
				elif line.startswith('blastskip'):
					blastskip = int(line.split('=')[1].strip())
				elif line.startswith('cutadaptskip'):
					cutadaptskip = int(line.split('=')[1].strip())
				elif line.startswith('adapter'):
					a1 = line.split('=')[1].strip().split(',')[0]
					a2 = line.split('=')[1].strip().split(',')[1]
					if blastskip == -1:
						blastskip = 0
					if cutadaptskip == -1:
						cutadaptskip = 0
					job_queue.put((lib, fq1, fq2, blastskip, cutadaptskip, a1, a2))
					blastskip = -1
					cutadaptskip = -1

def create_processes(nproc, job_queue, db, nthread, *tmp_dir):
	''' initiate processes '''
	for _ in range(nproc):
		process = multiprocessing.Process(target=clean_fastq,
					args=(job_queue, db, nthread, tmp_dir))
		process.daemon = True
		process.start()

def check_file_if_exists(file):
	''' check if file exists '''
	if os.path.exists(file):
		return 1
	else:
		return 0

def make_dirs_if_necessary(*dirs):
	''' create any number of directories if necessary '''
	for dir in dirs:
		if not os.path.exists(dir):
			os.makedirs(dir)

def rpgc_preqc(args):
	command = args.command
	target = args.target
	nproc = args.nproc
	nthread = args.nthread
	keep = args.keep
	db = args.db
	configs = args.configs
	print command
	print target
	print nproc
	print nthread
	print keep
	print db
	print configs
	tmp_fasta_dir = os.path.join(target, "tmp_fasta")
	tmp_decontam_dir = os.path.join(target, "tmp_decontam")
	tmp_rmadapter_dir = os.path.join(target, "tmp_rmadapter")
	make_dirs_if_necessary(target, tmp_fasta_dir, tmp_decontam_dir, tmp_rmadapter_dir)

	job_queue = multiprocessing.JoinableQueue()
	create_processes(nproc/nthread, job_queue, db, nthread, tmp_fasta_dir, tmp_decontam_dir, tmp_rmadapter_dir)
	add_jobs(configs, job_queue)
	try:
		job_queue.join()
	except KeyboardInterrupt:
		sys.stdout.write("Terminated by Keyboard\n")
		sys.exit()
	else:
		sys.stdout.write("Quit normally\n")

	# remove all the intermediate files
	if not keep:
		try:
			shutil.rmtree(tmp_decontam_dir)
			shutil.rmtree(tmp_rmadapter_dir)
			shutil.rmtree(tmp_fasta_dir)
		except OSError, e:
			if e.errno != 2:		# code 2 - no such file or directory
				raise
