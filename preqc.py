import sys
import os
import re
import subprocess
import shlex
import multiprocessing
import argparse
import shutil

sys.path.insert(0, os.path.join(os.getcwd(), "lib"))
from io import make_dir_if_necessary
from io import make_dirs_if_necessary
from io import check_if_files_exist

class Preqc:
	def __init__(self, args):
		self.target = os.path.realpath(args.target)
		make_dir_if_necessary(self.target)
		self.nproc = args.nproc
		self.nthread = args.nthread
		self.keep = args.keep
		self.db = args.db
		self.min_eval = args.min_eval
		self.config = args.config
		if self.db is not None:
			check_if_files_exist("Preqc", self.config, self.db)
		else:
			check_if_files_exist("Preqc", self.config)
		self.min_overlap = args.min_overlap
		self.min_readlen = args.min_readlen

	def _print_cmd_setting(self):
		sys.stdout.write("[Preqc] Config file(s): %s\n" %(self.config))
		sys.stdout.write("[Preqc] Target directory: %s\n" %(self.target))
		sys.stdout.write("[Preqc] Number of Preqc jobs running simultaneously: %d\n" %(self.nproc))
		if self.db is not None:
			sys.stdout.write("[Preqc] BLAST database: %s\n" %(self.db))
			sys.stdout.write("[Preqc] Number of threads running BLAST: %d\n" %(self.nthread))
		sys.stdout.write("[Preqc] Minimum overlap between adapter and read: %d\n" %(self.min_overlap))
		sys.stdout.write("[Preqc] Minimum read length to kepp after Cutadapt: %d\n" %(self.min_readlen))

	def _fastq2fasta(self, fq2fa_dict):
		for fq, fa in sorted(fq2fa_dict.iteritems()):
			sys.stdout.write("[Preqc] %s\tconverting %s to %s ...\n" %(multiprocessing.current_process().name, os.path.basename(fq), os.path.basename(fa)))
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

	def _validate_fasta(fq, fa):
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

	def _run_blast(self, fa2blast_dict):
		returncode = 0
		for fa, out_blast in sorted(fa2blast_dict.iteritems()):
			run_blast_cmd = "blastn -db %s -query %s -outfmt 6 -evalue %s -num_threads %d -out %s " %(self.db, fa, self.min_eval, self.nthread, out_blast)
			sys.stdout.write("%s\n" %(run_blast_cmd))
			p = subprocess.Popen(shlex.split(_run_blast_cmd))
			p.wait()
			if p.returncode != 0:
				returncode = 1
		return returncode

	def _remove_adapters(self, adapter, info_file, sum_file, tmp_decontam_fastq, tmp_rmadapter_fastq):
		sys.stdout.write("[Preqc] %s removing adapters from %s" %(multiprocessing.current_process().name, tmp_decontam_fastq))
		cutadapt_cmd = " cutadapt -b %s -O %d -m %d --info-file %s -o %s %s " %(adapter, self.min_overlap, self.min_readlen, info_file, tmp_rmadapter_fastq, tmp_decontam_fastq)
		sys.stdout.write(multiprocessing.current_process().name + "\t" + cutadapt_cmd + "\n")
		p = subprocess.Popen(shlex.split(cutadapt_cmd), stdout=open(sum_file, 'w'))
		p.wait()
		return p.returncode

	def _remove_contamination(self, blast_file, fastq_file, out_decontam_fastq):
		blast_hits = self._get_blast_hits(blast_file)
		sys.stdout.write("[Preqc] %s removing %d contaminations from %s" %(multiprocessing.current_process().name, len(blast_hits), fastq_file))
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
		sys.stdout.flush()

	def _get_blast_hits(self, blast_file):
		cmd_1 = "awk \'{print $1}\' %s" %(blast_file)
		cmd_2 = "uniq -"

		p1 = subprocess.Popen(shlex.split(cmd_1), stdout = subprocess.PIPE)
		p2 = subprocess.Popen(shlex.split(cmd_2), stdin = p1.stdout, stdout = subprocess.PIPE)
		hits = {}
		for hit in p2.communicate()[0].splitlines():
			hits[hit] = 1

		return hits			# return a dictionary of reads IDs who have blast hits (these reads probably come from contamination DNA)

	def _fix_pairing_brutal(self, fq1, fq2, oq1, oq2, ose):
		''' fix reads pairing '''
		sys.stdout.write("[Preqc] %s\tfixing reads order for %s %s\n" %(multiprocessing.current_process().name, os.path.basename(fq1), os.path.basename(fq2)))
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

	def _outsource_tasks(self, task_q):
		with open(self.config, 'r') as fCONFIG:
			blastskip = -1
			cutadaptskip = -1
			decontamskip = -1
			for line in fCONFIG:
				if not line.startswith('#') or not line.strip():		# skip comments
					if line.startswith('['):
						lib = line.strip().lstrip('[').rstrip(']')
					elif line.startswith("dir"):
						dir = line.split('=')[1].strip()
					elif line.startswith('fastq'):
						fq1 = os.path.join(dir, line.split('=')[1].strip().split(',')[0])
						fq2 = os.path.join(dir, line.split('=')[1].strip().split(',')[1])
						check_if_files_exist("Preqc", fq1, fq2)
					elif line.startswith('decontamskip'):
						decontamskip = int(line.split('=')[1].strip())
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
						if decontamskip == -1:
							decontamskip = 0
						task_q.put((lib, fq1, fq2, blastskip, decontamskip, cutadaptskip, a1, a2))
						blastskip = -1
						cutadaptskip = -1
						decontamskip = -1

	def _init_processes(self, task_q, *tmp_dir):
		''' initiate processes '''
		for _ in range(self.nproc):
			process = multiprocessing.Process(target=self._preqc,
						args=(task_q, tmp_dir))
			process.daemon = True
			process.start()

	def _preqc(self, task_q, tmp_dir):
		tmp_fasta_dir, tmp_decontam_dir, tmp_rmadapter_dir = tmp_dir
		while True:
			try:
				lib, fq1, fq2, blastskip, decontamskip, cutadaptskip, a1, a2 = task_q.get()
				basename1 = os.path.basename(fq1)[:os.path.basename(fq1).index('.')]
				basename2 = os.path.basename(fq2)[:os.path.basename(fq2).index('.')]
				tmp_fasta_1 = os.path.join(tmp_fasta_dir, "%s.fasta" %(basename1))
				tmp_fasta_2 = os.path.join(tmp_fasta_dir, "%s.fasta" %(basename2))

				'''
					check the existence of fasta file
					if exists, validate the file
					if not, convert the one needs to be converted
				'''
				blast_success = -1
				tmp_blast_1 = os.path.join(tmp_decontam_dir, "%s.decontam.blast" %(basename1))
				tmp_blast_2 = os.path.join(tmp_decontam_dir, "%s.decontam.blast" %(basename2))
				if not blastskip:
					exist_1 = os.path.exists(tmp_fasta_1)
					exist_2 = os.path.exists(tmp_fasta_2)
					if exist_1 and exist_2:
						validate_1 = self._validate_fasta(fq1, tmp_fasta_1)
						validate_2 = self._validate_fasta(fq2, tmp_fasta_2)
						if not validate_1 and not validate_2:
							self._fastq2fasta({fq1:tmp_fasta_1, fq2:tmp_fasta_2})
						else:
							if not validate_1:
								self._fastq2fasta({fq1:tmp_fasta_1})
							if not validate_2:
								self._fastq2fasta({fq2:tmp_fasta_2})
					elif not exist_1 and not exist_2:
						self._fastq2fasta({fq1:tmp_fasta_1, fq2:tmp_fasta_2})
					elif not exist_1:
							self._fastq2fasta({fq1:tmp_fasta_1})
					elif not exist_2:
							self._fastq2fasta({fq2:tmp_fasta_2})
					blast_success = self._run_blast({tmp_fasta_1:tmp_blast_1, tmp_fasta_2:tmp_blast_2})
				else:
					sys.stdout.write("[Preqc] Turn off BLAST for %s\n" %(lib))
					if os.path.exists(tmp_blast_1) and os.path.exists(tmp_blast_2):
						blast_success = 1			# this check is weak!
					else:
						blast_success = 0
				if not blast_success:
					sys.stderr.write("[Preqc] Error: %s BLAST failed on lib %s\n" %(multiprocessing.current_process().name, lib))
				else:
					decontam_fastq_1 = os.path.join(tmp_decontam_dir, "%s.decontam.fastq" %(basename1))
					decontam_fastq_2 = os.path.join(tmp_decontam_dir, "%s.decontam.fastq" %(basename2))
					if not decontamskip:
						self._remove_contamination(tmp_blast_1, fq1, decontam_fastq_1)
						self._remove_contamination(tmp_blast_2, fq2, decontam_fastq_2)
					else:
						sys.stdout.write("[Preqc] Turn off remove contamination for %s\n" %(lib))
					rmadapter_fastq_1 = os.path.join(tmp_rmadapter_dir, "%s.decontam.rmadapter.fastq" %(basename1))
					rmadapter_fastq_2 = os.path.join(tmp_rmadapter_dir, "%s.decontam.rmadapter.fastq" %(basename2))
					info_1 = os.path.join(tmp_rmadapter_dir, "%s.cutadapt.info" %(basename1))
					info_2 = os.path.join(tmp_rmadapter_dir, "%s.cutadapt.info" %(basename2))
					sum_1 = os.path.join(tmp_rmadapter_dir, "%s.cutadapt.summary" %(basename1))
					sum_2 = os.path.join(tmp_rmadapter_dir, "%s.cutadapt.summary" %(basename2))
					if not cutadaptskip:
						self._remove_adapters(a1, info_1, sum_1, decontam_fastq_1, rmadapter_fastq_1)
						self._remove_adapters(a2, info_2, sum_2, decontam_fastq_2, rmadapter_fastq_2)
					else:
						sys.stdout.write("[Preqc] Turn off Cutadapt for %s\n" %(lib))
					final_fastq_1 = os.path.join(os.path.join(tmp_rmadapter_dir, os.pardir), "%s.cleaned.fastq" %(basename1))
					final_fastq_2 = os.path.join(os.path.join(tmp_rmadapter_dir, os.pardir), "%s.cleaned.fastq" %(basename2))
					final_fastq_se = os.path.join(os.path.join(tmp_rmadapter_dir, os.pardir), "%s_se.cleaned.fastq" %(lib))
					self._fix_pairing_brutal(rmadapter_fastq_1, rmadapter_fastq_2, final_fastq_1, final_fastq_2, final_fastq_se)
			finally:
				task_q.task_done()

	def start(self):
		self._run()

	def _run(self):
		self._print_cmd_setting()
		tmp_fasta_dir = os.path.join(self.target, "tmp_fasta")
		tmp_decontam_dir = os.path.join(self.target, "tmp_decontam")
		tmp_rmadapter_dir = os.path.join(self.target, "tmp_rmadapter")
		make_dirs_if_necessary(tmp_fasta_dir, tmp_decontam_dir, tmp_rmadapter_dir)

		task_q = multiprocessing.JoinableQueue()
		self._init_processes(task_q, tmp_fasta_dir,
							 tmp_decontam_dir, tmp_rmadapter_dir)
		self._outsource_tasks(task_q)
		try:
			task_q.join()
		except KeyboardInterrupt:
			sys.stdout.write("[Preqc] Terminated by Keyboard\n")
			sys.exit()
		else:
			# remove all the intermediate files
			if not self.keep:
				try:
					shutil.rmtree(tmp_decontam_dir)
					shutil.rmtree(tmp_rmadapter_dir)
					shutil.rmtree(tmp_fasta_dir)
				except OSError, e:
					if e.errno != 2:		# code 2 - no such file or directory
						raise
			sys.stdout.write("[Preqc] Quit normally\n")

