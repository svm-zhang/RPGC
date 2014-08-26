import os
import sys
import datetime

def check_if_file_exists(file, prg):
	"""check if file exists"""
	if not os.path.exists(file):
		sys.stderr.write("[%s] Error: cannot find %s\n" %(prg, file))
		sys.exit(1)

def check_if_files_exist(prg, *files):
	for file in files:
		check_if_file_exists(file, prg)

def make_dir_if_necessary(dir):
	"""make a set of directories
	   (including those intermediate directories)
	   if not existing
	"""
	if not os.path.exists(dir) :
		os.makedirs(dir)

def make_dirs_if_necessary(*dirs):
	for dir in dirs:
		make_dir_if_necessary(dir)

def time_stamper():
	"""Returns dash-separated timestamp string"""
	return datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
