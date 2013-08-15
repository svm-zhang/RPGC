from os import makedirs
from os.path import exists
from datetime import datetime

def check_file_existence(logger, verbosity, *filenames) :
	""" check the existence of a set of files """
	for filename in filenames :
		if not exists(filename) :
			logger.error("[IO] cannot find the file: %s" %(filename))
		else :
			if verbosity :
				logger.debug("[IO] %s ... ... [checked]" %(filename))

def make_dir_if_not_exist(*dirnames) :
	""" make a set of directories (including those intermediate directories) if not existing """
	for dirname in dirnames :
		if not exists(dirname) :
			makedirs(dirname)

def time_stamper():
	"""Returns dash-separated timestamp string"""
	return datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
