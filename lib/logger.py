import logging
import sys

class HandlerFilter(logging.Filter) :
	def __init__(self, level) :
		self.level = level
	def filter(self, message) :
		return message.levelno == self.level

def logger_init(external_log_file = None) :
	logger = logging.getLogger("RPGC")
	logger.setLevel(logging.DEBUG)

	stdoutHandler = logging.StreamHandler(sys.stdout)
	stdoutHandler.setFormatter(logging.Formatter("%(asctime)s - %(name)s_STDOUT - INFO - %(message)s", datefmt="%Y-%m-%d %I:%m:%S %p"))
	stdoutHandler.addFilter(HandlerFilter(logging.DEBUG))

	stderrHandler = logging.StreamHandler(sys.stderr)
	stderrHandler.setFormatter(logging.Formatter("%(asctime)s - %(name)s_STDERR - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %I:%m:%S %p"))
	stderrHandler.addFilter(HandlerFilter(logging.ERROR))

	if external_log_file is not None :
		externalHandler = logging.FileHandler(external_log_file, 'w')
		externalHandler.addFilter(HandlerFilter(logging.INFO))
		externalHandler.setFormatter(logging.Formatter("%(asctime)s - EXTERNAL_PROG - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %I:%m:%S %p"))
		logger.addHandler(externalHandler)

	logger.addHandler(stdoutHandler)
	logger.addHandler(stderrHandler)

	return logger
