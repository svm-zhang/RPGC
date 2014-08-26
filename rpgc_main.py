# this script is a first attempt to use argparse subparser

import argparse
import os
import sys

#from preqc import rpgc_preqc

class SubcommandHelpFormatter(argparse.RawDescriptionHelpFormatter):
	def _format_action(self, action):
		parts = super(argparse.RawDescriptionHelpFormatter, self)._format_action(action)
		if action.nargs == argparse.PARSER:
			sub_cmd = "\n"
			for i, part in enumerate(parts.split("\n")):
				if i == 0:
					continue
				else:
					sub_cmd += part + "\n"
			return sub_cmd
		return parts

def handle_cmd():
	main_parser = argparse.ArgumentParser(description="try argparse subparser", formatter_class=SubcommandHelpFormatter)
	sub_parsers = main_parser.add_subparsers(title="Commands", metavar="[command]", dest="command")

	preqc_parser = sub_parsers.add_parser("preqc", help="housekeeping the data like crazy")
	preqc_parser.add_argument("-target", metavar="DIR", dest="target", required=True, help="specify directory where all the outputs will be placed")
	preqc_parser.add_argument("-nproc", metavar="INT", type=int, dest="nproc", help="number of cores to request")
	preqc_parser.add_argument("-nokeep", dest="keep", action="store_false", help="specify whether or not to keep all intermediate files")
	preqc_parser.add_argument("config", metavar="CONFIG", help="specify configuration file where all data information is stored")
	blast_group = preqc_parser.add_argument_group("BLAST", "Arguments for runing BLAST")
	blast_group.add_argument("-e", metavar="NUM", dest="min_eval", default=10, help="specify miminum E Value for BLAST")
	blast_group.add_argument("-nthread", metavar="INT", type=int, dest="nthread", default=1, help="number of threads to use for blast.")
	blast_group.add_argument("-db", metavar="FILE", dest="db", help="specify a database for blast")
	cutadapt_group = preqc_parser.add_argument_group("Cutadapt", "Arguments for running Cutadapt. Currently only two common arguments are supported")
	cutadapt_group.add_argument("-O", metavar="INT", dest="min_overlap", type=int, default=3, help="specify the minimum length of overlap for cutadapt [3]")
	cutadapt_group.add_argument("-m", metavar="INT", dest="min_readlen", type=int, default=0, help="specify the minimum length of reads to keep after cutadapt [0]")
	preqc_parser.set_defaults(func=rpgc_preqc)

	mapping_parser = sub_parsers.add_parser("map", help="mapping reads from multiple libraries against given reference")
	post_mapping_parser = sub_parsers.add_parser("prepbam", help="preparing BAM files for variants calling")

	return main_parser.parse_args()

def rpgc_preqc(args):
	from preqc import Preqc
	Preqc(args).start()

def main():
	args = handle_cmd()
	args.func(args)

if __name__ == "__main__":
	main()
