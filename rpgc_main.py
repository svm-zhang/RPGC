# this script is a first attempt to use argparse subparser

import argparse
import os
import sys

#sys.path.insert(0, os.path.join(os.getcwd(), os.pardir))
from preqc import rpgc_preqc

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
	preqc_parser.add_argument("-nthread", metavar="INT", type=int, dest="nthread", help="number of threads to use for blast.")
	preqc_parser.add_argument("-nokeep", dest="keep", action="store_false", help="specify whether or not to keep all intermediate files")
	preqc_parser.add_argument("-db", metavar="FILE", dest="db", required=True, help="specify a database for blast")
	preqc_parser.add_argument("configs", metavar="CONFIG", nargs='+', type=argparse.FileType('r'), help="specify configuration file where all data information is stored")
	preqc_parser.set_defaults(func=rpgc_preqc)

	mapping_parser = sub_parsers.add_parser("map", help="mapping reads from multiple libraries against given reference")
	post_mapping_parser = sub_parsers.add_parser("prepBAM", help="preparing BAM files for variants calling")

	return main_parser.parse_args()



def main():
	args = handle_cmd()
	args.func(args)
	if args.command == "preqc":
		args.func(args)
	elif args.command == "map":
		print "map"
	elif args.command == "s2b":
		print "sam2bam"

if __name__ == "__main__":
	main()
