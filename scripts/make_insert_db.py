#!/usr/bin/python

"""Creates a FASTA file for the insert.

This file contains the 3' and 5' ends of the insert.
"""

import argparse

from Bio import SeqIO


def Main():
	parser = argparse.ArgumentParser(description='Filter reads.')
	parser.add_argument("insert_filename",
      					help="Path to FASTA file containing insert sequence.")
	parser.add_argument("out_filename",
					    help="Output filename.")
	parser.add_argument("-l", "--end_length", default=20, type=int,
					    help="The number of 3' and 5' bases to store.")
	args = parser.parse_args()

	insert_seq = SeqIO.read(args.insert_filename, 'fasta')
	insert_5p = insert_seq[:args.end_length]
	insert_3p = insert_seq[-args.end_length:]

	insert_5p.name += '_5p'
	insert_5p.id += '_5p'
	insert_3p.name += '_3p'
	insert_3p.id += '_3p'

	SeqIO.write([insert_5p, insert_3p], args.out_filename, 'fasta')


if __name__ == '__main__':
	Main()

