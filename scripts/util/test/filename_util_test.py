#!/usr/bin/python

import unittest

from scripts.util import filename_util


class FilenameUtilTests(unittest.TestCase):
	"""Unit tests for filename_util.py."""

	def testMakeFnameBasic(self):
		out_ext = 'fq'
		test_cases = (('foo.fa.gz', 'foo.fq'),
					  ('foo.fa', 'foo.fq'),
					  ('asdkmals_askdm12.bz2', 'asdkmals_askdm12.fq'),
					  ('/home/ban/maggots/s123_sd.mmq', 's123_sd.fq'))
		for input_fname, expected_out in test_cases:
			out_fname = filename_util.MakeFname(input_fname, out_ext)
			self.assertEquals(expected_out, out_fname)

	def testMakeFnameDestDir(self):
		out_ext = 'fq'
		dest_dir = '_fake_dir/tmp'
		test_cases = (('foo.fa', '%s/foo.fq' % dest_dir),
					  ('asdkmals_askdm12.bz2',
					   '%s/asdkmals_askdm12.fq' % dest_dir),
					  ('/home/ban/maggots/s123_sd.mmq',
					   '%s/s123_sd.fq' % dest_dir))
		for input_fname, expected_out in test_cases:
			out_fname = filename_util.MakeFname(
				input_fname, out_ext, dest_dir=dest_dir)
			self.assertEquals(expected_out, out_fname)

	def testMakeFnamePostfix(self):
		out_ext = 'fq'
		postfix = 'postfix_asd'
		test_cases = (('foo.fa', 'foo_%s.fq' % postfix),
					  ('asdkmals_askdm12.bz2',
					   'asdkmals_askdm12_%s.fq' % postfix),
					  ('/home/ban/maggots/s123_sd.mmq',
					   's123_sd_%s.fq' % postfix))
		for input_fname, expected_out in test_cases:
			out_fname = filename_util.MakeFname(
				input_fname, out_ext, postfix=postfix)
			self.assertEquals(expected_out, out_fname)

	def testMakeFnamePostfixDestDir(self):
		out_ext = 'fq'
		dest_dir = '_fake_dir/tmp'
		postfix = 'postfix_asd'
		test_cases = (('foo.fa', '%s/foo_%s.fq' % (dest_dir, postfix)),
					  ('asdkmals_askdm12.bz2',
					   '%s/asdkmals_askdm12_%s.fq' % (dest_dir, postfix)),
					  ('/home/ban/maggots/s123_sd.mmq',
					   '%s/s123_sd_%s.fq' % (dest_dir, postfix)))
		for input_fname, expected_out in test_cases:
			out_fname = filename_util.MakeFname(
				input_fname, out_ext, dest_dir=dest_dir,
				postfix=postfix)
			self.assertEquals(expected_out, out_fname)


if __name__ == '__main__':
	unittest.main()
