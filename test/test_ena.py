#!/usr/bin/env python3

#=======================================================================
# Authors: Ben Woodcroft
#
# Unit tests.
#
# Copyright
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

import unittest
import os.path
import sys

import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
kingfisher = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','bin','kingfisher')

from in_tempdir import in_tempdir

class Tests(unittest.TestCase):
    maxDiff = None
    
    def test_fastq_gz(self):
        cmd_stub = '{} get -r SRR12118866 -m'.format(kingfisher)
        with in_tempdir():
            extern.run("{} {}".format(cmd_stub,'ena-ascp'))
            self.assertTrue(os.path.getsize('SRR12118866_1.fastq.gz')==4117481)
            self.assertTrue(os.path.getsize('SRR12118866_2.fastq.gz')==4945891)

    def test_fasta_gz_via_ena_ftp(self):
        with in_tempdir():
            extern.run('{} get -r SRR12118866 -m ena-ftp --output-format-possibilities fasta.gz'.format(
                kingfisher))
            self.assertTrue(os.path.getsize('SRR12118866_1.fasta.gz')==746749)
            self.assertTrue(os.path.getsize('SRR12118866_2.fasta.gz')==899862)

    def test_fastq_via_ena_ascp(self):
        with in_tempdir():
            extern.run('{} get -r SRR12118866 -m ena-ascp --output-format-possibilities fastq'.format(
                kingfisher))
            self.assertTrue(os.path.getsize('SRR12118866_1.fastq')==19930812)
            self.assertTrue(os.path.getsize('SRR12118866_2.fastq')==19930812)

    def test_fasta_via_ena_ascp(self):
        with in_tempdir():
            extern.run('{} get -r SRR12118866 -m ena-ascp --output-format-possibilities fasta.gz fasta'.format(
                kingfisher))
            self.assertTrue(os.path.getsize('SRR12118866_1.fasta')==10391100)
            self.assertTrue(os.path.getsize('SRR12118866_2.fasta')==10391100)

    def test_fasta_gz_via_ena_ascp(self):
        with in_tempdir():
            extern.run('{} get -r SRR12118866 -m ena-ascp --output-format-possibilities fasta.gz'.format(
                kingfisher))
            self.assertTrue(os.path.getsize('SRR12118866_1.fasta.gz')==746749)
            self.assertTrue(os.path.getsize('SRR12118866_2.fasta.gz')==899862)

if __name__ == "__main__":
    unittest.main()
