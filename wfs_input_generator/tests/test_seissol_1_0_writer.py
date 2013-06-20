#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test suite for the Seissol 1.0 writer.

:copyright:
    Marek Simon (marek.simon@geophysik.uni-muenchen.de), 2013
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
from wfs_input_generator import InputFileGenerator

import glob
import inspect
import numpy as np
from obspy.core import UTCDateTime
from obspy.core.event import readEvents
import os
import unittest


class Seissol_1_0_WriterTestCase(unittest.TestCase):
    
    """
    Test case for the Seissol 1.0 writer.
    """
    
    def setUp(self):
        # Most generic way to get the actual data directory.
        self.data_dir = os.path.join(os.path.dirname(os.path.abspath(
            inspect.getfile(inspect.currentframe()))), "data")

    def test_real_world_example(self):
        
        """
        Test that compares the created input files to those from a real world
        example.

        The only artificial thing is the source-time function but that is
        trivial to verify.

        This is a fairly comprehensive tests but should be used in comparision
        with other unit tests.
        """
        
        gen = InputFileGenerator()

        seissol_example_path = os.path.join(self.data_dir, "seissol_example")
        gen.add_stations([os.path.join(self.data_dir, "dataless.seed.BW_FURT"),\
			  os.path.join(self.data_dir, "dataless.seed.BW_RJOB")])
        gen.add_events(readEvents(os.path.join(self.data_dir, "event2.xml")))
        
        # Configure it.
        gen.config.mesh = 'most_simple_tet'
        gen.config.model = 'PREM'
        gen.config.working_directory = seissol_example_path
        gen.config.max_time = 1000.0
        gen.config.number_of_processors = 16
        # Write the input files to a dictionary.
        
        input_files = gen.write(format = 'seissol_1_0', output_dir = seissol_example_path)

        # The rest is only for asserting the produced files.
        for filename in glob.glob(os.path.join(seissol_example_path, "*_example")):
            with open(filename, "rt") as open_file:
                real_file = open_file.read()
            filename = os.path.basename(filename[:-8])

            if filename not in input_files:
                msg = "File '%s' has not been generated" % filename
                raise AssertionError(msg)

            lines = real_file.splitlines()
            new_lines = input_files[filename].splitlines()

            if len(lines) != len(new_lines):
                msg = ("File '%s' does not have the same number of lines "
                    "for the real (%i lines) and generated (%i lines) "
                    "input file") % (filename, len(lines), len(new_lines))
                raise AssertionError(msg)

            for line, new_line in zip(lines, new_lines):
                if line != new_line:
                    msg = "Line differs in file '%s'.\n" % filename
                    msg += "Expected: \"%s\"\n" % line
                    msg += "Got:      \"%s\"\n" % new_line
                    raise AssertionError(msg)



def suite():
    return unittest.makeSuite(Seissol_1_0_WriterTestCase, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
