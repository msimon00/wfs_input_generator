#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DESCRIPTION

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2013
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
"""
from wfs_input_generator import InputFileGenerator

import inspect
import os
import unittest


class InputFileGeneratorTestCase(unittest.TestCase):
    """
    Test case for the general InputFileGenerator.
    """
    def setUp(self):
        # Most generic way to get the actual data directory.
        self.data_dir = os.path.join(os.path.dirname(os.path.abspath(
            inspect.getfile(inspect.currentframe()))), "data")

    def test_readingSEEDFiles(self):
        """
        Tests the reading of SEED files.
        """
        seed_file_1 = os.path.join(self.data_dir, "dataless.seed.BW_FURT")
        seed_file_2 = os.path.join(self.data_dir, "dataless.seed.BW_RJOB")

        gen = InputFileGenerator()
        gen.add_stations([seed_file_1, seed_file_2])

        self.assertEqual([
            {"id": "BW.FURT",
             "latitude": 48.162899,
             "longitude": 11.2752,
             "elevation_in_m": 565.0,
             "local_depth_in_m": 0.0},
            {"id": "BW.RJOB",
             "latitude": 47.737167,
             "longitude": 12.795714,
             "elevation_in_m": 860.0,
             "local_depth_in_m": 0.0}], gen._stations)


def suite():
    return unittest.makeSuite(InputFileGeneratorTestCase, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
