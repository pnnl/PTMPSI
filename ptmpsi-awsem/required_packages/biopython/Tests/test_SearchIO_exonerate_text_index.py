# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO exonerate-text indexing."""

import os
import unittest

from search_tests_common import CheckIndex


class ExonerateTextIndexCases(CheckIndex):
    fmt = "exonerate-text"

    def test_exn_22_m_est2genome(self):
        """Test exonerate-text indexing, single."""
        filename = os.path.join("Exonerate", "exn_22_m_est2genome.exn")
        self.check_index(filename, self.fmt)

    def test_exn_22_q_multiple(self):
        """Test exonerate-text indexing, single."""
        filename = os.path.join("Exonerate", "exn_22_q_multiple.exn")
        self.check_index(filename, self.fmt)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
