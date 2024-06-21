"""
Regression tests for 2D acoustics with adjoint flagging.

Note that only the forward problem is being tested, but first
the code in adjoint/ must be run.
Regression data for that code is in adjoint/regression_data but is not
being tested against.
"""

import sys
import os
import subprocess
import unittest

thisfile = os.path.realpath(__file__)
testdir = os.path.split(thisfile)[0]

import clawpack.amrclaw.test as test


class Acoustics2DAdjointTest(test.AMRClawRegressionTest):
    r"""Basic test for a 2D acoustics adjoint-flagging forward problem test case"""


    def runTest(self, save=False):
        
        # Build and run adjoint code
        adjoint_path = os.path.join(self.test_path, "adjoint")
        self.stdout.write("Running adjoint sub\n")
        self.stdout.write(f"  adjoint path: {adjoint_path}\n")
        self.stdout.write(f"  temp path: {str(self.temp_path)}\n")
        exe_cmd = ['make', '-C', adjoint_path, '.exe']
        subprocess.check_call(['make', '-C', adjoint_path, '.exe'],
                                                    stdout=self.stdout,
                                                    stderr=self.stderr)
        exe_cmd = ['make', '-C', adjoint_path, '.output']
        subprocess.check_call(['make', '-C', adjoint_path, '.output'],
                                                    stdout=self.stdout,
                                                    stderr=self.stderr)

        # Write out data files
        self.load_rundata()
        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=1)
        self.check_gauges(save=save, gauge_id=2)

        self.success = True



if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Acoustics2DAdjointTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
