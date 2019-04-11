# -*- coding: utf-8 -*-
# test_stats.py
"""
Test stats tests implementation

Created on Tue Jul 17 09:24:38 2018

@author: wenyi
"""

import unittest
#import json
import cebspy.stats as stats
#import cebspy.tests.examples as examples

class TestCEBSStats(unittest.TestCase):
    """
    For unit tests of cebspy.stats package
    The default location of the file is at 'cebspy/tests' folder
    To run the test from commanline at parent folder:
    
    Example
    -------
    C:.../cebspy>python -m unittest -v tests/test_stats.py
    test_examplestats (tests.test_stats.TestCEBSStats) ... ok
    test_ttest (tests.test_stats.TestCEBSStats) ... ok
    
    ----------------------------------------------------------------------
    Ran 2 tests in 0.001s
    
    OK
    # Altetnative command using 'discover' sub-command
    C:.../cebspy>python -m unittest discover tests/ -v
    """
    def test_examplestats(self):
        sample = [362.8, 337.9, 341.4, 338.8, 285.1, 336.8, 343.0, 340.0,
                  339.4, 324.2]
        results = stats.samplestats.sample_stats(sample)
        self.assertEqual(results.get('size'), len(sample))
        self.assertEqual(results.get('mean'), sum(sample) / len(sample))
    
    def test_ttest(self):
        sample1 = [362.8, 337.9, 341.4, 338.8, 285.1, 336.8, 343.0, 340.0]
        sample2 = [422.2, 454.3, 429.8, 376.8, 408.7, 488.9, 405.9, 444.2, 369.1, 441.7]
        results = stats.ttest.t_test(sample1, sample2) # assume var_equal = False (default)
        self.assertTrue(results['has_output'])
        self.assertFalse(results['has_errors'])
        self.assertIn('is_finished', results)
        self.assertEqual(list(results.keys()), ['method', 'has_output', 'has_errors', 'is_finished', 'output'])
        
    def test_williams(self):
        sample1 = [362.8, 337.9, 341.4, 338.8, 285.1, 336.8, 343.0, 340.0]
        sample2 = [422.2, 454.3, 429.8, 376.8, 408.7, 488.9, 405.9, 444.2, 369.1, 441.7]
        results = stats.ttest.t_test(sample1, sample2) # assume var_equal = False (default)
        self.assertTrue(results['has_output'])
        self.assertFalse(results['has_errors'])
        self.assertIn('is_finished', results)
        self.assertEqual(list(results.keys()), ['method', 'has_output', 'has_errors', 'is_finished', 'output'])
     
if __name__ == '__main__':
    unittest.main()
