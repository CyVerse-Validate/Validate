import unittest
import os
import sys
sys.path.append(os.getcwd()[:os.getcwd().index('DevelopingRelease')])
from DevelopingRelease.WinnowDevel import HMeasure


class HMeasureTest(unittest.TestCase):
    def test_h_measure(self):
        true_class = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
        scores = [1.8, 0.9, 1.6, 0.5, 1, 0.1, 0.2, 2.6, -0.4, -0.1]
        self.assertAlmostEqual(HMeasure.h_measure(true_class, scores), 0.696296296296)

    def test_auc_solver(self):
        out_scores = self.test_get_score_distributions()
        self.assertAlmostEquals(0.2, HMeasure.auc_solver(out_scores['s0'], out_scores['f1'], out_scores['s1']))

    def test_get_score_distributions(self):
        true_class = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
        scores = [1.8, 0.9, 1.6, 0.5, 1, 0.1, 0.2, 2.6, -0.4, -0.1]
        n1 = float(sum(true_class))
        n0 = float(len(scores) - n1)
        expected = {'f0': [0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0, 1.0, 1.0],
                    'f1': [0.0, 0.2, 0.4, 0.6000000000000001, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 1.0, 1.0],
                    's1': [0.0, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0],
                    's0': [0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                    's_len': 12}
        self.assertEquals(expected, HMeasure.get_score_distributions(true_class, scores, n1, n0))
        return expected

    def test_cum_sum(self):
        for x, y in zip(HMeasure.cum_sum([1, 2, 3, 4, 5]), [1.0, 3.0, 6.0, 10.0, 15.0]):
            self.assertEquals(x, y)

    def test_misclass_counts(self):
        true_class = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
        scores = [1.8, 0.9, 1.6, 0.5, 1, 0.1, 0.2, 2.6, -0.4, -0.1]
        expected = {'fp': 4.0, 'tp': 1.0, 'fn': 4.0, 'er': 0.8, 'youden': -0.6, 'f': 0.2, 'recall': 0.2, 'prec': 0.2,
                    'tn': 1.0, 'tpr': 0.2, 'fpr': 0.8, 'sens': 0.2, 'spec': 0.2}
        self.assertEquals(expected, HMeasure.misclass_counts(scores, true_class, 0.5))


def get_test_suite():
    """
    Returns a test suite with all tests for the HMeasure functions

    :return:
    """
    return unittest.TestLoader().loadTestsFromTestCase(HMeasureTest)

if __name__ == "__main__":
    unittest.main()
