import unittest
import os
import sys

sys.path.append(os.getcwd()[:os.getcwd().index('DevelopingRelease')])
from DevelopingRelease.WinnowDevel import performetrics


class PerformetricsTest(unittest.TestCase):

    def test_rmse(self):
        """
        Tests method rmse (Returns root mean squared error given the known truth and effects)
        Checks for equality of returned value and expected root mean squared error value
        """
        beta_col = [1, 2, 3, 4, 5, 6]
        beta_tf = [1, 0, 1, 1, 0, 0]
        expected = 13.0
        self.assertEqual(performetrics.rmse(beta_col, beta_tf), expected)

    def test_mae(self):
        """
        Tests method mae (Returns the mean absolute error of a data set)
        Checks for equality of returned value and expected mean absolute error value
        """
        beta_col = [1, 2, 3, 4, 5, 6]
        beta_tf = [1, 0, 1, 1, 0, 0]
        expected = 3.0
        self.assertEqual(performetrics.mae(beta_col, beta_tf), expected)

    def test_mattcorr(self):
        """
        Tests method mattcorr (Returns Matthews correlation coefficient value)
        Checks for equality of returned value and expected correlation coefficient value
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 0.1690308509457033
        self.assertEqual(performetrics.mattcorr(snp_tf, threshold, score), expected)

    def test_auc(self):
        """
        Tests method auc (Returns the area under curve for binary classification)
        Checks for equality of returned value and expected area under curve value
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 0.56944444444444442
        self.assertEqual(performetrics.auc(snp_tf, score), expected)

    def test_tp(self):
        """
        Tests method tp (Returns total number of SNPs correctly identified as significant)
        Checks for equality of returned number and expected number of true positives
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 4
        self.assertEqual(performetrics.tp(snp_tf, threshold, score), expected)

    def test_fp(self):
        """
        Tests method fp (Returns the number of SNPs incorrectly identified as significant)
        Checks for equality of returned number and expected number of false positives
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 3
        self.assertEqual(performetrics.fp(snp_tf, threshold, score), expected)

    def test_tn(self):
        """
        Tests method tn (Returns the number of SNPs correctly identified as not significant)
        Checks for equality of returned number and expected number of true negatives
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 3
        self.assertEqual(performetrics.tn(snp_tf, threshold, score), expected)

    def test_fn(self):
        """
        Tests method fn (Returns the number of SNPs incorrectly identified as not significant)
        Checks for equality of returned number and expected number of false negatives
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 2
        self.assertEqual(performetrics.fn(snp_tf, threshold, score), expected)

    def test_tpr(self):
        """
        Tests method tpr (Returns the proportion of true positives identified from the entire data set)
        Checks for equality of returned value and expected proportion value
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 0.6666666666666666
        self.assertEqual(performetrics.tpr(snp_tf, threshold, score), expected)

    def test_fpr(self):
        """
        Tests method fpr (Returns proportion of false positives identified from the data set)
        Checks for equality of returned value with expected proportion value
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 0.5
        self.assertEqual(performetrics.fpr(snp_tf, threshold, score), expected)

    def test_prevalence(self):
        """
        Tests method prevalence (Returns the prevalence value - total number of positives divided by sample size)
        Checks for equality of returned value and expected prevalence value
        """
        snpTF = [True,False,True,True,True,False,False,True,False,False,True,False]
        threshold = 0.05
        score = [0.003,0.65,0.004,0.006,0.078,0.003,0.0001,0.513,0.421,0.0081,0.043,0.98]
        expected = 0.5
        self.assertEqual(performetrics.prevalence(snpTF, threshold, score), expected)

    def test_error(self):
        """
        Tests method error (Returns error value of the analysis - total false identifications divided by sample size)
        Checks for equality of returned value and expected error value
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 0.4166666666666667
        self.assertEqual(performetrics.error(snp_tf, threshold, score), expected)
    
    def test_accuracy(self):
        """
        Tests method accuracy (Returns accuracy value of analysis - total true identifications divided by sample size)
        Checks for equality of returned value and expected accuracy value
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 0.58333333
        self.assertEqual(format_float(performetrics.accuracy(snp_tf, threshold, score)), expected)

    def test_sens(self):
        """
        Tests method sens (Returns sensitivity value of analysis - number of correctly identified positives
        divided by number of known-truth positives)
        Checks for equality of returned value and expected sensitivity value
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 0.6666666666666666
        self.assertEqual(performetrics.sens(snp_tf, threshold, score), expected)

    def test_spec(self):
        """
        Tests method spec (Returns specificity value - number of correctly identified negatives
        divided by number of known-truth negatives)
        Checks for equality of returned value and expected specificity value
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 0.5
        self.assertEqual(performetrics.spec(snp_tf, threshold, score), expected)

    def test_precision(self):
        """
        Tests method precision (Returns precision - number of correctly identified positives divided by
        total identified positives)
        Checks for equality of returned value and expected precision value
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 0.5714285714285714
        self.assertEqual(performetrics.precision(snp_tf, threshold, score), expected)
    
    def test_fdr(self):
        """
        Tests method fdr (Returns false discovery rate - number of false positives divided by total number of positives)
        Checks for equality of returned value and expected fdr value
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 0.42857143
        self.assertEqual(format_float(performetrics.fdr(snp_tf, threshold, score)), expected)

    def test_youden(self):
        """
        Tests method youden (Returns the Youden statistic for the data: Sensitivity + Specificity - 1
        Checks for equality of returned value and expected Youden statistic value
        """
        snp_tf = [True, False, True, True, True, False, False, True, False, False, True, False]
        threshold = 0.05
        score = [0.003, 0.65, 0.004, 0.006, 0.078, 0.003, 0.0001, 0.513, 0.421, 0.0081, 0.043, 0.98]
        expected = 0.16666666666666652
        self.assertEqual(performetrics.youden(snp_tf, threshold, score), expected)

    def test_avgcovarweight(self):
        """
        Tests method avgcovarweight (Returns average covariate weight based on the model produced from the data set)
        Checks for equality of returned value and expected average covariate weight value
        """
        covar = [0.08410, 0.00161, 0.25200, 0.00161, 0.24000, 0.25000, 0.00161, 0.25000, 0.00161, 0.25200, 0.25000,
                 0.00161, -0.02440, 0.05360]
        expected = 0.11538214285714285
        self.assertAlmostEquals(performetrics.avgcovarweight(covar), expected)


def format_float(x):
    """
    Truncates floats to 8 decimal places

    :param float_list:
    :return: a list of float truncated to 8 decimal places
    """
    return float('%.8f' % x)


def get_test_suite():
    """
    Returns a test suite with all tests
    """
    return unittest.TestLoader().loadTestsFromTestCase(PerformetricsTest)

if __name__ == "__main__":
    unittest.main()