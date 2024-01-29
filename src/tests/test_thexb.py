import unittest
import statistics
from pathlib import Path
# Dependency Imports
from pyfaidx import Fasta
# THExBuilder Imports
from thexb.STAGE_minifastas import get_seq
from thexb.STAGE_pairwise_estimator import coverage_and_median
from thexb.STAGE_pairwise_estimator import p_distance
from thexb.STAGE_pairwise_filter import Mean, Median, StandardDeviation
from thexb.STAGE_iqtree import remove_heterotachy_info as rhi
from thexb.STAGE_iqtree_external import remove_heterotachy_info as rhie

class TestTHExBuilder(unittest.TestCase):
    ########## Fasta Windower ##########
    def test_minifasta_get_seq(self):
        # -- Test inputs --
        test_seq_dict = {"Sample1": list(Fasta("src/tests/data/seq.fasta")["Sample1"][:].seq)}
        test_start_pos1 = 0
        test_end_pos1 = 20
        test_start_pos2 = 35
        test_end_pos2 = 50
        # -- Test Results --
        r1 = ["AGTTGTGTGTTGTGTCTGTT"]
        r2 = ["CTGTCTGTTCTGTCA"]
        # -- Module Results --
        er1 = get_seq(test_seq_dict, "Sample1", test_start_pos1, test_end_pos1)
        er2 = get_seq(test_seq_dict, "Sample1", test_start_pos2, test_end_pos2)
        # -- Assert results are valid --
        self.assertEqual(r1, er1)
        self.assertEqual(r2, er2)
        self.assertEqual(len(er1[0]), 20)
        self.assertEqual(len(er2[0]), 15)
        return
    
    
    ########## Pairwise Estimator ##########
    def test_pwe_coverage_and_median(self):
        # -- Test inputs --
        seq_files = [f for f in Path("src/tests/data/pwe/").iterdir() if ".fai" not in f.name]
        # -- Test Results --
        test_avg_coverage = {"Sample1": 1.0, 'Sample2': 0.5125}
        test_cov_median = {"Sample1": 1.0, 'Sample2': 0.5375}
        # -- Module Results --
        avg_coverage, cov_median = coverage_and_median(seq_files)
        # -- Assert results are valid --
        self.assertDictEqual(avg_coverage, test_avg_coverage)
        self.assertDictEqual(cov_median, test_cov_median)
        return


    def test_pwe_p_distance(self):
        # -- Test inputs --
        seq_files = [f for f in Path("src/tests/data/pwe/").iterdir() if ".fai" not in f.name]
        # -- Test Results --
        test_avg_p_distance = {'Sample2': round(statistics.mean([(0/80), (80/80), (37/80)]), 4)}
        test_pdist_median = {'Sample2': round(statistics.median([(0/80), (80/80), (37/80)]), 4)}
        # -- Module Results --
        avg_p_distance, pdist_median = p_distance(seq_files, 'Sample1')
        # -- Assert results are valid --
        self.assertDictEqual(avg_p_distance, test_avg_p_distance)
        self.assertDictEqual(pdist_median, test_pdist_median)

    
    ########## Pairwise Filter ##########
    def test_pwf_Mean(self):
        # -- Test inputs --
        test_list1 = [5, 10, 15, 20]
        test_list2 = [1, 3, 15, 25]
        # -- Test Results --
        test_mean1 = 12.5
        test_mean2 = 11
        # -- Module Results --
        mean1 = Mean(test_list1)
        mean2 = Mean(test_list2)
        # -- Assert results are valid --
        self.assertEqual(mean1, test_mean1)
        self.assertEqual(mean2, test_mean2)


    def test_pwf_Median(self):
        # -- Test inputs --
        test_list1 = [5, 15, 15, 20]
        test_list2 = [1, 3, 15, 25]
        # -- Test Results --
        test_median1 = 15
        test_median2 = 9
        # -- Module Results --
        median1 = Median(test_list1)
        median2 = Median(test_list2)
        # -- Assert results are valid --
        self.assertEqual(median1, test_median1)
        self.assertEqual(median2, test_median2)
    
    
    def test_pwf_StandardDeviation(self):
        # -- Test inputs --
        test_list1 = [5, 15, 15, 20]
        test_list2 = [1, 3, 15, 25]
        # -- Test Results --
        test_stdev1 = 6.292
        test_stdev2 = 11.20
        # -- Module Results --
        stdev1 = round(StandardDeviation(test_list1), 3)
        stdev2 = round(StandardDeviation(test_list2), 2)
        # -- Assert results are valid --
        self.assertAlmostEqual(stdev1, test_stdev1)
        self.assertAlmostEqual(stdev2, test_stdev2)

    ########## Pairwise Filter ##########
    def test_iqtree_remove_heterotachy_info(self):
        # -- Test inputs --
        test_tree1 = "(A[0.0000346/0.0000106/0.0837614/0.0739605]:1,(B[0.0000346/0.0000106/0.0837614/0.0739605]:1,C[0.0000346/0.0000106/0.0837614/0.0739605]:1));"
        test_tree2 = "(A,(B,C));"
        test_tree3 = "(A:1,(B:1,C:1):1):1;"
        # -- Test Results --
        test_tree1_clean = "(A:1,(B:1,C:1));"
        test_tree2_clean = "(A,(B,C));"
        test_tree3_clean = "(A:1,(B:1,C:1):1):1;"
        # -- Module Results --
        tree1 = rhi(test_tree1)
        tree2 = rhi(test_tree2)
        tree3 = rhi(test_tree3)
        # -- Assert results are valid --
        self.assertEqual(test_tree1_clean, tree1)
        self.assertEqual(test_tree2_clean, tree2)
        self.assertEqual(test_tree3_clean, tree3)
    

    def test_iqtree_external_remove_heterotachy_info(self):
        # -- Test inputs --
        test_tree1 = "(A[0.0000346/0.0000106/0.0837614/0.0739605]:1,(B[0.0000346/0.0000106/0.0837614/0.0739605]:1,C[0.0000346/0.0000106/0.0837614/0.0739605]:1));"
        test_tree2 = "(A,(B,C));"
        test_tree3 = "(A:1,(B:1,C:1):1):1;"
        # -- Test Results --
        test_tree1_clean = "(A:1,(B:1,C:1));"
        test_tree2_clean = "(A,(B,C));"
        test_tree3_clean = "(A:1,(B:1,C:1):1):1;"
        # -- Module Results --
        tree1 = rhie(test_tree1)
        tree2 = rhie(test_tree2)
        tree3 = rhie(test_tree3)
        # -- Assert results are valid --
        self.assertEqual(test_tree1_clean, tree1)
        self.assertEqual(test_tree2_clean, tree2)
        self.assertEqual(test_tree3_clean, tree3)

if __name__ == '__main__':
    unittest.main()
