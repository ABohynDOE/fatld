import random

import numpy as np

import fatld
import fatld.main


def test_version():
    assert fatld.__version__ == "0.1.0"


class TestBasicFactorMatrix:
    k = random.randint(2, 6)
    mat = fatld.main.basic_factor_matrix(k)
    coded_mat = fatld.main.basic_factor_matrix(k, False)

    def test_sum(self):
        row_sum = np.sum(self.mat, axis=0)
        assert all(row_sum == (2**self.k // 2))

    def test_diag(self):
        mult_mat = np.matmul(self.coded_mat.T, self.coded_mat) // 2**self.k
        diag_mat = np.eye(self.k)
        equal_mat = mult_mat == diag_mat
        assert equal_mat.all()

    def test_unique(self):
        vals = np.unique(self.mat)
        assert all(vals == (0, 1))

    def test_unique_coded(self):
        vals = np.unique(self.coded_mat)
        assert all(vals == (-1, 1))


def test_power2_decomposition():
    assert fatld.main.power2_decomposition(11) == [1, 1, 0, 1]


# TODO: test for the design class
