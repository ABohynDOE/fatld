import random

import numpy as np
import oapackage as oa  # type: ignore
import pytest  # type: ignore

from fatld.design import Design
from fatld.main import (
    basic_factor_matrix,
    power2_decomposition,
    twlp,
    alpha_wlp,
    beta_star_wlp,
)


class TestBasicFactorMatrix:
    k = random.randint(2, 6)
    mat = basic_factor_matrix(k)
    coded_mat = basic_factor_matrix(k, False)

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
    assert power2_decomposition(11) == [1, 1, 0, 1]


def test_power2_decomposition_warnings():
    with pytest.warns(UserWarning):
        power2_decomposition(11, length=2)


class TestTWLP:
    D = Design(32, 1, [29, 26, 22])
    ar = oa.array_link(D.array)

    def test_twlp(self):
        t = twlp(self.ar, max_length=5)
        assert t == [[0, 0], [1, 4], [0, 2]]

    def test_twlp_type1(self):
        t = twlp(self.ar, max_length=5, type_0=False)
        assert t == [[0, 0], [4, 1], [2, 0]]

    def test_twlp_length_None(self):
        t = twlp(self.ar, max_length=None)
        assert len(t) == self.D.m + self.D.n - 2

    def test_twlp_length_wrong(self):
        with pytest.warns(UserWarning):
            twlp(self.ar, max_length=2)


class TestSubDesign:
    D = Design(
        32,
        2,
        [
            5,
            6,
            7,
            9,
            10,
            11,
            13,
            14,
            15,
            17,
            18,
            19,
            20,
            21,
            22,
            23,
            24,
            25,
            26,
        ],
    )
    ar = oa.array_link(D.array).deleteColumn(0)
    sub_D = Design(
        32,
        1,
        [12, 20, 24, 5, 9, 6, 18, 11, 19, 21, 25, 14, 26, 15, 23, 29, 30],
    )

    def test_wlp(self):
        t = twlp(self.ar, max_length=5)
        assert t == self.sub_D.twlp(max_length=5)

    def test_alpha_wlp(self):
        t = alpha_wlp(self.ar)
        assert t == [round(i, 2) for i in self.sub_D.alpha_wlp()]

    def test_beta_wlp(self):
        b_vec, _ = beta_star_wlp(self.ar, max_length=7)
        assert b_vec == self.sub_D.beta_star_wlp(max_length=7)[0]
