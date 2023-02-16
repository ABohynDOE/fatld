import random

import numpy as np
import oapackage as oa  # type: ignore
import pytest  # type: ignore

import fatld
from fatld import Design
from fatld.main import basic_factor_matrix, gen2num, num2gen, power2_decomposition, twlp


def test_version():
    assert fatld.__version__ == "0.1.0"


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


class TestNum2Gen:
    def test_num2gen(self):
        assert num2gen(7) == "abc"

    def test_num2gen_mvalue1(self):
        assert num2gen(7, m=1) == "A3c"

    def test_num2gen_mvalue2(self):
        assert num2gen(13, m=2) == "A1B3"

    def test_num2gen_typeerror(self):
        with pytest.raises(TypeError):
            num2gen("a")

    def test_num2gen_valueerror(self):
        with pytest.raises(ValueError):
            num2gen(7, m=4)


def test_gen2num():
    assert gen2num("acd") == 13


def test_gen2num_type_error():
    with pytest.raises(TypeError):
        gen2num(7)


def test_gen2num_value_error():
    with pytest.raises(ValueError):
        gen2num("Abc")
