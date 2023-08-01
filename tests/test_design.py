"""
Test for the design class

Created on: 09/02/2023
Author: Alexandre Bohyn
"""
import numpy as np
import pandas as pd  # type: ignore
import pytest  # type: ignore

from collections import Counter

import fatld
import random


class TestDesignErrors:
    def test_wrong_init(self):
        assert True


class TestDesignMethods:
    # Design 6 from Table 3 of Wu (1993)
    design = fatld.Design(runsize=32, m=1, cols=[29, 26, 22])
    relation = fatld.relation.Relation(["bcef", "bdeg", "acdeh"], m=1)

    def test_twlp(self):
        twlp = self.design.twlp(max_length=5)
        assert twlp == [[0, 0], [1, 4], [0, 2]]

    def test_twlp_type1(self):
        twlp = self.design.twlp(type_0=False, max_length=5)
        assert twlp == [[0, 0], [4, 1], [2, 0]]

    def test_wlp(self):
        wlp = self.design.wlp(max_length=5)
        assert wlp == [0, 5, 2]

    def test_flatten(self):
        flat_array = self.design.flatten()
        assert all(np.unique(flat_array) == [0, 1])

    def test_flatten_zero_coding(self):
        flat_array = self.design.flatten(zero_coding=False)
        assert all(np.unique(flat_array) == [-1, 1])

    def test_relation(self):
        rel = self.design.defining_relation()
        assert rel.words == self.relation.words


class TestCols:
    design = fatld.Design(runsize=32, m=1, cols=[29, 27])

    def test_add_cols(self):
        design_more = self.design.add_factor(22)
        assert design_more.cols == [4, 8, 16, 22, 27, 29]

    def test_add_wrong_col(self):
        with pytest.raises(ValueError):
            self.design.add_factor(72)

    def test_add_already_used_cols(self):
        with pytest.raises(ValueError):
            self.design.add_factor(4)

    def test_remove_cols(self):
        design_less = self.design.remove_factor(29)
        assert design_less.cols == [4, 8, 16, 27]

    def test_remove_wrong_cols(self):
        with pytest.raises(ValueError):
            self.design.remove_factor(17)


class TestClarity:
    design = fatld.Design(32, 1, [27, 29])
    df = pd.DataFrame(
        {
            "4-4 clear": [0, 15, 10, 25],
            "4-2 clear": [0, 15, 4, 19],
            "2-2 clear": [0, 9, 10, 19],
            "Totally clear": [0, 9, 4, 13],
        },
        dtype=int,
    )
    df = df.set_index(pd.Index(["4-4", "4-2", "2-2", "Any type"]))

    design2 = fatld.Design(runsize=32, m=2, cols=[17, 21, 25, 29, 31])
    df2 = pd.DataFrame(
        {
            "4-4 clear": [9, 36, 8, 53],
            "4-2 clear": [9, 4, 15, 28],
            "2-2 clear": [2, 36, 9, 47],
            "Totally clear": [2, 4, 2, 8],
        },
        dtype=int,
    )
    df2 = df2.set_index(pd.Index(["4-4", "4-2", "2-2", "Any type"]))

    def test_clarity(self):
        tfi = self.design.clarity()
        assert tfi.equals(self.df)

    def test_clarity_m2(self):
        tfi = self.design2.clarity()
        assert tfi.equals(self.df2)

    def test_clear(self):
        assertion_list = []
        for row in ["4-4", "4-2", "2-2"]:
            if row == "Any type":
                row_index = "all"
            else:
                row_index = row
            for col in [
                "4-4 clear",
                "4-2 clear",
                "2-2 clear",
                "Totally clear",
            ]:
                if col == "Totally clear":
                    col_index = "all"
                else:
                    col_index = col.replace(" clear", "")
                assertion_list.append(
                    self.design.clear(row_index, col_index)
                    == self.df.loc[row, col]
                )
        assert all(assertion_list)

    def test_clear_error(self):
        with pytest.raises(ValueError):
            self.design.clear("4-2", "4-3")


def test_from_array():
    D = fatld.Design(32, 1, [17, 21, 29, 31])
    mat = D.array
    newD = fatld.design.from_array(mat)
    assert D.cols == newD.cols


def test_beta_wlp():
    D = fatld.Design(32, 2, [17, 18, 20, 24, 19, 28, 5, 9, 6, 10, 7, 11, 13])
    vec, perm = D.beta_wlp(max_length=5)
    assert vec == [8.0, 37.0, 107.0]


def test_beta_aberration_low_n():
    D = fatld.Design(
        32, 2, [31]
    )  # Here n = 2 which is lower than the max_length set
    qvec, perm = D.beta_wlp(max_length=6)
    assert qvec == [0.0, 0.0, 0.0, 1.0]


def test_number_interactions():
    assert fatld.design.nbr_interactions(7, 5) == 119


def test_tfi_model_matrix():
    tfi_model_matrix = fatld.design.build_tfi_model_matrix(n=7, max_length=5)
    int_length = np.sum(tfi_model_matrix, axis=0).tolist()
    c = Counter(int_length)
    d = dict(c)
    assert d == {1: 7, 2: 21, 3: 35, 4: 35, 5: 21}


def test_beta_star_m3():
    """Test the computation of B* vector for m=3"""
    D = fatld.Design(64, 3, [22, 27, 41, 46])
    b, p = D.beta_star_wlp(max_length=6)
    assert b == [[0.0, 0], [0.1, 0], [3.2, 0], [2.9, 0]]


def test_beta_star_sum():
    """
    Test that for any design, the sum of A_i^0  words is equal to the sum of
    the beta* values
    """
    available_cols = [
        i for i in range(1, 32) if i not in [1, 2, 3, 4, 8, 12, 16]
    ]
    cols = random.sample(available_cols, 5)
    D = fatld.Design(32, 2, cols)
    ai_values = [sum(i[1:]) for i in D.twlp(type_0=True)]
    beta_values = [b[0] for b in D.beta_star_wlp()[0]]
    assert np.round(sum(ai_values)) == np.round(sum(beta_values))


def test_alpha_sum():
    """
    Test that for any design, the sum of A_3  words is equal to the sum of
    the omega_2 and omega_4 values, and that the sum of the A-4 words is equal
    to the sum of the remaining omega values.
    """
    available_cols = [
        i for i in range(1, 32) if i not in [1, 2, 3, 4, 8, 12, 16]
    ]
    cols = random.sample(available_cols, 5)
    D = fatld.Design(32, 2, cols)
    a3 = D.wlp()[0]
    a4 = D.wlp()[1]
    omega_3 = sum(D.alpha_wlp()[0:2])
    omega_4 = sum(D.alpha_wlp()[2:])
    assert a3 == np.round(omega_3) and a4 == np.round(omega_4)


class TestSmallDesign:
    design = fatld.Design(16, 2, [10])

    def test_wlp(self):
        assert self.design.wlp() == [1]

    def test_alpha_wlp(self):
        assert self.design.alpha_wlp() == [0.667, 0.333, 0, 0, 0]

    def test_w2_userwarning(self):
        """A warning should be generated since the WLP stops at A_2 so it
        cannot be combined into a W2 WLP. Instead only [0,0,0] is returned
        """
        with pytest.warns(UserWarning):
            self.design.w2_wlp()
