"""
Test for the design class

Created on: 09/02/2023
Author: Alexandre Bohyn
"""
import numpy as np
import pandas as pd  # type: ignore
import pytest  # type: ignore
import itertools

from collections import Counter

import fatld


class TestDesignErrors:
    def test_wrong_init(self):
        assert True


class TestDesignMethods:
    # Design 6 from Table 3 of Wu (1993)
    design = fatld.Design(runsize=32, m=1, cols=[29, 26, 22])
    relation = fatld.relation.Relation(["bcef", "bdeg", "acdeh"], m=1)  # type: ignore

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
            for col in ["4-4 clear", "4-2 clear", "2-2 clear", "Totally clear"]:
                if col == "Totally clear":
                    col_index = "all"
                else:
                    col_index = col.replace(" clear", "")
                assertion_list.append(
                    self.design.clear(row_index, col_index) == self.df.loc[row, col]
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


def test_beta_aberration():
    D = fatld.Design(32, 1, [12, 20, 24, 29])

    # Compute the 12 permutations of the four levels
    perms = []
    for p in itertools.permutations(range(3)):
        for i in range(4):
            ordering = list(p)
            ordering.insert(i, 3)
            perms.append(ordering)
    perms = perms[:12]

    # Compute the unique qWLP
    unique_qwlp = []

    for i, p in enumerate(itertools.product(range(12), repeat=D.m)):
        perm_list = [perms[i] for i in p]
        beta_vector = fatld.design.beta_wlp(
            D,
            permutation_list=perm_list,
        )
        if beta_vector not in unique_qwlp:
            unique_qwlp.append(beta_vector)

    final_list = [
        [0.0, 0.0, 3.0, 4.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 2.4, 3.2, 0.6, 0.8, 0.8, 0.0, 0.2],
        [0.0, 0.6, 0.8, 2.4, 3.2, 0.2, 0.0, 0.8],
    ]

    assert unique_qwlp == final_list


def test_qwlp():
    D = fatld.Design(32, 2, [17, 18, 20, 24, 19, 28, 5, 9, 6, 10, 7, 11, 13])
    vec, perm = D.qwlp(max_length=5)
    assert vec == [8.0, 37.0, 107.0]


def test_beta_aberration_low_n():
    D = fatld.Design(32, 2, [31])  # Here n = 2 which is lower than the max_length set
    qvec, perm = D.qwlp(max_length=6)
    assert qvec == [0.0, 0.0, 0.0, 1.0]


def test_number_interactions():
    assert fatld.design.nbr_interactions(7, 5) == 119


def test_tfi_model_matrix():
    tfi_model_matrix = fatld.design.build_tfi_model_matrix(n=7, max_length=5)
    int_length = np.sum(tfi_model_matrix, axis=0).tolist()
    c = Counter(int_length)
    d = dict(c)
    assert d == {1: 7, 2: 21, 3: 35, 4: 35, 5: 21}
