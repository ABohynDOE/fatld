"""
Test for the design class

Created on: 09/02/2023
Author: Alexandre Bohyn
"""
import numpy as np
import pytest

import fatld


class TestDesignErrors:
    def test_wrong_init(self):
        assert True


class TestDesignMethods:
    # Design 6 from Table 3 of Wu (1993)
    design = fatld.Design(runsize=32, m=1, cols=[29, 26, 22])

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
