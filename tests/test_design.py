"""
Test for the design class

Created on: 09/02/2023
Author: Alexandre Bohyn
"""
import fatld


class TestDesign:
    # Design 6 from Table 3 of Wu (1993)
    Design = fatld.Design(runsize=32, m=1, cols=[29, 26, 22])

    def test_twlp(self):
        twlp = self.Design.twlp(max_length=5)
        assert twlp == [[0, 0], [1, 4], [0, 2]]

    def test_twlp_type1(self):
        twlp = self.Design.twlp(type_0=False, max_length=5)
        assert twlp == [[0, 0], [4, 1], [2, 0]]

    def test_wlp(self):
        wlp = self.Design.wlp(max_length=5)
        assert wlp == [0, 5, 2]
