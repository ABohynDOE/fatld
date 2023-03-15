import pytest  # type: ignore

from fatld.relation import (
    gen2num,
    num2gen,
    relabel_word,
    word_length,
    word_type,
    Relation,
)


@pytest.mark.parametrize(
    "num,mval,expected",
    [
        (7, None, "abc"),
        (7, 1, "A3c"),
        (13, 2, "A1B3"),
    ],
)  # type: ignore
def test_num2gen(num, mval, expected):
    assert num2gen(num, m=mval) == expected


def test_num2gen_typeerror():
    with pytest.raises(TypeError):
        num2gen("a")  # type: ignore


def test_num2gen_valueerror():
    with pytest.raises(ValueError):
        num2gen(7, m=4)


def test_gen2num():
    assert gen2num("acd") == 13


def test_gen2num_type_error():
    with pytest.raises(TypeError):
        gen2num(7)  # type: ignore


def test_gen2num_value_error():
    with pytest.raises(ValueError):
        gen2num("Abc")


def test_relabel():
    assert relabel_word(word="abcde", m=2) == "A3B3e"


def test_word_length():
    assert word_length(word="A1cde") == 4


def test_word_type():
    assert word_type(word="A1B3e") == 2


class TestRelation:
    rel_2lvl = Relation(words=["abcdf", "aceg"])
    rel_4lvl = Relation(words=["abcdf", "aceg"], m=1)

    def test_wrong_m(self):
        with pytest.raises(ValueError):
            Relation(["abcde"], m=4)

    def test_wrong_word(self):
        with pytest.raises(ValueError):
            Relation(["abde", "Acdf"])

    def test_single_word(self):
        with pytest.raises(TypeError):
            Relation("abcde")  # type: ignore

    def test_wrong_word_type(self):
        with pytest.raises(TypeError):
            Relation([21, 29, 27], m=1)  # type: ignore

    def test_expand(self):
        assert self.rel_2lvl.expand() == ["abcdf", "aceg", "bdefg"]

    def test_expand_relabel(self):
        assert self.rel_4lvl.expand(relabel=True) == ["A3cdf", "A1ceg", "A2defg"]

    def test_wlp(self):
        assert self.rel_2lvl.word_length_pattern() == [0, 1, 2]

    def test_wlp_relabel(self):
        assert self.rel_4lvl.word_length_pattern() == [[0, 0], [0, 2], [0, 1]]
