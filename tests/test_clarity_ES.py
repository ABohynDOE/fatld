"""
Compare the output of the `clarity` function with the data generated on 4^1 2^9
64-run resolution IV designs from ES.
"""
import fatld
import pytest  # type: ignore
import pandas as pd  # type: ignore


cols_df = pd.read_table("tests/data/columns.txt")
cols = cols_df["columns"].to_list()[:10]
ref_df = pd.read_table("tests/data/reference.txt", sep="\t")

names = "cols,t22_clear22,t22_clear42,tc22,tc42"
values = []
for i, c in enumerate(cols):
    r = ref_df.iloc[i].to_list()
    values.append((c, r[0], r[1], r[2], r[3]))


@pytest.mark.parametrize(names, values)  # type: ignore
def test_clarity_tc22(cols, t22_clear22, t22_clear42, tc42, tc22):
    cols_list = list(map(int, cols.split(",")))
    D = fatld.Design(64, 1, cols_list)
    assert D.clear("2-2", "all") == tc22


@pytest.mark.parametrize(names, values)  # type: ignore
def test_clarity_tc42(cols, t22_clear22, t22_clear42, tc42, tc22):
    cols_list = list(map(int, cols.split(",")))
    D = fatld.Design(64, 1, cols_list)
    assert D.clear("4-2", "all") == tc42


@pytest.mark.parametrize(names, values)  # type: ignore
def test_clarity_t22_clear22(cols, t22_clear22, t22_clear42, tc42, tc22):
    cols_list = list(map(int, cols.split(",")))
    D = fatld.Design(64, 1, cols_list)
    assert D.clear("2-2", "2-2") == t22_clear22


@pytest.mark.parametrize(names, values)  # type: ignore
def test_clarity_t42_clear22(cols, t22_clear22, t22_clear42, tc42, tc22):
    cols_list = list(map(int, cols.split(",")))
    D = fatld.Design(64, 1, cols_list)
    assert D.clear("2-2", "4-2") == t22_clear42
