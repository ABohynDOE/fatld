"""
Functions to create a FATL design

Author: Alexandre Bohyn
"""
# %% Packages
import warnings
import numpy as np

from typing import List


def basic_factor_matrix(k: int, zero_coding: bool = True) -> np.ndarray:
    """
    Create a matrix containing `k` basic factors.

    Parameters
    ----------
    k : int
        Number of basic factors
        zero_coding : bool, optional
            The matrix is in 0/1 coding instead of -1/1 coding, by default True

        Returns
    -------
    np.ndarray
        Basic factor matrix
    """
    mat = np.zeros((2**k, k), dtype=int)
    for i in range(k):
        unit = [0] * (2**k // 2 ** (i + 1)) + [1] * (2**k // 2 ** (i + 1))
        mat[:, i] = unit * 2**i
    if not zero_coding:
        mat = (mat * 2) - 1
    return mat


def power2_decomposition(n: int, length: int = None) -> List[int]:
    """
    Decompose a number into powers powers of 2 and returns the corresponding indices.

    Parameters
    ----------
    n : int
        Number to decompose
    length : int, optional
        Wanted length for the list of indices, by default None, returns the list when
        the last power of two `x` for which `x` < `n` is reached

    Returns
    -------
    List[int]
        List of the indices of the powers of two obtained in the decomposition

    Example
    -------
    >>> # 7 can be decomposed into 1*2 + 1*2 + 1*4
    >>> power2_decomposition(7)
    [1,1,1]
    >>> # 11 can be decomposed into 1*1 + 1*2 + 0*4 + 1*8, we want a list of length 7
    >>> power2_decomposition(11, length=7)
    [1,1,0,1,0,0,0]
    """
    powers = []
    i = 1
    while i <= n:
        if i & n:
            powers.append(1)
        else:
            powers.append(0)
        i <<= 1
    if length is not None and len(powers) < length:
        powers += [0] * (length - len(powers))
    elif length is not None and len(powers) > length:
        warnings.warn("Power list larger than supplied length", UserWarning)
    return powers


def custom_design(runsize: int, column_list: List[int], zero_coding=True) -> np.ndarray:
    """
    Create a custom design with a specific run size, based on the column numbers
    provided.

    Parameters
    ----------
    runsize : int
        Number of runs
         : List[int]
        Column numbers to be used for the design
    zero_coding : bool, optional
        he matrix is in 0/1 coding instead of -1/1 coding, by default True

    Returns
    -------
    np.ndarray
        Design matrix with the factors corresponding to the columns given (and
        in the same order)
    """
    k = int(np.log2(runsize))
    bf_mat = basic_factor_matrix(k, zero_coding=True)
    power_list = [power2_decomposition(n=col, length=k) for col in column_list]
    yates_matrix = np.array(power_list).T
    custom_mat = (bf_mat @ yates_matrix) % 2
    return custom_mat
