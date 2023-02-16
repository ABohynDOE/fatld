"""
Functions to create a FATL design

Author: Alexandre Bohyn
"""
# %% Packages
import warnings
from typing import List, Optional

import numpy as np
import oapackage as oa  # type: ignore


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
        A `2^k` by `k` matrix, containing the `k` basic factors
    """
    mat = np.zeros((2**k, k), dtype=int)
    for i in range(k):
        unit = [0] * (2**k // 2 ** (i + 1)) + [1] * (2**k // 2 ** (i + 1))
        mat[:, i] = unit * 2**i
    if not zero_coding:
        mat = (mat * 2) - 1
    return mat


def power2_decomposition(n: int, length: Optional[int] = None) -> List[int]:
    """
    Decompose a number into powers of 2 and returns the corresponding indices.

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


def custom_design(runsize: int, column_list: List[int]) -> np.ndarray:
    """
    Create a custom design with a specific run size, based on the column numbers
    provided.

    Parameters
    ----------
    runsize : int
        Number of runs
    column_list : List[int]
        Column numbers to be used for the design

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


def twlp(
    ar: oa.array_link, type_0: bool = True, max_length: Optional[int] = None
) -> List[List[int]]:
    """
    Compute the type-specific word length pattern of a design, starting with words
    of length 3.

    Parameters
    ----------
    ar : oa.array_link
        Design matrix
    type_0 : bool, optional
        Return the word length pattern of type 0, by default True.
        If False, return the word length pattern of type `m`.
    max_length : int, optional
        Max word length to display in the word length pattern, by default None,
        all the word lengths are displayed.
        If this number is more than the maximum word length of the design, ignores it.

    Returns
    -------
    List[List[int]]
        Type-specific word length pattern, starting with words of length 3.
        Each sublist give the number of words of that length, sorted by type.
    """
    run_size = ar.shape[0]
    # Distance distribution of Hamming distance between the rows
    distance_distribution = oa.distance_distribution_mixed(array=ar, verbose=0)
    # Compute McWilliams transform of the distance distribution based on array
    # class
    array_class = oa.arraylink2arraydata(ar)
    dst_dist_transformed = oa.macwilliams_transform_mixed(
        B=distance_distribution,
        N=run_size,
        factor_levels_for_groups=array_class.factor_levels_column_groups(),
        verbose=0,
    )
    # Only full aliasing is possible in regular designs
    dst_dist_transformed_mat = np.array(dst_dist_transformed).astype(int)
    # Transformed distance distribution gives the shape of the WLP matrix
    r, c = dst_dist_transformed_mat.shape
    # Columns are length and rows are type
    max_type = r
    max_word_length = r + c - 1
    wlp_matrix = np.zeros((max_type, max_word_length), dtype=int)
    # Sum the transformed matrix over diagonal rows
    for idx in range(max_word_length):
        if idx < r and idx < c:
            comb = [(i, idx - i) for i in range(idx + 1)]
        else:
            comb = [(i, idx - i) for i in range(idx + 1) if ((i < r) and (idx - i) < c)]
        for coord in comb:
            wlp_matrix[coord[0], idx] = dst_dist_transformed_mat[coord[0], coord[1]]
    # Trim the matrix to max_length if needed
    if max_length is not None and (max_length < 3 or max_length > max_word_length):
        max_length = None
        warnings.warn("Incorrect max_length. Full type-specific WLP is outputted")
    if max_length is None:
        wlp_matrix_trimmed = wlp_matrix[:, 3:]
    else:
        wlp_matrix_trimmed = wlp_matrix[:, 3 : (max_length + 1)]
    # Create a list of lists
    wlp_type_list = wlp_matrix_trimmed.T.tolist()
    # Reverse the sublists if it's not type 0 because they are ordered in type 0
    # by default
    if type_0 is False:
        return [x[::-1] for x in wlp_type_list]
    else:
        return wlp_type_list


def num2gen(n: int) -> str:
    """Return the generator corresponding to a given column number

    Parameters
    ----------
    n : int
        Column number

    Returns
    -------
    str
        Corresponding generator

    Examples
    ------
    >>> num2gen(7)
    'abc'

    See Also
    --------
    gen2num : Convert a generator into the corresponding column number
    """
    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    powers = power2_decomposition(n)
    word = "".join([chr(97 + i) for i, x in enumerate(powers) if x == 1])
    return word


def gen2num(g: str) -> int:
    """Convert a generator into its corresponding column number

    Generators can only contain lowercase letters from a to z.

    Parameters
    ----------
    g : str
        Generator

    Returns
    -------
    int
        Corresponding column number

    Examples
    --------
    >>> gen2num('acd')
    13

    See Also
    --------
    num2gen : Convert a column number into its corresponding generator
    """
    if not isinstance(g, str):
        raise TypeError("g must be a non-empty string")
    ascii_code = [ord(i) for i in g]
    if any([i > 122 or i < 97 for i in ascii_code]):
        raise ValueError("Only lowercase letters are allowed in the generator")
    return sum([2 ** (i - 97) for i in ascii_code])
