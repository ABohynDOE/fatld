import warnings
from itertools import chain
from typing import List

import numpy as np
import oapackage as oa

from .main import basic_factor_matrix, custom_design, twlp


class Design:
    """
    Regular design with four- and two-level factors.

    ...

    Attributes
    ----------
    runsize : int
        Number of runs
    k : int
        Number of basic factors. Equals to the log2 of the runsize.
    m : int
        Number of four-level factors
    pf : List[List[int]]
        List of the pseudo-factors triplets (as lists of integers), used to
        define the four-level factors.
    bf: List[int]
        List of the column numbers of the basic factors not used
        in the pseudo-factors `pf`.
    p : int
        Number of added factors
    af : List[int]
        List of the column numbers of the added factors
    cols: List[int]
        List of the column numbers of all the two-level factors of the design.
    n: int
        Number of two-level factors in the design
    """

    def __init__(self, runsize: int, m: int, cols: List[int]):
        # TODO: think about a `strict` keyword
        # it would bypass columns check and uses only the columns suplied
        """
        Parameters
        ----------
        runsize : int
            Number of runs
        m : int
            Number of four-level factors
        cols: List[int]
            List of the added factors of the design. Cannot contain columns used
            in four-level factors and basic factors.
        """
        # Run size value check
        if not isinstance(runsize, int) or runsize <= 0:
            raise TypeError("Runsize must be a positive integer")
        if np.log2(runsize) % 1 != 0:
            raise ValueError("Runsize must be a power of two")
        self.runsize = runsize

        # Define value of k
        self.k = int(np.log2(runsize))

        # m value check
        if not isinstance(m, int) or runsize <= 0:
            raise TypeError("Runsize must be a positive integer")
        if m not in [1, 2, 3]:
            raise ValueError("m can only have a value of 1, 2 or 3")
        if 2 * m > self.k:
            raise ValueError(
                """`m` must be lower of equal to `k/2`
                Four-level factors are made of basic factors."""
            )
        self.m = m

        # Define the pseudo-factors based on the value of m
        self.pf = [
            [2 ** (2 * m), 2 ** ((2 * m) + 1), 2 ** (2 * m) + 2 ** ((2 * m) + 1)]
            for m in range(self.m)
        ]

        # Available basic factors are the ones not used in the four-level factors
        all_basic_factors = [2**i for i in range(self.k)]
        self.bf = [i for i in all_basic_factors if i not in chain(*self.pf)]

        # Cols value check
        if not isinstance(cols, list):
            raise TypeError("`cols` must be a list of integer")
        if not all([isinstance(i, int) for i in cols]):
            raise TypeError("All column numbers must be positive integers")
        if not all([(i < self.runsize and i > 0) for i in cols]):
            raise ValueError(
                "Column numbers cannot be negative or greater than the runsize"
            )
        # Added factors cannot be used in pseudo-factors
        if any([i in chain(*self.pf) for i in cols]):
            raise ValueError(
                f"""Some column numbers are used to define pseudo-factors 
                for the four-level factor: {self.pf}"""
            )
        # User should know if columns are basic factors
        if any([i in self.bf for i in cols]):
            cols = [i for i in cols if i not in self.bf]
            warnings.warn(
                f"Some column numbers are basic factors: {self.bf}", UserWarning
            )
        # Number of added factors
        self.p = len(cols)
        # List of added factors
        self.af = cols
        # List of all two-level columns in the design
        self.cols = self.bf + self.af
        self.cols.sort()
        # Number of two-level factors
        self.n = len(self.cols)

        # Generate the design matrix
        # 4-level part
        coeff_matrix_4lvl = np.zeros((2 * self.m, self.m), dtype=int)
        bf_matrix = basic_factor_matrix(self.k)
        for i in range(self.m):
            idx = [2 * i, ((2 * i) + 1)]
            coeff_matrix_4lvl[idx, i] = [2, 1]
            four_lvl_part = bf_matrix[:, 0 : (2 * self.m)] @ coeff_matrix_4lvl
        # 2-level part
        two_lvl_part = custom_design(self.runsize, self.cols)
        # Assemble into one matrix
        self.array = np.concatenate([four_lvl_part, two_lvl_part], axis=1)

    def __str__(self):
        return f"4^{self.m} 2^{self.n}-{self.p} design in {self.runsize} runs"

    def __repr__(self):
        return f"Design(runsize={self.runsize}, m={self.m}, cols={self.af})"

    def twlp(self, type_0: bool = True, max_length: int = None) -> List[List[int]]:
        """Type-specific word length pattern

        Compute the type-specific word length pattern of a design, starting with words
        of length 3.

        Parameters
        ----------
        type_0 : bool, optional
            Return the word length pattern of type 0, by default True.
            If False, return the word length pattern of type `m`.
        max_length : int, optional
            Max word length to display in the word length pattern, by default None,
            all the word lengths are displayed.
            If this number is more than the maximum word length of the design, ignores
            it.

        Returns
        -------
        List[List[int]]
            Type-specific word length pattern, starting with words of length 3.
            Each sublist give the number of words of that length, sorted by type.
        """
        ar = oa.array_link(self.array)
        return twlp(ar, type_0, max_length)

    def wlp(self, max_length: int = None) -> List[int]:
        """Generalized word length pattern

        Compute the word length pattern, i.e., the number of words in the defining
        relation of the design, and arrange them by length.

        Parameters
        ----------
        max_length : int, optional
            Max word length to display in the word length pattern, by default None,
            all the word lengths are displayed.
            If this number is more than the maximum word length of the design, ignores
            it.

        Returns
        -------
        List[int]
            Word length pattern, starting with words of length 3.
        """
        ar = oa.array_link(self.array)
        wlp_list = list(map(int ,ar.GWLP()[3:]))
        if max_length is not None and (max_length <= 3 or max_length > len(wlp_list)):
            warnings.warn("Wrong max_length value, ignoring")
            return wlp_list
        elif max_length is None:
            return wlp_list
        else:
            return wlp_list[0 : (max_length - 2)]
