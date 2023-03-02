import warnings
from itertools import chain, combinations
from typing import Dict, List, Optional

import numpy as np
import oapackage as oa  # type: ignore
import pandas as pd  # type: ignore

from .main import basic_factor_matrix, custom_design, twlp
from .relation import Relation, num2gen


class Design:
    """
    Regular design with four- and two-level factors.

    Parameters
    ----------
    runsize : int
        Number of runs
    m : int
        Number of four-level factors
    cols: List[int]
        List of the added factors of the design. Cannot contain columns used
        in four-level factors and basic factors.

    Attributes
    ----------
    runsize : int
        Number of runs
    m : int
        Number of four-level factors
    cols: List[int]
        List of the column numbers of all the two-level factors of the design.
    k : int
        Number of basic factors. Equals to the log2 of the runsize.
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
    n: int
        Number of two-level factors in the design
    """

    # TODO: refactor the docstring of the Design class
    def __init__(self, runsize: int, m: int, cols: List[int]):
        # TODO: think about a `strict` keyword
        # it would bypass columns check and uses only the columns suplied
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
        if not all([(self.runsize > i > 0) for i in cols]):
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
        self.af.sort()
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

    def twlp(
        self, type_0: bool = True, max_length: Optional[int] = None
    ) -> List[List[int]]:
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

    def wlp(self, max_length: Optional[int] = None) -> List[int]:
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
        wlp_list = list(map(int, ar.GWLP()[3:]))
        if max_length is not None and (max_length <= 3 or max_length > len(wlp_list)):
            warnings.warn("Wrong max_length value, ignoring")
            return wlp_list
        elif max_length is None:
            return wlp_list
        else:
            return wlp_list[0 : (max_length - 2)]

    def flatten(self, zero_coding: bool = True) -> np.ndarray:
        """
        Flatten each four-level factor of the design into two independent two-level
        factors.

        Parameters
        ----------
        zero_coding : bool, optional
            The matrix is in 0/1 coding instead of -1/1 coding, by default True

        Returns
        -------
        np.ndarray
            The flattened design containing `2*m + n` two-level factors

        """
        flat_4lvl_part = np.zeros((self.runsize, 3 * self.m), dtype=int)
        for i in range(self.m):
            flat_4lvl_part[:, 3 * i] = self.array[:, i] > 1
            flat_4lvl_part[:, (3 * i + 1)] = self.array[:, i] % 2 == 0
            # Using !XOR so that -/- and +/+ both give +
            flat_4lvl_part[:, (3 * i + 2)] = np.logical_not(
                np.logical_xor(self.array[:, i] > 1, self.array[:, i] % 2 == 0)
            )
        mat = np.concatenate((flat_4lvl_part, self.array[:, self.m :]), axis=1)
        if zero_coding is False:
            return mat * 2 - 1
        else:
            return mat

    def tfi_clearance(self) -> Dict[str, Dict[str, object]]:
        """
        Compute clearance of all two-factor interactions in the design.

        There are three types of two-factor interactions in a design:

        - 4-4: interaction between two pseudo-factors of four-level factors
        - 4-2: interaction between a pseudo-factor of a four-level factor and
          a two-level factor.
        - 2-2: interaction between two two-level factors.

        An interaction is called clear when it is not aliased with any other interaction
        of a specific type. Thus, an interaction is '4-4 clear' when it is not aliased
        with any 4-4 interactions, '4-2 clear' when it is not aliased with any 4-2
        interactions, and '2-2 clear' when it is not aliased with any 2-2 interactions.
        When an interaction is 4-4 clear, 4-2 clear, and 2-2 clear, we call it a
        "totally clear (TC)" interaction.
        Defining how clear is an interaction, is called the *clearance* of an
        interaction.

        Returns
        -------

        A dictionary where each entry corresponds to an interaction of the design,
        where the key is the name of the interaction, and the value is another
        dictionary containing four entries:

            - '4-4': the interaction is '4-4' clear (True or False)
            - '4-2': the interaction is '4-2' clear (True or False)
            - '2-2': the interaction is '2-2' clear (True or False)
            - 'Type': the type of the interaction, either '4-4', '4-2', or '2-2'

        """
        # 4lvl PF are labeled as uppercase + number
        label_list_4lvl_pf = [
            [f"{chr(65 + i)}{j}" for j in [1, 2, 3]] for i in range(self.m)
        ]
        # 2lvl factors labeled as lowercase
        label_list_2lvl = [chr(97 + i) for i in range(self.n)]
        # Combinations of PF from the same four-level factors cannot be considered as
        # for any PF, p1 x p2 = p3
        if self.m > 1:
            label_list_44_tfi = list(
                chain(
                    *[
                        [f"{i}.{j}" for i in x for j in y]
                        for x, y in combinations(label_list_4lvl_pf, 2)
                    ]
                )
            )
        else:
            label_list_44_tfi = []
        label_list_42_tfi = [
            f"{i}.{j}"
            for i in list(chain(*label_list_4lvl_pf))
            for j in label_list_2lvl
        ]
        label_list_22_tfi = [f"{i}.{j}" for i, j in combinations(label_list_2lvl, 2)]
        # All combinations of labels are all the TFI
        label_list_tfi = label_list_44_tfi + label_list_42_tfi + label_list_22_tfi
        # tfi can be one of three types: 4-4, 4-2, or 2-2
        type_list_tfi = (
            ["4-4"] * len(label_list_44_tfi)
            + ["4-2"] * len(label_list_42_tfi)
            + ["2-2"] * len(label_list_22_tfi)
        )
        # All alias-related variables are NA and will be filled later on
        tfi = dict()
        for tfi_idx, label in enumerate(label_list_tfi):
            tfi_type = type_list_tfi[tfi_idx]
            tfi[label] = {"4-4": True, "4-2": True, "2-2": True, "Type": tfi_type}
        # Zero coding is needed to perform multiplication of factors
        flat_array_coded = self.flatten(zero_coding=False)
        four_level_pf_range = range(3 * self.m)
        two_level_fac_range = range(3 * self.m, (3 * self.m + self.n))
        matrix_42_tfi = np.hstack(
            [
                flat_array_coded[:, [i]] * flat_array_coded[:, [j]]
                for i in four_level_pf_range
                for j in two_level_fac_range
            ]
        )
        matrix_22_tfi = np.hstack(
            [
                flat_array_coded[:, [i]] * flat_array_coded[:, [j]]
                for i, j in combinations(two_level_fac_range, 2)
            ]
        )
        if self.m > 1:
            list_44_tfi = [
                [
                    flat_array_coded[:, [i]] * flat_array_coded[:, [j]]
                    for i in range(3 * x, (3 * x) + 3)
                    for j in range(3 * y, (3 * y) + 3)
                ]
                for x, y in combinations(range(self.m), 2)
            ]
            matrix_44_tfi = np.hstack(list(chain(*list_44_tfi)))
            tfi_mat = np.concatenate(
                (matrix_44_tfi, matrix_42_tfi, matrix_22_tfi), axis=1
            )
        else:
            tfi_mat = np.concatenate((matrix_42_tfi, matrix_22_tfi), axis=1)
        tfi_aliasing_mat = tfi_mat.T @ tfi_mat
        # All TFI are aliased with themselves
        np.fill_diagonal(tfi_aliasing_mat, 0)
        row, col = np.nonzero(tfi_aliasing_mat)
        for idx, row_num in enumerate(row):
            tfi_label = label_list_tfi[row_num]
            col_num = col[idx]
            alias_type = type_list_tfi[col_num]
            tfi[tfi_label][alias_type] = False
        return tfi

    def clarity(self) -> pd.DataFrame:
        """
        Generate the clarity matrix of the design.

        The clarity matrix is a table where:

        - the rows represent the type of the interactions (4-4, 4-2, 2-2,
            or any type)
        - the columns represent the type of clarity that the interactions can have (
            4-4 clear, 4-2 clear, 2-2 clear, or totally clear)

        Returns
        -------
        clarity_matrix
            A pandas dataframe holding the clarity matrix where the rows are
            ["4-4", "4-2", "2-2", "Any type"] and the columns are ["4-4 clear",
            "4-2 clear", "2-2 clear", "Totally clear"]
        """
        clarity_matrix = np.zeros((4, 4), dtype=int)
        clarity_df = pd.DataFrame(
            clarity_matrix,
            columns=["4-4 clear", "4-2 clear", "2-2 clear", "Totally clear"],
            index=["4-4", "4-2", "2-2", "Any type"],
        )
        tfi = self.tfi_clearance()
        for interaction in tfi.values():
            int_type = interaction["Type"]
            all_clear = 0
            for clear_type in ["4-4", "4-2", "2-2"]:
                if interaction[clear_type]:
                    column_label = f"{clear_type} clear"
                    clarity_df.loc[int_type, column_label] += 1
                    clarity_df.loc["Any type", column_label] += 1
                    all_clear += 1
                if all_clear == 3:
                    clarity_df.loc[int_type, "Totally clear"] += 1
                    clarity_df.loc["Any type", "Totally clear"] += 1
        return clarity_df

    def clear(self, interaction_type: str, clear_from: str) -> int:
        """
        Compute the number of interactions of type `interaction_type` that are clear
        from interactions of the `clear_from` type.

        Parameters
        ----------
        interaction_type : str
            Type of the interaction studied. Can either be "4-4", "4-2", "2-2", or "all"
            to count for all types of interaction.
        clear_from : str
            Type of clarity to count. Can either be "4-4", "4-2", "2-2", or "all" to
            count the number of totally clear interaction.

        """
        possible_types = ["4-4", "4-2", "2-2", "all"]
        if interaction_type not in possible_types or clear_from not in possible_types:
            raise ValueError(
                """
                Wrong type of interaction provided.
                Must be one of '4-4', '4-2', '2-2' or 'all'.
                """
            )
        tfi = self.tfi_clearance()
        count = 0
        for _, value in tfi.items():
            if interaction_type != "all" and value["Type"] != interaction_type:
                continue
            if clear_from == "all":
                if all([value[t] for t in ["4-4", "4-2", "2-2"]]):
                    count += 1
            elif value[clear_from]:
                count += 1
        return count

    def defining_relation(self, raw: bool = False):
        """
        Generate the defining relation of the design.

        For each added factor, generate the word representing the added factor.
        Each word contains the generator + the letter representing the added factor.
        The pseudo-factors used to generate the four-level factors are replaced by
        four-level factors labels in the final words.
        For example, if factor `A` is created using pseudo-factors `a`, `b` and `ab`,
        then the replacement scheme is: `a -> A1`, `b -> A2`, and `ab -> A3`.

        Parameters
        ----------
        raw : bool, optional
            Return the relation without replacing the pseudo-factors by their four-level
            factor label, by default False

        Returns
        -------
        Relation
            Defining relation as a ``Relation`` object.
        """
        relation = [
            f"{num2gen(x)}{chr(97 + self.k + i)}" for i, x in enumerate(self.af)
        ]
        return Relation(words=relation, m=self.m)

    def add_factor(self, number: int):
        """
        Add a two-level factor to the design.

        Generate a design with an added two-level factor defined by the number supplied.

        Parameters
        ----------
        number : int
            Column number of the factor to add. Cannot correspond to a column number
            already in use (as two-level factor or pseudo-factor) in the design.

        Returns
        -------
        Design
            A `Design` object containing the added two-level factor

        """
        # Factor cannot be used in the columns, the pf or be out of range
        if number < 0 or number >= self.runsize:
            raise ValueError(
                """Supplied column number must be greater than zero and less
                   than the runsize."""
            )
        if number in self.cols or number in list(chain(*self.pf)):
            raise ValueError(f"Column number {number} already used in the design")
        new_cols = self.af + [number]
        return Design(self.runsize, self.m, new_cols)

    def remove_factor(self, number: int):
        """
        Remove a two-level factor from the design.

        Parameters
        ----------
        number : int
            Column number of the factor to remove. Must correspond to one of the factor
            used in the design.

        Returns
        -------
        Design
            A `Design` object containing the added two-level factor

        """
        # Factor cannot be used in the columns, the pf or be out of range
        if number not in self.af:
            raise ValueError(
                f"Column number {number} is not an added factor in the design."
            )
        new_cols = [i for i in self.af]
        new_cols.remove(number)
        return Design(self.runsize, self.m, new_cols)
