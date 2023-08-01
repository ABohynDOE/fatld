from __future__ import annotations
import warnings
from itertools import chain, combinations, product, permutations
from math import comb

import numpy as np
import oapackage as oa  # type: ignore
import pandas as pd  # type: ignore

from .main import basic_factor_matrix, custom_design, twlp
from .relation import Relation, num2gen, gen2num

# Setup the contrast matrix for the beta aberration
contrast_matrix = np.array(
    [[-3.0, 1.0, -1.0], [-1.0, -1.0, 3.0], [1.0, -1.0, -3.0], [3.0, 1.0, 1.0]]
)
scaled_contrast_matrix = np.divide(
    np.multiply(2.0, contrast_matrix), np.linalg.norm(contrast_matrix, axis=0)
)


class Design:
    """
    Regular design with four-level and two-level factors.

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
    n : int
        Number of two-level factors in the design
    k : int
        Number of basic factors. Equals to the log2 of the runsize.
    p : int
        Number of added factors

    bf : List[int]
        List of the column numbers of the basic factors not used in the
        pseudo-factors `pf`.
    pf : List[List[int]]
        List of the pseudo-factors triplets (as lists of integers), used to
        define the four-level factors.
    af : List[int]
        List of the column numbers of the added factors
    cols : List[int]
        List of the column numbers of all the two-level factors of the design.

    array : np.ndarray
        Design matrix of the design with the four-level factors first and then
        the two-level factors. The two-level factors are ordered by column
        number and not kept in the original order.
    """

    # TODO: refactor the docstring of the Design class
    def __init__(self, runsize: int, m: int, cols: list[int]):
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
            [
                2 ** (2 * m),
                2 ** ((2 * m) + 1),
                2 ** (2 * m) + 2 ** ((2 * m) + 1),
            ]
            for m in range(self.m)
        ]

        # Available basic factors are the ones not used in the four-level
        # factors
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
                f"Some column numbers are basic factors: {self.bf}",
                UserWarning,
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
        end_four_lvl_part = 2 * self.m
        four_lvl_part = (
            bf_matrix[:, 0:end_four_lvl_part] @ coeff_matrix_4lvl
        )  # noqa: E203
        # 2-level part
        two_lvl_part = custom_design(self.runsize, self.cols)
        # Assemble into one matrix
        self.array = np.concatenate([four_lvl_part, two_lvl_part], axis=1)

    def __str__(self):
        return f"4^{self.m} 2^{self.n}-{self.p} design in {self.runsize} runs"

    def __repr__(self):
        return f"Design(runsize={self.runsize}, m={self.m}, cols={self.af})"

    def twlp(
        self, type_0: bool = True, max_length: int | None = None
    ) -> list[list[int]]:
        """Type-specific word length pattern

        Compute the type-specific word length pattern of a design, starting
        with words of length 3.

        Parameters
        ----------
        type_0 : bool, optional
            Return the word length pattern of type 0, by default True.
            If False, return the word length pattern of type `m`.
        max_length : int, optional
            Max word length to display in the word length pattern, by default
            None, all the word lengths are displayed.
            If this number is more than the maximum word length of the design,
            ignores it.

        Returns
        -------
        List[List[int]]
            Type-specific word length pattern, starting with words of length 3.
            Each sublist give the number of words of that length, sorted by
            type.
        """
        ar = oa.array_link(self.array)
        return twlp(ar, type_0, max_length)

    def wlp(self, max_length: int | None = None) -> list[int]:
        """Generalized word length pattern

        Compute the word length pattern, i.e., the number of words in the
        defining relation of the design, and arrange them by length.

        Parameters
        ----------
        max_length : int, optional
            Max word length to display in the word length pattern, by default
            None, all the word lengths are displayed.
            If this number is more than the maximum word length of the design,
            ignores it.

        Returns
        -------
        List[int]
            Word length pattern, starting with words of length 3.
        """
        ar = oa.array_link(self.array)
        wlp_list = list(map(int, ar.GWLP()[3:]))
        if max_length is not None and (
            max_length <= 3 or max_length > len(wlp_list)
        ):
            warnings.warn("Wrong max_length value, ignoring")
            return wlp_list
        elif max_length is None:
            return wlp_list
        else:
            return wlp_list[0 : (max_length - 2)]  # noqa: E203

    def beta_star_wlp(
        self, max_length: int | None = None
    ) -> "tuple[list[tuple[float]], list[list[int]]]":
        """Compute the 𝛽* vector for a design.
        It contains tuples with the 𝛽*_i value, and the A_i^0 value for values
        of i starting with 3.

        Parameters
        ----------
        design : fatld.Design
            A design object
        max_length : int | None
            Maximum value of i for which the 𝛽* is computed. Default is None,
            so all length are considered.

        Returns
        -------
        beta_star_vector : list[tuple[float | int]]
            The 𝛽* vector with tuples of 𝛽* and A_i^0 values
        perm : list[list[int]]
            Permutation of the four-level factors associated with the 𝛽* vector
        """
        # Compute the 12 basic permutations of the four levels
        perms = []
        for p in permutations(range(3)):
            for i in range(4):
                ordering = list(p)
                ordering.insert(i, 3)
                perms.append(ordering)
        perms = perms[:12]

        # Compute all possible permutations for m factors
        all_perms = []
        for _, p in enumerate(product(range(12), repeat=self.m)):
            all_perms.append([perms[i] for i in p])

        # Isolate four-level and two-level part of the design
        m = self.m
        fl_matrix = self.array[:, :m]
        tl_matrix = self.array[:, m:]  # noqa: E201

        # Build all interactions between all two-level factors
        n = self.n
        if max_length is None:
            max_n_value = n
        else:
            max_n_value = max_length - 1
        tfi_model_matrix = build_tfi_model_matrix(n=n, max_length=max_n_value)
        tl_interaction_matrix = np.matmul(tl_matrix, tfi_model_matrix) % 2
        tl_interaction_matrix = (2 * tl_interaction_matrix) - 1
        tl_interaction_length = np.sum(tfi_model_matrix, axis=0).tolist()

        # Compute beta wlp for all permutations
        beta_vector_list = []
        for perms in all_perms:
            beta_vector = compute_correlations(
                tl_interaction_matrix=tl_interaction_matrix,
                tl_interaction_length=tl_interaction_length,
                fl_matrix=fl_matrix,
                permutations_list=perms,
                max_length=max_length,
            )
            beta_vector_list.append(beta_vector)

        # Find the best Beta vector and the corresponding permutation
        best_vector = beta_vector_list[0]
        best_perm = all_perms[0]

        for i, v in enumerate(beta_vector_list[1:]):
            if best_vector > v:
                best_vector = v
                best_perm = all_perms[i]

        # Return a combination of the A_i^0 values and the beta values
        beta_star_vector = []
        twlp = self.twlp(type_0=True)
        for i, x in enumerate(best_vector):
            if i >= len(twlp):
                beta_star_vector.append([x, 0])
            else:
                beta_star_vector.append([x, twlp[i][0]])

        return beta_star_vector, best_perm

    def beta_wlp(
        self, max_length: int | None = None
    ) -> tuple[list[float], list[int]]:
        """
        Find the permutation of the levels of the m four-level factors that
        minimize the qWLP for the design.

        Parameters
        ----------
        max_length : int, optional
            Max word length to consider in the qWLP, by default None,
            all the word lengths are computed (careful as it is computationally
            intensive for large values of n).

        Returns
        -------
        best_wlp, permutations : tuple[list[float], list[int]]
            Returns the minimal qWLP and the corresponding permutations of the
            factor levels. The qWLP starts with words of length 3.
        """
        beta_vector, perm = self.beta_star_wlp(max_length=max_length)
        return [x[0] for x in beta_vector], perm

    def w2_wlp(self) -> tuple[list[int], str]:
        """
        Compute the W_2 word length pattern of the design.
        The structure of the W_2 vector is:
        A3.0, A2.1, A4.0, A5.0, A3.1, A6.0, A7.0, ...

        Returns
        -------
        W2_vector: list[int]
            The vector containing the A_x.i values of the W_2 word length
            pattern
        factor: str
            A string indicating which factor must be used to obtain W_2 optimal
            blocking. Can either be 'A', 'B', or 'C'.

        """
        design_oa = oa.array_link(self.array)
        design_wlp = list(map(int, design_oa.GWLP()))

        best_factor = "A"
        best_w2_vector = []
        for factor_index in range(self.m):
            # Remove the blocking factor and evaluate the WLP of the design
            trmt_wlp = list(
                map(int, design_oa.deleteColumn(factor_index).GWLP())
            )
            block_wlp = [design_wlp[i] - x for i, x in enumerate(trmt_wlp)]
            w2_vector = build_w2_vector(block_wlp, trmt_wlp)
            # Nothing is returned if the vectors cannot be combined
            if w2_vector is None:
                continue
            # Set as default if it is the first factor we evaluate
            if factor_index == 0 or best_w2_vector == []:
                best_w2_vector = w2_vector
            else:
                if w2_vector < best_w2_vector:
                    best_w2_vector = w2_vector
                    best_factor = chr(65 + factor_index)
        return best_w2_vector, best_factor

    def alpha_wlp(self, rounding: bool = True) -> list[float]:
        """
        Compute the alpha word length pattern of the design.
        The alpha wlp contains 5 values ordered in the following way:
        ω4, ω2, ω42, ω22, ω44

        Parameters
        ----------
        rounding : bool, optional
            Round the ω values to 2 decimals, by default True

        Returns
        -------
        list[float]
            The alpha word length pattern
        """
        twlp = self.twlp(max_length=4)
        # For cases with designs with not enough words
        if len(twlp) == 0:
            twlp = [[0, 0, 0], [0, 0, 0, 0]]
        elif len(twlp) == 1:
            twlp.append([0, 0, 0, 0])
        a3_vector = twlp[0] + [0] * (3 - len(twlp[0]))  # Right pad with 0
        a4_vector = twlp[1] + [0] * (4 - len(twlp[1]))
        omega_weights = {
            "w2": [3 / 3, 2 / 3, 1 / 3],
            "w4": [0 / 3, 1 / 3, 2 / 3],
            "w22": [6 / 6, 3 / 6, 1 / 6, 0 / 6],
            "w42": [0 / 6, 3 / 6, 4 / 6, 3 / 6],
            "w44": [0 / 6, 0 / 6, 1 / 6, 3 / 6],
        }
        omega_values = dict()
        for name, weight_vector in omega_weights.items():
            if len(name) == 2:
                value = [w * a for w, a in zip(weight_vector, a3_vector)]
            else:
                value = [w * a for w, a in zip(weight_vector, a4_vector)]
            omega_values[name] = sum(value)
        alpha_wlp = [
            omega_values[i] for i in ["w4", "w2", "w42", "w22", "w44"]
        ]
        if rounding:
            alpha_wlp = [round(i, 3) for i in alpha_wlp]
        return alpha_wlp

    def resolution(self) -> int | None:
        """
        Compute the resolution of the design.
        """
        wlp = self.wlp()
        res = next((i + 3 for i, x in enumerate(wlp) if x), None)
        return res

    def flatten(self, zero_coding: bool = True) -> np.ndarray:
        """
        Flatten each four-level factor of the design into three two-level
        pseudo-factors.

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
        mat = np.concatenate(
            (flat_4lvl_part, self.array[:, self.m :]), axis=1  # noqa: E203
        )
        if zero_coding is False:
            return mat * 2 - 1
        else:
            return mat

    def tfi_clearance(self) -> dict[str, dict[str, object]]:
        """
        Compute clarity of all two-factor interactions in the design.

        There are three types of two-factor interactions in a design:

        - 4-4: interaction between two pseudo-factors of four-level factors
        - 4-2: interaction between a pseudo-factor of a four-level factor and
               a two-level factor.
        - 2-2: interaction between two two-level factors.

        An interaction is called clear when it is not aliased with any other
        interaction of a specific type. Thus, an interaction is '4-4 clear'
        when it is not aliased with any 4-4 interactions, '4-2 clear' when it
        is not aliased with any 4-2 interactions, and '2-2 clear' when it is
        not aliased with any 2-2 interactions.
        When an interaction is 4-4 clear, 4-2 clear, and 2-2 clear, we call it
        a totally clear (TC)" interaction.
        Defining how clear is an interaction, is called the *clarity* of an
        interaction.

        Returns
        -------

        A dictionary where each entry corresponds to an interaction of the
        design, where the key is the name of the interaction, and the value is
        another dictionary containing four entries:

            - 'clear': the interaction is clear from main effects
            - '4-4': the interaction is '4-4' clear (True or False)
            - '4-2': the interaction is '4-2' clear (True or False)
            - '2-2': the interaction is '2-2' clear (True or False)
            - 'Type': the type of the interaction, either '4-4', '4-2',
            or '2-2'

        """
        # LABELS
        # MAIN EFFECTS

        # 4lvl PF are labeled as uppercase + index i = 1,2,3
        label_list_4lvl_pf = [
            [f"{chr(65 + i)}{j}" for j in [1, 2, 3]] for i in range(self.m)
        ]
        # 2lvl factors labeled as lowercase
        label_list_2lvl = [chr(97 + i) for i in range(self.n)]

        # TWO-FACTOR-INTERACTIONS

        # Combinations of PF from the same four-level factors cannot be
        # considered as for any PF, p1 x p2 = p3 so combine BETWEEN the
        # pseudo-factors triplets and not within
        if self.m > 1:
            label_list_44_tfi = list(
                chain(
                    *[
                        [f"{i}.{j}" for i in x for j in y]
                        for x, y in combinations(label_list_4lvl_pf, 2)
                    ]
                )
            )
        # There are no 4-4 interactions for m=1
        else:
            label_list_44_tfi = []
        # A 4-2 interaction is any combination of a PF and a two-level factor
        label_list_42_tfi = [
            f"{i}.{j}"
            for i in list(chain(*label_list_4lvl_pf))
            for j in label_list_2lvl
        ]
        label_list_22_tfi = [
            f"{i}.{j}" for i, j in combinations(label_list_2lvl, 2)
        ]

        # All combinations of labels are all the TFI
        label_list_tfi = (
            label_list_44_tfi + label_list_42_tfi + label_list_22_tfi
        )

        # MATRIX

        # Zero coding (0/1) is needed to perform multiplication of factors
        flat_array_coded = self.flatten(zero_coding=False)

        # MAIN EFFECTS

        two_level_fac_range = range(3 * self.m, (3 * self.m + self.n))
        matrix_2_me = flat_array_coded[:, two_level_fac_range]
        four_level_pf_range = range(3 * self.m)
        matrix_4_me = flat_array_coded[:, four_level_pf_range]

        # TWO-FACTOR INTERACTIONS

        matrix_22_tfi = np.hstack(
            [
                flat_array_coded[:, [i]] * flat_array_coded[:, [j]]
                for i, j in combinations(two_level_fac_range, 2)
            ]
        )
        matrix_42_tfi = np.hstack(
            [
                flat_array_coded[:, [i]] * flat_array_coded[:, [j]]
                for i in four_level_pf_range
                for j in two_level_fac_range
            ]
        )

        # We only compute 44 interactions if m > 1
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
            effect_mat = np.concatenate(
                (
                    matrix_44_tfi,
                    matrix_42_tfi,
                    matrix_22_tfi,
                    matrix_2_me,
                    matrix_4_me,
                ),
                axis=1,
            )
            # We don't need to define clarity for the ME
            tfi_mat = np.concatenate(
                (matrix_44_tfi, matrix_42_tfi, matrix_22_tfi),
                axis=1,
            )
        else:
            effect_mat = np.concatenate(
                (matrix_42_tfi, matrix_22_tfi, matrix_2_me, matrix_4_me),
                axis=1,
            )
            # We don't need to define clarity for the ME
            tfi_mat = np.concatenate(
                (matrix_42_tfi, matrix_22_tfi),
                axis=1,
            )

        # ALIASING

        # tfi can be one of three types: 4-4, 4-2, or 2-2
        type_list_tfi = (
            ["4-4"] * len(label_list_44_tfi)
            + ["4-2"] * len(label_list_42_tfi)  # noqa : W503
            + ["2-2"] * len(label_list_22_tfi)  # noqa : W503
            + ["clear"]  # noqa : W503
            * (
                len(label_list_2lvl) + 3 * len(label_list_4lvl_pf)
            )  # noqa : W503
        )

        # All alias-related variables are NA and will be filled later on
        tfi = dict()
        for tfi_idx, label in enumerate(label_list_tfi):
            tfi_type = type_list_tfi[tfi_idx]
            tfi[label] = {
                "clear": True,
                "4-4": True,
                "4-2": True,
                "2-2": True,
                "Type": tfi_type,
            }

        # Rows are the TFI we are investigating
        tfi_aliasing_mat = tfi_mat.T @ effect_mat
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
        - the columns represent the type of clarity that the interactions can
            have (4-4 clear, 4-2 clear, 2-2 clear, or totally clear)

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

    def clear(self, interaction_type: str, clear_from: str = "all") -> int:
        """
        Compute the number of interactions of type `interaction_type` that are
        clear from interactions of the `clear_from` type.

        Parameters
        ----------
        interaction_type : str
            Type of the interaction studied. Can either be "4-4", "4-2", "2-2",
            or "all" to count for all types of interaction.
        clear_from : str
            Type of clarity to count. Can either be "4-4", "4-2", "2-2", or
            "all" to count the number of totally clear interaction.

        """
        possible_types = ["4-4", "4-2", "2-2", "all"]
        if (
            interaction_type not in possible_types
            or clear_from not in possible_types
        ):
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
        Each word contains the generator + the letter representing the added
        factor.
        The pseudo-factors used to generate the four-level factors are replaced
        by four-level factors labels in the final words.
        For example, if factor `A` is created using pseudo-factors `a`, `b`
        and `ab`, then the replacement scheme is: `a -> A1`, `b -> A2`,
        and `ab -> A3`.

        Parameters
        ----------
        raw : bool, optional
            Return the relation without replacing the pseudo-factors by their
            four-level factor label, by default False

        Returns
        -------
        Relation
            Defining relation as a ``Relation`` object.
        """
        relation = [
            f"{num2gen(x)}{chr(97 + self.k + i)}"
            for i, x in enumerate(self.af)
        ]
        return Relation(words=relation, m=self.m)

    def add_factor(self, number: int):
        """
        Add a two-level factor to the design.

        Generate a design with an added two-level factor defined by the number
        supplied.

        Parameters
        ----------
        number : int
            Column number of the factor to add. Cannot correspond to a column
            number already in use (as two-level factor or pseudo-factor) in the
            design.

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
            raise ValueError(
                f"Column number {number} already used in the design"
            )
        new_cols = self.af + [number]
        return Design(self.runsize, self.m, new_cols)

    def remove_factor(self, number: int):
        """
        Remove a two-level factor from the design.

        Parameters
        ----------
        number : int
            Column number of the factor to remove. Must correspond to one of
            the factor used in the design.

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


def from_array(mat: np.ndarray, zero_coded: bool = True) -> Design:
    """
    Create a Design object from a matrix.

    Parameters
    ----------
    mat : np.ndarray
        Design matrix
    zero_coded : bool, optional
        The two-level factors in the design matrix are coded in 0/1, by
        default True.

    Returns
    -------
    Design
        A Design object with the factors detected in the matrix.
    """
    # Define runsize
    run_size, n_fac = mat.shape
    # Define number of four-level factors
    levels = [len(np.unique(mat[:, i])) for i in range(n_fac)]
    n_flvl = sum([i == 4 for i in levels])
    # Sort matrix and extract two-level part
    sorting_list = [mat[:, i] for i in range(n_flvl)]
    sorted_mat = mat[np.lexsort(sorting_list)]
    tlvl_mat = sorted_mat[:, [i == 2 for i in levels]]
    if not zero_coded:
        tlvl_mat = (tlvl_mat + 1) / 2
    # Compute all possible columns for that runsize
    k = int(np.log2(run_size))
    bf_matrix = basic_factor_matrix(k=k, zero_coding=True)
    full_generator_matrix = np.zeros((k, 2**k - 1), dtype=int)
    labels_list = []
    for i in range(2**k - 1):
        bin_string = f"{(i + 1):0{k}b}"[::-1]
        full_generator_matrix[:, i] = [int(j) for j in bin_string]
        gen = [chr(i + 97) for i, x in enumerate(bin_string) if int(x) == 1]
        labels_list.append("".join(gen))
    full_design_matrix = np.matmul(bf_matrix, full_generator_matrix) % 2
    cols_mat = np.matmul(full_design_matrix.T * 2 - 1, tlvl_mat * 2 - 1)
    cols_id, _ = np.nonzero(cols_mat)
    cols_gen = [labels_list[i] for i in cols_id]
    cols_num = [gen2num(i) for i in cols_gen]
    added_cols = [i for i in cols_num if np.log2(i) % 1 != 0]
    return Design(run_size, n_flvl, added_cols)


def fl_to_contrast(
    fl_matrix: np.ndarray,
    scaled_contrast_matrix: np.ndarray,
    permutations_list: "list[list[int]]",
) -> np.ndarray:
    """Decompose each four-level factor of `fl_matrix` into its three
    polynomial components: linear (L), quadratic (Q), and cubic (C) and apply a
    specific permutation to each factor.

    Parameters
    ----------
    fl_matrix : np.ndarray
        2D array of size N by m containing the m four-level factors with levels
        0, 1, 2, and 3.
    scaled_contrast_matrix : np.ndarray
        2D array of size 4 by 3 containing the scaled polynomial contrasts.
        The four rows correspond to the four levels of a four-level factor,
        while the three columns correspond to the 3 polynomial contrasts:
        linear (L), quadratic (Q), and cubic (C).
    permutations_list : list[list[int]]
        List of length m containing the permutations for the levels of each
        four-level factor.
        Each permutation must be a list of length four containing the levels 0,
        1, 2, and 3 in a specific order.

    Returns
    -------
    np.ndarray
        2D array of size N by 3*m containing the LQC-decomposed four-level
        factors with the permutations applied.

    Raises
    ------
    ValueError
        permutations_list must have the same length as the number of columns in
        the fl_matrix.
    :meta private:
    """
    # There must be one permutation for each 4LF
    N, m = fl_matrix.shape
    if len(permutations_list) != m:
        raise ValueError("permutations_list must have length %i" % m)

    # Matrix holding all contrasts for m 4LFs
    fl_contrast_matrix = np.zeros((N, 3 * m), dtype=float)

    # Apply each permutation to the contrast matrix sequentially
    for fl_index, permutation in enumerate(permutations_list):
        perm_contrast_matrix = scaled_contrast_matrix[permutation, :]
        # Iterate over the permutated contrast matrix
        for index, contrast_value in np.ndenumerate(perm_contrast_matrix):
            i, j = index  # Unfold index for easier array indexing
            # Fnd rows in 4lf equal to i
            row_indices = fl_matrix[:, fl_index] == i
            # Set these rows equal to the contrast value
            fl_contrast_matrix[row_indices, 3 * fl_index + j] = contrast_value

    return fl_contrast_matrix


def contrast_interactions(
    fl_contrast_matrix: np.ndarray, m: int
) -> "tuple[np.ndarray, list[int]]":
    """Compute all interactions of order 1 up to `m` for `m` four-level factors
    decomposed into their LQC contrasts.

    Parameters
    ----------
    fl_contrast_matrix : np.ndarray
        2D array of size N by 3*m containing the LQC-decomposition of the m
        four-level factors.
    m : int
        Number of four-level factors.

    Returns
    -------
    fl_full_contrast_matrix : np.ndarray
        2D array containing all interactions between the LQC contrasts of m
        four-level factors.
    fl_contrast_length :  list[int]
        List containing the order (L is 1, Q is 2, and C is 3) of all the
        interactions generated

    Raises
    ------
    ValueError
        The number of four-level factors `m` must be either 1, 2, or 3.
    """
    # Number of 4LF is either 1, 2, or 3
    if m not in [1, 2, 3]:
        raise ValueError("`m` must be 1, 2, or 3.")

    # Contrast length is 1, 2, 3 for L, Q, C
    fl_contrast_length = [i + 1 for _ in range(m) for i in range(3)]

    # Number of runs to define empty matrices
    N = fl_contrast_matrix.shape[0]

    # No interaction can be created with one four-level factor
    if m == 1:
        fl_full_contrast_matrix = fl_contrast_matrix

    # Main effects for the two four-level factors and the 3^2 interactions
    # between the 6 terms.
    elif m == 2:
        fl_contrast_interaction_matrix = np.zeros((N, 9), dtype=float)
        for fac_idx, pair in enumerate(product([0, 1, 2], [3, 4, 5])):
            fl_contrast_interaction_matrix[:, fac_idx] = np.prod(
                fl_contrast_matrix[:, pair], axis=1
            )
            normalized_indices = [i % 3 for i in pair]
            fl_contrast_length.append(sum(normalized_indices) + 2)
        fl_full_contrast_matrix: np.ndarray = np.hstack(
            (fl_contrast_matrix, fl_contrast_interaction_matrix)
        )

    # Three different terms:
    # - Main effects for the 3 four-level factors (3*3=9 terms)
    # - Interactions between 2 of the 3 four-level factors (3 ways of
    # creating 9)
    # - Interactions between the 3 four-level factors (27 terms)
    elif m == 3:
        # Interactions between 2 of the 3 four-level factors
        fl_contrast_2f_interaction_matrix = np.zeros((N, 27), dtype=float)
        for comb_idx, combi in enumerate([(0, 1), (0, 2), (1, 2)]):
            range1 = [3 * combi[0], 3 * combi[0] + 1, 3 * combi[0] + 2]
            range2 = [3 * combi[1], 3 * combi[1] + 1, 3 * combi[1] + 2]
            for pair_idx, pair in enumerate(product(range1, range2)):
                index = 9 * comb_idx + pair_idx
                fl_contrast_2f_interaction_matrix[:, index] = np.prod(
                    fl_contrast_matrix[:, pair], axis=1
                )
                normalized_indices = [i % 3 for i in pair]
                fl_contrast_length.append(sum(normalized_indices) + 2)

        # Interactions between 3 four-level factors
        fl_contrast_3f_interaction_matrix = np.zeros((N, 27), dtype=float)
        for triplet_idx, triplet in enumerate(
            product([0, 1, 2], [3, 4, 5], [6, 7, 8])
        ):
            fl_contrast_3f_interaction_matrix[:, triplet_idx] = np.prod(
                fl_contrast_matrix[:, triplet], axis=1
            )
            normalized_indices = [i % 3 for i in triplet]
            fl_contrast_length.append(sum(normalized_indices) + 3)

        # Join ME, 2FI, and 3FI into a single matrix
        fl_full_contrast_matrix: np.ndarray = np.hstack(
            (
                fl_contrast_matrix,
                fl_contrast_2f_interaction_matrix,
                fl_contrast_3f_interaction_matrix,
            )
        )
    return fl_full_contrast_matrix, fl_contrast_length  # type: ignore


def compute_correlations(
    tl_interaction_matrix: np.ndarray,
    tl_interaction_length: "list[int]",
    fl_matrix: np.ndarray,
    permutations_list: "list[list[int]]",
    max_length: int | None = None,
) -> "list[float]":
    """Compute the squared correlations between the two-level
    factors, and the LQC decomposition of the four-level factors.
    Summarise the correlations by summing them by interaction order.

    For example, an interaction between a xyz interaction (order 3) and a Q
    component (order 2), has order 5.

    Parameters
    ----------
    tl_interaction_matrix : np.ndarray
        2D array containing all interactions between two-level factors.
    tl_interaction_length : list[int]
        A list containing the order of all interactions between two-level
        factors. It must have the same size as the number of columns in
        `tl_interaction_matrix`.
    fl_matrix : np.ndarray
        2D array of size N by m containing the four-level factors.
    permutations_list : list[list[int]]
        List of length m containing the permutations for the levels of each
        four-level factor.
        Each permutation must be a list of length four containing the levels 0,
        1, 2, and 3 in a specific order.
    max_length : int | None, optional
        Maximum order to consider when summing the correlations. Default is
        None so all orders are considered.

    Returns
    -------
    beta_vector : list[float]
        A list of the summed squared correlations, starting with order 3 and
        going up to order `max_length`.
    """
    # Declare global variables
    global scaled_contrast_matrix

    # Compute runsize and number of four-level factors
    N, m = fl_matrix.shape

    # Unfold and permute levels of 4LFs
    fl_contrast_matrix = fl_to_contrast(
        fl_matrix=fl_matrix,
        scaled_contrast_matrix=scaled_contrast_matrix,
        permutations_list=permutations_list,
    )

    # Length of the contrasts unfolded from the 4LF
    fl_full_contrast_matrix, fl_contrast_length = contrast_interactions(
        fl_contrast_matrix=fl_contrast_matrix, m=m
    )

    # Build word length vector
    word_length = (
        np.array(fl_contrast_length)[:, None].T
        + np.array(tl_interaction_length)[:, None]
    )

    # Compute squared correlations
    correlation = (tl_interaction_matrix.T @ fl_full_contrast_matrix) / N
    correlation_sq = np.round(correlation**2, 2)

    # Compute Beta values
    max_possible_length = max(fl_contrast_length) + max(tl_interaction_length)
    if max_length is None:
        max_Aq_length = max_possible_length
    else:
        max_Aq_length = max_length + 1
    beta_vector = []
    for length in range(3, max_Aq_length):
        beta_value = np.round(np.sum(correlation_sq[word_length == length]), 2)
        beta_vector.append(beta_value)

    return beta_vector


def nbr_interactions(n: int, max_length: int) -> int:
    """
    Compute the total number of interactions between i factors among n, where i
    ranges between 1 and max_length.
    """
    val = 0
    for i in range(max_length):
        val += comb(n, (i + 1))
    return val


def build_tfi_model_matrix(n: int, max_length: int) -> np.ndarray:
    """
    Build the model matrix to compute all interactions between i factors among
    n with i ranging from 1 to n.
    In this model matrix, a 1 represents an active factor in the interaction,
    while a 0 zero represents an inactive factor.

    Parameters
    ----------
    n : int
        Number of two-level factors
    max_length : int
        Maximum number of factors to consider in the interactions.
        For example, for a two-factor interaction, max_length should be equal
        to 2.

    Returns
    -------
    np.ndarray
        Model matrix of size n-by-s, where s is the total number of
        interactions.
    """
    # Initialize matrix
    tfi_model_matrix = np.zeros(
        (n, nbr_interactions(n, max_length)), dtype=int
    )
    # Ones represent an active factor in the interaction
    index = 0
    for i in range(max_length):
        for c in combinations(range(n), i + 1):
            tfi_model_matrix[c, index] = 1
            index += 1
    return tfi_model_matrix


def build_w2_vector(
    blocking_wlp: list[int], treatment_wlp: list[int]
) -> list[int] | None:
    """
    Compute the W_2 word length pattern (WLP) of a design, given its treatment
    WLP and blocking WLP.
    Both WLP must start with words of length 0 (so that list indices matches
    word lengths).
    """
    if len(blocking_wlp) != len(treatment_wlp):
        raise ValueError("Both WLP must have the same length")
    if len(blocking_wlp) <= 3:
        if any([i != 0 for i in blocking_wlp]) or any(
            [i != 0 for i in treatment_wlp]
        ):
            warnings.warn(
                "A vector contains a non-zero A_i value for i < 3",
                UserWarning,
            )
            return None
        else:
            warnings.warn(
                "WLPs too short to combine. Returning a zero vector",
                UserWarning,
            )
            return [0, 0, 0]
    w2_vector = [treatment_wlp[3], blocking_wlp[3]]
    n_words = len(blocking_wlp) - 1
    max_i = 1 + (n_words - 1) // 2
    for i in range(2, max_i):
        w2_vector.append(treatment_wlp[2 * i])
        w2_vector.append(treatment_wlp[2 * i + 1])
        w2_vector.append(blocking_wlp[i + 2])
    return w2_vector
