from __future__ import annotations
import warnings
import numpy as np
import oapackage as oa  # type: ignore

from itertools import product, permutations, combinations
from math import comb

# Setup the contrast matrix for the beta aberration
contrast_matrix = np.array(
    [[-3.0, 1.0, -1.0], [-1.0, -1.0, 3.0], [1.0, -1.0, -3.0], [3.0, 1.0, 1.0]]
)
scaled_contrast_matrix = np.divide(
    np.multiply(2.0, contrast_matrix), np.linalg.norm(contrast_matrix, axis=0)
)


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


def power2_decomposition(n: int, length: int | None = None) -> list[int]:
    """
    Decompose a number into powers of 2 and returns the corresponding indices.

    Parameters
    ----------
    n : int
        Number to decompose
    length : int, optional
        Wanted length for the list of indices, by default None, returns the
        list when the last power of two `x` for which `x` < `n` is reached

    Returns
    -------
    List[int]
        List of the indices of the powers of two obtained in the decomposition

    Example
    -------
    >>> # 7 can be decomposed into 1*2 + 1*2 + 1*4
    >>> power2_decomposition(7)
    [1, 1, 1]
    >>> # 11 can be decomposed into 1*1 + 1*2 + 0*4 + 1*8, we want a list of
    >>> # length 7
    >>> power2_decomposition(11, length=7)
    [1, 1, 0, 1, 0, 0, 0]
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


def custom_design(runsize: int, column_list: list[int]) -> np.ndarray:
    """
    Create a custom design with a specific run size, based on the column
    numbers provided.

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
    bf_mat = basic_factor_matrix(k=k, zero_coding=True)
    power_list = [power2_decomposition(n=col, length=k) for col in column_list]
    yates_matrix = np.array(power_list).T
    custom_mat = (bf_mat @ yates_matrix) % 2
    return custom_mat


def twlp(
    ar: oa.array_link, type_0: bool = True, max_length: int | None = None
) -> list[list[int]]:
    """
    Compute the type-specific word length pattern of a design, starting with
    words of length 3.

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
        If this number is more than the maximum word length of the design,
        ignores it.

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
            comb = [
                (i, idx - i)
                for i in range(idx + 1)
                if ((i < r) and (idx - i) < c)
            ]
        for coord in comb:
            wlp_matrix[coord[0], idx] = dst_dist_transformed_mat[
                coord[0], coord[1]
            ]
    # Trim the matrix to max_length if needed
    if max_length is not None and (
        max_length < 3 or max_length > max_word_length
    ):
        max_length = None
        warnings.warn(
            "Incorrect max_length. Full type-specific WLP is outputted"
        )
    if max_length is None:
        wlp_matrix_trimmed = wlp_matrix[:, 3:]
    else:
        end = max_length + 1
        wlp_matrix_trimmed = wlp_matrix[:, 3:end]
    # Create a list of lists
    wlp_type_list = wlp_matrix_trimmed.T.tolist()
    # Reverse the sublists if it's not type 0 because they are ordered in type
    # 0 by default
    if type_0 is False:
        return [x[::-1] for x in wlp_type_list]
    else:
        return wlp_type_list


def read_OA(filename: str) -> np.ndarray:
    """
    Read a `.oa` file, generated by the `OApackage` package and return all the
    matrices it contains as a 3-D numpy array, where the third dimension is the
    index of the designs.
    """
    with open(filename, "r") as f:
        info = f.readline()
        n_fac, run_size, n_des = [int(i) for i in info.strip().split(" ")]
        designs = np.zeros((run_size, n_fac, n_des), dtype=int)
        row_index = 0
        design_index = 0
        for line in f:
            # Matrices lines have more than three characters
            if len(line) < 3:
                design_index = int(line.strip())
                row_index = 0
            else:
                row_vector = [int(i) for i in line.strip().split(" ")]
                if row_vector == [-1]:
                    break
                designs[row_index, :, design_index] = row_vector
                row_index += 1
    return designs


def alpha_wlp(ar: oa.array_link) -> list[float]:
    """Compute the α word length pattern of an array.
    The alpha wlp contains 5 values ordered in the following way: ω4, ω2, ω42,
    ω22, ω44

    Parameters
    ----------
    ar : oa.array_link
        Matrix of the FATL design

    Returns
    -------
    list[float]
        The alpha word length pattern
    """
    twlp_vector = twlp(ar, max_length=4)
    # For cases with designs with not enough words
    if len(twlp_vector) == 0:
        twlp_vector = [[0, 0, 0], [0, 0, 0, 0]]
    elif len(twlp_vector) == 1:
        twlp_vector.append([0, 0, 0, 0])
    a3_vector = twlp_vector[0] + [0] * (
        3 - len(twlp_vector[0])
    )  # Right pad with 0
    a4_vector = twlp_vector[1] + [0] * (4 - len(twlp_vector[1]))
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
    alpha_wlp = [omega_values[i] for i in ["w4", "w2", "w42", "w22", "w44"]]
    return [round(i, 2) for i in alpha_wlp]


def beta_star_wlp(
    ar: oa.array_link, max_length: int | None = None
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

    # Compute the number of four-level factors
    mat = ar.getarray()
    levels = [len(np.unique(mat[:, i])) for i in range(mat.shape[1])]
    m = sum([i == 4 for i in levels])

    # Compute all possible permutations for m factors
    all_perms = []
    for _, p in enumerate(product(range(12), repeat=m)):
        all_perms.append([perms[i] for i in p])

    # Isolate four-level and two-level part of the design
    fl_matrix = mat[:, :m]
    tl_matrix = mat[:, m:]  # noqa: E201

    # Build all interactions between all two-level factors
    n = mat.shape[1] - m
    if max_length is None or max_length > n:
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
    twlp_vector = twlp(ar)
    for i, x in enumerate(best_vector):
        if i >= len(twlp_vector):
            beta_star_vector.append([x, 0])
        else:
            beta_star_vector.append([x, twlp_vector[i][0]])

    return beta_star_vector, best_perm


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
