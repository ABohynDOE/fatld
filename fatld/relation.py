import re
from collections import Counter
from itertools import combinations, chain
from typing import List, Optional

from .main import power2_decomposition


def num2gen(n: int, m: Optional[int] = None) -> str:
    """
    Return the generator corresponding to a given column number

    If a value is provided for the number of four-level factors, the
    pseudo-factors used to generate the four-level factors will be replaced by
    their corresponding labels in the generator.
    The labels are:

    - `A1=a`, `A2=b`, `A3=ab` for the first factor
    - `B1=c`, `B2=d`, `B3=cd` for the second factor
    - `C1=e`, `C2=f`, `C3=ef` for the third factor

    Parameters
    ----------
    n : int
        Column number

    m : Optional[int]
        Number of four-level factors. If no value is provided, the generator
        will be considered to come from a two-level design. By default, no
        value is provided.


    Examples
    ------
    >>> num2gen(7)
    'abc'
    >>> num2gen(7, m=1)
    'A1c'

    See Also
    --------
    gen2num : Convert a generator into the corresponding column number
    """
    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if m is not None and m not in [1, 2, 3]:
        raise ValueError("Only accepted values for m are 1, 2, 3")
    powers = power2_decomposition(n)
    word = "".join([chr(97 + i) for i, x in enumerate(powers) if x == 1])
    if m is None:
        return word
    else:
        for i in range(m):
            pseudo_factors = [
                chr(97 + 2 * i),
                chr(97 + 2 * i + 1),
                f"{chr(97+2*i)}{chr(97+2*i+1)}",
            ]
            pf_labels = [f"{chr(65+i)}{x}" for x in [1, 2, 3]]
            word = (
                word.replace(pseudo_factors[2], pf_labels[2])
                .replace(pseudo_factors[1], pf_labels[1])
                .replace(pseudo_factors[0], pf_labels[0])
            )
        return word


def gen2num(g: str) -> int:
    """
    Convert a generator into its corresponding column number

    Generators can only contain lowercase letters from a to z.

    Parameters
    ----------
    g : str
        Generator

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


def relabel_word(word: str, m: int) -> str:
    """
    Relabel basic factors in a word to their corresponding pseudo-factors.

    Relabel a word containing only two-level factors by replacing the basic
    factors used as pseudo-factors by their corresponding pseudo-factor labels.
    The pseudo-factor labels are A_i, B_i and C_i with i = 1,2,3, for the
    first, second, and third four-level factor, respectively.

    Parameters
    ----------
    word : str
        Word to relabel
    m : int
        Number of four-level factors. Define up to which factor the relabeling
        occurs.

    Examples
    --------
    >>> # The value of `m` define how many factors are relabeled
    >>> relabel_word(word='abcde', m=1)
    'A3cde'
    >>> relabel_word(word='abcde', m=3)
    'A3B3C1'

    """
    for i in range(m):
        pf_factors = [
            chr(97 + 2 * i),
            chr(97 + 2 * i + 1),
            f"{chr(97 + 2 * i)}{chr(97 + 2 * i + 1)}",
        ]
        pf_labels = [f"{chr(65 + i)}{x}" for x in [1, 2, 3]]
        # We start with p1p2 to avoid replacing p1/p2 first and not the
        # interaction
        s = (
            word.replace(pf_factors[2], pf_labels[2])
            .replace(pf_factors[1], pf_labels[1])
            .replace(pf_factors[0], pf_labels[0])
        )
        word = s
    return word


def word_length(word: str) -> int:
    """
    Compute the length of a word that can contain pseudo-factors.

    Examples
    --------
    >>> word_length('abcd')
    4
    >>> word_length('A1cdC1')
    4
    """
    return len(re.findall(r"([A-Ca-z]\d?)", word))


def word_type(word: str) -> int:
    """
    Compute the type of a word that can contain pseudo-factors.

    Examples
    --------
    >>> word_type('abcde')
    0
    >>> word_type('A1cde')
    1
    """
    return len(set(re.findall("[A-C]", word)))


class Relation:
    # TODO: document the relation class
    def __init__(self, words: List[str], m: int = 0):
        # `m` can only be 0, 1,2 or 3 since there can't be more four-level
        # factors, 0 means there are no four-level factors
        if m not in [0, 1, 2, 3]:
            raise ValueError("m can only take value of 0, 1, 2, and 3")
        # Words must be strings and only contain basic factors
        if not isinstance(words, List) or any(
            [not isinstance(w, str) for w in words]
        ):
            raise TypeError("Words must be given as a list of strings")
        if any([(ord(i) > 122 or ord(i) < 97) for w in words for i in w]):
            raise ValueError(
                "Words can only contain lowercase letters (representing "
                "basic factors)"
            )
        self.m = m
        self.words = words
        self.relabel_words = [
            relabel_word(word=w, m=self.m) for w in self.words
        ]

    def __repr__(self):
        return f"{self.relabel_words}"

    def expand(self, relabel: bool = False) -> List[str]:
        """
        Expand a defining subgroup into the complete defining relation.

        For a defining subgroup with `p` members, expand it into all `2^p -1`
        interactions between the members.

        Parameters
        ----------
        relabel: bool
            Relabel the pseudo-factors in the words of the defining relation.
            Default is False

        Examples
        --------
        >>> r = Relation(words=['abcdf', 'aceg'], m=1)
        >>> r.expand()
        ['abcdf', 'aceg', 'bdefg']
        >>> # You can also relabel the pseudo-factors in the words
        >>> r.expand(relabel=True)
        ['A3cdf', 'A1ceg', 'A2defg']
        """
        full_relation = [w for w in self.words]
        for i in range(1, len(self.words)):
            word_combination = combinations(self.words, i + 1)
            for comb in word_combination:
                word_freq = Counter("".join(comb))
                remaining_factors = [
                    factor
                    for factor, freq in word_freq.items()
                    if freq % 2 == 1
                ]
                remaining_factors.sort()
                new_word = "".join(remaining_factors)
                full_relation.append(new_word)
        if relabel:
            return [
                relabel_word(word=word, m=self.m) for word in full_relation
            ]
        else:
            return full_relation

    def word_length_pattern(self) -> List[List[int]]:
        """
        Compute the type-specific word length pattern of the defining relation.

        The type-specific word length pattern is a list of sublists, where the
        sublist at index `i` is the number of word of length `3+i`, ordered by
        increasing type.

        Examples
        -------
        >>> # Create a relation with two words of type 1 and length 4
        >>> r = Relation(words=['abcdf', 'aceg'], m=1)
        >>> # The third word is bdefg and has length 5
        >>> r.word_length_pattern()
        [[0, 0], [0, 2], [0, 1]]
        >>> # If m=0, then output a normal word length pattern
        >>> Relation(words=['abcdf', 'aceg']).word_length_pattern()
        [0, 1, 2]
        """
        full_relation = self.expand(relabel=True)
        lengths = [word_length(w) for w in full_relation]
        types = [word_type(w) for w in full_relation]
        # Each sublist must be initiated as an individual object (can't do
        # shallow copy)
        wlp = [[0 for _ in range(self.m + 1)] for _ in range(max(lengths) + 1)]
        for i, x in enumerate(lengths):
            t = types[i]
            wlp[x][t] += 1
        if self.m == 0:
            return list(chain(*wlp[3:]))  # type: ignore
        return wlp[3:]
