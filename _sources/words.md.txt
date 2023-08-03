# Words

Added factors are created as products of basic factors represented by generators, but they can be also be represented by lowercase letters.
The combination of the generator used to create a added factor and the letter representing it is called a **word**.

For example, in a 32-run design, there are $log_2(32)=5$ basic factors $(a,b,c,d,e)$.
Let us consider the first added factor, represented by the letter $f$.
If it is created using with generator $abcde$, then, the word representing that factor is $abcdef$.
This word has length 6 and type 0 since it contains no pseudo-factors.

(word-length)=

## Length

The length of a word is the number of factors that composes it.
A word of length means that a main effect is aliased with a two-factor interaction, while a word of length four means that two two-factor interactions are aliased together.

The length of any word can be computed using the function {func}`fatld.relation.word_length`.
This function works for words containing both two-level factors and pseudo-factors.

```python
>>> from fatld.relation import word_length
>>> word_length('abcdef')
6
>>> word_length('A1cdg') # A1 counts a a single pseudo-factor that has length one !
4
```

```{important}
Do not forget to add the letter corresponding to the added factor at the end of the generator, else the length will be off be one !
```

(word-type)=

## Type

The type of a word is the number of pseudo-factors from *different* four-level factors that it contains.
The pseudo-factor of a four-level factor can always be presented as the product of the two other pseudo-factors.
Therefore, a word cannot contain two pseudo-factors from the same four-level factor.

For example, factor $\bf{A}$ is created using the three pseudo-factors $(A_1=a, A_2=b, A_3=ab)$, so that $A_3=A_1A_2$.
Therefore, if a word contains both $A_1$ and $A_2$, it can simply be replaced by $A_3$.

The type of a word can be computed using the function {func}`fatld.relation.word_type`.
This function works for words containing both two-level factors and pseudo-factors.

```python
>>> from fatld.relation import word_type
>>> word_type('abcef')
0
>>> word_type('A1cdg')
1
>>> # acdeg relabeled (with 2 four-level factors) is A1C3eg
>>> word_type(relabel_word('acdeg', m=2))
2
```

```{warning}
Remember that pseudo-factors are always represented as *uppercase* letters with an index.
Using another notation will not yield a correct result.
```
