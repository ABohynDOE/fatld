# Factors

## Two-level factors

### Basic factors

A two-level factor with $N$ runs has $k=log_2(N)$ **basic factors**.
These factors specify all the treatment combinations of a $2^k$ full factorial design.

### Added factors

A design with $2^k$ runs can only have $k$ basic factors.
Assuming that the factor levels are -1 and +1, additional factors are generated as products of these basic factors.
These factors are called **added factors**.
If we represent the basic factors with individual lowercase letters starting with $a$, then a product of basic factors can be written by combining the letters of the basic factors.
Such products are called **generators**.

For example, a 32-run design has 5 basic factors.
An added factor created as the product of the first, third and fourth basic factors would have $acd$ as generator.

### Generators as numbers

When a large number of basic factors is involved, the letter notation of the generators can be cumbersome.
The solution is to store them as numbers instead of words.
Since every number can be uniquely decomposed into powers of 2, each generator can be associated with a combination of powers of 2.
To do this, each basic factor is uniquely assigned to a specific power of 2.

For example, $a$ is represented by 1, $b$ by 2, and $c$ by 4.
The sum of these powers creates the number defining the generator $acd$.
The figure below show that process for the generator $acd$ .

```{image} _static/column_numbering.png
:align: center
:width: 300
```

In the `fatld` package, the function {func}`fatld.relation.gen2num` can be used to generate the number corresponding to a generator.
Inversely, the function {func}`fatld.relation.num2gen` can be used to generate the generator corresponding to a column number.

```python
>>> from fatld.relation import num2gen, gen2num
>>> column = 23
>>> num2gen(column)
abce
>>> gen2num('abce')
23
```

````{tip}
To convert several numbers to their corresponding word, you can use list comprehension:

```{code-block} python
>>> [num2gen(n) for n in [21, 23, 17]]
['ace', 'abce', 'ae']
```
````

(four-level-factors)=

## Four-level factors

Four-level factors are are represented by uppercase letters, and are constructed from pairs of two-level factors using the grouping scheme of {cite:t}`wu1993minimum`.
In this scheme, the four levels of a factor are created from the combination of the levels of two two-level factors.
The two two-level factors and their interaction are called **pseudo-factors** and are represented as $A_1$, $A_2$, and $A_3$.

In this package, four-level factors are always created from pairs of basic factors.
Therefore, the pseudo-factors used to generate the four-level factors are always the same:

- $(a, b, ab)$, or column numbers (1, 2, 3), for factor $A$
- $(c, d, cd)$, or column numbers (4, 8, 12), for factor $B$
- $(e, f, ef)$, or column numbers (16, 32, 48), for factor $C$

Since some two-level factors are used as pseudo-factors, some generators need to be relabeled to represent that.
For example, in a design with one four-level factor defined as $\mathbf{A} =(A_1=a, A_2=b, A_3=ab)$, the generator $abce$ would be relabeled to $A_3ce$.
The function {func}`fatld.relation.relabel_word` can be used to relabel a generator, given the number of four-level factors that are considered.

```python
>>> from fatld.relation import relabel_word
>>> new_word = relabel_word(word='abce', m=1)
>>> new_word
A3ce
```
