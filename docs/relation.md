# Defining relation

The aliasing of all the factors of a design can be summarized by a **defining relation**.
The defining relation of a design contains all the words used to generate the added factors of the design.
Therefore, it uniquely defines any regular design.

For example, consider the 32-run design with one four-level factor $A$, and two added two-level factors 27 and 30.
Since the design involves 32 runs, there are 5 basic factors ($a$ to $e$), so that the added factors are labeled $f$ and $g$.
The column numbers, 27 and 30, correspond to the generators $abde$ and $bcde$, respectively.
After relabeling the pseudo-factors, the generators become $A_3de$ and $A_2cde$.
Therefore, the defining relation of the design can be written as
$$
(A_3def, \; A_2cdeg)
$$

In the ``fatld`` package, we can create the defining relation of the design using the {class}`fatld.relation.Relation` class that create a ``Relation`` object, which holds a defining relation.

```python
>>> from fatld.relation import Relation, num2gen
>>> added_factors = [27, 30]
>>> letters = ['f', 'g']
>>> words = [num2gen(x) + letters[i] for i,x in enumerate(added_factors)]
>>> words
['abdef', 'bcdeg']
>>> # We specify `m` so that the words can be relabeled
>>> r = Relation(words, m=1)
>>> r
['A3def', 'A2cdeg']
```

````{tip}
If you are using a ``Design`` object, created using {class}`fatld.design.Design`, you can obtain the defining relation of the design using the {meth}`fatld.design.Design.defining_relation` method, instead of creating it by hand.
````

## Expanding

Two or more words multiplied together can generate another word.
This means that a defining relation with $k$ words can be expanded to $2^k-1$ words in total, by multliplying all the words together.
All the aliasing between the factors of a design can be infered from such an expanded relation.

```python
>>> full_relation = r.expand(relabel=True)
>>> full_relation
['A3def', 'A2cdeg', 'A1cfg']
```

When the full expanded relation is too long you can summarise it using the word length pattern.
This pattern counts the number of words of each type and each length.

```{code} python
>>> [word_length(w) for w in full_relation]
[4, 5, 4]
>>> [word_type(w) for w in full_relation]
[1, 1, 1]
>>> r.word_length_pattern()
[[0, 0], [0, 2], [0, 1]]
```
