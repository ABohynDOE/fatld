# Defining relation

The aliasing pattern of a design can be summarized by a **defining relation** that contains all the words associated with the added factors of the design.
Defining relations have two interesting features:

- **Expansion**: expand the relation containing $p$ words to the full defining relation containing $2^{p}-1$ words.
- **Word length pattern**: compute the word length pattern of the design based on the full defining relation.

For example, consider the 64-run design with two four-level factors, $A$ and $B$, and three added two-level factors: 21, 42, and 53.
Since the design involves 64 runs, there are 6 basic factors ($a$ to $f$), so that the four added factors are labeled $g$, $h$, and $i$.
Using the four generators and the four labels, we can create the defining relation of the design using the {class}`fatld.relation.Relation` class that create a ``Relation`` object, which holds a defining relation.

```python
>>> from fatld.relation import Relation, num2gen
>>> added_factors = [21, 42, 53]
>>> letters = [chr(97+6+i) for i in range(4)]
>>> subgroup = [f"{num2gen(x)}{letters[i]}" for i,x in enumerate(added_factors)]
>>> subgroup
['aceg', 'bdfh', 'acefi']
>>> r = Relation(subgroup, m=2)
>>> r
['A1B1eg', 'A2B2fh', 'A1B1efi']
```

This class has two methods {meth}`fatld.relation.Relation.expand` and {meth}`fatld.relation.Relation.word_length_pattern` that allow you to expand the defining relation, and to compute its word length pattern, respectively.

```python
>>> full_relation = r.expand(relabel=True)
>>> full_relation
['A1B1eg', 'A2B2fh', 'A1B1efi', 'A3B3efgh', 'fgi', 'A3B3ehi', 'A2B2ghi']
>>> [word_length(w) for w in full_relation]
[4, 4, 5, 6, 3, 5, 5]
>>> [word_type(w) for w in full_relation]
[2, 2, 2, 2, 0, 2, 2]
>>> r.word_length_pattern()
[[1, 0, 0], [0, 0, 2], [0, 0, 3], [0, 0, 1]]
```
