# Defining relation

The aliasing pattern of a design can be summarized by a **defining relation** that contains all the words associated with the added factors of the design.
Defining relations have two interesting features:

- **Expansion**: expand the relation containing `p` words to the full defining relation containing $2^{p}-1$ words.
- **Word length pattern**: compute the word length pattern of the design based on the full defining relation.

For example, consider the 128-run design with two four-level factors, $A$ and $B$, and four added two-level factors: 21, 42, 103, and 121.
Since the design involves 128 runs, there are 7 basic factors (`a` to `g`), so that the four added factors are labeled `h`, `i`, `j`, and `k`.
Using the four generators and the four labels, we can create the defining relation of the design using the {class}`fatld.relation.Relation` class that create a `Relation` object, which holds a defining relation.

```python
>>> from fatld.relation import Relation
>>> added_factors = [21, 42, 103, 121]
>>> letters = [chr(97+7+i) for i in range(4)]
>>> subgroup = [f"{num2gen(x)}{letters[i]}" for i,x in enumerate(added_factors)]
>>> subgroup
['aceh', 'bdfi', 'abcfgj', 'adefgk']
>>> r = Relation(subgroup, m=2)
>>> r
['A1B1eh', 'A2B2fi', 'A3B1fgj', 'A1B2efgk']
```

This class has two methods {meth}`fatld.relation.Relation.expand` and {meth}`fatld.relation.Relation.word_length_pattern` that allow you to expand the defining relation, and to compute its word length pattern, respectively.

```python
>>> full_relation = r.expand(relabel=True)
>>> full_relation
['A1B1eh', 'A2B2fi', 'A3B1fgj', 'A1B2efgk', 'A3B3efhi', 'A2efghj', 'B3fghk', 'A1B3gij', 'A3egik', 'A2B3ejk', 'B2eghij', 'A2B1ghik', 'A3B2hjk', 'B1efijk', 'A1fhijk']
>>> [word_length(w) for w in full_relation]
[4, 4, 5, 6, 6, 6, 5, 5, 5, 5, 6, 6, 5, 6, 6]
>>> [word_type(w) for w in full_relation]
[2, 2, 2, 2, 2, 1, 1, 2, 1, 2, 1, 2, 2, 1, 1]
>>> r.word_length_pattern()
[[0, 0, 0], [0, 0, 2], [0, 2, 4], [0, 4, 3]]
```
