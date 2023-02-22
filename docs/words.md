# Words

Added factors are also represented by lowercase letters.
The combination of the generator used to create a added factor and the letter representing it is called a **word**.
Words have two important properties:

- **length**: the number of factors present in the word
- **type**: the number of pseudo-factors, coming from *different* four-level factors, present in the word

For example, in a 32-run design, there are $log_2(32)=5$ basic factors $(a,b,c,d,e)$.
If we consider the added factor $f$, with generator $abcde$, then the word representing that factor is $abcdef$.
This word has length 6 and type 0 since it contains no pseudo-factors.

## Length

The length of a word can be computed using the function {func}`fatld.relation.word_length`.
This function works for words containing only two-level factors and words containing pseudo-factors.

```python
>>> from fatld.relation import word_length
>>> word_length('abcdef')
6
>>> word_length('A1cdg')
4
```

```{tip}
Do not forget to add the letter corresponding to the added factor at the end of the generator, else the length will be off be one !
```

## Type

The type of a word can be computed using the function {func}`fatld.relation.word_type`.
This function works for words containing only two-level factors and words containing pseudo-factors.

```python
>>> from fatld.relation import word_type
>>> word_type('abcef')
0
>>> word_type('A1cdg')
1
>>> word_type(relabel_word('acdeg', m=2)) # relabeled word is A1C3eg
2
```
