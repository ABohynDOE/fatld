# Word length pattern

## Generalized

The word length pattern (WLP) is a summary of the length of all the words in the defining relation of a design.
It is a good indicator of the quality of a design.
Indeed, the length of a word gives information about the aliasing between factors:

- length 3: aliasing between a main effect (ME) and a two-factor interaction (TFI)
- length 4: aliasing between two TFI

Even though the word length pattern can be computed from a defining relation, as seen in section {ref}`word-length`, it is easier to compute it from a ``Design`` object using the {meth}`fatld.design.Design.wlp` method.
The WLP always goes from the words of length 3 up to words of length $n+m$, but it can truncated using the `max_length` keyword.

```{code-block} python
>>> import fatld
>>> D = fatld.Design(runsize=32, m=1, cols=[21, 27, 29])
>>> D.wlp()
[1, 3, 3, 0, 0]
>>> D.wlp(max_length=5)
[1, 3, 3]
```

## Type-specific

As explained in section {ref}`word-type`, words can have different types, so it is interesting to differentiate words of the same length by their type.
Indeed, the aliasing between a ME and a TFI is not as important if the TFI is between two two-level factors (type 0) or between a two-level factor and a pseudo-factor (type 1).

The type-specific word length pattern (tWLP) is a good summary of the lengths and types of all the words in the defining relation of a design.

The tWLP can easily be obtained from a ``Design`` object using the {meth}`fatld.design.Design.twlp` method.
It starts with words of length 3 and ordered from type 0 to type $m$.
You can obtain the type-$m$ WLP using the `type_0` keyword.
When set to `False`, the function will sort the words of the same length by descending type.

```{code-block} python
>>> import fatld
>>> D = fatld.Design(runsize=32, m=1, cols=[21, 27, 29])
>>> # Words of type 0 are given first, then words of type 1
>>> D.twlp()
[[1, 0], [0, 3], [0, 3], [0, 0], [0, 0]]
>>> # This behavior can be reverted
>>> D.wlp(max_length=5, type_0=False)
[[0, 1], [3, 0], [3, 0]]
```
