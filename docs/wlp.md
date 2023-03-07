# Word length pattern

The word length pattern (WLP) can only be obtained from a ``Design`` object.
The generalized WLP is obtained using the {meth}`fatld.design.Design.wlp` and it starts with words of length 3.
You can truncate it using the `max_length` keyword.

The type-specific WLP is obtained using the {meth}`fatld.design.Design.twlp` and it also starts with words of length 3 and ordered from type 0 to type $m$.
You can obtain the type-$m$ WLP using the `type_zero` keyword.
