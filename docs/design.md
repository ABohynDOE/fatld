# Four-and-two-level designs

When constructing a design with four-level factors and two-level factors, the four-level factors are always constructed from pairs of basic factors, (see section on {ref}`four-level-factors`).
The two-level factors in the design include the remaining basic factors (if any) and the added factors.

Consider the example of a 64-run design with 2 four-level factors and 7 two-level factors.
There are $log_2(64)=6$ basic factors, among which $2*m=4$ are used to generate four-level factors.
Therefore two basic two-level factors are left and only five added two-level factors still need to be defined.
In this example, we use the factors 5, 13, 17, 27.
The figure below details this process.

```{image} _static/design_generation.png
:align: center
:width: 300
```

To create the design using the package, we use the {func}`fatld.design.Design` class.
To create a design, we only need to specify:

- the number of runs as ``runsize``
- the number of four-level factors as ``m``
- the column numbers of the added factors as ``cols``

The basic factors that weren't used to generate four-levels factors will automatically be added to the columns of the design.
The full list of columns can be obtained using the ``Design.cols`` attribute.

```{warning}
Columns are always sorted by numbers so basic factors will not always be first in the set of columns of the design.
```

```{code-block} python
>>> from fatld.design import Design
>>> D = Design(runsize=64, m=2, cols=[5, 13, 17, 27])
>>> D.cols
[5, 13, 16, 17, 27, 32]
```

## Characteristics

Several information about a design are available as attributes of the ``Design`` object.
The list includes:

- ``runsize``: the runsize
- ``m``: the number of *four-level* factors
- ``n``: the number of *two-level* factors
- ``k``: the number of *basic* two-level factors
- ``p``: the number of *added* two-level factors

The full list of available attributes is available in the [documentation](design-documentation)

```{code-block} python
>>> from fatld.design import Design
>>> D = Design(runsize=64, m=2, cols=[5, 13, 17, 27])
>>> D.runsize
64
>>> D.n
7
>>> D.p
5
>>> D.n == (D.p + D.k - 2 * D.m)
True
```

## Design matrix

bla bla
