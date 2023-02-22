# Four-and-two-level designs

When constructing a design with four-level factors and two-level factors, the four-level factors are always constructed from pairs of basic factors, (see section on {ref}`four-level-factors`).
The two-level factors in the design include the remaining basic factors (if any) and the added factors.

Consider the example of a 64-run design with 2 four-level factors and 7 two-level factors.
There are $log_2(64)=6$ basic factors, among which $2*m=4$ are used to generate four-level factors.
Therefore 2 basic two-level factors are left and 5 added two-level factors need to be generated.
In this example, we add the factors 5, 13, 17, 27.
The figure below details this process.

```{image} _static/design_generation.png
:align: center
:width: 300
```

To create the design using the package, we use the {func}`fatld.design.Design` class.
