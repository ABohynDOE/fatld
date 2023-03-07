# Two-factor interactions

<!-- TFI explanation:

- clear
- non-clear -->
In designs with four-level factors and two-level factors,
Four-level factors are composed of three pseudo-factors.
There are two types of letters: factors and pseudo-factors thus, there are three types of interactions:

- 4-4: interaction between two four-level factors
- 4-2: interaction between a four-level factor and a two-level factor
- 2-2: interaction between two two-level factors

Interactions are clear (CSW 93) when they are not aliased with any main effect and with any other two-factor interaction.
Since there are three type of interactions, an interaction can be clear in three different ways:

- 4-4 clear: not aliased with any other M.E. or 4-4 interaction
- 4-2 clear: not aliased with any other M.E. or 4-2 interaction
- 2-2 clear: not aliased with any other M.E. or 2-2 interaction

If an interaction is 4-4 clear, 4-2 clear, and 2-2 clear, we say that it is **totally clear**.
For simplicity, we define the *clairty* of an interaction as the different ways it can be clear.

The clarity of the interactions of a design is a good indicator of how good the design can accomodate TFI in a model.

To vizualise this information easily, we use the {meth}`fatld.design.Design.clarity_matrix` method of a design.

<!-- add example in code block-->

And the individual numbers in a clarity_matrix can be obtained using the {meth}`fatld.design.Design.clear_from`.
