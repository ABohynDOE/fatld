# Four And Two Level Designs

[![PyPI version](https://badge.fury.io/py/fatld.svg)](https://badge.fury.io/py/fatld)
[![CI](https://github.com/ABohynDOE/fatld/actions/workflows/CI.yml/badge.svg)](https://github.com/ABohynDOE/fatld/actions/workflows/CI.yml)

The `fatld` package contains functionality to generate and characterize designs with four-and-two-level (FATL) factors.
Design characteristics include word length pattern, defining relation, and number of clear interactions.
For more information about the package see the documentation at [https://abohyndoe.github.io/fatld/](https://abohyndoe.github.io/fatld/).
A large collection of FATL designs can be explored interactively using a web app at [https://abohyndoe.shinyapps.io/fatldesign-selection-tool/](https://abohyndoe.shinyapps.io/fatldesign-selection-tool/).

## Usage

The package can be used from Python:

```python
>>> import fatld
>>> D = fatld.Design(runsize=32, m=1, cols=[21, 27, 29])
>>> D.wlp()
[1, 3, 3, 0, 0]
>>> D.defining_relation()
['A1cef', 'A3deg', 'A1cdeh']
>>> print("There are %s 2-2 interactions clear from any main effect or other two-factor interaction." % D.clear('2-2'))
There are 6 2-2 interactions clear from any main effect or other two-factor interaction.
>>> print("The design contains %s four-level factors and %s two-level factors" % (D.m, D.n))
The design contains 1 four-level factors and 6 two-level factors
```

For more examples see the documentation.

## Installation

The Python interface to the package is available on pypi. Installation can be done using the following command:

```bash
pip install fatld
```

## Development

All development is done in a virtual environment based on `poetry`, activate it using:

```shell
poetry shell
```

### Code style

* Try to follow the [PEP 8](https://www.python.org/dev/peps/pep-0008/) style guide. A usefull tool for automated formatting is [black](https://pypi.python.org/pypi/black). We do allow lines upto 120 characters.
* Document functions using the [Numpy docstring](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html#example-numpy) convention

* Linting is based on `ruff`, configuration is found in the [pyproject.toml](pyproject.toml) file.

* Tests are ran using `pytest` and a coverage report can be generated using `coverage` inside the virtual environment:

    ```shell
    coverage run -m pytest tests
    coverage report -m
    ```

### Submitting code

If you would like to contribute, please submit a pull request. (See the [Github Hello World](https://guides.github.com/activities/hello-world/) example, if you are new to Github).

By contributing to the repository you state you own the copyright to those contributions and agree to include your contributions as part of this project under the BSD license.

### Bugs reports and feature requests

To submit a bug report or feature request use the [Github issue tracker](https://github.com/ABohynDOE/fatld/issues).
Search for existing and closed issues first. If your problem or idea is not yet addressed, please open a new issue.

### Contact

For further information please contact alexandre dot bohyn at kuleuven dot be
