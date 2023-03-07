# `fatld`: four-and-two-level designs

[![PyPI version](https://badge.fury.io/py/fatld.svg)](https://badge.fury.io/py/fatld)
[![CI](https://github.com/ABohynDOE/fatld/actions/workflows/CI.yml/badge.svg)](https://github.com/ABohynDOE/fatld/actions/workflows/CI.yml)

Generate and characterize designs with four-and-two-level (FATL) factors.

## Development

Virtual environment based on `poetry`, activate it using:

```shell
poetry shell
```

Linting is based on `ruff`, configuration is found in the [pyproject.toml](pyproject.toml) file.

Tests are ran using `pytest` and a coverage report can be generated using `coverage` inside the virtual environment:

```shell
coverage run -m pytest tests
coverage report -m
```
