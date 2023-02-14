# `fatld`: four-and-two-level designs

Generate and characterize designs with four-and-two-level (FATL) factors.

## Development

Virtual environement based on `poetry`, activate it using:

```shell
poetry shell
```

Linting is based on `ruff`, configuration is found in the [pyproject.toml](pyproject.toml) file.

Tests are ran using `pytest` and a coverage report can be generated using `coverage` inside the virtual environement:

```shell
coverage run -m pytest tests
coverage report -m
```
