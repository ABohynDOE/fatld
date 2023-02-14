# fatld

Virtual environement based on `poetry`, activate it using:

```shell
poetry shell
```

Tests are ran using `pytest`.
Linting based on `ruff`.
Coverage report can be generated using `coverage` inside the virtual environement:

```shell
coverage run -m pytest tests
coverage report -m
```
