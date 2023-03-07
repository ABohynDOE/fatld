# Welcome

## Building regular four-and-two-level designs in python

This package provides several tools to work with **four-and-two-level designs** (FATLD), more specifically, you can:

- Generate designs using the {class}`fatld.design.Design` class
- Characterize these designs using the methods documented in the [documentation](design-documentation).

The {class}`fatld.design.Design` class and its methods are detailed in the following sections.

```{eval-rst}
.. toctree::
   :caption: Design
   :maxdepth: 2

   design
   wlp
   tfi

```

Since this package focuses on *regular* designs, aliasing between the factors is an important focus too.
You can explore aliasing relationships using *words* and *defining relations*, using several functions, presented in the following sections.

```{eval-rst}
.. toctree::
   :caption: Aliasing
   :maxdepth: 2

   factors
   words
   relation

```

The complete module documentation is also available here.

```{eval-rst}
.. toctree::
   :caption: Module documentation
   :maxdepth: 2

   references
   documentation
```

Along with {ref}`genindex`: a list of all the functions, classes and methods of the package.

(contact)=

## Usage

The package can be downloaded from [pypi](https://pypi.org/project/fatld/) and installed on your machine using the following command:

```{code} bash
pip install fatld
```

## Contact & Contribution

For any questions/suggestions related to this project, you can contact me at [alexandre dot bohyn at kuleuven dot be](mailto:alexandre.bohyn@kuleuven.be).

This is an on-going project so any help is welcome.
If you want to contribute, you can either create an [issue](https://github.com/ABohynDOE/fatld/issues/new) on Github or send me an email.

## License

This online work is licensed under a MIT license. Visit [here](https://opensource.org/license/mit/) for more information about the license.
