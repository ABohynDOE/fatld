# Welcome

## Building regular four-and-two-level designs in python

This package provides several tools to work with **four-and-two-level designs** (FATLD).
These designs only have factors with two level or four levele.
They are also *regular*, meaning that the aliasing between the factors is either complete or null.

To explore the relationships between factors in a design, you can use *words* and *defining relations*.
These concepts are explained in the following sections, along with example using functions from of this package:

```{eval-rst}
.. toctree::
   :maxdepth: 1

   factors
   words
   relation

```

However, the core of this package is the {class}`fatld.design.Design` class that allows you to generate FATL designs, and characterize them using several criteria.
The class and its related methods are detailed, using examples, in these following sections:

```{eval-rst}
.. toctree::
   :maxdepth: 1

   design
   wlp
   tfi

```

The complete module documentation is also available here.

```{eval-rst}
.. toctree::
   :maxdepth: 2

   references
   documentation
```

Along with {ref}`genindex` a list of all the functions, classes and methods of the package.

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
