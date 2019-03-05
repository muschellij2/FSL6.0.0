armawrap
========


A wrapper around `Armadillo <http://arma.sourceforge.net/>`_ which provides a
`Newmat <http://www.robertnz.net/nm_intro.htm>`_ style API.

The current version of ``armawrap`` is written against:

- newmat10D
- Armadillo 5.200


Please cite:

    | Conrad Sanderson.
    | Armadillo: An Open Source C++ Linear Algebra Library for
    | Fast Prototyping and Computationally Intensive Experiments.
    | Technical Report, NICTA, 2010.


Usage
-----

Like Armadillo, ``armawrap`` is a header library, so no compilation is
required.  Just put ``path/to/armawrap`` on your include path, then ``#include
"newmat.h"``, and code against the Newmat API.


Integration into FSL
--------------------


1. Compile the dummy ``newmat`` wrapper. In order to avoid having to modify
   loads of Makefiles which have ``'-lnewmat'``, you can simply create a dummy
   ``'libnewmat.a'`` which contains no symbols (``dummy_newmat/newmat.cpp`` is
   literally an empty file).

2. Modify ``$FSLDIR/config/<your_machine_type>/externallibs.mk``, changing
   these lines::

     LIB_NEWMAT = ${FSLEXTLIB}
     INC_NEWMAT = ${FSLEXTINC}/newmat

   to::

     LIB_NEWMAT = <location_of_dummy_libnewmat.a_file> -llapack -lblas
     INC_NEWMAT = <location_of_armawrap_directory>     -DARMA_USE_LAPACK -DARMA_USE_BLAS

   The Armadillo README file suggests that, under OSX, you should replace
   ``'-llapack -lblas'`` with ``'-framework Accelerate'``.


Running the tests
-----------------


.. image:: https://git.fmrib.ox.ac.uk/fsl/armawrap/badges/master/coverage.svg


To run the tests, simply call the ``run_tests.sh`` script. This script assumes
that you have a C++ development environment configured to compile Armadillo
code::


  path/to/armawrap/tests/run_tests.sh


The ``run_tests.sh`` script will return an exit code of zero if all tests
passed, or non-zero if any tests failed.


Differences between ``armawrap`` and ``newmat``
-----------------------------------------------


The following behavioural differences between ``armawrap`` and ``newmat`` are
known to exist.


Matrix assignment
^^^^^^^^^^^^^^^^^


With ``armawrap``, it is possible to perform lossy assignments between
matrices of different types:

```
Matrix m(3, 3);
m << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9;
DiagonalMatrix d;
d = m;
cout << d << endl;
```

This will result in the diagonal elements being filled:

```
1 0 0
0 5 0
0 0 9
```

Under ``newmat`` this will cause a runtime exception to be raised, decause
``newmat`` only allows assignments which do not lose information.


Insertion of vectors into diagonal/triangular/symmetric matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


In ``armawrap``, it is possible to insert a vector into a diagonal (or
triangular, or symmetric) matrix, with preditable results:

```
DiagonalMatrix d;
RowVector r(5);
r << 5 << 4 << 3 << 2 << 1;
d << r;
cout << d << endl;
```

This will result in a matrix with the values of ``r`` along the diagonal:
```
5
  4
    3
      2
        1
```

This is not possible in ``newmat`` in the above code, ``d`` will be set to a
``(1, 1)`` matrix containing ``5``.

Column vectors produce slightly different behaviour:

```
DiagonalMatrix d;
ColumnVector c(5);
c << 5 << 4 << 3 << 2 << 1;
d << r;
```

In this code, ``d`` will be set to a ``(5, 5)`` matrix containing ``5`` as its
first element, and zeros everywhere else.


Calling ``.Storage`` on a sub-view
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


In ``armawrap``, you can call ``.Storage`` on a subview, e.g.:


```
SymmetricMatrix m(5);
cout << m.Row(4).Storage() << endl;
```


This will output ``4`` (as symmetric matrices are internally represented as
lower-triangular matrices).


This is not possible in ``newmat``.


.. note:: The `.Storage` method is not implemented on sub-views of the band
          matrix types.


``armawrap`` development
------------------------


The ``armawrap`` source code is located at:

  https://git.fmrib.ox.ac.uk/fsl/armawrap


All changes to the source must occur via merge requests to the master
branch. New releases are denoted by a tag off the master branch, with the tag
name equal to the new version number.


``armawrap`` follows `Semantic Versioning <https://www.semver.org>`_
conventions. The current version of ``armawrap`` can be found in
``armawrap.hpp``.


When you wish to make a new release, ensure that the following criteria have
been met before tagging the new release:


 - Version number updated in ``armawrap.hpp``
 - Change log updated


License
-------


``armawrap`` is released under the Apache License, Version 2.0. ``armawrap`
includes a copy of the Armadillo source, which is released under the Mozilla
Public License, Version 2.0.


The ``armawrap/tests`` directory includes a copy of the ``newmat`` 10D source
code, which is released under a liberal open source license - see the `newmat
<http://robertnz.net/nm10.htm#use>`_ web site for details.
