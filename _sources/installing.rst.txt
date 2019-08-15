.. _installing:

================================
Installation and Getting Started
================================

If you are familiar with Python and installing python libraries, feel free to skip
below to the installation section.


``opcsim`` requires python2.7+ or python3.3+ to be installed on your computer. If
you do not have an existing installation, I would recommend checking out the
`Anaconda Python distribution <https://www.continuum.io/downloads>`_. If you
encounter issues installing Anaconda, check out `StackOverflow <https://stackoverflow.com/search?q=anaconda>`_ or a simple
Google search.

Once you have python installed, you can go ahead and install the
``opcsim`` python library.


**Dependencies**

+ Python2.7 or Python3.3+
+ `numpy <http://www.numpy.org/>`_
+ `scipy <https://www.scipy.org/>`_
+ `pandas <http://pandas.pydata.org/>`_
+ `matplotlib <http://matplotlib.org/>`_
+ `seaborn <http://seaborn.pydata.org/api.html>`_


---------------------
Installing ``opcsim``
---------------------

There are several ways to install the package from source depending on your
local development environment and familiarity with installing python libraries. If you
are new to python (and don't have `git` installed), install using the `Install
from Source` option.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Install directly from pypi (best option)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installing directly from pypi using pip is the easiest way forward and will
automagically install any dependencies that aren't already installed.

.. code-block:: shell

    $ pip install opcsim [--upgrade]


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Install from GitHub using `pip`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note: must have git installed

.. code-block:: shell

    $ pip install --upgrade git+git://github.com/dhhagan/opcsim.git


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Clone Repository and Install from GitHub
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you plan on contributing to ``opcsim``, you will probably want to fork the
library, and then clone it to your local environment. You can then install from
source directly.

.. code-block:: shell

    $ git clone https://github.com/dhhagan/opcsim.git
    $ cd opcsim/
    $ python3 setup.py install


~~~~~~~~~~~~~~~~~~~
Install from Source
~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    $ wget https://github.com/dhhagan/opcsim/archive/master.zip
    $ unzip master.zip
    $ cd opcsim-master/
    $ python3 setup.py install


-------
Testing
-------

Testing is automated using `unittests`. To run the unittests with coverage
reporting, run the following commands from the main directory:

.. code-block:: shell

    $ coverage run --source opcsim setup.py test
    $ coverage report -m

Unittests are also run automatically through continuous integration via TravisCI
upon every pull request and code coverage is tracked online with `Code Climate <https://codeclimate.com/>`_.

-------------------------------
Reporting Bugs and other Issues
-------------------------------

Please report any bugs or issues you find through the `GitHub issues tracker
<https://github.com/dhhagan/opcsim/issues/new>`_. Please provide as much
information as possible that will make it easier to solve/fix the problem. Useful
information to include would be the operating system, python version, and version
of the ``opcsim`` library as well as any dependencies. If there are issues with
graphics, screenshots are very helpful!
