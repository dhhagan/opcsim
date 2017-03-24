.. _installing:

Getting Started
---------------

If you are familiar with Python and Python installations, feel free to skip
below to the installation section.

``opcsim`` requires python2.7+ or python3.3+ to be installed on your computer. If
you do not have an existing installation, I would recommend checking out the
`Anaconda Python distribution <https://www.continuum.io/downloads>`_. If you
encounter issues installing Anaconda, check out StackOverflow or a simple
Google search.

Once you have python installed, you can go ahead and install the
``opcsim`` python library.


Dependencies
~~~~~~~~~~~~

+ Python2.7 or Python3.3+
+ `numpy <http://www.numpy.org/>`_
+ `scipy <https://www.scipy.org/>`_
+ `pandas <http://pandas.pydata.org/>`_
+ `matplotlib <http://matplotlib.org/>`_
+ `seaborn <http://seaborn.pydata.org/api.html>`_

Installation
------------

There are several ways to install the package from source depending on your
local dev environment and familiarity with installing python packages. If you
are new to python (and don't have `git` installed), install using the `Install
from Source` option.

Install using `pip`
~~~~~~~~~~~~~~~~~~~

Note: must have git installed

.. code-block:: shell

    $ pip install --upgrade git+git://github.com/dhhagan/opcsim.git


Clone Repository and Install from GitHub
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    $ git clone https://github.com/dhhagan/opcsim.git
    $ cd opcsim
    $ python3 setup.py install


Install from Source
~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    $ wget https://github.com/dhhagan/opcsim/archive/master.zip
    $ unzip master.zip
    $ cd opcsim-master
    $ python3 setup.py install


Testing
~~~~~~~

Testing is automated using `unittests`. To run the unittests with coverage
reporting, run the following commands from the main directory:

.. code-block:: shell

    $ coverage run --source opcsim setup.py test
    $ coverage report -m


Bugs and Issues
~~~~~~~~~~~~~~~

Please report any bugs or issues you find through the `GitHub issues tracker
<https://github.com/dhhagan/opcsim/issues/new>`_. Please provide as much
information as possible that will make it easier to solve/fix the problem.
