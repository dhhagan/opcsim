[![PyPI version](https://badge.fury.io/py/opcsim.svg)](https://badge.fury.io/py/opcsim)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://github.com/dhhagan/opcsim/blob/master/LICENSE)
[![Build Status](https://travis-ci.org/dhhagan/opcsim.svg?branch=master)](https://travis-ci.org/dhhagan/opcsim)
[![Test Coverage](https://api.codeclimate.com/v1/badges/62e396e65ce4ade478fc/test_coverage)](https://codeclimate.com/github/dhhagan/opcsim/test_coverage)
[![Maintainability](https://api.codeclimate.com/v1/badges/62e396e65ce4ade478fc/maintainability)](https://codeclimate.com/github/dhhagan/opcsim/maintainability)

# opcsim
opcsim is a Python library for simulating low-cost Optical Particle Counters and
their response to various aerosol distributions.

## Documentation

Full online documentation can be found [here][1].

The docs include a [tutorial][2], an [example gallery][3], and an [API Reference][4].

In addition, documentation can be built locally for development purposes. To do so, please check out the complete details in the *contributing to opcsim* section of the documentation.

## Dependencies

Opcsim supports Python 2.7 and 3.4+.

Installation requires [scipy][5], [numpy][6], [pandas][7], [matplotlib][8],
and [seaborn][9].


## Installation

To install (or upgrade to) the latest stable release:

    $ pip install opcsim [--upgrade] --process-dependency-links

To install the development version directly from GitHub using pip:

    $ pip install git+https://github.com/dhhagan/opcsim.git

In addition, you can either clone the repository and install from source or download/unzip the zip file and install from source:

    $ git clone https://github.com/dhhagan/opcsim.git
    $ cd /opcsim
    $ python setup.py install


## Testing

All tests are run via Travis.ci pre-merge. For results of these tests, please click on the link in the above travis badge. In addition, you can run tests locally.

To run tests locally:

    $ python setup.py test

To run the tests locally using coverage (adds `coverage` as a dependency):

    $ coverage run --source opcsim setup.py test

    # View the report
    $ coverage report -m


## Development

**opcsim** development takes place on GitHub. Issues and bugs can be submitted and tracked via the [GitHub Issue Tracker][10] for this repository.


[1]: https://dhhagan.github.io/opcsim/
[2]: https://dhhagan.github.io/opcsim/tutorial.html
[3]: https://dhhagan.github.io/opcsim/examples/index.html
[4]: https://dhhagan.github.io/opcsim/api.html
[5]: https://www.scipy.org/
[6]: http://www.numpy.org/
[7]: http://pandas.pydata.org/
[8]: http://matplotlib.org/
[9]: https://seaborn.pydata.org/
[10]: https://github.com/dhhagan/opcsim/issues
