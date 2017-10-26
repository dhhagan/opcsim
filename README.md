[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://github.com/dhhagan/opcsim/blob/master/LICENSE)
[![Build Status](https://travis-ci.com/dhhagan/opcsim.svg?token=h93DryaPQEsWsPpUCPym&branch=master)](https://travis-ci.com/dhhagan/opcsim)

# opcsim
opcsim is a Python library for simulating low-cost Optical Particle Counters and
their response to various aerosol distributions.

## Documentation

Full online documentation can be found [here][1].

The docs include a [tutorial][2], an [example gallery][3], and an [API Reference][4].

## Dependencies

Opcsim supports Python 2.7 and 3.4+.

Installation requires [scipy][5], [numpy][6], [pandas][7], [matplotlib][8],
and [seaborn][9].


## Installation

To install the latest stable release:

    $ pip install opcsim

To install the development version directly from GitHub using pip:

    $ pip install git+https://github.com/dhhagan/opcsim.git


## Testing

To run tests locally:

    $ python setup.py test

To run the tests locally using coverage:

    $ coverage run --source opcsim setup.py test

    # View the report
    $ coverage report -m

## Development

Opcsim development takes place on GitHub. Issues and bugs can be submitted and
tracked via the [GitHub Issue Tracker][10] for this repository.


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
