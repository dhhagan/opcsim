try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

__version__ = '0.1.0'

AUTHOR          = 'David H Hagan'
AUTHOR_EMAIL    = 'dhagan@mit.edu'

setup(
    name = 'opcsim',
    version = __version__,
    packages = ['opcsim'],
    description = "Python package for ...",
    author = AUTHOR,
    author_email = AUTHOR_EMAIL,
    maintainer = AUTHOR,
    maintainer_email = AUTHOR_EMAIL,
    license = 'MIT',
    url = 'https://github.com/dhhagan/opcsim',
    keywords = ['atmospheric chemistry'],
    test_suite = 'tests',
    install_requires = [
        'pandas',
        'numpy',
        'seaborn',
        'scipy'
    ],
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
		'Intended Audience :: Education',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
		'Topic :: Scientific/Engineering :: Atmospheric Science',
		'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ]
)
