.. _contributing:

Contributing to *opcsim*
========================

All development, issue-tracking, and conversation surrounding `opcsim` takes place on GitHub. If you aren't familiar with GitHub or, more generally Git, we recommend you check out some tutorials that will catch you up to speed. Some of our favorites are `Codecademy <https://www.codecademy.com/learn/learn-git>`_ and the `GitHub Guide <https://guides.github.com/activities/hello-world/>`_.

Language Syntax
---------------

The documentation for ``opcsim`` is built with `Sphinx <http://www.sphinx-doc.org/en/master/rest.html>`_ and `sphinx-autodoc <http://www.sphinx-doc.org/en/stable/ext/autodoc.html>`_ using reStructuredText (reST). There are three primary parts that are all pulled together to build the final documentation; they include the API documentation which is pulled directly from the source-code docstrings, the tutorials which are written as jupyter notebooks, and the gallery/examples which are individual python files. More details on each can be found below.

Editing the Main (Overview) Text (.rst files)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are several files that contain some text/information that is not pulled from docstrings (including this file) which can be found in the ``/docs`` directory. They are:

1. ``introduction.rst``
2. ``installing.rst``
3. ``contributing.rst``
4. ``api.rst``
5. ``tutorial.rst``

You can edit the individual files, paying close attention to some that bring in text through the automatic generation process (take a look at ``api.rst`` for example). For the most part, you can simply edit these files using the sphinx rst syntax.

Editing docstrings
^^^^^^^^^^^^^^^^^^

Each function and/or class within the source code has an accompanying docstring. Using sphinx-autodoc, we pull these docstrings and format them in a human-readable way to generate nice, clean documentation! Thus, it is imperative to use the correct formatting of the docstrings, so autodoc can do its job.

Docstrings are written using reST, which allows for very rich documentation, including things like equations, images, and examples. For examples, simply check out any of the files in the source code, or click on "source" on one of the documentation pages. For a primer on how and why we document python code, check out `this post <http://docs.python-guide.org/en/latest/writing/documentation/>`_.

Editing and/or Adding Tutorials
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tutorials are written as jupyter notebooks and are located in the ``/docs/tutorial`` folder. If editing an existing tutorial, one can simply fire up the jupyter notebook as one typically would (``$ jupyter notebook``) from terminal/command line, edit, and save. Changes will take effect once the changes are pushed to the master branch of the repository.

If interested in adding a new tutorial, you can create a new notebook file in the same directory and then edit/save as you did above. It is important to note that the first/top cell in the notebook needs to follow a special format, which can be copied and pasted from an existing tutorial notebook file by clicking within the cell, and copying the raw text. Once you are finished with the new tutorial, you must add the name of the file to the ``Makefile`` within the tutorial directory so that it is built when ``$ make notebooks`` is executed. To do so, simply add a new line below the last entry in the file with form ``tools/nb_to_doc.py <name-of-file-without-extension>``.


Editing and/or Adding Examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each example is a single python executable file that will output a figure, along with the formatted and highlighted code snippet. They are located within the ``/docs/examples`` folder and you can easily add new ones by simply adding a new file. We recommend you start by copying another file (``$ cp old-file.py new-file.py``) and then editing it. Make sure to save the file with a name that is fairly descriptive! Other notes:

* Where it says "_thumb" in the docstring at the top of the file, the following two numbers describe the center point of the thumbnail image that will be created. You can play around with these numbers (always between 0-1) to see what looks best for your example.


Building a Local Version of the Documentation
---------------------------------------------

Once you have edited the documentation you are working on, you can render it by building a local copy. To do so, follow these instructions:

First, navigate to the *docs* directory::

    $ cd /docs

Next, clean out the old files::

    $ make clean

Next, build the notebooks::

    $ make notebooks

Finally, build the html::

    $ make html


Once these steps are completed, you can open up the *index.html* file located in ``/docs/_build/html/`` in your favorite browser. Once you are satisfied with the local version, feel free to send a Pull Request to the repository and it will be integrated into the official documentation.
