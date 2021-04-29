.. gemgis documentation master file, created by
   sphinx-quickstart on Mon Nov  2 22:04:17 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Contributing
===========================================================

Thank you for considering to contribute to GemGIS!

As GemGIS is an open-source and a community-driven project, we need people like you to make the package even better, simpler to use for beginners, and to extend its functionality for more advanced users. There are several ways to contribute to GemGIS:

- Submitting bug report or requests for new features
- Adding new or expanding existing tutorials and examples
- Fixing typos and improving documentation
- Adding new functionality to the package

How to get started with GemGIS?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Work through the examples and tutorials or work with your own data.
- If you run into a problem and you think it is not related to GemGIS directly but rather to the installation or your data or something just does not work, have a look at our `Github Discussion <https://github.com/cgre-aachen/gemgis/discussions>`_ to ask your question. Maybe someone already had the same issue before and the answer is already available.
- If you run into a problem and you think it is related to GemGIS directly please `open a new issue on Github <https://github.com/cgre-aachen/gemgis/issues/new?assignees=&labels=&template=bug_report.md&title=>`_. Our template will guide you through the necessary information.
- If you work with a dataset and you think something is missing in GemGIS please `open a new feature request on Github <https://github.com/cgre-aachen/gemgis/issues/new?assignees=&labels=&template=feature_request.md&title=>`_. Our template will guide you through the necessary information. We will then try to implement the requested feature or support you in contributing yourself to the package.
- If you already have code for a new feature/tutorial/example, we would appreciate it if you open a pull request and contribute this way. Our team will then support you in terms of testing your code, suggesting improvements and adding tutorials or examples for your implemented feature.
- If you find typos or other error in our README or the documentation open a new issue or fix the error right away by opening a pull request.

How to get in touch with the GemGIS developers?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can reach the developers of GemGIS either by opening a new issue on Github (feature request or bug report), open a pull request to contribute directly to the package, ask a question in the Github Discussions forum or via the `Slack workspace of the Software Underground <https://softwareunderground.org/slack>`_.

How to contribute code?
~~~~~~~~~~~~~~~~~~~~~~~~

**Is this your first contribution?**

Please take a look at these resources to learn about git and pull requests but do not hesitate to ask questions if you get stuck:

* `How to Contribute to Open Source <https://opensource.guide/how-to-contribute/>`_.
* `Aaron Meurer's tutorial on the git workflow <http://www.asmeurer.com/git-workflow/>`_.
* `How to Contribute to an Open Source Project on GitHub <https://egghead.io/courses/how-to-contribute-to-an-open-source-project-on-github>`_.

General guidelines
__________________

We follow the `git pull request workflow <http://www.asmeurer.com/git-workflow/>`_ to
make changes to our codebase.
Every change made goes through a pull request, even our own, so that our
`continuous integration <https://en.wikipedia.org/wiki/Continuous_integration>`_ services (Github Actions)
have a change to check that the code is up to standards and passes all our tests.
This way, the *main* branch is always stable.

General guidelines for pull requests (PRs)
___________________________________________

* **Open an issue first** describing what you want to do. If there is already an issue
  that matches your PR, leave a comment there instead to let us know what you plan to
  do.
* Each pull request should consist of a **small** and logical collection of changes.
* Larger changes should be broken down into smaller components and integrated
  separately.
* Bug fixes should be submitted in separate PRs.
* Describe what your PR changes and *why* this is a good thing. Be as specific as you
  can. The PR description is how we keep track of the changes made to the project over
  time.
* Do not commit changes to files that are irrelevant to your feature or bugfix (eg:
  `.gitignore`, IDE project files, etc).
* Write descriptive commit messages. Chris Beams has written a
  `guide <https://chris.beams.io/posts/git-commit/>`_ on how to write good commit
  messages.
* Be willing to accept criticism and work on improving your code; we don't want to break
  other users' code, so care must be taken not to introduce bugs.
* Be aware that the pull request review process is not immediate, and is generally
  proportional to the size of the pull request.

Setting up your environment
___________________________

We highly recommend using `Anaconda <https://www.anaconda.com/download/>`_ and the `conda`
package manager to install and manage your Python packages.
It will make your life a lot easier!

The repository includes a conda environment file `environment_dev.yml` with the
specification for all development requirements to build and test the project.
Once you have forked and clone the repository to your local machine, you use this file
to create an isolated environment on which you can work.

Docstrings
__________

**All docstrings** should follow the
`numpy style guide <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.
All functions/classes/methods should have docstrings with a full description of all
arguments and return values.

To play nicely with Jupyter and IPython, **keep docstrings
limited to 79 characters** per line. We don't have a good way of enforcing this
automatically yet, so please do your best.

Testing your code
_________________

Automated testing helps ensure that our code is as free of bugs as it can be.
It also lets us know immediately if a change we make breaks any other part of the code.

All of our test code and data are stored in the `tests` subpackage.
We use the `pytest <https://pytest.org/>`_ framework to run the test suite.

Please write tests for your code so that we can be sure that it won't break any of the
existing functionality.
Tests also help us be confident that we won't break your code in the future.

If you're **new to testing**, see existing test files for examples of things to do.
**Don't let the tests keep you from submitting your contribution!**
If you're not sure how to do this or are having trouble, submit your pull request
anyway.
We will help you create the tests and sort out any kind of problem during code review.

Documentation
_____________

Most documentation sources are in the `docs` folder.
We use `sphinx <http://www.sphinx-doc.org/>`_ to build the web pages from these sources.
To build the HTML files::

   cd docs
   make all


This will build the HTML files in `docs/_build/html`.
Open `docs/_build/html/index.html` in your browser to view the pages.

The API reference is manually assembled in `docs/api_reference/index.rst`.
The *autodoc* sphinx extension will automatically create pages for each
function/class/module listed there.

Code Review
___________

After you've submitted a pull request, you should expect to hear at least a comment
within a couple of days.
We may suggest some changes or improvements or alternatives.

Some things that will increase the chance that your pull request is accepted quickly:

* Write a good and detailed description of what the PR does.
* Write tests for the code you wrote/modified.
* Readable code is better than clever code (even with comments).
* Write documentation for your code (docstrings) and leave comments explaining the
  *reason* behind non-obvious things.
* Include an example of new features in the gallery or tutorials.
* Follow the `PEP8 <http://pep8.org>`_ style guide for code and the
  `numpy guide <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
  for documentation.

Pull requests will automatically have tests run by Github Actions.
Github will show the status of these checks on the pull request.
Try to get them all passing (green).
If you have any trouble, leave a comment in the PR or get in touch with us


Attribution
===========
This contributing document is largely based upon the work by the Fatiando a Terra project.