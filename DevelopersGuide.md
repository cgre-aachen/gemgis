Before a release
----------------
Set the new version number (major.minor.patch) in ``setup.py``, also in the ``conf.py`` file of the documentation and ``__init__.py`` file of the package.

Github release
--------------
    # add new tag
    $ git tag X.X -m "Add X.X tag for PyPI"
    # push git tag
    $ git push --tags origin main

PyPi release
------------

A PyPi release it automatically being created using Github Actions when a new Github release is published.


Contributing to GemGIS
----------------------

### Type of commits:

ENH: Enhancement, new functionality
BUG: Bug fix
DOC: Additions/updates to documentation
TST: Additions/updates to tests
BLD: Updates to the build process/scripts
PERF: Performance improvement
CLN: Code cleanup


### Doc strings:

- Use sphinx reference a lot
- Use decorator for arguments


Sphinx Gallery Building
-----------------------

Example and tutorial notebooks are stored in the root folder under examples. In order to use Sphinx Gallery, these notebooks have to be converted to ``.py`` files. This can be done using the Anaconda prompt. Navigate to th ``gemgis`` folder containing the ``.py`` file and execute

    # python ipynb_to_py.py 
    
The Python files are then stored in the respective folder and will be accessed by the Sphinx Gallery to build the documentation. Navigate to the docs folder and execute

    # make html