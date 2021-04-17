Before a release.
----------------
# set version number in setup.py, also in the config file of the documentation and init of the package
[ ] setup.py
[ ] gemgis.__init__
# ~[ ] in the config for sphinx~ this is taken from gemgis.__init__

Github release
--------------
    # add new tag
    $ git tag X.X -m "Add X.X tag for PyPI"
    # push git tag
    $ git push --tags origin main

PyPi release
------------
    #  First create the dist - create dist folder in the main directory containing the compressed package data
    python setup.py sdist bdist_wheel

    # Second upload the distributions - upload all contents of the dist folder
    twine upload dist/*


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