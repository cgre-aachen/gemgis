from setuptools import setup, find_packages
from os import path

version = '1.0.0'

# Loading Readme for Description on PyPi
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='gemgis',
    version=version,
    packages=find_packages(exclude=('test', 'data', 'notebooks')),
    include_package_data=True,
    install_requires=[],
    url='https://github.com/cgre-aachen/gemgis',
    license='LGPL v3',
    author='Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann, Marius Pischke',
    author_email='alexander.juestel@rwth-aachen.de',
    description="""GemGIS is a Python-based, open-source spatial data processing library.
                It is capable of preprocessing spatial data such as vector data 
                raster data, data obtained from online services and many more data formats.
                GemGIS wraps and extends the functionality of packages known to the geo-community
                such as GeoPandas, Rasterio, OWSLib, Shapely, PyVista, Pandas, and NumPy.""",
    keywords=['geology', 'geographic', 'structural geology', 'GIS', 'spatial data'],
    long_description=long_description,
    long_description_content_type='text/markdown'
)
