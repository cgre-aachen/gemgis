from setuptools import setup, find_packages

version = '0.1.2'

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='gemgis',
    version=version,
    packages=find_packages(exclude=('test', 'data', 'notebooks')),
    include_package_data=True,
    install_requires=[
        "numpy",
        "scooby",
        "pandas",
        # "rasterio",
        # "geopandas",
        "typing",
        "shapely",
        "matplotlib",
        "scikit-image",
        "xmltodict",
        "pyvista",
        "mplstereonet",
        "owslib",
        "requests",
        "sklearn",
        "descartes"
        "PyPDF2"
        "tqdm"
    ],
    url='https://github.com/cgre-aachen/gemgis',
    license='LGPL v3',
    author='Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann',
    author_email='alexander.juestel@rwth-aachen.de',
    description='GemGIS is a Python-based, open-source geographic information processing library.' \
                'It is capable of preprocessing spatial data such as vector data (shape files, '
                'geojson files,' \
                'geopackages),' \
                'raster data, data obtained from WMS services or XML/KML files.' \
                'Preprocessed data can be stored in a dedicated Data Class to be passed to the '
                'geomodeling package GemPy' \
                'in order to accelerate to model building process.',
    keywords=['geology', 'geographic', 'structural geology', 'GIS']
)
