"""
Contributors: Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann, Marius Pischke

GemGIS is a Python-based, open-source spatial data processing library.
It is capable of preprocessing spatial data such as vector data
raster data, data obtained from online services and many more data formats.
GemGIS wraps and extends the functionality of packages known to the geo-community
such as GeoPandas, Rasterio, OWSLib, Shapely, PyVista, Pandas, and NumPy.

GemGIS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GemGIS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (LICENSE) for more details.

This file is modified from https://github.com/cgre-aachen/gempy/blob/master/scripts/convert_jupyter_to_py.py

"""

import os
import json
import pypandoc as pdoc


def convert_ipynb_to_gallery(nb, new_file):
    python_file = ""

    nb_dict = json.load(open(nb, encoding="utf8", errors='ignore'))
    cells = nb_dict['cells']

    for i, cell in enumerate(cells):
        if i == 0:
            if cell['cell_type'] != 'markdown':
                rst_source = os.path.basename(file_name[:-5])
                rst_source = bytes(rst_source, 'utf-8').decode('utf-8', 'ignore')
                python_file = '"""\n' + rst_source + '\n"""'
                source = ''.join(cell['source'])
                python_file = python_file + '\n' * 2 + source

            else:
                b = cell['source']
                print(b)
                a = bytes(cell['source'][0], 'utf-8').decode('utf-8', 'ignore')
                print(a)
                md_source = ''.join(a)
                rst_source = pdoc.convert_text(md_source, 'rst', 'md')
                print(rst_source)
                rst_source = bytes(rst_source, 'utf-8').decode('utf-8', 'ignore')
                python_file = '"""\n' + rst_source + '\n"""'
        else:
            if cell['cell_type'] == 'markdown':
                md_source = ''.join(cell['source'])
                rst_source = pdoc.convert_text(md_source, 'rst', 'md')
                rst_source = rst_source.encode().decode('utf-8', 'ignore')
                commented_source = '\n'.join(['# ' + x for x in
                                              rst_source.split('\n')])

                python_file = python_file + '\n\n\n' + '# %%' + '\n' + commented_source

            elif cell['cell_type'] == 'code':
                source = ''.join(cell['source'])
                python_file = python_file + '\n' * 2 + '# %% \n' + source

    python_file = python_file.replace("\n%", "\n# %")
    open(new_file, 'w', newline='', errors='ignore').write(python_file)


directory = '../examples/notebooks'

for root, dirs, files in os.walk(directory):
    for file in files:
        if file.endswith('.ipynb'):
            if file.find('checkpoint') != -1:
                continue

            file_name = root + '/' + file

            nb = root + '/' + file

            new_dir = '../examples/models'

            new_file = new_dir + '/' + file.replace('.ipynb', '.py')
            print(new_file)
            try:
                os.makedirs(new_dir)
            except:
                pass

            convert_ipynb_to_gallery(nb, new_file)