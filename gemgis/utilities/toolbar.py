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

"""

import os
import ipywidgets as widgets
from ipyleaflet import WidgetControl
from ipyfilechooser import FileChooser
from IPython.display import display


def toolbar(m):
    """Function creating the main toolbar for the interactive ipyleaflet map

    Parameters
    __________

        m : ipyleaflet.Map
            ipyleaflet Map

    """

    # Defining widget width and padding
    padding = "0px 0px 0px 5px"  # upper, right, bottom, left

    # Creating toolbar button
    toolbar_button = widgets.ToggleButton(
        value=False,
        tooltip="Toolbar",
        icon="wrench",
        layout=widgets.Layout(width="28px", height="28px", padding=padding),
    )

    # Creating close button
    close_button = widgets.ToggleButton(
        value=False,
        tooltip="Close the tool",
        icon="times",
        button_style="primary",
        layout=widgets.Layout(height="28px", width="28px", padding=padding),
    )

    # Creating horizontal box
    toolbar = widgets.HBox([toolbar_button])

    # Defining function for toggling closing button
    def close_click(change):
        if change["new"]:
            toolbar_button.close()
            close_button.close()
            filechooser_widget.close()
            toolbar.close()

    close_button.observe(close_click, "value")

    # Creating grid
    rows = 2
    cols = 2
    grid = widgets.GridspecLayout(rows, cols, grid_gap="0px", layout=widgets.Layout(width="62px"))

    # Defining buttons
    icons = ["folder-open", "map", "gears", "question"]

    # Filling grid with buttons
    for i in range(rows):
        for j in range(cols):
            grid[i, j] = widgets.Button(description="", button_style="primary", icon=icons[i * rows + j],
                                        layout=widgets.Layout(width="28px", padding="0px"))

    # Creating vertical box from toolbar button
    toolbar = widgets.VBox([toolbar_button])

    # Defining function for toggling toolbar button
    def toolbar_click(change):
        if change["new"]:
            toolbar.children = [widgets.HBox([close_button, toolbar_button]), grid]
        else:
            toolbar.children = [toolbar_button]

    toolbar_button.observe(toolbar_click, "value")

    # Creating widget control for toolbar
    toolbar_ctrl = WidgetControl(widget=toolbar, position="topright")

    # Adding toolbar to map
    m.add_control(toolbar_ctrl)

    output = widgets.Output()
    output_ctrl = WidgetControl(widget=output, position="topright")

    buttons = widgets.ToggleButtons(
        value=None,
        options=["Open", "Reset", "Close"],
        tooltips=["Open", "Reset", "Close"],
        button_style="primary",
    )
    buttons.style.button_width = "80px"

    data_dir = os.path.abspath('.')

    fc = FileChooser(data_dir)
    fc.use_dir_icons = True
    fc.filter_pattern = ['*.shp', '*.geojson']

    filechooser_widget = widgets.VBox([fc, buttons])

    def button_click(change):
        if change["new"] == "Open" and fc.selected is not None:
            if fc.selected.endswith(".shp"):
                m.add_shapefile(fc.selected, layer_name=fc.selected_filename)
            elif fc.selected.endswith(".geojson"):
                m.add_geojson(fc.selected, layer_name=fc.selected_filename)
        elif change["new"] == "Reset":
            fc.reset()
        elif change["new"] == "Close":
            fc.reset()
            m.remove_control(output_ctrl)

    buttons.observe(button_click, "value")

    def tool_click(b):
        with output:
            output.clear_output()
            if b.icon == "folder-open":
                display(filechooser_widget)
                m.add_control(output_ctrl)
            elif b.icon == "gears":
                import whiteboxgui.whiteboxgui as wbt

                if hasattr(m, "whitebox") and m.whitebox is not None:
                    if m.whitebox in m.controls:
                        m.remove_control(m.whitebox)

                tools_dict = wbt.get_wbt_dict()
                wbt_toolbox = wbt.build_toolbox(
                    tools_dict, max_width="800px", max_height="500px"
                )

                wbt_control = WidgetControl(
                    widget=wbt_toolbox, position="bottomright"
                )

                m.whitebox = wbt_control
                m.add_control(wbt_control)

    for i in range(rows):
        for j in range(cols):
            tool = grid[i, j]
            tool.on_click(tool_click)
