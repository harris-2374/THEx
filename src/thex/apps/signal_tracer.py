"""
Author: Andrew Harris
Date: 1/28/2021
Python Version 3.8
This Dash tool visualizes p-distance data created by p-distance-calculator.py
It shows data by single or multiple-samples, single or multi-chromosomes.
"""
import base64
import io

import numpy as np
import pandas as pd
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt
import dash_daq as daq
import plotly.express as px
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from thex.app import app
from thex.apps.utils import tree_utils, signal_utils, graph_options, data_utils
from thex.apps.docs import docs_pdist_contents

############################### Templates ###############################
GRAPH_TEMPLATES = graph_options.graph_templates()
SNAPSHOT_FILE_OPTIONS = graph_options.snapshot_file_type()
COLOR_SWATCHES = graph_options.color_swatches()
SCALE_OPTIONS = graph_options.figure_output_scales()
FONT_FAMILIES = graph_options.font_families()

colorscales = px.colors.named_colorscales()


docs_sidebar = html.Div(
    [
        html.Div("Contents", className="pdist-docs-title"),
        html.Hr(),
        dbc.Nav(
            [
                dbc.Button("Usage", id="pdist-docs-overview-button", className="pdist-docs-buttons", color="info", outline=True),
                dbc.Button("Input File Structure", id="pdist-docs-input-structure-button", className="pdist-docs-buttons", color="info", outline=True),
                dbc.Button("Graph Customization + Export Options", id="pdist-docs-graphing-options-button", className="pdist-docs-buttons", color="info", outline=True),
            ],
            vertical=True,
            pills=True,
            className="pdist-docs-buttons",
        ),
    ],
    # style={"background-color": "#575151"}
)

docs_pdist_df = pd.DataFrame(
    {"Chromosome": ["Chr1", "Chr1", "Chr2", "Chr2"],
     "Window": [100, 100, 200, 200],
     "Sample": ["Tiger1", "Tiger2", "Tiger1", "Tiger2"],
     "Value": [0.03, 0.05, 0.02, 0.05]}
)

tv_video = html.Div(children = [
        html.Video(
            controls = True,
            id = 'movie_player',
            src = "https://www.w3schools.com/html/mov_bbb.mp4",
            autoPlay=True,
            style={"height": "100%", "width": "100%"},
        ),
    ])


docs_input_structure_contents = html.Div(
    children=[
        dbc.Col([
            html.H2(
                children=["Input File Structure"],
            ),
            dcc.Markdown(
                children=[
                """
                The input for the Signal Tracer is a flat, tab-delimited file consisting of four columns: 
                Chromosome, Window, Sample, and pDistance. Ensure the headers of your input file match the example headers!
                """
                ],
            ),
            html.Div([
                dt.DataTable(
                    columns=[{"name": i, "id": i} for i in docs_pdist_df.columns],
                    data=docs_pdist_df.to_dict('records'),
                    style_header={
                        'fontWeight': 'bold',
                        'backgroundColor': 'rgb(230, 230, 230)',
                        'textAlign': 'center'
                    },
                    style_cell={
                        'color': 'black',
                        'textAlign': 'center'
                    },
                    style_data_conditional=[                
                        {
                            "if": {"state": "selected"},
                            "backgroundColor": "inherit !important",
                            "border": "inherit !important",
                        }
                    ],
                ),
            ]),
            html.Hr(style={"background-color": "white"}),
        ], width={"size": 6, "offset": 3}),
    ],
    className="pdist-docs-content-scroll-div"
)

docs_graphing_options_contents = html.Div(
    children=[
        dbc.Col([
            html.H2(
                children=["Graph Customization"],
            ),
            dcc.Markdown(
                children=[
                """
                Graph Customization allow you to customize the look of the graphs in an easy-to-use way. Dropdowns provide lists of themes, color palettes, gridline toggles, and several output 
                formats to choose from. Making changes to the graph customizations will update all graphs currently loaded on the browser. Through these customizations, one can create publication ready graphs that 
                can be dropped directly into your in-progress manuscript. One note to mention about the output dimension option, we have allowed users the choice to specify 
                the output dimensions in pixels, however leaving the inputs blank will revert each graph to default dimensions that are optimized for a word document or pdf file.
                """
                ],
            ),
            html.Hr(style={"background-color": "white"}),
            html.H2(
                children=["Export Options"],
            ),
            dcc.Markdown(
                children=[
                """
                Signal Tracer offers a single export option that allows you to extract local information from your input file. Just like in Tree Viewer, zooming into a region of a chromosome in single 
                chromosome view and then clicking “File -> Current View” will bring up a download prompt where you can choose what to name the new file and where to download it to. 
                """
                ],
            ),
        ], width={"size": 6, "offset": 3}),
    ],
)

docs_content = html.Div(id="pdist-docs-content", className="pdist-docs-body")

navbar = dbc.Navbar(
    [
        dbc.Col(dbc.Button(
            id="home-button",
            children=[
                html.Img(src=app.get_asset_url("p-distance-logo-V7-71pt.svg"), height='40px'),
            ],
            className="logo-png",
            href="/",
            outline=True,
        ), width='auto'),
        dbc.Col(dbc.NavbarBrand("Signal Tracer", className="m-2"), width='auto'),
        dbc.Col([
            dbc.DropdownMenu(
                children=[
                    dbc.DropdownMenuItem("Import Options", header=True, className="file-menu-header"),
                    dbc.DropdownMenuItem("New Session", id="upload-file-button"),
                    dbc.DropdownMenuItem("GFF/GTF", id="pdist-gff-init-upload-button"),
                    dbc.DropdownMenuItem(divider=True),
                    dbc.DropdownMenuItem("Export Options", header=True, className="file-menu-header"),
                    dbc.DropdownMenuItem("Current View", id="pdist-curr-view-button"),
                ],
                className="nav-button",
                label="File",
                nav=True,
                in_navbar=True,
            ),
            # Toolbar
            dbc.NavItem(
                dbc.NavLink(
                    f"Toolbar",
                    id=f"pdist-toggle-toolbar",
                    className="nav-button",
                ),
                className="nav-button",
            ),
            # Graph Options
            dbc.NavItem(
                dbc.NavLink(
                    f"Graph Customization",
                    id=f"pdist-toggle-plot-options",
                    className="nav-button",
                ),
                className="nav-button",
            ),
            # Docs
            dbc.NavItem(
                dbc.NavLink(
                    f"Docs",
                    id=f"pdist-toggle-docs",
                    className="nav-button",
                ),
                className="nav-button",
            ),
        ]),
    ],
    color="orange",
    expand=True,
)

gff_input_button = dbc.InputGroup(
    [
        dbc.InputGroupAddon(
            children=[
                dcc.Upload(
                    id="gff-upload-data",
                    children=[
                        dbc.Button(
                            id="pdist-gff-upload-file-button",
                            children=["GFF/GTF file"], 
                            className="tv-input-buttons",
                            color="info"
                        )
                    ]
                ),
            ],
            addon_type="prepend",
        ),
        html.Div(id="pdist-gff-filename-div", className="tv-filename"),
    ],
    size='lg',
)

main_input_button = dbc.InputGroup(
    [
        dbc.InputGroupAddon(
            children=[
                dcc.Upload(
                    id='pd-upload-data',
                    children=[
                        dbc.Button(
                            id='pd-upload-file-button',
                            children=['Signal Tracer File'],
                            className="tv-input-buttons",
                            color='info'
                        ),
                    ],
                ),
            ],
            addon_type="prepend",
        ),
        html.Div(id="input-filename", children=["No file chosen..."], className="tv-filename"),
    ],
    size='lg',
)

view_toggle = daq.ToggleSwitch(
    id="view-option",
    disabled=True,
    value=False,
    label=[{'label': 'View | Chromosome ', 'style': {"font-size": '12pt'}}, {'label': 'Whole Genome', 'style': {"font-size": '12pt'}}],
    color="#375a7f",
    style={"color": "black"},
    size=25,
    theme='dark',
)


color_step_buttons = html.Div([
    dbc.Button("<", id='st-step-down-button', size='sm'),
    html.Div(["Shift"], style={'margin': '5px 5px 5px 5px', 'color': 'black'}),
    dbc.Button(">", id='st-step-up-button', size='sm'),
], style={'display': 'inline-flex'})

######################################  Layout  #####################################
def signalTracer_layout():
    layout = dbc.Container(
        id='container',
        children=[
            # Current View Download Component
            dcc.Download(id="pdist-current-view-download"),
            dcc.Download(id="pdist-gff-current-view-download"),
            dcc.Store(id="color-palette"),
            # Tool tips
            dbc.Tooltip(
                "Return to homepage",
                target="home-button",
                delay={"show": 500, "hide": 250}
            ),
            # Main Input Modal
            dbc.Modal(
                [
                    dbc.ModalHeader(
                        children=[
                            dbc.Row([
                                dbc.Label("New Session Upload", className="modal-label"),
                            ], justify="center", align="center")
                        ],
                        className='modal-header'
                    ),
                    dbc.ModalBody([dbc.Row(main_input_button, no_gutters=True)]),
                    dbc.ModalFooter(
                        children=[
                            dcc.Loading([html.Div(id='pdist-input-df', className="hidden-div")], className="toolbar-loading"),
                            dbc.Button("Submit", id="signal-submit-button",
                                       className="submit-buttons", color="info"),
                            dbc.Button("Close", id="signal-close-button",
                                       className="submit-buttons", color="info"),
                        ],
                    ),
                ],
                id="pdist-new-session-modal",
                size="lg",
                keyboard=False,
            ),
            # File Type Error Modal
            dbc.Modal(
                [
                    dbc.ModalHeader(
                        children=[
                            dbc.Label("Error: Invalid input file type", className='modal-label'),
                        ],
                        className='modal-header'
                    ),
                    dbc.ModalBody(
                        children=[
                            dbc.FormText(
                                R"Input file is not a valid file type for this tool. Please ensure the input file type is a csv, tsv, or xlsx file.", 
                                color="white",
                                className='modal-info-text'
                            ),
                        ],
                    ),
                    dbc.ModalFooter(
                        dbc.Button("Close", id="pdist-file-error-close-button", className="submit-buttons")
                    ),
                ],
                id="st-file-type-error-modal",
                size="xl",
                keyboard=False,
            ),
            # Header Error Modal
            dbc.Modal(
                [
                    dbc.ModalHeader(
                        children=[
                            # dbc.Label("Error: Invalid column headers"),
                            html.H2("Error: Invalid file headers!")
                        ],
                        className='modal-header'
                    ),
                    dbc.ModalBody(
                        children=[
                            dbc.FormText(
                                "Required column headers = Chromosome | Window | Sample | Value", 
                                color="white",
                                className='modal-info-text'
                            ),
                        ],
                    ),
                    dbc.ModalFooter(
                        dbc.Button("Close", id="pdist-file-error-close-button", className="submit-buttons", color="info")
                    ),
                ],
                id="st-header-error-modal",
                size="xl",
                keyboard=False,
            ),
            # Value Column Error Modal
            dbc.Modal(
                [
                    dbc.ModalHeader(
                        children=[
                            # dbc.Label("Error: Invalid column headers"),
                            html.H2("Error: Invalid data in value column")
                        ],
                        className='modal-header'
                    ),
                    dbc.ModalBody(
                        children=[
                            dbc.FormText(
                                "Ensure data in Value column is an integer or float value", 
                                color="white",
                                className='modal-info-text'
                            ),
                        ],
                    ),
                    dbc.ModalFooter(
                        dbc.Button("Close", id="pdist-file-error-close-button", className="submit-buttons", color="info")
                    ),
                ],
                id="st-value-col-error-modal",
                size="xl",
                keyboard=False,
            ),
            # Alternative Input Files
            dbc.Modal(
                children=[
                    dbc.ModalHeader(
                        children=[
                            dbc.Row([
                                dbc.Label("Alternative Input Files", className="modal-label"),
                            ], justify="center", align="center"),
                        ],
                        className="modal-header"
                    ),
                    dbc.ModalBody([dbc.Row(gff_input_button, no_gutters=True)]),
                    dbc.ModalFooter(
                        children=[
                            dcc.Loading([html.Div(id="pdist-gff-data-upload", className="hidden-div")], className="toolbar-loading"),
                            dbc.Button("Submit", id="pdist-gff-submit-button",
                                    className="submit-buttons", color="secondary")
                        ],
                    ),
                ],
                id="pdist-alt-input-modal",
                size="lg",
            ),
            # Navbar
            dbc.Row(
                children=[navbar],
                className="navbar-div",
            ),
            dbc.Alert(
                "Current View export is not available for whole genome view",
                id="alert-wg-cv",
                is_open=False,
                duration=4000,
                color='danger'
            ),
            # Toolbar
            dbc.Row(
                children=[
                    # Docs
                    dbc.Card(
                        children=[
                            dbc.Collapse(
                                dbc.CardBody(
                                    children=[
                                        dbc.Row(
                                            children=[
                                                dbc.Col(
                                                    children=[docs_content],
                                                    style={"height": "50vh", "overflow-y": "auto"}
                                                ),
                                            ],
                                            no_gutters=True,
                                        ),
                                    ],
                                    className="narbar-pdist-collapse-toolbar",
                                ),
                                id="pdist-collapse-documentation",
                                className="narbar-pdist-collapse-toolbar",
                            ),
                        ],
                        id="pdist-collapse-toolbar-card",
                        className="narbar-card"
                    ),
                    # Graph Options
                    dbc.Card(
                        children=[
                            dbc.Collapse(
                                dbc.CardBody(
                                    children=[
                                        dbc.Row(
                                            children=[
                                                # Graph Theme
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            children=['Color Theme:'],
                                                            className="text-divs",
                                                        ),
                                                        dcc.Dropdown(
                                                            id="template-option",
                                                            className='dropdown-style',
                                                            options=[{'label': i, 'value': i}
                                                                    for i in GRAPH_TEMPLATES],
                                                            value="plotly_dark",
                                                            style={
                                                                "color": "black"
                                                            },
                                                            clearable=False,
                                                        ),
                                                    ],
                                                    className='sm-margin',
                                                    width=2,
                                                ),
                                                # --- Font Family ---
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            "Font Family:",
                                                            style={'color': 'black'},
                                                        ),
                                                        dcc.Dropdown(
                                                            id="st-font-family-option",
                                                            options=[{"label": i, "value": i} for i in FONT_FAMILIES],
                                                            value="Arial",
                                                            clearable=False,
                                                            className="dropdown-style",
                                                        ),
                                                    ],
                                                    className='sm-margin',
                                                    width=1,
                                                ),
                                                # Snapshot fileType
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            children=['Snapshot File Type:'],
                                                            className="text-divs",
                                                        ),
                                                        dcc.Dropdown(
                                                            id="snapshot-file-option",
                                                            className='dropdown-style',
                                                            options=[{'label': i, 'value': i}
                                                                    for i in SNAPSHOT_FILE_OPTIONS],
                                                            value="svg",
                                                            style={
                                                                "color": "black"
                                                            },
                                                            clearable=False,
                                                        ),
                                                    ],
                                                    className='sm-margin',
                                                    width=1,
                                                ),
                                                # --- Axis line width ---
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            children=[
                                                                "Axis line width:"],
                                                            className="text-divs",
                                                        ),
                                                        dcc.Dropdown(
                                                            id="pdist-axis-line-width",
                                                            options=[{"label": i, "value": i} for i in range(1, 6)],
                                                            value=1,
                                                            clearable=False,
                                                            className="dropdown-style",
                                                        ),
                                                    ],
                                                    className="sm-margin",
                                                    width=1,
                                                ),
                                                # --- Pixel output size ---
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            "Output Dimensions: (px)",
                                                            className="text-divs",
                                                        ),
                                                        html.Div([
                                                            dcc.Input(id="tv-pixel-width", style={'width': '50%'}, type="number", placeholder="Width"),
                                                            html.Div(["x"], style={"margin": "0px 5px 0px 5px"}),
                                                            dcc.Input(id="tv-pixel-height", style={'width': '50%'}, type="number", placeholder="Height"),
                                                            
                                                        ], className="pixel-input-style",),
                                                    ],
                                                    className="sm-margin",
                                                    width=2,
                                                ),
                                            ],
                                            style={'background': 'orange', 'border-radius': '5px 5px 0px 0px'},
                                            no_gutters=True,
                                        ),
                                        dbc.Row(
                                            children=[
                                                # Line Colors
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            children=['Line Colors:'],
                                                            className="text-divs",
                                                        ),
                                                        html.Div([
                                                            dcc.Dropdown(
                                                                id="line-color-option",
                                                                className='dropdown-style',
                                                                options=[{'label': f"{i} - {len(COLOR_SWATCHES[i])} colors", 'value': i}
                                                                        for i in COLOR_SWATCHES],
                                                                value="Plotly",
                                                                style={'color': 'black', 'width': '98%'},
                                                                clearable=False,
                                                            ),
                                                            color_step_buttons,
                                                        ], style={'display': 'inline-flex', 'width': '100%'},)
                                                    ],
                                                    className='sm-margin',
                                                    width=2,
                                                ),
                                                # Font Size
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            children=['Font Size:'],
                                                            className="text-divs",
                                                        ),
                                                        dcc.Dropdown(
                                                            id='font_size_options',
                                                            className='dropdown-style',
                                                            options=graph_options.font_size_options(),
                                                            value=12,
                                                            style={'color': 'black'},
                                                            clearable=False,
                                                        ),
                                                    ],
                                                    className='sm-margin',
                                                    width=1,
                                                ),
                                                # Snapshot Scale
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            children=['Snapshot Scale:'],
                                                            className="text-divs",
                                                        ),
                                                        dcc.Dropdown(
                                                            id="snapshot-scale",
                                                            className='dropdown-style',
                                                            options=[{'label': i, 'value': i}
                                                                    for i in SCALE_OPTIONS],
                                                            value=1,
                                                            style={
                                                                "color": "black"
                                                            },
                                                            clearable=False,
                                                        ),
                                                    ],
                                                    className='sm-margin',
                                                    width=1,
                                                ),
                                                # Marker Width
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            children=['Marker Width:'],
                                                            className="text-divs",
                                                            style={
                                                                "color": "black",
                                                            },
                                                        ),
                                                        dcc.Input(
                                                            id="marker_width_options",
                                                            type="number", value=1.0,
                                                            style={'width': '100%', 'height': '50%'},
                                                        ),
                                                    ],
                                                    className='sm-margin',
                                                    width=1,
                                                ),
                                                # --- Axis Grid Lines ---
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            children=["Grid Lines:"],
                                                            className="text-divs",
                                                        ),
                                                        dbc.Checklist(
                                                            id="axis-gridlines",
                                                            options=[{'label': 'X-axis', 'value': 'xaxis'},
                                                                        {'label': 'Y-axis', 'value': 'yaxis'}],
                                                            value=["xaxis"],
                                                            switch=True,
                                                            inline=True,
                                                            className="topo-organization",
                                                        ),
                                                    ],
                                                    className="sm-margin",
                                                    width=2,
                                                ),
                                                # --- Edit Graphs ---
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            children=["Editable Graphs:"],
                                                            className="text-divs",
                                                        ),
                                                        html.Div(
                                                            children=[
                                                                html.Div(["Off"]),
                                                                dbc.Checklist(
                                                                    id="editable-graph-option",
                                                                    options=[
                                                                        {"label": "", "value": "editable_graphs_switch"},
                                                                    ],
                                                                    switch=True,
                                                                    className="topo-organization",
                                                                ),
                                                                html.Div(["On"]),
                                                            ],
                                                            className="text-divs ilf",
                                                        ),
                                                    ],
                                                    className="sm-margin",
                                                    width=1,
                                                    align='center',
                                                ),
                                            ],
                                            style={'background': 'orange', 'border-radius': '0px 0px 5px 5px', 'margin-bottom': '10px'},
                                            no_gutters=True,
                                        ),
                                    ],
                                    className="narbar-pdist-collapse-toolbar",
                                ),
                                id=f"pdist-collapse-plot-options",
                                className="narbar-collapse-toolbar",
                            ),
                        ],
                        id="pdist-collapse-toolbar-card",
                        className="narbar-card"
                    ),
                    # Toolbar
                    dbc.Card(
                        children=[
                            dbc.Collapse(
                                dbc.CardBody(
                                    dbc.Row(
                                        children=[
                                            # Graph type
                                            dbc.Col(
                                                style={'padding': '5px', 'display': 'inline-flex'},
                                                align='center',
                                                children=[
                                                    # Graph type
                                                    html.Div(
                                                        id='graph-type-div',
                                                        children=[
                                                            html.H6(
                                                                children=["Graph Type | "],
                                                                style={'color': 'black', "padding-bottom": "2px"},
                                                            ),
                                                            dcc.Dropdown(
                                                                id='graph-type-dropdown',
                                                                options=[
                                                                    {'label': 'Line', 'value': 'Line'},
                                                                    {'label': 'Scatter', 'value': 'Scatter'},
                                                                ],
                                                                value='Line',
                                                                className='dropdown-style',
                                                                clearable=False,
                                                            ),
                                                        ],
                                                        style={'padding': '5px', 'width': '10%'},
                                                    ),
                                                    # Chromosome choice
                                                    html.Div(
                                                        id='chromosome-choice-div',
                                                        children=[
                                                            html.H6(
                                                                children=[view_toggle],
                                                                style={'color': 'black', "padding-bottom": "2px"},
                                                            ),
                                                            dcc.Dropdown(
                                                                id='pdist-chromosome-options',
                                                                className='dropdown-style',
                                                            ),
                                                        ],
                                                        style={'padding': '5px'},
                                                    ),
                                                    # Graph Switches
                                                    html.Div(
                                                        children=[
                                                            html.H6("Graph Switches | ", style={'color': 'black', "padding-bottom": "3px"}),
                                                            dbc.FormGroup(
                                                                children=[
                                                                    dbc.Checklist(
                                                                        id="pdist-alt-graph-switches",
                                                                        options=[
                                                                            {"label": "GFF/GTF", "value": "gff_switch"},
                                                                        ],
                                                                        switch=True,
                                                                        className="checklist-style",
                                                                        inline=True,
                                                                    ),
                                                                ],
                                                                inline=True,
                                                            )
                                                        ],
                                                        style={'padding': '5px'},
                                                    ),
                                                ],
                                            ),
                                        ],
                                        style={'background': 'orange', 'border-radius': '5px', 'margin-bottom': '10px'},
                                        no_gutters=True,
                                    ),
                                    className="narbar-pdist-collapse-toolbar",
                                ),
                                id=f"pdist-collapse-toolbar",
                                className="narbar-collapse-toolbar",
                            ),
                        ],
                        id="pdist-collapse-toolbar-card",
                        className="narbar-card"
                    ),
                ],
                className="narbar-pdist-collapse-toolbar",
                no_gutters=True,
            ),
            # Graph
            dbc.Row(
                # id='body',
                children=[
                    # Graphs
                    dbc.Col(
                        id='graph-column',
                        width=12,
                        children=[
                            dcc.Loading(
                                id="loading-1",
                                children=[
                                    dcc.Graph(
                                        id="pdist-graph",
                                        figure=(tree_utils.init_data_graph('plotly_dark')),
                                        config=dict(
                                            doubleClick="autosize",
                                            displaylogo=False,
                                            displayModeBar=False,
                                            toImageButtonOptions=dict(
                                                format='svg',
                                                filename="Graph_Name",
                                            ),  
                                        ),
                                        style={"padding": "2px"}
                                    )
                                ],
                            ),
                        ],
                    ),
                    dbc.Col(
                        id='gff-col',
                        className="hidden-div",
                        width=12,
                        children=[
                            dcc.Loading(
                                children=[
                                    html.Div(id='pdist-gff-graph')
                                ],
                            ),
                        ],
                    ),
                ],
                className='pdist-row-divs',
                no_gutters=True
            ),
        ],
        fluid=True,
    )
    return layout


##################################### Callbacks #####################################
@app.callback(
    [Output(f"pdist-collapse-toolbar", "is_open"),
     Output(f"pdist-collapse-documentation", "is_open"),
     Output(f"pdist-collapse-plot-options", "is_open"),],
    [Input(f"pdist-toggle-toolbar", "n_clicks"),
     Input(f"pdist-toggle-docs", "n_clicks"),
     Input(f"pdist-toggle-plot-options", "n_clicks"),],
    [State(f"pdist-collapse-toolbar", "is_open"),
     State(f"pdist-collapse-documentation", "is_open"),
     State(f"pdist-collapse-plot-options", "is_open"),],
)
def toggle_toolbar(n1, n2, n3, toolbar_is_open, docs_is_open, plot_options_is_open):
    ctx = dash.callback_context

    if not ctx.triggered:
        return True, False, False
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if button_id == "pdist-toggle-toolbar" and n1:
        return not toolbar_is_open, False, False
    elif button_id == "pdist-toggle-docs" and n2:
        return False, not docs_is_open, False
    elif button_id == "pdist-toggle-plot-options" and n3:
        return False, False, not plot_options_is_open
    else:
        return False, False, False


@app.callback(
    Output("pdist-docs-content", "children"),
    [Input("pdist-collapse-documentation", "is_open"),],)
def render_page_content(is_open):
    """Opens/closes Docs toolbar when menu button is pressed"""
    if is_open:
        page_content = html.P(
            children=[docs_pdist_contents],
        )
        return page_content
    else:
        raise PreventUpdate

# ------------------------------------------------------------------------------------
################################# Set display options ################################
@app.callback(
    Output("gff-col", "className"),
    [Input("pdist-alt-graph-switches", "value"),],
)
def hide_gff_div(altSwitches):
    if not altSwitches:
        return "hidden-div"
    elif "gff_switch" in altSwitches:
        return "visible-div"
    else:
        return "hidden-div"


@app.callback(
    [Output("pdist-alt-graph-switches", "options"),
     Output("pdist-alt-graph-switches", "value"),],
    [Input("pdist-gff-data-upload", "children"),
     Input("pdist-input-df", "children"),
     Input('pdist-chromosome-options', 'value'),]
)
def main_graph_button_options(gffData, pdistData, chrom_choice):
    if not gffData:
        gffOption = {"label": "GFF/GTF - No file chosen...", "value": "gff_switch", "disabled": True}
        return [gffOption], []
    elif not pdistData:
        gffOption = {"label": "GFF/GTF", "value": "gff_switch", "disabled": True}
        return [gffOption], []
    elif chrom_choice == "Whole Genome View":
        gffOption = {"label": "GFF/GTF", "value": "gff_switch", "disabled": True}
        return [gffOption], []
    else:
        gffOption = {"label": "GFF/GTF", "value": "gff_switch"}
        return [gffOption], dash.no_update
 

@app.callback(
    [Output("pdist-chromosome-options", 'disabled'),
     Output("view-option", "disabled"),
    ],
    [Input("pdist-input-df", "children"),],
)
def disable_components(topo_data):
    if not topo_data:
        chromOptions = True
        viewOption = True
        return chromOptions, viewOption
    else:
        chromOptions = False
        viewOption = False
        return chromOptions, viewOption

# -----------------------------------------------------------------------------------
################################### File Upload #####################################
@app.callback(
    [Output("pdist-alt-input-modal", "is_open"),
     Output("pdist-gff-data-upload", "children"),
     Output("pdist-gff-filename-div", "children"),
     Output("pdist-gff-upload-file-button", "color"),],
    [Input("pdist-gff-init-upload-button", "n_clicks"),
     Input("pdist-gff-submit-button", "n_clicks"),
     Input("signal-submit-button", "n_clicks"),
     Input("gff-upload-data", "contents"), ],
    [State("gff-upload-data", "filename")])
def upload_gff_file(gffBtn, gffSubmit, pdistBtn, gffContents, gffFilename):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "No clicks yet"
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if button_id == "pdist-gff-init-upload-button":
        if not gffFilename:
            gffFilename = "No file chosen..."
        return True, None, gffFilename, "primary"
    elif button_id == "pdist-gff-upload-file-button":
        return True, None, None, "primary"
    elif button_id == "signal-submit-button":
        return False, None, None, "primary"
    elif button_id == "pdist-gff-submit-button":
        if not gffContents:
            return True, None, gffFilename, "primary"
        _, content_string = gffContents.split(",")
        tvDecoded = base64.b64decode(content_string)
        if ("gff" in gffFilename) or ("gtf" in gffFilename):
            gffDF = pd.read_csv(
                io.StringIO(tvDecoded.decode("utf-8")),
                sep="\t",
                names=["chromosome", "source", "feature", "start",
                       "end", "score", "strand", "frame", "attribute"],
                comment="#",
            )
            gffDF["attribute"] = gffDF["attribute"].apply(lambda x: data_utils.get_gene_name(x))
            return False, gffDF.to_json(), gffFilename, "success"
        else:
            return True, None, "Please provide a valid GFF/GTF file"
    elif gffFilename:
        if tree_utils.validate_gff_gtf_filename(gffFilename):
            return True, None, gffFilename, "success"
        else:
            return True, None, f"File: {gffFilename}", "danger"
    else:
        raise PreventUpdate


# Upload input file to DataFrame
@app.callback([Output('pdist-input-df', 'children'),
               Output('st-file-type-error-modal', 'is_open'),
               Output("pdist-new-session-modal", "is_open"),
               Output("st-header-error-modal", "is_open"),
               Output("st-value-col-error-modal", "is_open"),
               Output("input-filename", 'children'),
               Output("pd-upload-file-button", "color"),
            ],
              [Input('pdist-file-error-close-button', 'n_clicks'),
               Input("upload-file-button", "n_clicks"),
               Input("signal-submit-button", "n_clicks"),
               Input("signal-close-button", "n_clicks"),
               Input('pd-upload-data', 'contents'),],
              [State('pd-upload-data', 'filename'),],
)
def load_input_files(btn1, btn2, btn3, btn4, fileContents, pdistFilename):
    ctx = dash.callback_context
    # Set button_id
    if not ctx.triggered:
        button_id = 'No clicks yet'
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    # Set button color
    if not pdistFilename:
        pdist_color = "primary"
    elif 'txt' in pdistFilename:
        pdist_color = "success"
    elif 'csv' in pdistFilename:
        pdist_color = "success"
    elif 'tsv' in pdistFilename:
        pdist_color = "success"
    elif 'xls' in pdistFilename:
        pdist_color = "success"
    # Handle input files
    if button_id == 'pdist-file-error-close-button':
        return None, False, True, False, False, [pdistFilename], "primary"
    elif button_id == 'upload-file-button':
        if not pdistFilename:
            pdistFilename = "No file chosen..."
        return dash.no_update, False, True, False, False, [pdistFilename], pdist_color
    elif button_id == 'signal-close-button':
        return dash.no_update, False, False, False, False, [pdistFilename], pdist_color
    elif not fileContents:
        if not pdistFilename:
            pdistFilename = "No file chosen..."
        return None, False, True, False, False, [pdistFilename], pdist_color
    elif button_id == "signal-submit-button":
        _, tvcontent_string = fileContents.split(',')
        tvDecoded = base64.b64decode(tvcontent_string)
        if '.csv' in pdistFilename:
            # Open csv file
            csv_df = pd.read_csv(io.StringIO(tvDecoded.decode('utf-8')))
            validate_header = signal_utils.validate_signal_tracer_headers(csv_df)
            validate_values = signal_utils.validate_signal_tracer_values(csv_df)
            if not validate_header:
                return None, False, False, True, False, [pdistFilename], "danger"
            elif not validate_values:
                return None, False, False, False, True, [pdistFilename], "danger"
            else:
                return csv_df.to_json(), False, False, False, False, [pdistFilename], "success"
        elif '.tsv' in pdistFilename:
            # Open tsv file
            tsv_df = pd.read_csv(io.StringIO(tvDecoded.decode('utf-8')), sep='\t')
            validate_header = signal_utils.validate_signal_tracer_headers(tsv_df)
            validate_values = signal_utils.validate_signal_tracer_values(tsv_df)
            if not validate_header:
                return None, False, False, True, False, [pdistFilename], "danger"
            elif not validate_values:
                return None, False, False, False, True, [pdistFilename], "danger"
            else:
                return tsv_df.to_json(), False, False, False, False, [pdistFilename], "success"
        elif '.txt' in pdistFilename:
            # Open txt file
            txt_df = pd.read_csv(io.StringIO(tvDecoded.decode('utf-8')), sep='\t')
            validate_header = signal_utils.validate_signal_tracer_headers(txt_df)
            validate_values = signal_utils.validate_signal_tracer_values(txt_df)
            if not validate_header:
                return None, False, False, True, False, [pdistFilename], "danger"
            elif not validate_values:
                return None, False, False, False, True, [pdistFilename], "danger"
            else:
                return txt_df.to_json(), False, False, False, False, [pdistFilename], "success"
        elif 'xls' in pdistFilename:
            # Open excel file
            xlsx_df = pd.read_excel(io.BytesIO(tvDecoded))
            validate_header = signal_utils.validate_signal_tracer_headers(xlsx_df)
            validate_values = signal_utils.validate_signal_tracer_values(xlsx_df)
            if not validate_header:
                return None, False, False, True, False, [pdistFilename], "danger"
            elif not validate_values:
                return None, False, False, False, True, [pdistFilename], "danger"
            else:
                return xlsx_df.to_json(), False, False, False, False, [pdistFilename], "success"
        else:
            return None, True, False, False, False, [pdistFilename], "danger"
    else:
        valid_file_type = signal_utils.validate_file_type(pdistFilename)
        if valid_file_type:
            return dash.no_update, False, dash.no_update, False, False, [pdistFilename], "success"
        else:
            return dash.no_update, False, dash.no_update, False, False, [pdistFilename], "danger"

# -----------------------------------------------------------------------------------
######################### Set option dropdown menu selections #######################
@app.callback(
    [Output('pdist-chromosome-options', 'options'),
     Output('pdist-chromosome-options', 'value'),
    ],
    [Input('view-option', 'value'),
     Input('pdist-input-df', 'children'),
     Input("pdist-new-session-modal", "is_open"),
    ]
)
def set_chromosome_choice(chrom_quantity, fileData, input_modal):
    if not fileData:
        raise PreventUpdate
    elif input_modal:
        raise PreventUpdate
    elif chrom_quantity:
        return [dict(label='Whole Genome View', value='Whole Genome View')], 'Whole Genome View'
    else:
        read_file = pd.read_json(fileData)
        options = [dict(label=chrom, value=chrom) for chrom in sorted(read_file['Chromosome'].unique())]
        return options, options[0]['value']


@app.callback(
    Output("color-palette", "data"),
    [Input('pdist-input-df', 'children'),
     Input("line-color-option", 'value'),
     Input("color-palette", "data"),
     Input("st-step-up-button", "n_clicks"),
     Input("st-step-down-button", "n_clicks"),
    ],
)
def set_line_colors(
    data,
    color_option,
    curr_colors,
    step_up_btn,
    step_down_btn,
):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if not data:
        return dash.no_update
    elif not curr_colors:
        color_set = COLOR_SWATCHES[color_option]
        return color_set
    elif button_id == 'st-step-up-button':
        step_color_set = curr_colors[1:] + [curr_colors[0]]
        return step_color_set
    elif button_id == 'st-step-down-button':
        step_color_set = [curr_colors[-1]] + curr_colors[:-1]
        return step_color_set
    else:
        color_set = COLOR_SWATCHES[color_option]
        return color_set


#----------------------------------------------------------------------------------------------------
### Graphing CB's ###
@app.callback(
    [Output('pdist-graph', 'figure'),
     Output('pdist-graph', 'config'),
     Output("marker_width_options", "value"),
    ],
    [Input("pdist-new-session-modal", "is_open"),
     Input('pdist-input-df', 'children'),
     Input('pdist-chromosome-options', 'value'),
     Input('template-option', 'value'),
     Input("snapshot-file-option", "value"),
     Input("marker_width_options", "value"),
     Input("color-palette", "data"),
     Input("tv-pixel-height", "value"),
     Input("tv-pixel-width", "value"),
     Input("snapshot-scale", "value"),
     Input("font_size_options", "value"),
     Input("axis-gridlines", "value"),
     Input("pdist-axis-line-width", "value"),
     Input("editable-graph-option", "value"),
     Input("view-option", "value"),
     Input("graph-type-dropdown", "value"),
     Input("st-font-family-option", "value"),
    ]
)
def update_main_graph(
    input_modal,
    currData,
    chromosome, 
    template,
    snapshot_file_type,
    marker_width,
    line_color,
    pixel_height,
    pixel_width,
    snapshot_scale,
    font_size,
    axis_gridlines,
    axis_marker_width,
    editable_graphs,
    view,
    graph_type,
    font_family,
):
    # Prevent update if True
    if not currData:
        raise PreventUpdate
    elif input_modal:
        raise PreventUpdate
    elif not line_color:
        raise PreventUpdate
    # Set editable value
    if not editable_graphs:
        editable_graphs = False
    else:
        editable_graphs = True
    # Load in data and clean data
    df = pd.read_json(currData)
    df.columns = ["Chromosome", "Window", "Sample", "Value"]
    df = df.sort_values(by=["Chromosome", "Window"])
    # Get sample and chromosome data
    chromosomes = [i for i in df['Chromosome'].unique()]
    samples = [i for i in df['Sample'].unique()]
    samples.sort(reverse=True)

    # Get overall max Signal value to set y-axis ranges for all plots
    y_max = float(df["Value"].max())
    x_max = df["Window"].max()

    # Add colors 
    try:
        assert len(line_color) > len(samples)
        colors = line_color
    except AssertionError:
        colors = line_color + line_color[:(len(samples)-len(line_color))]
    # Set gridline bools
    xaxis_gridlines, yaxis_gridlines = tree_utils.get_gridline_bools(axis_gridlines)

    if view:
        if graph_type == "Line":
            fig = signal_utils.whole_genome_line(
                df,
                chromosomes,
                samples,
                colors,
                marker_width,
                template,
                font_size,
                y_max,
                x_max,
                xaxis_gridlines,
                yaxis_gridlines,
                font_family,
            )
        elif graph_type == "Scatter":
            fig = signal_utils.whole_genome_scatter(
                df,
                chromosomes,
                samples,
                colors,
                marker_width,
                template,
                font_size,
                y_max,
                x_max,
                xaxis_gridlines,
                yaxis_gridlines,
                font_family,
            )
        # Config pixel size
        if (not pixel_height) or (not pixel_width):
            config = dict(
                doubleClick="autosize",
                displaylogo=False,
                editable=editable_graphs,
                edits=dict(
                    annotationText=True,
                    axisTitleText=False,
                    titleText=False,
                ),
                modeBarButtonsToRemove=[
                    "zoom2d",
                    "autoScale2d",
                    "resetScale2d",
                    "pan2d",
                    "lasso2d",
                    "select2d",
                ],
                toImageButtonOptions=dict(
                    format=snapshot_file_type,
                    filename="Graph_Name",
                ),
            )
        else:
            config = dict(
                doubleClick="autosize",
                displaylogo=False,
                editable=editable_graphs,
                edits=dict(
                    annotationText=True,
                    axisTitleText=False,
                    titleText=False,
                ),
                modeBarButtonsToRemove=[
                    "zoom2d",
                    "autoScale2d",
                    "resetScale2d",
                    "pan2d",
                    "lasso2d",
                    "select2d",
                ],
                toImageButtonOptions=dict(
                    format=snapshot_file_type,
                    filename="Graph_Name",
                    width=int(pixel_width),
                    height=int(pixel_height),
                ),
            )
        return fig, config, marker_width
    else:
        if graph_type == "Line":
            fig = signal_utils.single_chromosome_graph_line(
                df,
                chromosome,
                template,
                marker_width,
                colors,
                font_size,
                xaxis_gridlines,
                yaxis_gridlines,
                font_family,
                samples,
            )
        elif graph_type == "Scatter":
            fig = signal_utils.single_chromosome_graph_scatter(
                df,
                chromosome,
                template,
                marker_width,
                colors,
                font_size,
                xaxis_gridlines,
                yaxis_gridlines,
                font_family,
                samples,
            )
        # Config pixel size
        if (not pixel_height) or (not pixel_width):
            pixel_width, pixel_height = 750, 400
            
        config = dict(
            doubleClick="autosize",
            displaylogo=False,
            editable=editable_graphs,
            edits=dict(
                annotationText=True,
                titleText=False,
            ),
            modeBarButtonsToRemove=[
                "zoom2d",
                "autoScale2d",
                "resetScale2d",
                "pan2d",
                "lasso2d",
                "select2d",
            ],
            toImageButtonOptions=dict(
                format=snapshot_file_type,
                filename="Graph_Name",
                width=int(pixel_width),
                height=int(pixel_height),
            ),
        )
        return fig, config, marker_width


@app.callback(
    Output("pdist-gff-graph", "children"),
    [Input("pdist-new-session-modal", "is_open"),
    Input('pdist-input-df', 'children'),
    Input("pdist-gff-data-upload", "children"),
    Input("pdist-alt-graph-switches", "value"),
    Input('pdist-chromosome-options', 'value'),
    Input('template-option', 'value'),
    Input("snapshot-file-option", "value"),
    Input("line-color-option", 'value'),
    Input("axis-gridlines", "value"),
    Input("pdist-axis-line-width", "value"),
    Input('pdist-graph', 'relayoutData'),
    Input("editable-graph-option", "value"),
    ]
)
def update_gff_graph(
    input_modal,
    currData,
    gff_data,
    altSwitches,
    chromosome, 
    template,
    snapshot_file_type,
    line_color,
    axis_gridlines,
    axis_marker_width,
    relayout_data,
    editable_graphs,
):
    # Prevent update if True
    if not currData:
        raise PreventUpdate
    elif input_modal:
        raise PreventUpdate
    # Set editable value
    if not editable_graphs:
        editable_graphs = False
    else:
        editable_graphs = True
    # Load in data and clean data
    df = pd.read_json(currData)
    # df = df.melt(id_vars=["Chromosome", "Window"])
    df.columns = ["Chromosome", "Window", "Sample", "Value"]
    df = df.sort_values(by=["Chromosome", "Window"])
    # Get sample and chromosome data
    samples = [i for i in df['Sample'].unique()]
    samples.sort(reverse=True)

    # Get overall max Signal value to set y-axis ranges for all plots
    x_max = df["Window"].max()

    # Set gridline bools
    xaxis_gridlines, yaxis_gridlines = tree_utils.get_gridline_bools(axis_gridlines)

    try:
        relayout_keys = [k for k in relayout_data.keys()]
        if "xaxis.range" in relayout_keys[0]:
            dataRange = [relayout_data['xaxis.range[0]'], relayout_data['xaxis.range[1]']]
            dataMin, dataMax = [relayout_data['xaxis.range[0]'], relayout_data['xaxis.range[1]']]
        elif "xaxis2.range" in relayout_keys[0]:
            dataRange = [relayout_data['xaxis2.range[0]'], relayout_data['xaxis2.range[1]']]
            dataMin, dataMax = [relayout_data['xaxis2.range[0]'], relayout_data['xaxis2.range[1]']]
        else:
            raise KeyError
    except KeyError:
        dataRange = [0, x_max]
        dataMin, dataMax = 0, x_max
        pass
    except TypeError:
        dataRange = [0, x_max]
        dataMin, dataMax = 0, x_max
        pass
    except IndexError:
        dataRange = [0, x_max]
        dataMin, dataMax = 0, x_max
        pass

    if not altSwitches:
        gff_graph = None
    elif "gff_switch" in altSwitches:
        try:
            gff_df = pd.read_json(gff_data)
        except ValueError:
            raise PreventUpdate
        gff_df = gff_df[(gff_df["start"] >= dataMin) & (gff_df["end"] <= dataMax) & (gff_df["chromosome"] == chromosome)]
        if abs(dataMin - dataMax) > 50000000:
            gff_graph = [
                html.Div(
                    children=[
                        dcc.Graph(
                            id="gff_graph",
                            figure=tree_utils.zoom_in_gff(template),
                            className="gff-style",
                            config=dict(
                                displayModeBar=False,
                            ),
                        ),
                    ],
                )
            ]
        elif len(gff_df) == 0:
            gff_graph = [
                html.Div(
                    children=[
                        dcc.Graph(
                            id="gff_graph",
                            figure=tree_utils.no_tree_data(template, "No data in region"),
                            className="gff-style",
                            config=dict(
                                displayModeBar=False,
                            ),
                        ),
                    ],
                )
            ]
        else:
            gff_figure = tree_utils.build_gff_figure(
                gff_df,
                dataRange,
                template,
                axis_marker_width,
                xaxis_gridlines,
                yaxis_gridlines,
                "Arial",
            )
            gff_graph = [
                html.Div(
                    children=[
                        dcc.Graph(
                            id="gff_graph",
                            figure=gff_figure,
                            config=dict(
                                displaylogo=False,
                                doubleClick="autosize",
                                editable=editable_graphs,
                                modeBarButtonsToRemove=[
                                    "resetScale2d",
                                    "lasso2d",
                                    "select2d",
                                ],
                                toImageButtonOptions=dict(
                                    format=snapshot_file_type,
                                    filename="Graph_Name",
                                    height=200,
                                    width=1500,
                                ),
                            ),
                            className="gff-style",
                        ),
                    ],
                )
            ]

    return gff_graph

#----------------------------------------------------------------------------------------------------
# Current View Export
@app.callback(
    [Output("pdist-current-view-download", "data"),
     Output("pdist-gff-current-view-download", "data"),
     Output("alert-wg-cv", "is_open")],
    [Input("pdist-curr-view-button", "n_clicks")],
    [State("pdist-input-df", "children"),
     State("pdist-gff-data-upload", "children"),
     State("pdist-chromosome-options", "value"),
     State('pdist-graph', 'relayoutData'),],
    prevent_initial_call=True,
)
def current_window_summary(n, data_json, gff_json, current_chromosome, relayout_data):
    if not data_json:
        raise PreventUpdate
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if current_chromosome == "Whole Genome View":
        return None, None, True
    elif not button_id:
        raise PreventUpdate
    elif button_id == "pdist-curr-view-button":
        df = pd.read_json(data_json)
        relayout_keys = [k for k in relayout_data.keys()]
        if "xaxis.range" in relayout_keys[0]:
            dataMin, dataMax = [relayout_data['xaxis.range[0]'], relayout_data['xaxis.range[1]']]
            df = df[(df["Window"] >= dataMin) & (df["Window"] <= dataMax) & (df["Chromosome"] == current_chromosome)]   
        else:
            df = df[df["Chromosome"] == current_chromosome]
        
        if not gff_json:
            gff_df = pd.DataFrame()
        else:
            gff_df = pd.read_json(gff_json)
            gff_df = gff_df[(gff_df["start"] >= dataMin) & (gff_df["start"] <= dataMax) & (gff_df["chromosome"] == current_chromosome)]
        
        # Return according to available data
        if (len(df) == 0) and (len(gff_df) > 1):
            return None, dcc.send_data_frame(gff_df.to_csv, f"current_view_{current_chromosome}_{int(dataMin)}_{int(dataMax)}.gff", sep="\t", index=False), False
        elif (len(df) > 1) and (len(gff_df) == 0):
            return dcc.send_data_frame(df.to_excel, f"current_view_{current_chromosome}_{int(dataMin)}_{int(dataMax)}.xlsx", index=False), None, False
        else:
            return dcc.send_data_frame(df.to_excel, f"current_view_{current_chromosome}_{int(dataMin)}_{int(dataMax)}.xlsx", index=False), dcc.send_data_frame(gff_df.to_csv, f"current_view_{current_chromosome}_{int(dataMin)}_{int(dataMax)}.gff", sep="\t", index=False, header=False), False
    else:
        raise PreventUpdate

