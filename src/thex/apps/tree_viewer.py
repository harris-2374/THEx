import base64
import io

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt
# import numpy as np
import pandas as pd
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from ete3 import Tree
from ete3.coretype.tree import TreeError

from thex.app import app
from thex.apps.utils import tree_utils, data_utils, graph_options
from thex.apps.docs import docs_treeviewer_contents

# --- Graph options ---
GRAPH_TEMPLATES = graph_options.graph_templates()
COLOR_SWATCHES = graph_options.color_swatches()
OUTPUT_PIXEL_SIZES = graph_options.pixel_sizes()
SCALE_OPTIONS = graph_options.figure_output_scales()
SNAPSHOT_FILE_OPTIONS = graph_options.snapshot_file_type()

############################### Documentation Components ###############################
    
docs_content = html.Div(id="tv-docs-content", className="tv-docs-body")

navbar = dbc.Navbar(
    [
        html.A(
            # Use row and col to control vertical alignment of logo / brand
            dbc.Row(
                [
                    dbc.Col(dbc.Button(
                        id="home-button",
                            children=[
                                html.Img(src=app.get_asset_url("tree5.svg"), height='40px'),
                            ],
                            className="logo-png",
                            href="/",
                            outline=True,
                        ), width='auto'),
                    dbc.Col(dbc.NavbarBrand("Tree Viewer", className="m-2")),
                    
                ],
                align="center",
                no_gutters=True,
                className="nav-toolbar",
            ),
        ),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("Import Options", header=True, className="file-menu-header"),
                dbc.DropdownMenuItem("New Session", id="upload-file-button"),
                dbc.DropdownMenuItem("GFF/GTF", id="gff-init-upload-button"),
                dbc.DropdownMenuItem(divider=True),
                dbc.DropdownMenuItem("Export Options", header=True, className="file-menu-header"),
                dbc.DropdownMenuItem("Current View", id="curr-view-export-button"),
                dbc.DropdownMenuItem("File Pruning", id="tree-prune-export-button"),
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
                id=f"toggle-toolbar",
                className="nav-button",
            ),
            className="nav-button",
        ),
        # Graph Options
        dbc.NavItem(
            dbc.NavLink(
                f"Graph Customization",
                id=f"toggle-plot-options",
                className="nav-button",
            ),
            className="nav-button",
        ),
        # Docs
        dbc.NavItem(
            dbc.NavLink(
                f"Docs",
                id=f"toggle-docs",
                className="nav-button",
            ),
            className="nav-button",
        ),
    ],
    color="orange",
    expand=True,
)

tree_viewer_file_input_button = dbc.InputGroup(
    [
        dbc.InputGroupAddon(
            children=[
                dcc.Upload(
                    id="upload-data",
                    children=[
                        dbc.Button(
                            id="upload-tv-file-button",
                            children=["Tree Viewer File"],
                            className="tv-input-buttons",
                            color="info",
                        )
                    ],
                ),
            ],
            addon_type="prepend",
        ),
        html.Div(id="tv-filename-div", className="tv-filename"),
    ],
    size='lg',
)

chromosome_length_input_button = dbc.InputGroup(
    [
        dbc.InputGroupAddon(
            children=[
                dcc.Upload(
                    id="chromFile-upload",
                    children=[
                        dbc.Button(
                            id="chromFile-button",
                            children=["Chromosome Lengths"],
                            className="tv-input-buttons",
                            color="info",
                        ),
                    ],
                ),
            ],
            addon_type="prepend",
        ),
        html.Div(id="chrom-filename-div", className="tv-filename"),
    ],
    size='lg',
)

gff_input_button = dbc.InputGroup(
    [
        dbc.InputGroupAddon(
            children=[
                dcc.Upload(
                    id="gff-upload-data",
                    children=[
                        dbc.Button(
                            id="gff-upload-file-button",
                            children=["GFF/GTF file"], 
                            className="tv-input-buttons",
                            color="info"
                        )
                    ]
                ),
            ],
            addon_type="prepend",
        ),
        html.Div(id="gff-filename-div", className="tv-filename"),
    ],
    size='lg',
)

######################################  Layout  #####################################
def tv_layout():
    layout = dbc.Container(
        children=[
            # --- Data Storage ---
            dcc.Store(id="topology-color-chart-store"),
            dcc.Store(id="topology-freq-order-store"),
            dcc.Store(id="current-data-range"),
            # Tooltips
            dbc.Tooltip(
                "Order topologies based on whole-genome frequencies or by current chromosome frequencies",
                target="topo-freq-order",
                placement="top",
            ),
            # Open/close tooltips
            dbc.Tooltip(
                "Open/Close",
                target="toggle-toolbar",
                placement="top",
                delay={"show": 500, "hide": 250}
            ),
            dbc.Tooltip(
                "Open/Close",
                target="toggle-plot-options",
                placement="top",
                delay={"show": 500, "hide": 250}
            ),
            dbc.Tooltip(
                "Open/Close",
                target="toggle-docs",
                placement="top",
                delay={"show": 500, "hide": 250}
            ),
            dbc.Tooltip(
                "Download data in current view set by range slider",
                target="curr-view-export-button",
                # placement="top",
                delay={"show": 500, "hide": 250}
            ),
            dbc.Tooltip(
                "Preview rug/tile plot - zooming+panning disabled",
                target="main-graph-switches",
                delay={"show": 500, "hide": 250}
            ),
            dbc.Tooltip(
                "Return to homepage",
                target="home-button",
                delay={"show": 500, "hide": 250}
            ),
            # --- Modals ---
            # TreeViewer + BED file input
            dbc.Modal(
                children=[
                    dbc.ModalHeader(
                        children=[
                            dbc.Row([
                                dbc.Label("New Session Upload", className="modal-label"),
                            ], justify="center", align="center")
                        ],
                        className="modal-header"
                    ),
                    dbc.ModalBody(
                        children=[
                            dbc.Row(
                                children=[
                                    tree_viewer_file_input_button,
                                    chromosome_length_input_button,
                                ],
                                no_gutters=True
                            ),
                        ],
                    ),
                    dbc.ModalFooter(
                        children=[
                            dbc.Col([html.H5(id="tv-modal-warning-label")]),
                            dcc.Loading(
                                [html.Div(id="input-data-upload", className="hidden-div")],
                                className="toolbar-loading",
                                # color="#df691a"
                            ),
                            dbc.Button("Submit", id="input-files-submit-button",
                                       className="submit-buttons", color="info"),
                            dbc.Button("Close", id="input-files-close-button",
                                       className="submit-buttons", color="info"),
                        ],
                    ),
                ],
                id="main-input-modal",
                size="lg",
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
                            # dbc.Label("Alternative Input Files", className="modal-label"), 
                        ],
                        className="modal-header"
                    ),
                    dbc.ModalBody(
                        children=[
                            dbc.Row(
                                gff_input_button,
                                no_gutters=True,
                            ),
                        ],
                    ),
                    dbc.ModalFooter(
                        dbc.Button("Submit", id="gff-submit-button",
                                   className="submit-buttons", color="info")
                    ),
                ],
                id="alt-input-modal",
                size="lg",
            ),
            # Tree Pruning Modal
            dbc.Modal(
                children=[
                    dbc.ModalHeader(
                        children=[
                            dbc.Row([
                                dbc.Label("Tree Viewer File Pruning", className="modal-label"),
                            ], justify="center", align="center")
                        ],
                        className="modal-header"
                    ),
                    dbc.ModalBody(
                        children=[
                            dbc.Row(
                                children=[
                                    dbc.Col(
                                        children=[
                                            html.Div(
                                                children=[
                                                    html.P(
                                                        """
                                                        The file pruning export option allows you to select taxa and create a new Tree Viewer input files with a pruned set of Newick trees. This process 
                                                        will also clear the TopologyID column, enabling you to re-bin the new tree topologies. There is an option to run the binning process before downloading 
                                                        the new file, but please note that re-binning large trees across large genomes (e.g, >2 Gbp) may require considerable run time and it will lock 
                                                        your session until it completes. In such circumstances we recommend running the --topobinner command in THExBuilder.
                                                        """
                                                    ),
                                                ],
                                            ),
                                            html.Hr(style={"background-color": "white"}),
                                            html.Div(
                                                children=[
                                                    html.H4("Select taxa to keep:"),
                                                    dcc.Dropdown(
                                                        className="dropdown-style",
                                                        id="tree-prune-export-dropdown",
                                                        multi=True,
                                                    ),
                                                ],
                                                className="flex",
                                            ),
                                            html.Hr(style={"background-color": "white"}),
                                            html.Div([
                                                dcc.Checklist(
                                                    id="topobinner-checkbox",
                                                    options=[{"label": "Run topobinner", "value": "run"}],
                                                    labelStyle={'font-size': '15pt', 'margin-left': '5px'}
                                                ),
                                            ], style={"display": "inline-block", "width": "100%"}),
                                        ],
                                        width=12
                                    ),
                                ],
                            ),
                        ],
                    ),
                    dbc.ModalFooter(
                        children=[
                            dcc.Loading([html.Div(id="topobinner-loading", className="hidden-div")], className="toolbar-loading"),
                            dbc.Button("Submit", id="tree-prune-submit-button",
                                    className="submit-buttons", color="info"),
                            dbc.Button("Close", id="tree-prune-close-button",
                                    className="submit-buttons", color="info")
                        ],
                    ),
                ],
                id="tree-input-modal",
                size="xl",
                keyboard=False,
            ),
            # Empty DIVs for data chache
            html.Div(id="chromFile-data-upload", className="hidden-div"),
            html.Div(id="gff-data-upload", className="hidden-div"),
            html.Div(id="window-size", className="hidden-div"),
            html.Div(id="input-place-holder", className="hidden-div"),
            # Download Components
            dcc.Download(id="tv-current-view-download"),
            dcc.Download(id="gff-current-view-download"),
            dcc.Download(id="tree-prune-download"),
            # Navbar
            dbc.Row(
                children=[navbar],
                className="navbar-div",
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
                                            style={"border-radius": "5px"}
                                        ),
                                    ],
                                ),
                                id="collapse-documentation",
                            ),
                        ],
                        id="collapse-docs-card",
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
                                                # --- Graph Color Theme ---
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            children=[
                                                                "Color Theme:"],
                                                            className="text-divs",
                                                        ),
                                                        dcc.Dropdown(
                                                            id="template-option",
                                                            options=[{"label": i, "value": i}
                                                                    for i in GRAPH_TEMPLATES],
                                                            value="plotly_dark",
                                                            clearable=False,
                                                            className="dropdown-style",
                                                        ),
                                                    ],
                                                    className="col-divs",
                                                    sm=2,
                                                    lg=2,
                                                ),
                                                # --- Tree Style ---
                                                # dbc.Col(
                                                #     children=[
                                                #         html.Div(
                                                #             children=[
                                                #                 "Tree Style:"],
                                                #             className="text-divs",
                                                #         ),
                                                #         dcc.Dropdown(
                                                #             id="tree-shape-option",
                                                #             options=[
                                                #                 {"label": "Rectangular",
                                                #                     "value": "box"},
                                                #                 {"label": "Angular",
                                                #                     "value": "angular"},
                                                #             ],
                                                #             value="box",
                                                #             clearable=False,
                                                #             className="dropdown-style",
                                                #         ),
                                                #     ],
                                                #     className="col-divs",
                                                #     width=2,
                                                # ),
                                                # --- Data Color Palette ---
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            children=["Color Palette:"],
                                                            className="text-divs",
                                                        ),
                                                        dcc.Dropdown(
                                                            id="line_color",
                                                            options=[
                                                                {"label": f"{i} - {len(COLOR_SWATCHES[i])} colors", "value": i} for i in COLOR_SWATCHES],
                                                            value="Plotly",
                                                            clearable=False,
                                                            className="dropdown-style",
                                                        ),
                                                    ],
                                                    className="col-divs",
                                                    sm=2,
                                                    lg=2,
                                                ),
                                                # --- Output File Format ---
                                                dbc.Col(
                                                    children=[
                                                        html.H6(
                                                            children=[
                                                                "Output File Format:"],
                                                            className="text-divs",
                                                        ),
                                                        dcc.Dropdown(
                                                            id="snapshot-file-option",
                                                            options=[{"label": i, "value": i} for i in SNAPSHOT_FILE_OPTIONS],
                                                            value="svg",
                                                            clearable=False,
                                                            className="dropdown-style",
                                                        ),
                                                    ],
                                                    className="col-divs",
                                                    sm=2,
                                                    lg=2,
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
                                                            id="axis-line-width",
                                                            options=[{"label": i, "value": i} for i in range(1, 6)],
                                                            value=1,
                                                            clearable=False,
                                                            className="dropdown-style",
                                                        ),
                                                    ],
                                                    className="col-divs",
                                                    sm=2,
                                                    lg=2,
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
                                                    className="col-divs",
                                                    sm=2,
                                                    lg=2,
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
                                                    className="col-divs",
                                                    width='auto',
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
                                                    className="col-divs",
                                                    width='auto',
                                                ),
                                            ],
                                            className="row-div",
                                            no_gutters=True,
                                        ),
                                    ],
                                    className="narbar-collapse-plot-options",
                                ),
                                id="collapse-plot-options",
                            ),
                        ],
                        id="collapse-toolbar-card",
                        className="narbar-card"
                    ),
                    # Toolbar
                    dbc.Card(
                        children=[
                            dbc.Collapse(
                                dbc.CardBody(
                                    children=[
                                        dbc.Tabs(
                                            children=[
                                                # Single Chromosome View
                                                dbc.Tab(
                                                    children=[
                                                        dbc.Row(
                                                            children=[
                                                                # Chromosomes
                                                                dbc.Col(
                                                                    children=[
                                                                        dbc.Col(
                                                                            children=[
                                                                                # Chromosome
                                                                                html.Div(
                                                                                    children=[
                                                                                        html.H6(
                                                                                            children=[
                                                                                                html.Div(["Chromosome | "], style={"margin-bottom": "5px"}),
                                                                                            ],
                                                                                            className="text-divs",
                                                                                        ),
                                                                                        dcc.Dropdown(
                                                                                            id="chromosome_options",
                                                                                            className="dropdown-style",
                                                                                            optionHeight=20,
                                                                                        ),
                                                                                    ],
                                                                                    className="home-tab-button-cols",
                                                                                ),
                                                                                # Topologies
                                                                                html.Div(
                                                                                    children=[
                                                                                        html.H6(
                                                                                            children=[
                                                                                                html.Div(["Topologies | "], style={"margin-right": "2px"}),
                                                                                                html.Div(["(Genome"]),
                                                                                                dbc.Checklist(
                                                                                                    id="topo-freq-order",
                                                                                                    options=[
                                                                                                        {"label": "",
                                                                                                            "value": "topo_freq_organization"},
                                                                                                    ],
                                                                                                    switch=True,
                                                                                                    className="topo-organization",
                                                                                                ),
                                                                                                html.Div(["Chromosome)"]),
                                                                                            ],
                                                                                            className="text-divs",
                                                                                        ),
                                                                                        dcc.Dropdown(
                                                                                            id="topology_options",
                                                                                            className="dropdown-style",
                                                                                            multi=True,
                                                                                            clearable=False,
                                                                                        ),
                                                                                    ],
                                                                                    className="home-tab-button-cols",
                                                                                ),
                                                                            ], width=6
                                                                        ),
                                                                        dbc.Col(
                                                                            children=[
                                                                                # Additional Data
                                                                                html.Div(
                                                                                    children=[
                                                                                        html.H6(
                                                                                            children=[
                                                                                                html.Div(["Additional Data | "], style={"margin-bottom": "5px"}),
                                                                                                # html.Div(["(Off"]),
                                                                                                dbc.Checklist(
                                                                                                    id="alt-data-switch",
                                                                                                    options=[{"label": " ", "value": "alt_data_switch_on", "disabled": True}],
                                                                                                    switch=True,
                                                                                                    className="topo-organization",
                                                                                                ),
                                                                                                # html.Div(["On)"]),
                                                                                            ],
                                                                                            className="text-divs",
                                                                                        ),
                                                                                        dcc.Dropdown(
                                                                                            id="alt_data_options",
                                                                                            className="dropdown-style",
                                                                                            multi=True,
                                                                                            clearable=False,
                                                                                        ),
                                                                                    ],
                                                                                    className="home-tab-button-cols",
                                                                                ),
                                                                                # Trees
                                                                                html.Div(
                                                                                    children=[
                                                                                        html.H6([
                                                                                            html.Div(["Trees | "], style={"margin-bottom": "5px"}),
                                                                                            # html.Div(["(Off"]),
                                                                                            dbc.Checklist(
                                                                                                id="tree-graph-switch",
                                                                                                options=[{"label": " ", "value": "tree_switch", "disabled": True}],
                                                                                                switch=True,
                                                                                                className="topo-organization",
                                                                                            ),
                                                                                            # html.Div(["On)"]),
                                                                                        ], className="text-divs"),
                                                                                        dcc.Dropdown(
                                                                                            id="tree-taxa-choices",
                                                                                            className="dropdown-style",
                                                                                            multi=True,
                                                                                        ),
                                                                                    ],
                                                                                    className="home-tab-button-cols",
                                                                                ),
                                                                            ], width=6
                                                                        ),
                                                                    ],
                                                                    className="home-tab-divs",
                                                                    # width="auto",
                                                                    width=8,
                                                                ),
                                                                # Graph Toggle Switches
                                                                dbc.Col(
                                                                    children=[
                                                                        html.Div(
                                                                            children=[
                                                                                html.H6("Graph Switches:", className="text-divs",),
                                                                                dbc.FormGroup(
                                                                                    children=[
                                                                                        dbc.Label("Main Distribution", className="text-divs-sm"),
                                                                                        dbc.RadioItems(
                                                                                            id="main-graph-switches",
                                                                                            options=[
                                                                                                {"label": "Rug+Tile", "value": "topo_both"},
                                                                                                {"label": "Rug Preview", "value": "topo_rug_only"},
                                                                                                {"label": "Tile Preview", "value": "topo_tile_only"},
                                                                                            ],
                                                                                            value="topo_both",
                                                                                            className="checklist-style",
                                                                                            inline=True,
                                                                                        ),
                                                                                        html.Br(),
                                                                                        dbc.Label("Standarized File Type Graphs", className="text-divs-sm"),
                                                                                        dbc.Checklist(
                                                                                            id="alt-graph-switches",
                                                                                            options=[
                                                                                                {"label": "GFF/GTF - No file chosen...", "value": "gff_switch", "disabled": True},
                                                                                            ],
                                                                                            className="checklist-style",
                                                                                            inline=True,
                                                                                        ),
                                                                                    ],
                                                                                    className="formgroup",
                                                                                    inline=True,
                                                                                ),
                                                                            ],
                                                                            className="home-tab-button-cols",
                                                                        ),
                                                                    ],
                                                                    className="home-tab-divs-no-border",
                                                                ),
                                                            ],
                                                            className="nav-row-div",
                                                            no_gutters=True,
                                                        ),
                                                    ],
                                                    id="home-tab",
                                                    tab_id="home-tab",
                                                    label="Single Chromosome",
                                                    label_style={
                                                        "color": "white",
                                                        "border": "2px orange solid",
                                                        "border-radius": "5px",
                                                        "margin-right": "2px",
                                                    },
                                                ),
                                                # Whole Genome View
                                                dbc.Tab(
                                                    id="whole-genome-tab",
                                                    children=[
                                                        dbc.Row(
                                                            className="row-div",
                                                            no_gutters=True,
                                                            children=[
                                                                html.Div(
                                                                    children=[
                                                                        html.Div([
                                                                            html.H6("Chromosome Frequency Plot Type:", className="text-divs",),
                                                                            dcc.Dropdown(
                                                                                id='stats-chrom-plot-type',
                                                                                className="dropdown-style",
                                                                                clearable=False,
                                                                                searchable=False,
                                                                                options=[
                                                                                    {'label': 'Rug', 'value': 'rug'},
                                                                                    {'label': 'Bar', 'value': 'bar'},
                                                                                    {'label': 'Pie', 'value': 'pie'},
                                                                                    {'label': 'Tile', 'value': 'tile'},
                                                                                ],
                                                                                value='rug'
                                                                            ),
                                                                        ]),
                                                                    ],
                                                                    className="stats-toolbar",
                                                                )
                                                            ],
                                                        ),
                                                    ],
                                                    tab_id="whole-genome-view",
                                                    label="Whole Genome",
                                                    label_style={
                                                        "color": "white",
                                                        "border": "2px orange solid",
                                                        "border-radius": "5px",
                                                    },
                                                ),
                                            ],
                                            active_tab="home-tab",
                                            className="tabs-style",
                                            id="control-tabs",
                                        ),
                                    ],
                                    className="narbar-collapse-toolbar",
                                ),
                                id="collapse-toolbar",
                            ),
                        ],
                        id="collapse-toolbar-card",
                        className="narbar-card"
                    ),
                ],
                className="narbar-collapse-toolbar",
                no_gutters=True,
            ),
            # Graphs
            dbc.Row(
                children=[
                    # -------- Main Distriubtion --------
                    dbc.Row(
                        id='main-graph-div',
                        children=[
                            # --- Main Topology Graph ---
                            dbc.Col(
                                children=[
                                    html.Div(
                                        children=[
                                            dcc.Loading(
                                                children=[
                                                    html.Div(
                                                        id="main-dist-plot",
                                                        children=[
                                                            dcc.Graph(
                                                                id="topology_graph",
                                                                figure=tree_utils.init_data_graph("plotly_dark"),
                                                                className="topograph-style",
                                                                config=dict(displaylogo=False, doubleClick="autosize", displayModeBar=False),
                                                            ),
                                                            html.Div(id="preview-graphs-div"),
                                                            html.Div(id="additional-graphs-div"),
                                                        ],
                                                    ),
                                                ],
                                                id="topology_graph_div",
                                                className="hidden-div",
                                            ),
                                        ],
                                    ),
                                ],
                                width=12,
                            ),
                        ],
                        className="visible-div",
                        no_gutters=True,
                    ),
                    # -------- Whole genome graphs --------
                    dbc.Row(
                        id="stats-toolbar-div",
                        children=[
                            # Basic info
                            dbc.Row(
                                children=[
                                    # Topology frequencies
                                    dbc.Col(
                                        id='stats-pie-chart-col',
                                        children=[
                                            dcc.Loading(
                                                id="loading-stats-pie-chart",
                                                children=[
                                                    dcc.Graph(
                                                        id="stats-topology-freq-pie-graph",
                                                        figure=tree_utils.init_data_graph("plotly_dark"),
                                                        className="stats-graph-style",
                                                    )
                                                ],
                                            ),
                                        ],
                                        className="hidden-div",
                                        width=3,
                                    ),
                                    # RF Distance Plot of Top 3 topologies
                                    dbc.Col(
                                        id='stats-r1-d2',
                                        children=[
                                            dcc.Loading(
                                                id="loading-stats-RF-chart",
                                                children=[
                                                    dcc.Graph(
                                                        id="whole-genome-rf-graph",
                                                        figure=tree_utils.init_data_graph("plotly_dark"),
                                                        className="stats-graph-style",
                                                    ),
                                                ],
                                            ),
                                        ],
                                        className="hidden-div",
                                        width=6,
                                    ),
                                    # Basic Stats
                                    dbc.Col(
                                        id='stats-r1-d3',
                                        children=[
                                            dcc.Loading(
                                                id="loading-stats-basic-info-chart",
                                                children=[
                                                    html.Div(
                                                        id='stats-table-div',
                                                        className="stats-datatable"),
                                                ],
                                            ),
                                        ],
                                        className="hidden-div",
                                        width=3,
                                    ),
                                ],
                                no_gutters=True,
                            ),
                            # Per-chromosome Frequency Plots
                            dcc.Loading([
                                dbc.Row(
                                    id="whole-genome-graphs",
                                    className="stats-chrom-frequency",
                                    no_gutters=True,
                                ),
                            ]),
                        ],
                        className="hidden-div",
                        no_gutters=True,
                    ),   
                ],
                className="div-boxes",
                no_gutters=True,
            ),   
        ],
        id="tree-viewer-container",
        className="container-div",
        fluid=True,
    )
    return layout

##################################### Callbacks #####################################

@app.callback(
    Output("tv-docs-content", "children"),
    [Input("collapse-documentation", "is_open"),],)
def render_docs_content(is_open):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if is_open:
        page_content = html.P(
            children=[docs_treeviewer_contents],
        )
        return page_content
    else:
        raise PreventUpdate

# ---------------------- Toolbar Toggling ---------------------
@app.callback(
    [Output("collapse-toolbar", "is_open"),
     Output("collapse-documentation", "is_open"),
     Output("collapse-plot-options", "is_open"),],
    [Input("toggle-toolbar", "n_clicks"),
     Input("toggle-docs", "n_clicks"),
     Input("toggle-plot-options", "n_clicks"),],
    [State("collapse-toolbar", "is_open"),
     State("collapse-documentation", "is_open"),
     State("collapse-plot-options", "is_open"),],
)
def toggle_toolbar(n1, n2, n3, toolbar_is_open, docs_is_open, plot_options_is_open):
    ctx = dash.callback_context

    if not ctx.triggered:
        return True, False, False
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if button_id == "toggle-toolbar" and n1:
        return not toolbar_is_open, False, False
    elif button_id == "toggle-docs" and n2:
        return False, not docs_is_open, False
    elif button_id == "toggle-plot-options" and n3:
        return False, False, not plot_options_is_open
    else:
        return False, False, False

# ---------------------- Hidden div handlers ---------------------
@app.callback([Output("topology_graph_div", "className"),
               Output("main-graph-div", "className"),
               Output("stats-toolbar-div", "className"),
               Output("stats-pie-chart-col", "className"),
               Output("stats-r1-d2", "className"),
               Output("stats-r1-d3", "className"),],
              [Input("control-tabs", "active_tab"),
               Input("input-data-upload", "children"),
               Input("chromFile-data-upload", "children"),]
              )
def hide_graph_divs(control_tab, tvdata, chromdata):
    """This callback controls the graphs display and styling"""
    if not tvdata:
        raise PreventUpdate
    elif not chromdata:
        raise PreventUpdate
    elif control_tab == 'home-tab':
        # Topology distribution
        topology_graph_style = "dist-graph-div"
        distribution_div = "visible-div"
        # Hidden divs
        stats_toolbar_style = "hidden-div"
        stats_pie_chart = "hidden-div"
        stats_r1_d2 = "hidden-div"
        stats_r1_d3 = "hidden-div"
        return (topology_graph_style, distribution_div, stats_toolbar_style, 
                stats_pie_chart, stats_r1_d2, stats_r1_d3)
    ## ---- Whole Genome Graphs ----
    elif control_tab == "whole-genome-view":
        # Hidden divs
        topology_graph_style = "hidden-div"
        distribution_div = "hidden-div"
        # On divs
        stats_toolbar_style = "stats-main-div"
        stats_pie_chart = "stats-topo-frequency"
        stats_r1_d2 = "stats-topo-frequency"
        stats_r1_d3 = "stats-topo-frequency-no-buttom-border"
        return (topology_graph_style, distribution_div, stats_toolbar_style, 
                stats_pie_chart, stats_r1_d2, stats_r1_d3)
    else:
        raise PreventUpdate

# ---------------------- Update Graph Toggle Options + Dropdowns ---------------------
@app.callback(
    [Output("main-graph-switches", "options"), 
     Output("alt-graph-switches", "options"),
     Output("tree-graph-switch", "options"),
     Output("alt-data-switch", "options"),],
    [Input("gff-data-upload", "children"),
     Input("alt_data_options", "value"),
     Input("input-data-upload", "children"),
     Input("alt_data_options", "options"),]
)
def main_graph_button_options(gffData, altData, topoData, additional_graphs):
    if not gffData:
        gffOption = {"label": "GFF/GTF - No file chosen...", "value": "gff_switch", "disabled": True}
    else:
        gffOption = {"label": "GFF/GTF", "value": "gff_switch"}
    if not topoData:
        altDataOption = {"label": " ", "value": "alt_data_switch_on", "disabled": True}
        treeOption = {"label": " ", "value": "tree_switch", "disabled": True}
        combinedPlotOption = {"label": "Rug+Tile", "value": "topo_both", "disabled": True}
        rugPlotOption = {"label": "Rug Preview", "value": "topo_rug_only", "disabled": True}
        tilePlotOption = {"label": "Tile Preview", "value": "topo_tile_only", "disabled": True}
        # normRFPlotOption = {"label": "NormRF-distance", "value": "rf_dist", "disabled": True}
        gffOption = {"label": "GFF/GTF - No file chosen...", "value": "gff_switch", "disabled": True}
    else:
        altDataOption = {"label": " ", "value": "alt_data_switch_on"}
        treeOption = {"label": " ", "value": "tree_switch"}
        combinedPlotOption = {"label": "Rug+Tile", "value": "topo_both"}
        rugPlotOption = {"label": "Rug Preview", "value": "topo_rug_only"}
        tilePlotOption = {"label": "Tile Preview", "value": "topo_tile_only"}
        # normRFPlotOption = {"label": "NormRF-distance", "value": "rf_dist"}
    # If no additional data loaded, disable additional data switch
    if not additional_graphs:
        altDataOption = {"label": " ", "value": "alt_data_switch_on", "disabled": True}
    elif additional_graphs[0]['value'] == "":
        altDataOption = {"label": " ", "value": "alt_data_switch_on", "disabled": True}

    main_graph_options = [
        combinedPlotOption,
        rugPlotOption,
        tilePlotOption,
    ]
    additional_graph_options = [
        gffOption,
    ]
    tree_graph_options = [
        treeOption,
    ]
    alt_data_options = [
        altDataOption,
    ]
    return (main_graph_options, additional_graph_options,
            tree_graph_options, alt_data_options)


@app.callback(
    [Output("alt_data_options", 'disabled'),
     Output("tree-taxa-choices", 'disabled'),
     Output("chromosome_options", "disabled"),
     Output("topology_options", "disabled"),],
    [Input("input-data-upload", "children"),
     Input("alt-data-switch", "value"),
     Input("tree-graph-switch", "value"),],
)
def disable_dropdowns(topo_data, alt_switch, tree_switch):
    if not topo_data:
        chromDropdown = True
        topoDropdown = True
        altDropdown = True
        treeDropdown = True
        return altDropdown, treeDropdown, chromDropdown, topoDropdown
    else:
        chromDropdown = False
        topoDropdown = False
    if not alt_switch:
        altDropdown = True
    else:
        altDropdown = False
    if not tree_switch:
        treeDropdown = True
    else:
        treeDropdown = False
    return altDropdown, treeDropdown, chromDropdown, topoDropdown

# ---------------------- Data Upload Callbacks ---------------------
@app.callback(
    [Output("alt-input-modal", "is_open"),
     Output("gff-data-upload", "children"),
     Output("gff-filename-div", "children"),
     Output("gff-upload-file-button", "color"),],
    [Input("gff-init-upload-button", "n_clicks"),
     Input("gff-submit-button", "n_clicks"),
     Input("input-files-submit-button", "n_clicks"),
     Input("gff-upload-data", "contents"), ],
    [State("gff-upload-data", "filename")])
def upload_gff_file(gffBtn, gffSubmit, tvBtn, gffContents, gffFilename):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "No clicks yet"
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if button_id == "gff-init-upload-button":
        if not gffFilename:
            gffFilename = "No file chosen..."
        return True, None, gffFilename, "primary"
    elif button_id == "gff-upload-file-button":
        return True, None, None, "primary"
    elif button_id == "input-files-submit-button":
        return False, None, None, "primary"
    elif button_id == "gff-submit-button":
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
                # skiprows=data_utils.get_gff_header_len(
                #     tvDecoded.decode("utf-8")),
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
            return True, None, gffFilename, "danger"
    else:
        raise PreventUpdate


# Upload Tree Viewer and Chromosome length files
@app.callback(
    [Output("input-data-upload", "children"),
     Output("chromFile-data-upload", "children"),
     Output("main-input-modal", "is_open"),
     Output("tv-filename-div", "children"),
     Output("chrom-filename-div", "children"),
     Output("window-size", "children"),
     Output("upload-tv-file-button", "color"),
     Output("chromFile-button", "color"),
     Output("tv-modal-warning-label", "children"),
     Output("tree-taxa-choices", "options"),
     Output("tree-taxa-choices", "value")],
    [Input("upload-file-button", "n_clicks"),
     Input("input-files-submit-button", "n_clicks"),
     Input("input-files-close-button", "n_clicks"),
     Input("upload-data", "contents"),
     Input("chromFile-upload", "contents"), ],
    [State("input-data-upload", "children"),
     State("chromFile-data-upload", "children"),
     State("chromFile-upload", "filename"),
     State("upload-data", "filename")],)
def load_input_files(input_btn, submit_btn, close_btn, tvFile_contents,
                    chromFile_contents, tv_data, chrom_data, chromFilename,
                    tvFilename
):
    ctx = dash.callback_context
    button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    # Set button colors
    if not tvFilename:
        tv_color = "primary"
    elif 'txt' in tvFilename:
        tv_color = "success"
    elif 'csv' in tvFilename:
        tv_color = "success"
    elif 'tsv' in tvFilename:
        tv_color = "success"
    elif 'xls' in tvFilename:
        tv_color = "success"
    else:
        tv_color = "danger"
    if not chromFilename:
        chrom_color = "primary"
    elif 'bed' in chromFilename:
        chrom_color = "success"
    else:
        chrom_color = "danger"
    # Check buttons + return
    if button_id == "upload-file-button":
        if not tvFilename:
            tvFilename = "No file chosen..."
        if not chromFilename:
            chromFilename = "No file chosen..."
        return (dash.no_update, dash.no_update, True, tvFilename, chromFilename, 
                dash.no_update, tv_color, chrom_color, None, dash.no_update, dash.no_update)
    elif button_id == "input-files-close-button":
        return (dash.no_update, dash.no_update, False, tvFilename, chromFilename, 
                dash.no_update, tv_color, chrom_color, None, dash.no_update, dash.no_update)
    elif button_id == "input-files-submit-button":
        if (not tvFile_contents) and (not chromFile_contents):
            return (dash.no_update, dash.no_update, True, "No file chosen!", 
                    "No file chosen!", dash.no_update, "danger", "danger", None, 
                    dash.no_update, dash.no_update)
        elif not tvFile_contents:
            return (dash.no_update, dash.no_update, True, "No file chosen!", 
                    chromFilename, dash.no_update, "danger", chrom_color, None, dash.no_update, 
                    dash.no_update)
        elif not chromFile_contents:
            return (dash.no_update, dash.no_update, True, tvFilename, "No file chosen!", 
                    dash.no_update, tv_color, "danger", None, dash.no_update, dash.no_update)
        else:
            _, tvcontent_string = tvFile_contents.split(",")
            tvDecoded = base64.b64decode(tvcontent_string)
            _, ChromContent_string = chromFile_contents.split(",")
            chromDecoded = base64.b64decode(ChromContent_string)
            tvDF = None
            chromDF = None
            if "csv" in tvFilename:
                # Assume that the user uploaded a CSV file
                tvDF = pd.read_csv(io.StringIO(tvDecoded.decode("utf-8")))
                validated_tvDF = tree_utils.tv_header_validation(tvDF)
                # If validated_DF returns False, raise error 
                # indicating issue with input file 
                if not validated_tvDF:
                    return (None, None, True, f"{tvFilename} is malformed!", chromFilename, 
                            None, "danger", "warning", 
                            "WARNING: Tree Viewer file appears to be malformed. First four headers must be ['Chromosome', 'Window', 'NewickTree', 'TopologyID']", 
                            [], [])
                if len(tvDF.columns) == 4:
                    tvDF["None"] = [pd.NA]*len(tvDF)
                tvDF["TopologyID"] = tvDF["TopologyID"].apply(
                    lambda x: "NODATA" if type(x) != str else x)
                pass
            elif "txt" in tvFilename:
                # Assume that the user uploaded a CSV file
                tvDF = pd.read_csv(io.StringIO(tvDecoded.decode("utf-8")), sep='\t')
                validated_tvDF = tree_utils.tv_header_validation(tvDF)
                # If validated_DF returns False, raise error 
                # indicating issue with input file 
                if not validated_tvDF:
                    return (None, None, True, f"{tvFilename} is malformed!", chromFilename, 
                            None, "danger", "warning", 
                            "WARNING: Tree Viewer file appears to be malformed. First four headers must be ['Chromosome', 'Window', 'NewickTree', 'TopologyID']", 
                            [], [])
                if len(tvDF.columns) == 4:
                    tvDF["None"] = [pd.NA]*len(tvDF)
                tvDF["TopologyID"] = tvDF["TopologyID"].apply(
                    lambda x: "NODATA" if type(x) != str else x)
                pass
            elif "tsv" in tvFilename:
                # Assume that the user uploaded a CSV file
                tvDF = pd.read_csv(io.StringIO(tvDecoded.decode("utf-8")), sep='\t')
                validated_tvDF = tree_utils.tv_header_validation(tvDF)
                # If validated_DF returns False, raise error 
                # indicating issue with input file 
                if not validated_tvDF:
                    return (None, None, True, f"{tvFilename} is malformed!", chromFilename, 
                            None, "danger", "warning", 
                            "WARNING: Tree Viewer file appears to be malformed. First four headers must be ['Chromosome', 'Window', 'NewickTree', 'TopologyID']", 
                            [], [])
                if len(tvDF.columns) == 4:
                    tvDF["None"] = [pd.NA]*len(tvDF)
                tvDF["TopologyID"] = tvDF["TopologyID"].apply(
                    lambda x: "NODATA" if type(x) != str else x)
                pass
            elif "xls" in tvFilename:
                # Assume that the user uploaded an excel file
                tvDF = pd.read_excel(io.BytesIO(tvDecoded), engine="openpyxl", comment="#")
                # validated_tvDF = tree_utils.validate_tree_viewer_input(tvDF)
                validated_tvDF = tree_utils.tv_header_validation(tvDF)
                # If validated_DF returns False, raise error 
                # indicating issue with input file 
                if not validated_tvDF:
                    return (None, None, True, f"{tvFilename} is malformed!", chromFilename, 
                            None, "danger", "warning", 
                            "WARNING: Tree Viewer file appears to be malformed. First four headers must be ['Chromosome', 'Window', 'NewickTree', 'TopologyID']", 
                            [], [])
                if len(tvDF.columns) == 4:
                    tvDF["None"] = [pd.NA]*len(tvDF)
                tvDF["TopologyID"] = tvDF["TopologyID"].apply(lambda x: "NODATA" if type(x) != str else x)
                tvDF[["TopologyID"]] = tvDF[["TopologyID"]].fillna(value="NoData")
                pass
            else:
                return None, None, True, f"Invalid file type - {tvFilename}", chromFilename, None, "warning", "primary", "WARNING: Please provide a valid Tree Viewer file", [], []
            if "bed" in chromFilename:
                # Assume that the user uploaded a CSV file
                chromDF = pd.read_csv(
                    io.StringIO(chromDecoded.decode("utf-8")),
                    sep="\t",
                    names=["Chromosome", "Start", "End"],
                )
                pass
            else:
                return (None, None, True, tvFilename, f"Invalid file type - {chromFilename}", None, "primary", 
                        "warning", "WARNING: Please provide a valid BED file", [], [])
            # If Chromosome data type is int -> add chr to make categorical
            if tvDF["Chromosome"].dtype != 'object':
                tvDF["Chromosome"] = tvDF["Chromosome"].apply(lambda x: f"chr{x}")
            if chromDF["Chromosome"].dtype != 'object':
                chromDF["Chromosome"] = chromDF["Chromosome"].apply(lambda x: f"chr{x}")
            
            # Validate chromosome length file
            msg, validate = tree_utils.validate_chrom_lengths(chromDF, tvDF)
            if not validate:
                # Raise modal telling user there is an issue and they need to fix it
                return (None, None, True, tvFilename, chromFilename, None, "warning", "warning", msg, [], [])
            # Calculate needed variables and return to hidden div
            window_size = data_utils.get_window_size(tvDF["Window"])
            # Round up chromosome lengths to nearest whole window
            tvDF[["TopologyID"]] = tvDF[["TopologyID"]].fillna(value="NoData")
            # chromDF["End"] = chromDF["End"].apply(lambda x: data_utils.roundUp(x, window_size))
            # Collect taxa names from first tree and build option list-dict
            tree_taxa_options = [{"label": i, "value": i} for i in tree_utils.get_taxa_from_tree(tvDF["NewickTree"][0])]
            tree_taxa_values = [i['value'] for i in tree_taxa_options]
            return (tvDF.to_json(orient="columns"), chromDF.to_json(), False, tvFilename, chromFilename, 
                    window_size, "success", "success", None, tree_taxa_options, tree_taxa_values)
    else:
        if not tvFilename:
            tvFilename = "No file chosen..."
        if not chromFilename:
            chromFilename = "No file chosen..."
        return (dash.no_update, dash.no_update, True, tvFilename, chromFilename, dash.no_update, tv_color, chrom_color, 
                None, dash.no_update, dash.no_update)

# ---------------------- Dropdown Option Callbacks ---------------------
# Set topology color mapping
@app.callback([Output("topology-color-chart-store", "data"),
               Output("whole-genome-tab", "disabled"),],
              [Input("input-data-upload", "children"),
               Input("line_color", "value"), ])
def set_topology_colors(data, line_color):
    if not data:
        return None, True
    else:
        color_set = COLOR_SWATCHES[line_color]
        topo_colors = tree_utils.set_topology_colors(data, color_set)
        return topo_colors, False


# Set chromosome options + value
@app.callback([Output("chromosome_options", "options"),
               Output("chromosome_options", "value")],
              [Input("input-data-upload", "children"), ])
def set_chromosome_options(data_json):
    if not data_json:
        raise PreventUpdate
    else:
        df = pd.read_json(data_json)
        unique_chroms = df["Chromosome"].unique()
        chromosome_options = sorted([{"label": i, "value": i} for i in sorted(
            unique_chroms)], key=lambda i: (i["label"], i["value"]))
        chrom_val = chromosome_options[0]["value"]
        return chromosome_options, chrom_val


# Set topology dropdown options + value
@app.callback([Output("topology_options", "options"),
               Output("topology_options", "value"),
               Output("topology-freq-order-store", "data"), ],
              [Input("input-data-upload", "children"),
               Input("chromosome_options", "value"),
               Input("topo-freq-order", "value")],
               [State("topology_options", "value"),],)
def set_topology_options(data_json, chromValue, freq_order, topo_vals):
    if not data_json:
        return [{"label": "None", "value": "None", "disabled": True}], None, None
    else:
        df = pd.read_json(data_json)
        if not topo_vals:
            if not freq_order:
                pass
            else:
                df = df[df["Chromosome"] == chromValue]
            value_counts = df["TopologyID"].value_counts()
            topology_options = [{"label": i, "value": i} for i in value_counts.index]
            topoOrder = list(value_counts.index)
            topoOrder.reverse()
            return topology_options, [c["value"] for c in topology_options[:3]], topoOrder
        else:
            if not freq_order:
                return dash.no_update, topo_vals, dash.no_update
            else:
                df = df[df["Chromosome"] == chromValue]
            value_counts = df["TopologyID"].value_counts()
            topology_options = [{"label": i, "value": i} for i in value_counts.index]
            topoOrder = list(value_counts.index)
            topoOrder.reverse()
            return topology_options, [c["value"] for c in topology_options[:3]], topoOrder
        


# Set alt data dropdown options + value
@app.callback([Output("alt_data_options", "options"),
               Output("alt_data_options", "value"), ],
              [Input("input-data-upload", "children")])
def set_alt_data_options(data_json):
    if not data_json:
        return [{"label": "None", "value": "None", "disabled": True}], None
    else:
        df = pd.read_json(data_json)
        alt_data_cols = [col for col in df.columns][4:]
        if not alt_data_cols:
            return None, None
        elif (len(alt_data_cols) == 1) and (alt_data_cols[0] == "None"):
            return [{"label": "", "value": ""}], None
        else:
            alt_data_options = [{"label": i, "value": i} for i in alt_data_cols]
            return alt_data_options, alt_data_options[0]["value"]

# ---------------------- Graph Building Callbacks ----------------------
@app.callback(
    # Outputs
    [Output("topology_graph", "figure"),
     Output("additional-graphs-div", "children"),
     Output("current-data-range", "data"),
     Output("topology_graph", "config"),
    ],
    # Inputs
    [Input("input-data-upload", "children"),
     Input("chromFile-data-upload", "children"),
     Input("gff-data-upload", "children"),
     Input("main-graph-switches", "value"),
     Input("alt-graph-switches", "value"),
     Input("alt_data_options", "value"),
     Input("chromosome_options", "value"),
     Input("topology_options", "value"),
     Input("template-option", "value"),
     Input("snapshot-file-option", "value"),
     Input("topology-color-chart-store", "data"),
     Input("topology-freq-order-store", "data"),
     Input("window-size", "children"),
     Input('topology_graph', 'relayoutData'),
     Input("main-input-modal", "is_open"),
     Input("alt-data-switch", "value"),
     Input("tree-graph-switch", "value"),
     Input("tree-taxa-choices", "value"),
     Input("axis-line-width", "value"),
     Input("tv-pixel-height", "value"),
     Input("tv-pixel-width", "value"),
     Input("axis-gridlines", "value"),
     Input("editable-graph-option", "value"),
     Input("control-tabs", "active_tab"),
    #  Input("tree-shape-option", "value"),
    ],
)
def home_tab_graphs(
    topo_data,
    chromosome_length_data,
    gff_data,
    main_graph_switches,
    alt_graph_switches,
    alt_dropdown_options,
    chromosome,
    current_topologies,
    template,
    snapshot_file_type,
    color_mapping,
    topoOrder,
    window_size,
    relayout_data,
    input_modal_open,
    alt_data_switch,
    tree_switch,
    tree_taxa_values,
    axis_line_width,
    pixel_height,
    pixel_width,
    axis_gridlines,
    editable_graphs,
    active_tab,
    # tree_shape,
):
    if input_modal_open:
        raise PreventUpdate
    elif active_tab == "whole-genome-view":
        return dash.no_update, None, None, None
    elif not topo_data:
        raise PreventUpdate
    elif not chromosome_length_data:
        raise PreventUpdate
    elif not current_topologies:
        # If not topologies are selected, return No Data Loaded graph
        return [], None, None
    else:
        topoOrder = [e for e in topoOrder if e in current_topologies]
        whole_topology_df = pd.read_json(topo_data)
        chromosome_df = pd.read_json(chromosome_length_data)
        chromosome_df["End"] = chromosome_df["End"].apply(lambda x: data_utils.roundUp(x, window_size))
        topology_df = whole_topology_df[whole_topology_df["Chromosome"] == chromosome]
        chromosome_df = chromosome_df[chromosome_df["Chromosome"] == chromosome]
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
            dataRange = [0, max(chromosome_df["End"])]
            dataMin, dataMax = 0, max(chromosome_df["End"])
            pass
        except TypeError:
            dataRange = [0, max(chromosome_df["End"])]
            dataMin, dataMax = 0, max(chromosome_df["End"])
            pass
        except IndexError:
            dataRange = [0, max(chromosome_df["End"])]
            dataMin, dataMax = 0, max(chromosome_df["End"])
            pass
        # Put alt data option into list if only one selected.
        if type(alt_dropdown_options) == str:
            alt_dropdown_options = [alt_dropdown_options]
        # Set gridline bools
        xaxis_gridlines, yaxis_gridlines = tree_utils.get_gridline_bools(axis_gridlines)

        # Filter out data not within selected range
        topology_df = topology_df[(topology_df["Window"] >= dataMin) & (topology_df["Window"] <= dataMax)]
        topology_df = topology_df.sort_values(by=["Chromosome", "Window", "TopologyID"])
        topology_df = topology_df.reset_index(drop=True)
        
        # -- Set editable value
        if not editable_graphs:
            editable_graphs = False
        else:
            editable_graphs = True
        
        # Graph
        topology_graphs = []
        if not current_topologies:
            pass
        else:
            # Only keep data for selected topologies
            topology_df_filtered = topology_df[topology_df["TopologyID"].isin(current_topologies)]
            # Create figure depending on toggle-toolbar, add graph to output list
            if "topo_rug_only" in main_graph_switches:
                # Build histogram + Heatmap Figures
                topo_dist_graph = tree_utils.build_rug_plot(
                    topology_df_filtered,
                    chromosome,
                    template,
                    current_topologies,
                    color_mapping,
                    dataRange,
                    topoOrder,
                    axis_line_width,
                    xaxis_gridlines,
                    yaxis_gridlines,
                )
            elif "topo_tile_only" in main_graph_switches:
                # Build histogram + Heatmap Figures
                topo_dist_graph = tree_utils.build_tile_plot(
                    topology_df_filtered,
                    chromosome_df,
                    template,
                    current_topologies,
                    color_mapping,
                    dataRange,
                    window_size,
                    axis_line_width,
                    xaxis_gridlines,
                    yaxis_gridlines,
                )
            else:
                # Build histogram + Heatmap Figures
                topo_dist_graph = tree_utils.build_histogram_with_rug_plot(
                    topology_df_filtered,
                    chromosome,
                    chromosome_df,
                    template,
                    current_topologies,
                    window_size,
                    color_mapping,
                    dataRange,
                    topoOrder,
                    axis_line_width,
                    xaxis_gridlines,
                    yaxis_gridlines,
                )

        if not alt_data_switch:
            pass
        elif "alt_data_switch_on" in alt_data_switch:
            # Collect y-max from additional int/float data
            y_maxes = tree_utils.get_y_max_list(alt_dropdown_options, topology_df)
            # Generate graphs
            for y_max, alt_data_options in zip(y_maxes, alt_dropdown_options):
                alt_graph = tree_utils.build_alt_data_graph(
                    alt_data_options,
                    chromosome_df,
                    color_mapping,
                    topology_df,
                    window_size,
                    template,
                    dataRange,
                    y_max,
                    axis_line_width,
                    xaxis_gridlines,
                    yaxis_gridlines,
                )
                topology_graphs.append(
                    html.Div(
                        dcc.Graph(
                            figure=alt_graph,
                            config=dict(
                                editable=True,
                                displaylogo=False,
                                doubleClick="autosize",
                                modeBarButtonsToRemove=[
                                    "zoomIn2d",
                                    "zoomOut2d",
                                    "autoScale2d",
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
                            className="alt-data-graph",
                        )
                    )
                )

        if not tree_switch:
            pass
        else:
            tree_divs = []
            if (type(current_topologies) == str) or (type(current_topologies) == int):
                current_topologies = [current_topologies]
            if 'dark' in template:
                tree_title_className = "tree-title-dark"
            else:
                tree_title_className = "tree-title-light"

            for topo in current_topologies:
                # Filter data
                curr_topo_df = topology_df[topology_df["TopologyID"] == topo]

                # Grab tree (Should only be one possible tree to use!!!)
                try:
                    first_tree = Tree(curr_topo_df["NewickTree"].unique()[0])
                    first_tree.prune(tree_taxa_values, preserve_branch_length=True)
                # No data for topology on current chromosome
                except IndexError:
                    tree_divs.append(
                        dbc.Col(
                            children=[
                                html.Div(f"{topo}", className=tree_title_className),
                                html.Div(
                                    children=[
                                        dcc.Graph(
                                            figure=tree_utils.no_tree_data(template, "Topology not present on chromosome"),
                                            style={"width": "100%", "height": "40vh", "background-color": "black"},
                                            config=dict(
                                                toImageButtonOptions=dict(
                                                    format=snapshot_file_type,
                                                    filename="Graph_Name",
                                                    height=400,
                                                    width=500,
                                                    scale=3,
                                                ),
                                            ),
                                        )
                                    ],
                                    className="tree-div"
                                ),
                            ],
                            width=4
                        ),
                    )
                    continue
                except TreeError:
                    tree_divs.append(
                        dbc.Col(
                            children=[
                                html.Div(f"{topo}", className=tree_title_className),
                                html.Div(
                                    children=[
                                        dcc.Graph(
                                            figure=tree_utils.no_tree_data(template, "No Taxa Selected"),
                                            style={"width": "100%", "height": "40vh", "background-color": "black"},
                                            config=dict(
                                                toImageButtonOptions=dict(
                                                    format=snapshot_file_type,
                                                    filename="Graph_Name",
                                                    height=400,
                                                    width=500,
                                                    scale=3,
                                                ),
                                            ),
                                        )
                                    ],
                                    className="tree-div"
                                ),
                            ],
                            width=4
                        ),
                    )
                    continue
                except ValueError:
                    # Assumes taxa in dropdown selection
                    # is not found in a particular topology/tree
                    # Solution is to check list and remove taxa
                    # not present in tree
                    tree_taxa = first_tree.get_leaf_names()
                    trimmed_taxa_list = [t for t in tree_taxa_values if t in tree_taxa]
                    first_tree.prune(trimmed_taxa_list, preserve_branch_length=True)
                     

                tree = tree_utils.DrawTree(io.StringIO(first_tree.write()), template, topo, color_mapping)

                # Generate Tree figrues
                fig = tree.create_square_tree()
                # if tree_shape == "box":
                #     fig = tree.create_square_tree()
                # elif tree_shape == "angular":
                #     fig = tree.create_angular_tree()

                if len(tree_taxa_values) > 20:
                    tree_height = "80vh"
                elif len(tree_taxa_values) > 80:
                    tree_height = "120vh"
                else:
                    tree_height = "40vh"
                tree_style = {"width": "100%", "height": tree_height, "background-color": "black"}

                tree_divs.append(
                    dbc.Col(
                        children=[
                            html.Div(f"{topo}", className=tree_title_className),
                            html.Div(
                                children=[
                                    dcc.Graph(
                                        figure=fig,
                                        style=tree_style,
                                        config=dict(
                                            editable=editable_graphs,
                                            toImageButtonOptions=dict(
                                                format=snapshot_file_type,
                                                filename="Graph_Name",
                                                height=400,
                                                width=500,
                                                scale=3,
                                            ),
                                            modeBarButtonsToRemove=[
                                                "lasso2d",
                                                "select2d",
                                            ],
                                            displaylogo=False,
                                        ),
                                    )
                                ],
                                className="tree-div"
                            ),
                        ],
                        width=4,
                    ),
                )
                continue
            topology_graphs.append(
                dbc.Row(
                    children=tree_divs,
                    no_gutters=True
                )
            )

        if not alt_graph_switches:
            pass
        else:
            for graph in alt_graph_switches:
                if graph == "gff_switch":
                    try:
                        gff_df = pd.read_json(gff_data)
                    except ValueError:
                        raise PreventUpdate
                    gff_df = gff_df[(gff_df["start"] >= dataMin) & (gff_df["end"] <= dataMax) & (gff_df["chromosome"] == chromosome)]
                    if abs(dataMin - dataMax) > 50000000:
                        topology_graphs.append(
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
                        )
                    elif len(gff_df) == 0:
                        topology_graphs.append(
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
                        )
                    else:
                        gff_figure = tree_utils.build_gff_figure(
                            gff_df,
                            dataRange,
                            template,
                            axis_line_width,
                            xaxis_gridlines,
                            yaxis_gridlines,
                        )
                        topology_graphs.append(
                            html.Div(
                                children=[
                                    dcc.Graph(
                                        id="gff_graph",
                                        figure=gff_figure,
                                        config=dict(
                                            editable=editable_graphs,
                                            displaylogo=False,
                                            doubleClick="autosize",
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
                        )

        # Config
        # -- Set pixel default
        if (not pixel_height) or (not pixel_width):
            pixel_width = 1500
            pixel_height = 400

        tv_main_dist_config = dict(
            doubleClick="autosize",
            displaylogo=False,
            modeBarButtonsToRemove=[
                "resetScale2d",
                "pan2d",
                "lasso2d",
                "select2d",
                "autoScale2d",
            ],
            toImageButtonOptions=dict(
                format=snapshot_file_type,
                filename="Graph_Name",
                width=int(pixel_width),
                height=int(pixel_height),
            ),
        )
        return topo_dist_graph, topology_graphs, dataRange, tv_main_dist_config


@app.callback(
    [Output("stats-topology-freq-pie-graph", "figure"),
     Output("whole-genome-rf-graph", "figure"),
     Output("stats-topology-freq-pie-graph", "config"),
     Output("whole-genome-rf-graph", "config"),
     Output("stats-table-div", "children"),
    ],
    [Input("input-data-upload", "children"),
     Input("chromFile-data-upload", "children"),
     Input("template-option", "value"),
     Input("topology-color-chart-store", "data"),
     Input("control-tabs", "active_tab"),
     Input("window-size", "children"),
     Input("snapshot-file-option", "value"),
     Input("axis-line-width", "value"),
     Input("tv-pixel-height", "value"),
     Input("tv-pixel-width", "value"),
     Input("editable-graph-option", "value"),],
    [
        State("topology_options", "value"),
    ],
)
def whole_genome_overview_plots(
    data_json,
    chromlength_data,
    template,
    color_mapping,
    active_tab,
    window_size,
    snapshot_file_type,
    axis_line_width,
    pixel_height, 
    pixel_width,
    editable_graphs,
    topologies,
):

    def get_topo_trees(df, top_3_topos):
        trees = []
        for t in top_3_topos:
            topo_df = df[df["TopologyID"] == t]
            topo_df.reset_index(drop=True, inplace=True)
            trees.append(topo_df["NewickTree"][0])
        return trees

    def calc_rf_dist(ref, alt):
        ref = Tree(ref)
        alt = Tree(alt)
        try:
            compare_results = ref.compare(alt)
        except TreeError:
            compare_results = ref.compare(alt, unrooted=True)
        return compare_results["norm_rf"]

    if active_tab == "home-tab":
        return dash.no_update, dash.no_update, None, None, None
    elif not data_json:
        raise PreventUpdate
    elif not chromlength_data:
        raise PreventUpdate
    else:
        # Set editable value
        if not editable_graphs:
            editable_graphs = False
        else:
            editable_graphs = True
        # Load session data
        df = pd.read_json(data_json)
        chrom_lengths = pd.read_json(chromlength_data)
        # Create pie chart
        topo_freq_df = pd.DataFrame(df["TopologyID"].value_counts()/len(df))
        topo_freq_df = topo_freq_df.reset_index()
        topo_freq_df.columns = ['TopologyID', 'Frequency']
        topo_freq_df.sort_values(by="Frequency")
        fig_pie = tree_utils.build_topology_frequency_pie_chart(topo_freq_df, template, color_mapping)

        # Create RF graph
        topos = topo_freq_df["TopologyID"][0:len(topologies)]
        topo_trees = get_topo_trees(df, topos)
        rf_df = pd.DataFrame(
            {"TopologyID": topos[1:], "normRF-Distance": [calc_rf_dist(topo_trees[0], t) for t in topo_trees[1:]]}
        )
        rf_fig = tree_utils.build_rf_graph(rf_df, topos[0], template, color_mapping, axis_line_width)

        # Create stats datatable
        basic_info_df = tree_utils.make_basic_info_df(df, topo_freq_df, chrom_lengths, window_size)
        stats_table = dbc.Table.from_dataframe(basic_info_df, striped=True, bordered=True, hover=False, size='md')
        # Configs
        if (not pixel_height) or (not pixel_width):
            pixel_width = 750
            pixel_height = 400
        config = dict(
            doubleClick="autosize",
            displaylogo=False,
            editable=editable_graphs,
            modeBarButtonsToRemove=[
                "zoomIn2d",
                "zoomOut2d",
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
        # return fig_pie, rf_fig, config, config, columns, data
        return fig_pie, rf_fig, config, config, stats_table


@app.callback(
    Output("whole-genome-graphs", "children"),
    [Input("input-data-upload", "children"),
     Input("chromFile-data-upload", "children"),
     Input("template-option", "value"),
     Input("topology-color-chart-store", "data"),
     Input("stats-chrom-plot-type", 'value'),
     Input("topology-freq-order-store", "data"),
     Input("topology_options", "value"),
     Input("window-size", "children"),
     Input("control-tabs", "active_tab"),
     Input("snapshot-file-option", "value"),
     Input("axis-line-width", "value"),
     Input("tv-pixel-height", "value"),
     Input("tv-pixel-width", "value"),
     Input("axis-gridlines", "value"),
     Input("editable-graph-option", "value"),
    ],
)
def whole_genome_bottom_plots(
    data_json,
    chromosome_lengths,
    template,
    color_mapping,
    chrom_plot_type,
    topoOrder,
    current_topologies,
    window_size,
    active_tab,
    snapshot_file_type,
    axis_line_width,
    pixel_height,
    pixel_width,
    axis_gridlines,
    editable_graphs,
):
    if active_tab == "home-tab":
        return None
    elif not data_json:
        raise PreventUpdate
    elif not chromosome_lengths:
        raise PreventUpdate
    else:
        # Set editable value
        if not editable_graphs:
            editable_graphs = False
        else:
            editable_graphs = True
        # Read in json data
        df = pd.read_json(data_json)
        chromosome_df = pd.read_json(chromosome_lengths)
        # Set whole genome graph collection list
        whole_genome_graphs = []
        # Set up graph config
        if (not pixel_height) or (not pixel_width):
            pixel_width, pixel_height = 1500, 1123
        config = dict(
            doubleClick="autosize",
            displaylogo=False,
            editable=editable_graphs,
            modeBarButtonsToRemove=[
                "zoomin2D"
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
        # Set gridline bools
        xaxis_gridlines, yaxis_gridlines = tree_utils.get_gridline_bools(axis_gridlines)
        
        # Check if there are more than 100,000 rows, 
        # if so return too many data point graph
        if len(df) > 100000:
            fig = tree_utils.no_tree_data(template, "Too many data points to plot")
            whole_genome_graphs.append(
                dbc.Col([
                    dcc.Graph(
                        figure=fig,
                        config=config,
                        className="stats-2-graph-style",
                    )
                ], width=12)
            )
            return whole_genome_graphs
        df["GenomeFrequency"] = df["TopologyID"].map(df["TopologyID"].value_counts()/len(df))
        df_grouped = df.groupby(by="Chromosome")
        # Round up chromosome lengths for graph range
        chromosome_df["End"] = chromosome_df["End"].apply(lambda x: data_utils.roundUp(x, window_size))
        # Remove topologies not currently chosen
        topoOrder = [e for e in topoOrder if e in current_topologies]
                
        # Group Chromosomes into groups of 20
        chrom_groups = list(tree_utils.mygrouper(20, df["Chromosome"].unique()))
        if chrom_plot_type == 'bar':
            topoFreqTable = tree_utils.make_topo_freq_table(df_grouped)
            # Iterate through chromosome groups + append graph
            for chromGroup in chrom_groups:
                fig = tree_utils.build_stats_chromosome_freq_bar_plot(
                    topoFreqTable,
                    template,
                    color_mapping,
                    current_topologies,
                    axis_line_width,
                    chromGroup,
                    xaxis_gridlines,
                    yaxis_gridlines,
                )
                whole_genome_graphs.append(
                    dbc.Col([
                        dcc.Graph(
                            figure=fig,
                            config=config,
                            className="stats-2-graph-style",
                        )
                    ], width=12)
                )
        elif chrom_plot_type == 'pie':
            topoFreqTable = tree_utils.make_topo_freq_table(df_grouped)
            # Iterate through chromosome groups + append graph
            for chromGroup in chrom_groups:
                fig = tree_utils.build_stats_chromosome_freq_pie_plot(
                    topoFreqTable,
                    template,
                    color_mapping,
                    chromGroup,
                )
                whole_genome_graphs.append(
                    dbc.Col([
                        dcc.Graph(
                            figure=fig,
                            config=config,
                            className="stats-2-graph-style",
                        )
                    ], width=12)
                )
        elif chrom_plot_type == 'rug':
            # Iterate through chromosome groups + append graph
            for chromGroup in chrom_groups:
                fig = tree_utils.build_stats_chromosome_freq_rug_plot(
                    df,
                    chromosome_df,
                    chromGroup,
                    template,
                    color_mapping,
                    current_topologies,
                    topoOrder,
                    window_size,
                    axis_line_width,
                    xaxis_gridlines,
                    yaxis_gridlines,
                )
                whole_genome_graphs.append(
                    dbc.Col([
                        dcc.Graph(
                            figure=fig,
                            config=config,
                            className="stats-2-graph-style",
                        )
                    ], width=12)
                )
        elif chrom_plot_type == 'tile':
            # Iterate through chromosome groups + append graph
            for chromGroup in chrom_groups:
                fig = tree_utils.build_stats_chromosome_freq_tile_plot(
                    df,
                    chromosome_df,
                    template,
                    color_mapping,
                    current_topologies,
                    topoOrder,
                    window_size,
                    axis_line_width,
                    chromGroup,
                    xaxis_gridlines,
                    yaxis_gridlines,
                )
                whole_genome_graphs.append(
                    dbc.Col([
                        dcc.Graph(
                            figure=fig,
                            config=config,
                            className="stats-2-graph-style",
                        )
                    ], width=12)
                )
        return whole_genome_graphs


# # ---------------------- Export Option Callbacks ---------------------
@app.callback(
    [Output("tv-current-view-download", "data"),
     Output("gff-current-view-download", "data"),],
    [Input("curr-view-export-button", "n_clicks"),
     Input("input-data-upload", "children"),
     Input("chromosome_options", "value"),
     Input("gff-data-upload", "children"),],
    [State('current-data-range', 'data'),],
    prevent_initial_call=True,
)
def current_window_summary(n, data_json, curr_chromosome, gff_json, relayout_data):
    if not data_json:
        raise PreventUpdate
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if not button_id:
        raise PreventUpdate
    elif button_id == "curr-view-export-button":
        dataMin, dataMax = relayout_data
        df = pd.read_json(data_json)
        df = df[(df["Window"] >= dataMin) & (df["Window"] <= dataMax) & (df["Chromosome"] == curr_chromosome)]
        if not gff_json:
            gff_df = pd.DataFrame()
        else:
            gff_df = pd.read_json(gff_json)
            gff_df = gff_df[(gff_df["start"] >= dataMin) & (gff_df["start"] <= dataMax) & (gff_df["chromosome"] == curr_chromosome)]
        
        # Return according to available data
        if (len(df) == 0) and (len(gff_df) > 1):
            return None, dcc.send_data_frame(gff_df.to_csv, "current_view_annotations.gff", sep="\t", index=False)
        elif (len(df) > 1) and (len(gff_df) == 0):
            return dcc.send_data_frame(df.to_excel, "current_view_data.xlsx", index=False), None
        else:
            return dcc.send_data_frame(df.to_excel, "current_view_data.xlsx", index=False), dcc.send_data_frame(gff_df.to_csv, "current_view_annotations.gff", sep="\t", index=False, header=False)
    else:
        raise PreventUpdate


@app.callback(
    [Output("tree-prune-download", "data"),
     Output("tree-input-modal", "is_open"),
     Output("tree-prune-export-dropdown", "options"),
     Output("topobinner-loading", "children"),],
    [Input("tree-prune-export-button", "n_clicks"),
     Input("tree-prune-submit-button", "n_clicks"),
     Input("tree-prune-close-button", "n_clicks"),
     Input("input-data-upload", "children"),
     Input("tree-taxa-choices", "options"),
     Input("tree-prune-export-dropdown", "value"),
     Input("topobinner-checkbox", "value")],
    prevent_initial_call=True,
)
def tree_pruning_export(
    btn,
    btn2,
    btn3,
    data_json,
    taxa_options,
    prune_taxa_choices,
    topobinner,
):
    if not data_json:
        raise PreventUpdate
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if not button_id:
        raise PreventUpdate
    elif button_id == "tree-prune-export-button":
        return None, True, taxa_options, None
    elif button_id == "tree-prune-close-button":
        return None, False, taxa_options, None
    elif button_id == "tree-prune-submit-button":
        df = pd.read_json(data_json)
        # Clear TopologyID column
        df["TopologyID"] = [pd.NA]*len(df)
        # Prune trees with selected taxa
        df["NewickTree"] = df["NewickTree"].apply(lambda x: tree_utils.prune_tree(x, prune_taxa_choices))
        # Add comment row on top to track what was done -- Not implemented yet
        # taxa_join = "-".join(prune_taxa_choices)
        # df.columns = [f"# TreePrune: {taxa_join}"] + [pd.NA]*(len(df.columns)-1)
        # Run Topobinner if checked
        if topobinner == ["run"]:
            df = tree_utils.tv_topobinner(df)

        return (dcc.send_data_frame(df.to_excel, "pruned_TreeViewer_input.xlsx", index=False), 
                False, dash.no_update, None)
    else:
        raise PreventUpdate


# ---------------------- Real Time Interactive Callbacks ---------------------
@app.callback(
    Output("rfdist-div", "children"),
    [Input("topology_graph", "hoverData"),],
)
def calc_NormRF_distance(
    hoverdata,
):
    print("HOVER WORKING", hoverdata)
    return None

# @app.callback(
#     Output("rfdist-div", "children"),
#     [Input("alt-graph-switches", "value"),
#      Input("topology_graph", "hoverData"),
#      Input("chromosome_options", "value"),
#      Input("template-option", "value"),
#      Input("snapshot-file-option", "value"),
#      Input("topology-color-chart-store", "data"),
#      Input("topology_options", "value"),
#      Input("input-data-upload", "children"),
#     ],
# )
# def calc_NormRF_distanc(
#     alt_graph_options,
#     hoverdata,
#     curr_chrom,
#     template,
#     snapshot_file_type,
#     color_mapping,
#     topoology_options,
#     data,
# ):
#     rf_plot = []
#     if not data:
#         raise PreventUpdate
#     elif not alt_graph_options:
#         return None
#     elif "rf_dist" not in alt_graph_options:
#         return None
#     elif not hoverdata:
#         rf_plot.append(
#             dbc.Row(
#                 children=[
#                     dbc.Col(
#                         children=[
#                             html.Div(
#                                 children=[
#                                     dcc.Graph(
#                                         figure=tree_utils.init_RF_graph(
#                                             template),
#                                         style={
#                                             "height": "25vh", "width": "100%", "background-color": "black"},
#                                         config=dict(
#                                             toImageButtonOptions=dict(
#                                                 format=snapshot_file_type,
#                                                 filename="Graph_Name",
#                                             ),
#                                         ),
#                                     )
#                                 ],
#                                 style={
#                                     "width": "100%",
#                                 }
#                             ),
#                         ],
#                         width=12
#                     ),
#                 ],
#                 no_gutters=True
#             )
#         )
#         return rf_plot
#     else:
#         try:
#             input_data = pd.read_json(data)
#             input_data = input_data.drop(
#                 columns=input_data.columns.to_list()[4:])
#             x_pos = tree_utils.get_RFxpos(hoverdata, input_data)
#             chrom_info = input_data[(input_data["Chromosome"] == curr_chrom)]
#             chrom_info = chrom_info.sort_values(by="Window")
#             # chrom_info = chrom_info.reset_index(drop=True)
#             poi_info = chrom_info[chrom_info["Window"] == x_pos]
#             poi = poi_info.index.to_list()[0]
#             rfdist_df = chrom_info[(chrom_info.index >= poi - 10) &
#                                    (chrom_info.index <= poi + 10) & (chrom_info.index != poi)]
#             rfdist_df = pd.concat([poi_info, rfdist_df])
#             rfdist_df = rfdist_df.sort_values(by="Window")
#             rfdist_df = rfdist_df.reset_index(drop=True)
#             t1 = poi_info["NewickTree"].to_list()[0]
#             rf_dist_collect = []
#             for t2 in rfdist_df["NewickTree"]:
#                 try:
#                     rfd = tree_utils.RFDistance(t1, t2).NormRF()
#                 except TreeError:
#                     raise PreventUpdate

#                 if type(rfd) == str:
#                     rfd = np.nan
#                 rf_dist_collect.append(rfd)
#                 continue
#             rfdist_df["NormRF"] = rf_dist_collect
#             heatmap = tree_utils.build_RF_heatmap(
#                 rfdist_df,
#                 template,
#                 color_mapping,
#             )
#             rf_plot.append(
#                 html.Div(
#                     children=[
#                         dcc.Graph(
#                             figure=heatmap,
#                             style={"height": "25vh", "width": "100%",
#                                    "background-color": "black"},
#                             config=dict(
#                                 toImageButtonOptions=dict(
#                                     format=snapshot_file_type,
#                                     filename="Graph_Name",
#                                 ),
#                             ),
#                         )
#                     ],
#                     # style={
#                     #     "width": "100%",
#                     # }
#                 )
#             )
#             return rf_plot
#         except IndexError:
#             return None



