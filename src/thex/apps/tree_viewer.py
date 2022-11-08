import base64
import io

from pathlib import Path

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt
import dash_daq as daq
# import numpy as np
import pandas as pd
import scipy
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash_table.Format import Format, Scheme
from ete3 import Tree
from ete3.coretype.tree import TreeError
from ete3.parser.newick import NewickError

from configparser import MissingSectionHeaderError

from thex.app import app
from thex.apps.docs import docs_treeviewer_contents
from thex.apps.utils import data_utils, graph_options
from thex.apps.utils import tree_utils


############################### Graph Options ###############################
GRAPH_TEMPLATES = graph_options.graph_templates()
COLOR_SWATCHES = graph_options.color_swatches()
SCALE_OPTIONS = graph_options.figure_output_scales()
SNAPSHOT_FILE_OPTIONS = graph_options.snapshot_file_type()
FONT_FAMILIES = graph_options.font_families()

############################### Documentation Components ###############################
    
docs_content = html.Div(id="tv-docs-content", className="tv-docs-body")

navbar = dbc.Navbar(
    [
        dbc.Col(
            children=[
                dbc.Button(
                    id="home-button",
                    children=[
                        html.Img(src=app.get_asset_url("tree5.svg"), height='40px'),
                    ],
                    className="logo-png",
                    href="/",
                    outline=True,
                ), 
            ],
            width='auto'
        ),
        dbc.Col(dbc.NavbarBrand("Tree Viewer", className="m-2"), width='auto'), 
        dbc.Col([
            dbc.DropdownMenu(
                children=[
                    dbc.DropdownMenuItem("Import Options", header=True, className="file-menu-header"),
                    dbc.DropdownMenuItem("Load File/Project", id="new-session-button"),
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
            # Stats
            dbc.NavItem(
                dbc.NavLink(
                    f"Summary/Stats",
                    id=f"toggle-statsbar",
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
        ]),
    ],
    color="orange",
    expand=True,
    sticky='top',
    style={'width': '110%'}
)

tree_viewer_file_input_button = dbc.InputGroup(
    [
        dbc.InputGroupAddon(
            children=[
                dcc.Upload(
                    id="upload-data",
                    children=[
                        dbc.Button(
                            id="tv-button",
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
                            id="gff-new-session-button",
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

project_input_button = dbc.InputGroup(
    [
        dbc.InputGroupAddon(
            children=[
                dcc.Upload(
                    id="new-session-project",
                    children=[
                        dbc.Button(
                            id="new-session-project-button",
                            children=["Load Project"], 
                            className="tv-input-buttons",
                            color="info"
                        )
                    ]
                ),
            ],
            addon_type="prepend",
        ),
        html.Div(id="project-name-div", className="tv-filename"),
    ],
    size='lg',
)

correlation_input_form = dbc.Form(
    [
        # Significance test
        dbc.FormGroup(
            [
                dbc.Label("Significance Test:", className="mr-2"),
                dcc.Dropdown(
                    id='correlation-test-type',
                    options=[
                        {'label': 'Kruskal-Wallis H-test', 'value': 'Kruskal-Wallis H-test'},
                        {'label': 'One-Way ANOVA', 'value': 'One-Way ANOVA'}
                    ],
                    value="Kruskal-Wallis H-test",
                    className="dropdown-style",
                ),
            ],
        ),
        # Post hoc test
        dbc.FormGroup(
            [
                dbc.Label("Post Hoc Test:", className="mr-2"),
                dcc.Dropdown(
                    id='correlation-posthoc-type',
                    options=[
                        {'label': "Dunn's test", 'value': "Dunn's test"},
                        {'label': 'TukeyHSD', 'value': 'TukeyHSD'},
                        # {'label': 'Mann-Whitney rank test', 'value': 'Mann-Whitney rank test'},
                    ],
                    value="Dunn's test",
                    className="dropdown-style",
                ),
            ],
        ),
        # p-value adjustment
        dbc.FormGroup(
            [
                dbc.Label("p-value Adjustment:", className="mr-2"),
                dcc.Dropdown(
                    id='correlation-pval-adjustment-type',
                    options=[
                        {'label': 'Bonferroni', 'value': 'bonferroni'},
                        {'label': 'Benjamini/Hochberg', 'value': 'fdr_bh'},
                        # {'label': 'One step correction - Sidak', 'value': 'sidak'},
                        # {'label': 'Step-down method using Sidak adjustments', 'value': 'holm-sidak'},
                        # {'label': 'Step-down method using Bonferroni adjustments', 'value': 'holm'},
                        # {'label': 'Step-up method (independent) - simes-hochberg', 'value': 'simes-hochberg'},
                        # {'label': 'Closed method based on Simes tests (non-negative)', 'value': 'hommel'},
                        # {'label': 'Benjamini/Yekutieli (negative)', 'value': 'fdr_by'},
                        # {'label': 'Two stage fdr correction (non-negative) - fdr_tsbh', 'value': 'fdr_tsbh'},
                        # {'label': 'Two stage fdr correction (non-negative) - fdr_tsbky', 'value': 'fdr_tsbky'},
                    ],
                    value="bonferroni",
                    className="dropdown-style",
                ),
            ],
        ),
        # Additional Data
        dbc.FormGroup(
            [
                dbc.Label("Additional Data:", className="mr-2"),
                dcc.Dropdown(
                    id='correlation-add-data-type',
                    className="dropdown-style",
                ),
            ],
        ),
        # Topologies
        dbc.FormGroup(
            [
                dbc.Label("Topologies:", className="mr-2"),
                dcc.Dropdown(
                    id='correlation-topologies',
                    className="dropdown-style",
                    multi=True,
                ),
            ],
        ),
        # Range
        dbc.FormGroup(
            [
                dbc.Label("Range:"),
                dcc.Dropdown(
                    id='correlation-range-selection',
                    className="dropdown-style",
                    options=[
                        {'label': 'Current Chromosome', 'value': 'Current Chromosome'},
                        {'label': 'Current View', 'value': 'Current View'},
                        {'label': 'Whole Genome', 'value': 'Whole Genome'},
                    ],
                    value='Current Chromosome'
                ),
            ],
        ),
        # Alpha
        dbc.FormGroup(
            [
                dbc.Label("Alpha:"),
                dcc.Input(
                    id='alpha-input',
                    type='number',
                    value=0.05,
                    min=0.0,
                    max=1.0,
                ),
            ],
        ),
        dbc.Button("Run", color="primary", id="correlation-run-btn", style={'width': "100%"}),
    ],
)

# -- Toggle Switches --
chrom_genome_toggle = daq.ToggleSwitch(
    id="toggle-chrom-whole-genome",
    value=False,
    label=[{'label': 'View | Genome', 'style': {"font-size": '12pt'}}, {'label': 'Chromosome', 'style': {'font-size': '12pt'}}],
    color="#375a7f",
    style={"color": "#black", "padding-bottom": "2px"},
    size=30,
    theme='dark',
)

additional_data_toggle = daq.ToggleSwitch(
    id="alt-data-switch",
    disabled=False,
    value=False,
    label=[{'label': 'Additional Data | Off ', 'style': {"font-size": '12pt'}}, {'label': 'On', 'style': {"font-size": '12pt'}}],
    color="#375a7f",
    style={"color": "#black", "padding-bottom": "2px"},
    size=30,
    theme='dark',
)

tree_graph_toggle = daq.ToggleSwitch(
    id="tree-graph-switch",
    disabled=False,
    value=False,
    label=[{'label': 'Trees | Off ', 'style': {"font-size": '12pt'}}, {'label': 'On', 'style': {"font-size": '12pt'}}],
    color="#375a7f",
    style={"color": "#black", "padding-bottom": "2px"},
    size=30,
)

tree_styles = dbc.RadioItems(
    id='tree-shape',
    className="btn-group",
    inputClassName="btn-check",
    labelClassName="btn-outline-secondary btn btn-sm",
    labelCheckedClassName="active",
    options=[
        {'label': 'Rectangular', 'value': 'rectangle'},
        {'label': 'Circular', 'value': 'circle'},
        # {'label': 'Radial', 'value': 'radial'},
    ],
    value='rectangle',
    style={'paddingBottom': '5px', 'color': 'black'},
)

topology_order_toggle = daq.ToggleSwitch(
    id="topo-freq-order",
    disabled=False,
    value=False,
    label=[{'label': 'Topologies | Genome ', 'style': {"font-size": '12pt'}}, {'label': 'Local', 'style': {"font-size": '12pt'}}],
    color="#375a7f",
    style={"color": "black", "padding-bottom": "2px"},
    size=30,
    theme='dark',
)

quantile_graph_toggle = daq.ToggleSwitch(
    id="quantile-toggle",
    disabled=True,
    value=False,
    label=[{'label': 'Quantile Graph | Off', 'style': {"font-size": '12pt'}}, {'label': 'On', 'style': {"font-size": '12pt'}}],
    color="#375a7f",
    style={"color": "black", 'padding-bottom': '2px'},
    size=30,
    theme='dark',
)

# -- Toolbar graph switch col --
graph_switch_col1 = dbc.Form(
    children=[
        dbc.Label("Main Distribution Graph Options:", className="text-divs-sm"),
        dbc.Row([
            dbc.Label("Type:", style={'color': 'black', 'paddingLeft': '15px'}, width='auto'),
            dbc.Col([
                dbc.RadioItems(
                    id="main-graph-switches",
                    className="btn-group",
                    inputClassName="btn-check",
                    labelClassName="btn btn-outline-secondary",
                    labelCheckedClassName="active",
                    options=[
                        {'label': 'Rug+Tile', 'value': 'topo_both', 'disabled': True},
                        {'label': 'Rug', 'value': 'topo_rug_only', 'disabled': True},
                        {'label': 'Tile', 'value': 'topo_tile_only', 'disabled': True},
                    ],
                    value="topo_both",
                ),
            ], className="radio-group"),
        ], style={"paddingBottom": "5px"}),
        dbc.Label("Whole Genome Graph Options:", className="text-divs-sm"),
        dbc.Row([
            dbc.Label("Type:", style={'color': 'black', 'paddingLeft': '15px'}, width='auto'),
            dbc.Col([
                dbc.RadioItems(
                    id='whole-genome-graph-type',
                    className="btn-group",
                    inputClassName="btn-check",
                    labelClassName="btn btn-outline-secondary",
                    labelCheckedClassName="active",
                    options=[
                        {'label': 'Rug', 'value': 'rug', 'disabled': True},
                        {'label': 'Bar', 'value': 'bar', 'disabled': True},
                        {'label': 'Pie', 'value': 'pie', 'disabled': True},
                        {'label': 'Tile', 'value': 'tile', 'disabled': True},
                    ],
                    value='rug',
                ),
            ], className="radio-group"),
        ], style={"paddingBottom": "5px"}),
        dbc.Row([
            dbc.Label("Height:", style={'color': 'black', 'paddingLeft': '15px'}, width='auto'),
            dbc.Col([
                dbc.RadioItems(
                    id="wg-squish-expand",
                    className="btn-group",
                    inputClassName="btn-check",
                    labelClassName="btn btn-outline-secondary",
                    labelCheckedClassName="active",
                    options=[
                        {'label': 'Collapse', 'value': 'collapse', 'disabled': True},
                        {'label': 'Squish', 'value': 'squish', 'disabled': True},
                        {'label': 'Expand', 'value': 'expand', 'disabled': True},
                    ],
                    value='squish',
                ),
            ], className="radio-group"),
        ], style={"paddingBottom": "5px"}),
        dbc.Row([
            dbc.Label("# chr per graph:", style={'color': 'black', 'paddingLeft': '15px'}, width='auto'),
            dbc.Col([
                dcc.Input(
                    id='whole-genome-per-graph-chrom-count',
                    type="number",
                    style={'width': '50px'},
                ),
            ]),
        ], style={"paddingBottom": "5px"}),
    ],
)

# -- Topology-quantile components --
topology_quantile_input_form = dbc.Form(
    [
        # Additional Data
        dbc.FormGroup(
            [
                dbc.Label("Additional Data:", className="mr-2"),
                dcc.Dropdown(
                    id='topo-quantile-additional-data',
                    className="dropdown-style",
                ),
            ],
        ),
        # Topology
        dbc.FormGroup(
            [
                dbc.Label("Topologies:", className="mr-2"),
                dcc.Dropdown(
                    id='topo-quantile-topologies',
                    className="dropdown-style",
                    multi=True,
                ),
            ],
        ),
        # Chromosome selection
        dbc.FormGroup(
            [
                dbc.Label("Chromosome:", className="mr-2"),
                dcc.Dropdown(
                    id='topo-quantile-chrom-dropdown',
                    className="dropdown-style",
                ),
            ],
        ),
        # Sex chrom selection
        dbc.FormGroup(
            [
                dbc.Label("Sex Chromosome(s):", className="mr-2"),
                dcc.Dropdown(
                    id='sex-chrom-dropdown',
                    className="dropdown-style",
                    multi=True
                ),
            ],
        ),
        # n-Quantiles
        dbc.FormGroup(
            [
                dbc.Label("n-Quantiles:", className="mr-2"),
                dcc.Input(
                    id='topo-n-quantiles',
                    type='number',
                    placeholder='default: 4',
                ),
            ],
        ),
        # Range
        dbc.FormGroup(
            [
                dbc.Label("Range:"),
                dbc.RadioItems(
                        id='topo-quantile-range-selection',
                        options=[
                            {'label': 'Single Chromosome', 'value': 'per_chrom'},
                            {'label': 'Whole Genome', 'value': 'whole_genome'},
                            {'label': 'Autosomes vs Sex Chromosomes', 'value': 'a_vs_x'}
                        ],
                        value='per_chrom',
                    ),
            ],
        ),
        dbc.Button("Run", color="primary", id="topo-quantile-run-btn", style={'width': "100%"}),
    ],
)

tab_style = {
    'borderTop': '1px solid #d6d6d6',
    'backgroundColor': '#637589',
    'color': 'white',
}

tab_selected_style = {
    'borderTop': '2px solid orange',
    'backgroundColor': '#375a7f',
    'color': 'white',
}

color_step_buttons = html.Div([
    dbc.Button("<", id='step-down-button', size='sm'),
    html.Div(["Shift"], style={'margin': '5px 5px 5px 5px', 'color': 'black'}),
    dbc.Button(">", id='step-up-button', size='sm'),
], style={'display': 'inline-flex'})

######################################  Layout  #####################################
def tv_layout():
    layout = dbc.Container(
        children=[
            # --- Data Storage ---
            dcc.Store(id="topology-color-chart-store"),
            dcc.Store(id="current-color-list"),
            dcc.Store(id="topology-freq-order-store"),
            dcc.Store(id="current-data-range"),
            dcc.Store(id="range-init"),
            dcc.Store(id="project-dir", data=None),
            
            # --- Triggers ---
            html.Div(id='main-div-trigger'),
            html.Div(id='init-graph-trigger'),
            html.Div(id='break-init-graph-trigger'),
            # --- Data Handoffs ---
            html.Div(id='gff-data-handoff', className="hidden-div"),
            # Empty DIVs for data chache
            html.Div(id="chromFile-data-upload", className="hidden-div"),
            html.Div(id="gff-data-upload", className="hidden-div"),
            html.Div(id="window-size", className="hidden-div"),
            html.Div(id="input-place-holder", className="hidden-div"),
            # Download Components
            dcc.Download(id="tv-current-view-download"),
            dcc.Download(id="gff-current-view-download"),
            dcc.Download(id="tree-prune-download"),
            dcc.Download(id='project-ini-download'),
            dcc.Download(id='tree-viewer-download'),
            dcc.Download(id='chrom-length-download'),
            # --- Tooltips ---
            # Topology Freq order toggle
            dbc.Tooltip(
                "Order topologies by whole-genome frequency or current view",
                target="topo-freq-order",
                placement="top",
            ),
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
                "Return to homepage",
                target="home-button",
                delay={"show": 500, "hide": 250}
            ),
            # --- Modals ---
            # Input modal
            dbc.Modal(
                children=[
                    dbc.ModalBody(
                        children=[
                            dcc.Tabs(
                                id='input-modal-tabs',
                                value='individual-input-tab',
                                children=[
                                    # Single file input
                                    dcc.Tab(
                                        value='individual-input-tab',
                                        label="Individual File Upload",
                                        style=tab_style,
                                        selected_style = tab_selected_style,
                                        children=[
                                            dbc.Row([
                                                dbc.Col(html.Br(), width=12),
                                                dbc.Col([
                                                    dbc.Card(
                                                        [
                                                            dbc.CardBody(
                                                                html.Div([
                                                                    """
                                                                    Only upload a Tree Viewer input file and chromosome length bed file to quickly visualize data sets without tracking important information throughout the session. 
                                                                    Additional files like gene annotations (GFF/GTF) can still be added later during the session through the "File" dropdown menu.
                                                                    """,
                                                                ], style={'text-align': 'left', 'color': 'black', 'font-size': '1em'}),
                                                            ),
                                                            dbc.CardFooter([
                                                                dbc.Button("Tree Viewer Input File Template", id='tree-viewer-download-button', color='primary', style={'margin-right': '10px'}),
                                                                dbc.Button("Chromosome Length Bed File Template", id='chrom-length-download-button', color='primary'),
                                                            ]),
                                                        ],
                                                        color='light',
                                                    ),
                                                ], width=12),
                                                dbc.Col(html.Br(), width=12),
                                                dbc.Col([
                                                    tree_viewer_file_input_button,
                                                    chromosome_length_input_button,
                                                ], width=12),
                                                dbc.Col(html.Br(), width=12),
                                            ]),
                                        ],
                                    ),
                                    # Project
                                    dcc.Tab(
                                        value='project-input-tab',
                                        label="Project Upload",
                                        style=tab_style,
                                        selected_style = tab_selected_style,
                                        children=[
                                            dbc.Row([
                                                dbc.Col(html.Br(), width=12),
                                                dbc.Col([
                                                    dbc.Card(
                                                        [
                                                            dbc.CardBody([
                                                                html.P([
                                                                    """
                                                                    Choose a project ".ini" file to load all required and additional input files at once. This simplifies the input process
                                                                    and allows for fast switching between data sets. Projects also enable users to upload and visualize multiple GFF/GTF files 
                                                                    simultaneously by providing multiple file paths separated by ";". Click the button below to download a template ".ini" project file
                                                                    that follows the required file format.
                                                                    """,
                                                                ], style={'text-align': 'left', 'color': 'black', 'font-size': '1em'}),
                                                                # html.P([
                                                                #     """
                                                                #     ** In Development ** 
                                                                #     """,
                                                                # ], style={'text-align': 'left', 'color': 'black', 'font-size': '1em'}),
                                                                # html.P([
                                                                #     """
                                                                #     Provided a project directory is given within the project ".ini" file, changes made during the session will be logged to the project directory.
                                                                #     This includes information about pruned files, current view export coordinates, and other important changes. No user specific information is logged,
                                                                #     including but not limited to absolute file pathways, environment information, or user computing information. 
                                                                #     """,
                                                                # ], style={'text-align': 'left', 'color': 'black', 'font-size': '1em'}),
                                                            ]),
                                                            dbc.CardFooter([
                                                                dbc.Button("Download project.ini Template", id='project-ini-download-button', color='primary')
                                                            ]),
                                                        ],
                                                        color='light',
                                                    ),
                                                ], width=12),
                                                dbc.Col(html.Br(), width=12),
                                                dbc.Col([
                                                    project_input_button,
                                                ], width=12),
                                            ]),
                                        ],
                                    ),
                                ],
                            ),
                        ],
                    ),
                    dbc.ModalFooter(
                        children=[
                            dbc.Col([html.H5(id="tv-modal-warning-label")]),
                            dcc.Loading(
                                [html.Div(id="input-data-upload", className="hidden-div")],
                                className="toolbar-loading",
                            ),
                            dbc.Button("Submit", id="input-files-submit-button",
                                       className="submit-buttons", color="info"),
                            dbc.Button("Close", id="input-files-close-button",
                                       className="submit-buttons", color="info"),
                        ],
                    ),
                ],
                id="main-input-modal",
                size="xl",
                keyboard=False,
                is_open=True,
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
            # GFF/GTF malform error
            dbc.Modal(
                children=[
                    dbc.ModalHeader(
                        children=[
                            dbc.Row([
                                dbc.Label("GFF/GTF Error", className="modal-label"),
                            ], justify="center", align="center"),
                        ],
                        className="modal-header"
                    ),
                    dbc.ModalBody(
                        children=[
                            dbc.Row(
                                "GFF/GTF file appears to malformed - please verify input and rerun",
                                no_gutters=True,
                                justify='center',
                            ),
                        ],
                    ),
                ],
                id="gff-error-modal",
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
                                                        The file pruning export option allows you to create a subset of your Tree Viewer input file with each Newick tree pruned by a selection of taxa chosen in the dropdown below.
                                                        This process clears the current TopologyID column and prunes each tree, enabling you to re-bin the new tree topologies into novel bins.
                                                        There is an option to run the binning process before downloading the new file, but please note that re-binning large trees across large genomes (e.g., >2 Gbp) may require 
                                                        considerable run time and will lock your session until it completes. In such circumstances, we highly recommend running the --topobinner command in THExBuilder.
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
            # Too few topology error for stats tests
            dbc.Modal(
                children=[
                    dbc.ModalBody(html.H4(["Two or more topologies are required to run!"], style={"text-align": 'center'}))
                ],
                id="correlation-topo-alert",
                size="md",
                is_open=False,
            ),
            # --- Page Layout ---
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
                                                dbc.Col([
                                                    # --- Graph Color Theme ---
                                                    dbc.Col(
                                                        children=[
                                                            html.H6(
                                                                children=[
                                                                    "Color Theme:"],
                                                                style={'color': 'black'},
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
                                                        width=12,
                                                    ),
                                                    # --- Data Color Palette ---
                                                    dbc.Col(
                                                        children=[
                                                            html.H6(
                                                                children=["Color Palette:"],
                                                                style={'color': 'black'},
                                                            ),
                                                            html.Div([
                                                                dcc.Dropdown(
                                                                    id="line_color",
                                                                    options=[
                                                                        {"label": f"{i} - {len(COLOR_SWATCHES[i])} colors", "value": i} for i in COLOR_SWATCHES],
                                                                    value="Plotly",
                                                                    clearable=False,
                                                                    style={'color': 'black', 'width': '98%'}
                                                                ),
                                                                color_step_buttons,
                                                            ], style={'display': 'inline-flex', 'width': '100%'},)
                                                        ],
                                                        className="col-divs",
                                                        width=12,
                                                        align='center',
                                                    ),
                                                ], width=3),
                                                dbc.Col([
                                                    # --- Output File Format ---
                                                    dbc.Col(
                                                        children=[
                                                            html.H6(
                                                                children=[
                                                                    "Output File Format:"],
                                                                style={'color': 'black'},
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
                                                        width=12,
                                                    ),
                                                    # --- Axis line width ---
                                                    dbc.Col(
                                                        children=[
                                                            html.H6(
                                                                children=[
                                                                    "Axis line width:"],
                                                                style={'color': 'black'},
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
                                                        width=12,
                                                    ),
                                                ], width=3),
                                                dbc.Col([
                                                    # --- Font Family ---
                                                    dbc.Col(
                                                        children=[
                                                            html.H6(
                                                                "Font Family:",
                                                                style={'color': 'black'},
                                                            ),
                                                            dcc.Dropdown(
                                                                id="font-family-option",
                                                                options=[{"label": i, "value": i} for i in FONT_FAMILIES],
                                                                value="Arial",
                                                                clearable=False,
                                                                className="dropdown-style",
                                                            ),
                                                        ],
                                                        className="col-divs",
                                                        width=12
                                                    ),
                                                    # --- Pixel output size ---
                                                    dbc.Col(
                                                        children=[
                                                            html.H6(
                                                                "Output Dimensions: (px)",
                                                                style={'color': 'black'},
                                                            ),
                                                            html.Div([
                                                                dcc.Input(id="tv-pixel-width", style={'width': '50%'}, type="number", placeholder="Width"),
                                                                html.Div(["x"], style={"margin": "0px 5px 0px 5px"}),
                                                                dcc.Input(id="tv-pixel-height", style={'width': '50%'}, type="number", placeholder="Height"),  
                                                            ], className="pixel-input-style",),
                                                        ],
                                                        className="col-divs",
                                                        width=12,
                                                    ),
                                                ], width=3),
                                                dbc.Col([
                                                    # --- Axis Grid Lines ---
                                                    dbc.Col(
                                                        children=[
                                                            html.H6(
                                                                children=["Grid Lines:"],
                                                                style={'color': 'black'},
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
                                                        width=12,
                                                    ),
                                                    # --- Edit Graphs ---
                                                    dbc.Col(
                                                        children=[
                                                            html.H6(
                                                                children=["Editable Graphs:"],
                                                                style={'color': 'black'},
                                                            ),
                                                            html.Div(
                                                                children=[
                                                                    html.Div(["Off"], style={'paddingLeft': '5px'}),
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
                                                        width=12,
                                                    ),
                                                ], width=2),
                                            ],
                                            className="row-div",
                                            no_gutters=True,
                                        ),
                                        dbc.Row([
                                            dbc.Col([
                                                dbc.Button("Apply Updates", id='update-options', color='primary', block=True),
                                            ], width=12),
                                        ]),
                                    ],
                                    className="narbar-collapse-plot-options",
                                ),
                                id="collapse-plot-options",
                            ),
                        ],
                        className="narbar-card"
                    ),
                    # Summary/Stats
                    dbc.Card(
                        children=[
                            dbc.Collapse(
                                dbc.CardBody(
                                    id='basic-stats-cardbody',
                                    className="narbar-collapse-toolbar",
                                    children=[
                                        dbc.Tabs(
                                            children=[
                                                # Current View
                                                dbc.Tab(
                                                    id='basic-stats-tab',
                                                    label='Current View',
                                                    label_style={
                                                        "color": "white",
                                                        "border": "2px orange solid",
                                                        "border-radius": "5px",
                                                        "margin-right": "2px",
                                                    },
                                                    children=[
                                                        dbc.Row(
                                                            children=[
                                                                dbc.Col([
                                                                    dcc.Graph(
                                                                        id='basic-stats-pie-chart',
                                                                        figure=tree_utils.init_data_graph("plotly_dark"),
                                                                        config=dict(displaylogo=False),
                                                                        style={"margin": "0px 5px 0px 5px"},
                                                                    ),
                                                                ], width=6, align='center'),
                                                                dbc.Col([
                                                                    dcc.Graph(
                                                                        id='normRF-fig',
                                                                        figure=tree_utils.init_data_graph("plotly_dark"),
                                                                        config=dict(displaylogo=False),
                                                                        style={"margin": "0px 5px 0px 5px"},
                                                                    ),
                                                                ],width=6, align='center'),
                                                            ],
                                                            no_gutters=True,
                                                        ),
                                                        dbc.Row(
                                                            children=[
                                                                dbc.Col([
                                                                    html.Div([
                                                                        dt.DataTable(
                                                                            id='basic-stats',
                                                                            export_format='xlsx',
                                                                            style_as_list_view=True,
                                                                            style_header={'backgroundColor': 'rgb(30, 30, 30)',},
                                                                            style_cell={
                                                                                'backgroundColor': '#555758',
                                                                                'color': 'white',
                                                                                'textAlign': 'left',
                                                                            },
                                                                        ),
                                                                    ], style={"margin": "10px 10px 0px 0px"}),
                                                                ], width=6, align='top'),
                                                                dbc.Col([
                                                                    html.Div([
                                                                        dt.DataTable(
                                                                            id='basic-alt-data-stats',
                                                                            export_format='xlsx',
                                                                            merge_duplicate_headers=True,
                                                                            page_action='none',
                                                                            page_size=10,
                                                                            style_as_list_view=True,
                                                                            style_header={'backgroundColor': 'rgb(30, 30, 30)'},
                                                                            style_cell={
                                                                                'backgroundColor': '#555758',
                                                                                'color': 'white',
                                                                                'textAlign': 'center',
                                                                            },
                                                                            style_table={'height': '50%', 'overflowY': 'auto'},
                                                                        ),
                                                                    ], style={"margin": "10px 10px 0px 0px"}),
                                                                ], width=6, align='top'),
                                                            ],
                                                            no_gutters=True,
                                                        ),
                                                    ],
                                                ),
                                                # Whole Genome
                                                dbc.Tab(
                                                    id='average-stats-tab',
                                                    label='Whole Genome',
                                                    label_style={
                                                        "color": "white",
                                                        "border": "2px orange solid",
                                                        "border-radius": "5px",
                                                        "margin-right": "2px",
                                                    },
                                                    children=[
                                                        dbc.Row(
                                                            children=[
                                                                dbc.Col([
                                                                    dcc.Graph(
                                                                        id="stats-topology-freq-pie-graph",
                                                                        figure=tree_utils.init_data_graph("plotly_dark"),
                                                                        style={"margin-right": "10px"},
                                                                    ),
                                                                ], sm=6, md=6, lg=5, align='center'),
                                                                dbc.Col([
                                                                    html.Div([
                                                                        dt.DataTable(
                                                                            id='basic-whole-genome-stats',
                                                                            export_format='xlsx',
                                                                            merge_duplicate_headers=True,
                                                                            page_size=10,
                                                                            style_as_list_view=True,
                                                                            style_header={'backgroundColor': 'rgb(30, 30, 30)'},
                                                                            style_cell={
                                                                                'backgroundColor': '#555758',
                                                                                'color': 'white',
                                                                                'textAlign': 'center',
                                                                            },
                                                                            style_table={'overflowY': 'auto', 'height': '50%'},
                                                                            style_filter={'color': 'white'},
                                                                            filter_action="native",     # allow filtering of data by user ('native') or not ('none')
                                                                            sort_action="native",       # enables data to be sorted per-column by user or not ('none')
                                                                            tooltip_header={
                                                                                'TopologyID': ['Statistics of additional data and individual topologies across the entire genome', ''],
                                                                            },
                                                                            css=[{
                                                                                'selector': '.dash-table-tooltip',
                                                                                'rule': 'background-color: #555758; font-family: monospace;'
                                                                            }],
                                                                        ),
                                                                    ], style={"margin-right": "10px"}),
                                                                ], sm=6, md=6, lg=7, align='center'),
                                                            ],
                                                            no_gutters=True,
                                                        )
                                                    ],
                                                ),
                                                # Topology-Quantile
                                                dbc.Tab(
                                                    id='topology-quantile-stats-tab',
                                                    label='Topology-Quantile',
                                                    label_style={
                                                        "color": "white",
                                                        "border": "2px orange solid",
                                                        "border-radius": "5px",
                                                        "margin-right": "2px",
                                                    },
                                                    children=[
                                                        dbc.Row(
                                                            children=[
                                                                dbc.Col([
                                                                    html.Div([
                                                                        topology_quantile_input_form,
                                                                    ], style={"background-color": "#555758", 'padding': '5px', 'margin-right': '10px', 'border': 'solid black 1px'}),
                                                                ], width=3, align='top'),
                                                                dbc.Col([
                                                                    html.Div([
                                                                        dcc.Graph(id='topo-quantile-graph-1', figure=tree_utils.init_data_graph("plotly_dark")),
                                                                        dt.DataTable(
                                                                            id='topo-quantile-dt-1',
                                                                            export_format='xlsx',
                                                                            merge_duplicate_headers=True,
                                                                            page_action='none',
                                                                            page_size=11,
                                                                            style_header={'backgroundColor': 'rgb(30, 30, 30)'},
                                                                            style_cell={
                                                                                'backgroundColor': '#555758',
                                                                                'color': 'white',
                                                                                'textAlign': 'left',
                                                                            },
                                                                            style_table={'height': '50%', 'overflowY': 'auto'},
                                                                        )],
                                                                        style={'width': '50%', 'padding': '5px'}
                                                                    ),
                                                                    html.Div([
                                                                        dcc.Graph(id='topo-quantile-graph-2', figure=tree_utils.init_data_graph("plotly_dark")),
                                                                        dt.DataTable(
                                                                            id='topo-quantile-dt-2',
                                                                            export_format='xlsx',
                                                                            merge_duplicate_headers=True,
                                                                            page_action='none',
                                                                            page_size=11,
                                                                            style_header={'backgroundColor': 'rgb(30, 30, 30)'},
                                                                            style_cell={
                                                                                'backgroundColor': '#555758',
                                                                                'color': 'white',
                                                                                'textAlign': 'left',
                                                                            },
                                                                            style_table={'height': '50%', 'overflowY': 'auto'},
                                                                        )],
                                                                        style={'width': '50%', 'padding': '5px'}
                                                                    ),
                                                                ], align='center', width=9, style={'display': 'inline-flex'})
                                                            ],
                                                            no_gutters=True,
                                                        ),
                                                    ],
                                                ),
                                                # Frequency Distributions
                                                dbc.Tab(
                                                    id='data-freq-dist-tab',
                                                    label='Frequency Distributions',
                                                    label_style={
                                                        "color": "white",
                                                        "border": "2px orange solid",
                                                        "border-radius": "5px",
                                                        "margin-right": "2px",
                                                    },
                                                    children=[
                                                        dbc.Row(
                                                            children=[
                                                                dbc.Col([
                                                                    dcc.Dropdown(
                                                                        id='data-freq-dist-dropdown',
                                                                        className="dropdown-style",
                                                                        multi=True,
                                                                    ),
                                                                ], width=10, style={'padding': '5px'}),
                                                                dbc.Col([
                                                                    dbc.Button('Run', id='data-freq-dist-button', color='info', block=True),
                                                                ], width=1, style={'padding': '5px'}),
                                                                dbc.Col([
                                                                    dcc.Loading(id='freq-dist-loading'),
                                                                ], width=1, style={'padding': '5px'}, align='center'),
                                                            ],
                                                            no_gutters=True,
                                                            style={'borderBottom': '2px solid grey'}
                                                        ),
                                                        dbc.Row(
                                                            id='data-freq-dist-row',
                                                            no_gutters=True,
                                                            style={'padding': '5px'}
                                                        ),
                                                    ],
                                                ),
                                                # Statistical Testing
                                                dbc.Tab(
                                                    id='stats-tests-tab',
                                                    label='Exploratory Statistical Testing',
                                                    label_style={
                                                        "color": "white",
                                                        "border": "2px orange solid",
                                                        "border-radius": "5px",
                                                        "margin-right": "2px",
                                                    },
                                                    children=[
                                                        dbc.Row(
                                                            children=[
                                                                dbc.Col([
                                                                    html.Div([
                                                                        correlation_input_form,
                                                                    ], style={"background-color": "#555758", 'padding': '5px', 'margin-right': '10px', 'border': 'solid black 1px'}),
                                                                ], width=3, align='top'),
                                                                dbc.Col([
                                                                    html.Div([
                                                                        dt.DataTable(
                                                                            id='correlation-test-output',
                                                                            export_format='xlsx',
                                                                            merge_duplicate_headers=True,
                                                                            page_action='none',
                                                                            page_size=11,
                                                                            style_header={'backgroundColor': 'rgb(30, 30, 30)'},
                                                                            style_cell={
                                                                                'backgroundColor': '#555758',
                                                                                'color': 'white',
                                                                                'textAlign': 'left',
                                                                            },
                                                                            style_table={'height': '50%', 'overflowY': 'auto'},
                                                                        ),
                                                                    ], style={"margin-right": "10px"}),
                                                                    html.Br(),
                                                                    html.Div([
                                                                        dt.DataTable(
                                                                            id='correlation-stats',
                                                                            export_format='xlsx',
                                                                            merge_duplicate_headers=True,
                                                                            page_action='none',
                                                                            page_size=11,
                                                                            style_header={'backgroundColor': 'rgb(30, 30, 30)'},
                                                                            style_cell={
                                                                                'backgroundColor': '#555758',
                                                                                'color': 'white',
                                                                                'textAlign': 'left',
                                                                            },
                                                                            style_table={'height': '50%', 'overflowY': 'auto'},
                                                                        ),
                                                                    ], style={"margin-right": "10px"}),
                                                                    html.Br(),
                                                                    html.Div([
                                                                        dt.DataTable(
                                                                            id='mean-freq-stats',
                                                                            export_format='xlsx',
                                                                            merge_duplicate_headers=True,
                                                                            page_action='none',
                                                                            page_size=11,
                                                                            style_header={'backgroundColor': 'rgb(30, 30, 30)'},
                                                                            style_cell={
                                                                                'backgroundColor': '#555758',
                                                                                'color': 'white',
                                                                                'textAlign': 'left',
                                                                            },
                                                                            style_table={'height': '50%', 'overflowY': 'auto'},
                                                                        ),
                                                                    ], style={"margin-right": "10px"}),
                                                                ], width=5, align='top'),
                                                                dbc.Col([
                                                                    dcc.Graph(
                                                                        id="correlation-stats-heatmap",
                                                                        figure=tree_utils.init_stats_graph("plotly_dark"),
                                                                    ),
                                                                ], width=4, align='top'),
                                                            ],
                                                            no_gutters=True,
                                                        ),
                                                    ],
                                                ),
                                            ],
                                            className="tabs-style",
                                        ),
                                    ],
                                ),
                                id="collapse-statsbar",
                            ),
                        ],
                        className="narbar-card"
                    ),
                    # Toolbar
                    dbc.Card(
                        children=[
                            dbc.Collapse(
                                dbc.CardBody(
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
                                                                                chrom_genome_toggle,
                                                                            ],
                                                                            style={'color': 'black'},
                                                                        ),
                                                                        dcc.Dropdown(
                                                                            id="chromosome-options",
                                                                            className="dropdown-style",
                                                                            optionHeight=20,
                                                                            clearable=False,
                                                                        ),
                                                                    ],
                                                                    className="home-tab-button-cols",
                                                                ),
                                                                # Topologies
                                                                html.Div(
                                                                    children=[
                                                                        html.H6(
                                                                            children=[topology_order_toggle],
                                                                            style={'color': 'black'},
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
                                                                # Additional Data
                                                                html.Div(
                                                                    children=[
                                                                        html.H6(
                                                                            [additional_data_toggle],
                                                                            style={'color': 'black'},
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
                                                            ], width=6
                                                        ),
                                                        dbc.Col(
                                                            children=[
                                                                
                                                                # Trees
                                                                html.Div(
                                                                    children=[
                                                                        html.H6([
                                                                            tree_graph_toggle,
                                                                        ], className="text-divs"),
                                                                        html.Div(tree_styles, className="radio-group"),
                                                                        dcc.Dropdown(
                                                                            id="tree-taxa-choices",
                                                                            className="dropdown-style",
                                                                            multi=True,
                                                                        ),
                                                                    ],
                                                                    className="home-tab-button-cols",
                                                                ),
                                                                # Quantile Graph
                                                                html.Div(
                                                                    children=[
                                                                        html.H6(
                                                                            children=[quantile_graph_toggle],
                                                                            style={'color': 'black'},
                                                                        ),
                                                                        html.Div([
                                                                            dcc.Input(
                                                                                id='n-quantiles',
                                                                                disabled=True,
                                                                                placeholder="n-quantiles (default: 25)",
                                                                                type='number',
                                                                                min=2,
                                                                                step=1,
                                                                                style={'margin-top': '1px', 'margin-right': '2px'}
                                                                            ),
                                                                        ], style={'display': 'inline-flex', 'width': '100%'}),
                                                                    ], 
                                                                    className="home-tab-button-cols"
                                                                ),
                                                                # Standarized File
                                                                html.Div(
                                                                    children=[
                                                                        html.H6(
                                                                            children=["Standardized Files"],
                                                                            style={'padding-bottom': '2px', 'color': 'black'}
                                                                        ),
                                                                        dbc.Checklist(
                                                                            id="alt-graph-switches",
                                                                            options=[
                                                                                {"label": "GFF/GTF - No file chosen...", "value": "gff_switch", "disabled": True},
                                                                            ],
                                                                            className="checklist-style",
                                                                            inline=True,
                                                                        ),
                                                                    ], 
                                                                    className="home-tab-button-cols"
                                                                ),
                                                            ], width=6
                                                        ),
                                                    ],
                                                    className="home-tab-divs",
                                                    width=7,
                                                ),
                                                # Graph Toggle Switches
                                                dbc.Col(
                                                    children=[
                                                        dbc.Col([graph_switch_col1]),
                                                    ], 
                                                    style={'display': 'inline-flex'},
                                                    width=5,
                                                ),
                                            ],
                                            className="nav-row-div",
                                            no_gutters=True,
                                        ),
                                        dbc.Row([
                                            dbc.Col([
                                                dbc.Button("Apply Updates", id='update-options', color='primary', block=True),
                                            ], width=12),
                                        ]),
                                    ],
                                    className="narbar-collapse-toolbar",
                                ),
                                id="collapse-toolbar",
                            ),
                        ],
                        className="narbar-card"
                    ),
                ],
                className="narbar-collapse-toolbar",
                no_gutters=True,
            ),
            # Graphs
            dbc.Row(
                children=[
                    # --------- Main Distriubtion ---------
                    dbc.Row(
                        id='main-graph-div',
                        className="visible-div",
                        no_gutters=True,
                        children=[
                            # --- Main Topology Graph ---
                            dbc.Col(
                                children=[
                                    dcc.Graph(
                                        id="topologyGraph",
                                        figure=tree_utils.init_data_graph("plotly_dark"),
                                        className="topograph-style",
                                        config=dict(displaylogo=False, doubleClick="reset", displayModeBar=False),
                                    ),
                                    dcc.Loading(
                                        id="topology_graph_div",
                                        className="hidden-div",
                                        type="dot",
                                        children=[
                                            html.Div(
                                                id="mainDistPlot",
                                                children=[
                                                    html.Div(id="preview-graphs-div"),
                                                    html.Div(id="additional-graphs-div"),
                                                ],
                                            ),
                                        ],
                                    ),
                                    html.Div(id="tree-graphs-div"),
                                ],
                                width=12,
                            ),
                        ],
                    ),
                    # -------- Whole genome graphs --------
                    dbc.Row(
                        id="whole-genome-div",
                        children=[
                            # Per-chromosome Frequency Plots
                            dcc.Loading([
                                dbc.Row(
                                    children=[html.Div(id="whole-genome-graphs", style={'height': '100%', 'min-height': '200px', 'width': '100%'})],
                                    className="stats-chrom-frequency",
                                    no_gutters=True,
                                ),
                                dbc.Row(
                                    children=[html.Div(id="whole-genome-alt-graphs", style={'height': '100%', 'min-height': '0px', 'width': '100%'})],
                                    # className="stats-chrom-frequency",
                                    no_gutters=True,
                                ),
                            ]),
                        ],
                        className="hidden-div",
                        no_gutters=True,
                    ),
                    # -------------- Quantile graph ----------------
                    dbc.Row(
                        no_gutters=True,
                        style={'width': '100%'},
                        children=[
                            dbc.Col(
                                html.Div(id="quantile-graphs-div"),
                                width=12,
                            ),
                        ],
                    ),
                    # -------------- Trees ----------------
                    dbc.Row(
                        no_gutters=True,
                        style={'width': '100%'},
                        children=[
                            # --- Main Topology Graph ---
                            dbc.Col(
                                children=[
                                    html.Div(id="tree-graphs-div"),
                                ],
                                width=12,
                            ),
                        ],
                    ),
                ],
                className="div-boxes",
                no_gutters=True,
            ),
            html.Br(id="footer-spacing", style={'height': '5px'})
        ],
        id="tree-viewer-container",
        fluid=True,
    )
    return layout

##################################### Callbacks #####################################
@app.callback(
    Output("tv-docs-content", "children"),
    [Input("collapse-documentation", "is_open"),],)
def render_docs_content(
    is_open,
):
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
     Output("collapse-statsbar", "is_open"),
     Output("collapse-documentation", "is_open"),
     Output("collapse-plot-options", "is_open"),
    ],
    [Input("toggle-toolbar", "n_clicks"),
     Input("toggle-statsbar", "n_clicks"),
     Input("toggle-docs", "n_clicks"),
     Input("toggle-plot-options", "n_clicks"),
    ],
    [State("collapse-toolbar", "is_open"),
     State("collapse-statsbar", "is_open"),
     State("collapse-documentation", "is_open"),
     State("collapse-plot-options", "is_open"),
    ],
)
def toggle_toolbar(
    n1, n2, 
    n3, n4, 
    toolbar_is_open, 
    statsbar_is_open, 
    docs_is_open, 
    plot_options_is_open,
):
    """This function handles the opening/closing of the collapsible menus"""
    ctx = dash.callback_context

    if not ctx.triggered:
        return True, False, False, False
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if button_id == "toggle-toolbar" and n1:
        return not toolbar_is_open, False, False, False
    elif button_id == "toggle-statsbar" and n2:
        return False, not statsbar_is_open, False, False
    elif button_id == "toggle-docs" and n3:
        return False, False, not docs_is_open, False
    elif button_id == "toggle-plot-options" and n4:
        return False, False, False, not plot_options_is_open

    else:
        return False, False, False, False

# ---------------------- Hidden div handler ---------------------
@app.callback(
    [Output("topology_graph_div", "className"),
     Output("main-graph-div", "className"),
     Output("whole-genome-div", "className"),
    ],
    [Input("toggle-chrom-whole-genome", "value"),
     Input("input-data-upload", "children"),
     Input("chromFile-data-upload", "children"),
    ],
)
def hide_graph_divs(
    view_toggle,
    tvdata,
    chromdata,
):
    """This callback controls the graphs display and styling"""
    if not tvdata:
        raise PreventUpdate
    elif not chromdata:
        raise PreventUpdate
    elif view_toggle:
        # Topology distribution
        topology_graph_style = "dist-graph-div"
        distribution_div = "visible-div"
        # Hidden divs
        whole_genome_style = "hidden-div"
        return (topology_graph_style, distribution_div, whole_genome_style)
    ## ---- Whole Genome Graphs ----
    elif not view_toggle:
        # Hidden divs
        topology_graph_style = "hidden-div"
        distribution_div = "hidden-div"
        # On divs
        whole_genome_style = "stats-main-div"
        return (topology_graph_style, distribution_div, whole_genome_style)
    else:
        raise PreventUpdate

# ---------------------- Update Graph Toggle Options + Dropdowns ---------------------
@app.callback(
    [Output("main-graph-switches", "options"), 
     Output("alt-graph-switches", "options"),
     Output("tree-graph-switch", "disabled"),
     Output("alt-data-switch", "disabled"),
     Output("quantile-toggle", "disabled"),
     Output("n-quantiles", "disabled"),
     Output("whole-genome-graph-type", "options"),
     Output("wg-squish-expand", "options"),
    ],
    [Input("gff-data-upload", "children"),
     Input("input-data-upload", "children"),
     Input("alt_data_options", "options"),
    ]
)
def main_graph_button_options(
    gffData,
    topoData,
    additional_graphs,
):
    if not gffData:
        gffOption = {"label": "GFF/GTF - No file chosen...", "value": "gff_switch", "disabled": True}
    else:
        gffOption = {"label": "GFF/GTF", "value": "gff_switch"}
    if not topoData:
        altDataOption = True
        treeOption = True
        quantileButton = True
        quantileToggle = True
        quantileInput = True
        combinedPlotOption = {"label": "Rug+Tile", "value": "topo_both", "disabled": True}
        rugPlotOption = {"label": "Rug", "value": "topo_rug_only", "disabled": True}
        tilePlotOption = {"label": "Tile", "value": "topo_tile_only", "disabled": True}
        gffOption = {"label": "GFF/GTF - No file chosen...", "value": "gff_switch", "disabled": True}
        wholeGenomeGraphsoptions=[
            {'label': 'Rug', 'value': 'rug', 'disabled': True},
            {'label': 'Bar', 'value': 'bar', 'disabled': True},
            {'label': 'Pie', 'value': 'pie', 'disabled': True},
            {'label': 'Tile', 'value': 'tile', 'disabled': True},
        ]
        wholeGenomeGraphSizeoptions=[
            {'label': 'Collapse', 'value': 'collapse', 'disabled': True},
            {'label': 'Squish', 'value': 'squish', 'disabled': True},
            {'label': 'Expand', 'value': 'expand', 'disabled': True},
        ]
    else:
        altDataOption = False
        treeOption = False
        quantileButton = False
        quantileToggle = False
        quantileInput = False
        combinedPlotOption = {"label": "Rug+Tile", "value": "topo_both"}
        rugPlotOption = {"label": "Rug", "value": "topo_rug_only"}
        tilePlotOption = {"label": "Tile", "value": "topo_tile_only"}
        wholeGenomeGraphsoptions=[
            {'label': 'Rug', 'value': 'rug'},
            {'label': 'Bar', 'value': 'bar'},
            {'label': 'Pie', 'value': 'pie'},
            {'label': 'Tile', 'value': 'tile'},
        ]
        wholeGenomeGraphSizeoptions=[
            {'label': 'Collapse', 'value': 'collapse'},
            {'label': 'Squish', 'value': 'squish'},
            {'label': 'Expand', 'value': 'expand'},
        ]
    # If no additional data loaded, disable additional data switch
    if not additional_graphs:
        altDataOption = True
    elif len(additional_graphs) > 0:
        altDataOption = False
    
    main_graph_options = [
        combinedPlotOption,
        rugPlotOption,
        tilePlotOption,
    ]
    additional_graph_options = [
        gffOption,
    ]
    return (main_graph_options, additional_graph_options,
            treeOption, altDataOption, quantileToggle,
            quantileInput, wholeGenomeGraphsoptions, wholeGenomeGraphSizeoptions)


@app.callback(
    [Output("alt_data_options", 'disabled'),
     Output("tree-taxa-choices", 'disabled'),
     Output("chromosome-options", "disabled"),
     Output("topology_options", "disabled"),],
    [Input("input-data-upload", "children"),
     Input("alt-data-switch", "value"),
     Input("tree-graph-switch", "value"),
     Input("toggle-chrom-whole-genome", "value"),
    ],
)
def disable_dropdowns(
    tv_input_json,
    alt_switch,
    tree_switch,
    view_toggle,
):
    if not tv_input_json:
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
    if not view_toggle:
        chromDropdown = True
    return altDropdown, treeDropdown, chromDropdown, topoDropdown

# ---------------------- Data Upload Callbacks ---------------------
@app.callback(
    [Output("alt-input-modal", "is_open"),
     Output("gff-data-upload", "children"),
     Output("gff-filename-div", "children"),
     Output("gff-new-session-button", "color"),
     Output("gff-error-modal", "is_open"),
    ],
    [Input("gff-init-upload-button", "n_clicks"),
     Input("gff-submit-button", "n_clicks"),
     Input("input-files-submit-button", "n_clicks"),
     Input("gff-upload-data", "contents"),
     Input("gff-data-handoff", "children"),
     Input("gff-filename-div", "children"),
    ],
    [State("gff-upload-data", "filename"),
     State("gff-data-upload", "children"),
    ],
)
def upload_gff_file(
    gffBtn,
    gffSubmit,
    tvBtn,
    gffContents,
    gffHandoff,
    gffFilenamediv,
    gffFilename,
    gffData,
):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = "No clicks yet"
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if button_id == "gff-init-upload-button":
        if gffFilenamediv:
            gffDF = gffData
            gffFilename = gffFilenamediv
            buttonColor = 'success'
        elif not gffFilename:
            gffDF = None
            gffFilename = "No file chosen..."
            buttonColor = "primary"
        return True, gffDF, gffFilename, buttonColor, False
    elif button_id == "gff-new-session-button":
        return True, None, None, "primary", False
    elif button_id == "gff-submit-button":
        if not gffContents:
            return True, None, gffFilename, "primary", False
        _, content_string = gffContents.split(",")
        tvDecoded = base64.b64decode(content_string)
        if ("gff" in gffFilename) or ("gtf" in gffFilename):
            try:
                assert tree_utils.valid_gff_gtf(io.StringIO(tvDecoded.decode("utf-8")))
            except AssertionError:
                return False, None, gffFilename, "danger", True
            gffDF = pd.read_csv(
                io.StringIO(tvDecoded.decode("utf-8")),
                sep="\t",
                names=["chromosome", "source", "feature", "start",
                       "end", "score", "strand", "frame", "attribute"],
                comment="#",
            )
            gffDF["attribute"] = gffDF["attribute"].apply(lambda x: data_utils.get_gene_name(x))
            return False, [gffDF.to_json()], gffFilename, "success", False
        else:
            return True, None, gffFilename, "Please provide a valid GFF/GTF file", False
    elif button_id == "gff-upload-data":
        if tree_utils.validate_gff_gtf_filename(gffFilename):
            return True, None, gffFilename, "success", False
        else:
            return True, None, gffFilename, "danger", False
    elif gffHandoff:
        # Currently this feature only takes 1 GFF/GTF file as input.
        # Future work should focus on allowing multiple files
        gffDFs = list()
        for i in gffHandoff:
            if ("gff" in i) or ("gtf" in i):
                try:
                    assert tree_utils.valid_gff_gtf(i)
                except AssertionError:
                    return False, None, gffFilename, "danger", True
                file = Path(i).name
                gffDF = pd.read_csv(
                    i,
                    sep="\t",
                    names=["chromosome", "source", "feature", "start",
                        "end", "score", "strand", "frame", "attribute"],
                    comment="#",
                )
                gffDF["attribute"] = gffDF["attribute"].apply(lambda x: data_utils.get_gene_name(x))
                for i in gffDF["start"]:
                    try:
                        int(i)
                    except:
                        return False, None, gffFilename, None, True
                gffDFs.append(gffDF.to_json())
            else:
                file = None
        return False, gffDFs, file, "success", False
    elif gffFilename:
        if tree_utils.validate_gff_gtf_filename(gffFilename):
            return True, None, gffFilename, "success", False
        else:
            return True, None, gffFilename, "danger", False
    else:
        raise PreventUpdate


# Upload Tree Viewer and Chromosome length files
@app.callback(
    [Output("input-data-upload", "children"),
     Output("chromFile-data-upload", "children"),
     Output("main-input-modal", "is_open"),
     Output("tv-filename-div", "children"),
     Output("chrom-filename-div", "children"),
     Output("project-name-div", "children"),
     Output("window-size", "children"),
     Output("tv-button", "color"),
     Output("chromFile-button", "color"),
     Output("new-session-project-button", "color"),
     Output("tv-modal-warning-label", "children"),
     Output("tree-taxa-choices", "options"),
     Output("tree-taxa-choices", "value"),
     Output("gff-data-handoff", "children"),
     Output("project-dir", "data"),
    ],
    [Input("new-session-button", "n_clicks"),
     Input("input-files-submit-button", "n_clicks"),
     Input("input-files-close-button", "n_clicks"),
     Input("upload-data", "contents"),
     Input("chromFile-upload", "contents"), 
     Input("new-session-project", "contents"),
    ],
    [State("input-data-upload", "children"),
     State("chromFile-data-upload", "children"),
     State("chromFile-upload", "filename"),
     State("upload-data", "filename"),
     State("new-session-project", "filename"),
     State("tv-filename-div", "children"),
     State("chrom-filename-div", "children"),
     State("input-modal-tabs", "value"),
     State("project-dir", "data"),
    ],
)
def load_input_files(
    input_btn,
    submit_btn,
    close_btn,
    tvFile_contents,
    chromFile_contents,
    project_contents,
    tv_data,
    chrom_data,
    chromFilename,
    tvFilename,
    projectName,
    current_tvFileName,
    current_chromFileName,
    active_input_tab,
    project_dir,
):
    ctx = dash.callback_context
    button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    # Init output options
    tvDF = dash.no_update
    chromDF = dash.no_update
    modalOpen = True
    windowSize = dash.no_update
    tvButtonColor = dash.no_update
    chromButtonColor = dash.no_update
    projectButtonColor = dash.no_update
    tvWarningLabel = dash.no_update
    treeTaxaOptions = dash.no_update
    treeTaxaValues = dash.no_update
    gffData = dash.no_update
    # Input file extensions
    try:
        tvFileType = tvFilename.split(".")[-1]
    except AttributeError:
        pass
    try:
        chromFileType = chromFilename.split(".")[-1]
    except AttributeError:
        pass
    # TreeViewer input button color
    if not tvFilename:
        tvButtonColor = "primary"
    elif 'txt' in tvFileType:
        tvButtonColor = "success"
    elif 'csv' in tvFileType:
        tvButtonColor = "success"
    elif 'tsv' in tvFileType:
        tvButtonColor = "success"
    elif 'xls' in tvFileType:
        tvButtonColor = "success"
    else:
        tvButtonColor = "danger"
    # Chrom length button color
    if not chromFilename:
        chromButtonColor = "primary"
    elif 'bed' in chromFileType:
        chromButtonColor = "success"
    else:
        chromButtonColor = "danger"
    # Project button color
    if not projectName:
        projectButtonColor = "primary"
    elif ".ini" in projectName:
        projectButtonColor = "success"
    else:
        projectButtonColor = "danger"
    # Repond depending on button clicked
    # - Open new session modal
    if button_id == "new-session-button":
        if current_tvFileName:
            tvFilename = current_tvFileName
        elif not tvFilename:
            tvFilename = "No file chosen..."
        
        if current_chromFileName:
            chromFilename = current_chromFileName
        elif not chromFilename:
            chromFilename = "No file chosen..."
        
        if not projectName:
            projectName = "No project.ini file chosen..."        
        return [
            tvDF,
            chromDF,
            modalOpen,
            tvFilename,
            chromFilename,
            projectName,
            windowSize,
            tvButtonColor,
            chromButtonColor,
            projectButtonColor,
            tvWarningLabel, 
            treeTaxaOptions,
            treeTaxaValues,
            gffData,
            project_dir,
        ]
    # - Close modal and return current values
    elif button_id == "input-files-close-button":
        modalOpen = False
        return [
            tvDF,
            chromDF,
            modalOpen,
            tvFilename,
            chromFilename,
            projectName,
            windowSize,
            tvButtonColor,
            chromButtonColor,
            projectButtonColor,
            tvWarningLabel, 
            treeTaxaOptions,
            treeTaxaValues,
            gffData,
            project_dir,
        ]
    # - Submit current data to Div
    elif button_id == "input-files-submit-button":
        if active_input_tab == "project-input-tab":
            # Set button colors
            tvButtonColor = 'primary'
            chromButtonColor = 'primary'
            # Load TreeViewer file into DF
            if not project_contents:
                modalOpen=True
                projectName="No file chosen..."
                projectButtonColor='danger'
                return [
                    None,
                    None,
                    modalOpen,
                    tvFilename,
                    chromFilename,
                    projectName,
                    windowSize,
                    tvButtonColor,
                    chromButtonColor,
                    projectButtonColor,
                    tvWarningLabel,
                    treeTaxaOptions,
                    treeTaxaValues,
                    gffData,
                    None,
                ]
            _, projectContentString = project_contents.split(",")
            projectDecoded = base64.b64decode(projectContentString)
            import configparser
            projectConfig = configparser.ConfigParser()
            try:
                projectConfig.read_string(projectDecoded.decode("utf-8"))
            except MissingSectionHeaderError:
                tvButtonColor = 'primary'
                chromButtonColor = 'primary'
                projectButtonColor = "danger"
                tvWarningLabel = "Invalid project file - validate and run again"
                return [
                    tvDF,
                    chromDF,
                    modalOpen,
                    tvFilename,
                    chromFilename,
                    projectName,
                    windowSize,
                    tvButtonColor,
                    chromButtonColor,
                    projectButtonColor,
                    tvWarningLabel, 
                    treeTaxaOptions,
                    treeTaxaValues,
                    gffData,
                    None,
                ]
            except UnicodeDecodeError:
                tvButtonColor = 'primary'
                chromButtonColor = 'primary'
                projectButtonColor = "danger"
                tvWarningLabel = "Invalid project file - validate and run again"
                return [
                    tvDF,
                    chromDF,
                    modalOpen,
                    tvFilename,
                    chromFilename,
                    projectName,
                    windowSize,
                    tvButtonColor,
                    chromButtonColor,
                    projectButtonColor,
                    tvWarningLabel, 
                    treeTaxaOptions,
                    treeTaxaValues,
                    gffData,
                    None,
                ]
            try:
                tvFilename = projectConfig['MAIN']['TreeViewerFile']
                chromFilename = projectConfig['MAIN']['ChromLengths']
                projectDir = projectConfig['MAIN']['ProjectDir']  # Unused for the moment, but will need later on as Projects develops
            except KeyError:
                tvButtonColor = 'primary'
                chromButtonColor = 'primary'
                projectButtonColor = "danger"
                tvWarningLabel = "Input header 'MAIN' is invalid - validate and rerun"
                return [
                    tvDF,
                    chromDF,
                    modalOpen,
                    tvFilename,
                    chromFilename,
                    projectName,
                    windowSize,
                    tvButtonColor,
                    chromButtonColor,
                    projectButtonColor,
                    tvWarningLabel, 
                    treeTaxaOptions,
                    treeTaxaValues,
                    gffData,
                    None,
                ]
            try:
                gffFiles = projectConfig['ADDITIONAL']['GFF_GTF'].split(";")
            except KeyError:
                tvButtonColor = 'primary'
                chromButtonColor = 'primary'
                projectButtonColor = "danger"
                tvWarningLabel = "Input header 'ADDITIONAL' is invalid - validate and rerun"
                return [
                    tvDF,
                    chromDF,
                    modalOpen,
                    tvFilename,
                    chromFilename,
                    projectName,
                    windowSize,
                    tvButtonColor,
                    chromButtonColor,
                    projectButtonColor,
                    tvWarningLabel, 
                    treeTaxaOptions,
                    treeTaxaValues,
                    gffData,
                    None,
                ]
            # If no gff/gtf file provided, set to None object
            if gffFiles == ["None"]:
                gffFiles = None
            else:
                gffData = gffFiles
            # Check if files exist
            try:
                assert Path(tvFilename).exists()
            except AssertionError:
                tvButtonColor = 'danger'
                chromButtonColor = 'warning'
                tvWarningLabel = "Tree Viewer input file does not exist - check input and re-submit"
                return [
                    tvDF,
                    chromDF,
                    modalOpen,
                    tvFilename,
                    chromFilename,
                    projectName,
                    windowSize,
                    tvButtonColor,
                    chromButtonColor,
                    projectButtonColor,
                    tvWarningLabel, 
                    treeTaxaOptions,
                    treeTaxaValues,
                    None,
                    None,
                ]
            try:
                assert Path(chromFilename).exists()
            except AssertionError:
                tvButtonColor = 'warning'
                chromButtonColor = 'danger'
                tvWarningLabel = "Chromosome length BED file does not exist - check input and re-submit"
                return [
                    tvDF,
                    chromDF,
                    modalOpen,
                    tvFilename,
                    chromFilename,
                    projectName,
                    windowSize,
                    tvButtonColor,
                    chromButtonColor,
                    projectButtonColor,
                    tvWarningLabel, 
                    treeTaxaOptions,
                    treeTaxaValues,
                    None,
                    None,
                ]
            if gffFiles:
                for i in gffFiles:
                    try:
                        assert Path(i).exists()
                    except AssertionError:
                        tvButtonColor = 'warning'
                        chromButtonColor = 'warning'
                        tvWarningLabel = f"{i} does not exist - check input and re-submit"
                        return [
                            tvDF,
                            chromDF,
                            modalOpen,
                            tvFilename,
                            chromFilename,
                            projectName,
                            windowSize,
                            tvButtonColor,
                            chromButtonColor,
                            projectButtonColor,
                            tvWarningLabel, 
                            treeTaxaOptions,
                            treeTaxaValues,
                            None,
                            None,
                        ]
            # Load file by suffix
            if "csv" in tvFilename:
                # Assume that the user uploaded a CSV file
                tvDF = pd.read_csv(tvFilename, comments="#")
                tvDF.sort_values(by=["Chromosome", "Window"], inplace=True)
                validated_tvDF = tree_utils.tv_header_validation(tvDF)
                # If validated_DF returns False, raise error 
                # indicating issue with input file 
                if not validated_tvDF:
                    tvDF = dash.no_update
                    tvFilename = f"{tvFilename} is malformed!"
                    tvButtonColor = 'danger'
                    chromButtonColor = 'warning'
                    tvWarningLabel = "WARNING: Tree Viewer file appears to be malformed. First four headers must be ['Chromosome', 'Window', 'NewickTree', 'TopologyID']"
                    return [
                        tvDF,
                        chromDF,
                        modalOpen,
                        tvFilename,
                        chromFilename,
                        projectName,
                        windowSize,
                        tvButtonColor,
                        chromButtonColor,
                        projectButtonColor,
                        tvWarningLabel, 
                        treeTaxaOptions,
                        treeTaxaValues,
                        gffData,
                        None,
                    ]
                if len(tvDF.columns) == 4:
                    tvDF["None"] = [pd.NA]*len(tvDF)
                tvDF["TopologyID"] = tvDF["TopologyID"].apply(
                    lambda x: "NODATA" if type(x) != str else x)
                pass
            elif "txt" in tvFilename:
                # Assume that the user uploaded a CSV file
                tvDF = pd.read_csv(io.StringIO(tvFilename, sep='\t'), comments="#")
                tvDF.sort_values(by=["Chromosome", "Window"], inplace=True)
                validated_tvDF = tree_utils.tv_header_validation(tvDF)
                # If validated_DF returns False, raise error 
                # indicating issue with input file 
                if not validated_tvDF:
                    tvDF = dash.no_update
                    tvFilename = f"{tvFilename} is malformed!"
                    tvButtonColor = 'danger'
                    chromButtonColor = 'warning'
                    tvWarningLabel = "WARNING: Tree Viewer file appears to be malformed. First four headers must be ['Chromosome', 'Window', 'NewickTree', 'TopologyID']"
                    return [
                        tvDF,
                        chromDF,
                        modalOpen,
                        tvFilename,
                        chromFilename,
                        projectName,
                        windowSize,
                        tvButtonColor,
                        chromButtonColor,
                        projectButtonColor,
                        tvWarningLabel, 
                        treeTaxaOptions,
                        treeTaxaValues,
                        gffData,
                        None,
                    ]
                if len(tvDF.columns) == 4:
                    tvDF["None"] = [pd.NA]*len(tvDF)
                tvDF["TopologyID"] = tvDF["TopologyID"].apply(
                    lambda x: "NODATA" if type(x) != str else x)
                pass
            elif "tsv" in tvFilename:
                # Assume that the user uploaded a CSV file
                tvDF = pd.read_csv(tvFilename, sep='\t', comments="#")
                tvDF.sort_values(by=["Chromosome", "Window"], inplace=True)
                validated_tvDF = tree_utils.tv_header_validation(tvDF)
                # If validated_DF returns False, raise error 
                # indicating issue with input file 
                if not validated_tvDF:
                    tvDF = dash.no_update
                    tvFilename = f"{tvFilename} is malformed!"
                    tvButtonColor = 'danger'
                    chromButtonColor = 'warning'
                    tvWarningLabel = "WARNING: Tree Viewer file appears to be malformed. First four headers must be ['Chromosome', 'Window', 'NewickTree', 'TopologyID']"
                    return [
                        tvDF,
                        chromDF,
                        modalOpen,
                        tvFilename,
                        chromFilename,
                        projectName,
                        windowSize,
                        tvButtonColor,
                        chromButtonColor,
                        projectButtonColor,
                        tvWarningLabel, 
                        treeTaxaOptions,
                        treeTaxaValues,
                        gffData,
                        None,
                    ]
                if len(tvDF.columns) == 4:
                    tvDF["None"] = [pd.NA]*len(tvDF)
                tvDF["TopologyID"] = tvDF["TopologyID"].apply(
                    lambda x: "NODATA" if type(x) != str else x)
                pass
            elif "xls" in tvFilename:
                # Assume that the user uploaded an excel file
                tvDF = pd.read_excel(tvFilename, engine="openpyxl", comment="#")
                tvDF.sort_values(by=["Chromosome", "Window"], inplace=True)
                # validated_tvDF = tree_utils.validate_tree_viewer_input(tvDF)
                validated_tvDF = tree_utils.tv_header_validation(tvDF)
                # If validated_DF returns False, raise error 
                # indicating issue with input file 
                if not validated_tvDF:
                    tvDF = dash.no_update
                    tvFilename = f"{tvFilename} is malformed!"
                    tvButtonColor = 'danger'
                    chromButtonColor = 'warning'
                    tvWarningLabel = "WARNING: Tree Viewer file appears to be malformed. First four headers must be ['Chromosome', 'Window', 'NewickTree', 'TopologyID']"
                    return [
                        tvDF,
                        chromDF,
                        modalOpen,
                        tvFilename,
                        chromFilename,
                        projectName,
                        windowSize,
                        tvButtonColor,
                        chromButtonColor,
                        projectButtonColor,
                        tvWarningLabel, 
                        treeTaxaOptions,
                        treeTaxaValues,
                        gffData,
                        None,
                    ]
                if len(tvDF.columns) == 4:
                    tvDF["None"] = [pd.NA]*len(tvDF)
                tvDF["TopologyID"] = tvDF["TopologyID"].apply(lambda x: "NODATA" if type(x) != str else x)
                tvDF[["TopologyID"]] = tvDF[["TopologyID"]].fillna(value="NoData")
                pass
            else:
                tvDF = dash.no_update
                tvFilename = f"Invalid file type - {tvFilename}"
                tvButtonColor = 'warning'
                chromButtonColor = 'primary'
                tvWarningLabel = "WARNING: Please provide a valid Tree Viewer file"
                return [
                    tvDF,
                    chromDF,
                    modalOpen,
                    tvFilename,
                    chromFilename,
                    projectName,
                    windowSize,
                    tvButtonColor,
                    chromButtonColor,
                    projectButtonColor,
                    tvWarningLabel, 
                    treeTaxaOptions,
                    treeTaxaValues,
                    gffData,
                    None,
                ]
            if "bed" in chromFilename:
                # Assume that the user uploaded a CSV file
                chromDF = pd.read_csv(
                    chromFilename,
                    sep="\t",
                    names=["Chromosome", "Start", "End"],
                )
                pass
            else:
                chromFilename = f"Invalid file type - {chromFilename}"
                tvButtonColor = 'primary'
                chromButtonColor = 'warning'
                return [
                    tvDF.to_json(),
                    chromDF,
                    modalOpen,
                    tvFilename,
                    chromFilename,
                    projectName,
                    windowSize,
                    tvButtonColor,
                    chromButtonColor,
                    projectButtonColor,
                    tvWarningLabel, 
                    treeTaxaOptions,
                    treeTaxaValues,
                    gffData,
                    None,
                ]
            # Set file names to just file name and not whole path
            tvFilename = Path(tvFilename).name
            chromFilename = Path(chromFilename).name
            # If Chromosome data type is int -> add chr to make categorical
            if tvDF["Chromosome"].dtype != 'object':
                tvDF["Chromosome"] = tvDF["Chromosome"].apply(lambda x: f"chr{x}")
            if chromDF["Chromosome"].dtype != 'object':
                chromDF["Chromosome"] = chromDF["Chromosome"].apply(lambda x: f"chr{x}")
            # Validate chromosome length file
            msg, validate = tree_utils.validate_chrom_lengths(chromDF, tvDF)
            if not validate:
                # Raise modal telling user there is an issue and they need to fix it
                tvButtonColor = 'warning'
                chromButtonColor = 'warning'
                tvWarningLabel = msg
                return [
                    None,
                    None,
                    modalOpen,
                    tvFilename,
                    chromFilename,
                    projectName,
                    windowSize,
                    tvButtonColor,
                    chromButtonColor,
                    projectButtonColor,
                    tvWarningLabel,
                    treeTaxaOptions,
                    treeTaxaValues,
                    gffData,
                    None,
                ]
            # Calculate needed variables and return to hidden div
            windowSize = data_utils.get_window_size(tvDF["Window"])
            # Fix TopologyID NaN values + remove heterotachy information from NewickTrees
            tvDF[["TopologyID"]] = tvDF[["TopologyID"]].fillna(value="NoData")
            tvDF["NewickTree"] = tvDF["NewickTree"].apply(lambda x: tree_utils.remove_heterotachy_info(x))
            # Collect taxa names from first tree and build option list-dict
            init_tree = tree_utils.get_valid_init_tree(tvDF["NewickTree"])
            treeTaxaOptions = [{"label": i, "value": i} for i in tree_utils.get_taxa_from_tree(init_tree)]
            treeTaxaValues = [i['value'] for i in treeTaxaOptions]
            modalOpen = False
            return [
                tvDF.to_json(),
                chromDF.to_json(),
                modalOpen,
                tvFilename,
                chromFilename,
                projectName,
                windowSize,
                tvButtonColor,
                chromButtonColor,
                projectButtonColor,
                tvWarningLabel,
                treeTaxaOptions,
                treeTaxaValues,
                gffData,
                projectDir,
            ]
        elif (not tvFile_contents) and (not chromFile_contents):
            tvFilename = 'No file chosen!'
            chromFilename = 'No file chosen!'
            tvButtonColor = 'danger'
            chromButtonColor = 'danger'
            return [
                tvDF,
                chromDF,
                modalOpen,
                tvFilename,
                chromFilename,
                projectName,
                windowSize,
                tvButtonColor,
                chromButtonColor,
                projectButtonColor,
                tvWarningLabel, 
                treeTaxaOptions,
                treeTaxaValues,
                gffData,
                None,
            ]
        elif not tvFile_contents:
            tvFilename = 'No file chosen...'
            tvButtonColor = 'danger'
            return [
                tvDF,
                chromDF,
                modalOpen,
                tvFilename,
                chromFilename,
                projectName,
                windowSize,
                tvButtonColor,
                chromButtonColor,
                projectButtonColor,
                tvWarningLabel, 
                treeTaxaOptions,
                treeTaxaValues,
                gffData,
                None,
            ]
        elif not chromFile_contents:
            chromFilename = 'No file chosen...'
            chromButtonColor = 'danger'
            return [
                tvDF,
                chromDF,
                modalOpen,
                tvFilename,
                chromFilename,
                projectName,
                windowSize,
                tvButtonColor,
                chromButtonColor,
                projectButtonColor,
                tvWarningLabel, 
                treeTaxaOptions,
                treeTaxaValues,
                gffData,
                None,
            ]
        else:
            _, tvcontent_string = tvFile_contents.split(",")
            tvDecoded = base64.b64decode(tvcontent_string)
            _, ChromContent_string = chromFile_contents.split(",")
            chromDecoded = base64.b64decode(ChromContent_string)
            if "csv" in tvFileType:
                # Assume that the user uploaded a CSV file
                tvDF = pd.read_csv(io.StringIO(tvDecoded.decode("utf-8")))
                tvDF.sort_values(by=["Chromosome", "Window"], inplace=True)
                validated_tvDF = tree_utils.tv_header_validation(tvDF)
                # If validated_DF returns False, raise error 
                # indicating issue with input file 
                if not validated_tvDF:
                    tvFilename = f"{tvFilename} is malformed!"
                    tvButtonColor = 'danger'
                    chromButtonColor = 'warning'
                    tvWarningLabel = "WARNING: Tree Viewer file appears to be malformed. First four headers must be ['Chromosome', 'Window', 'NewickTree', 'TopologyID']"
                    return [
                        tvDF.to_json(),
                        chromDF,
                        modalOpen,
                        tvFilename,
                        chromFilename,
                        projectName,
                        windowSize,
                        tvButtonColor,
                        chromButtonColor,
                        projectButtonColor,
                        tvWarningLabel, 
                        treeTaxaOptions,
                        treeTaxaValues,
                        gffData,
                        None,
                    ]
                if len(tvDF.columns) == 4:
                    tvDF["None"] = [pd.NA]*len(tvDF)
                tvDF["TopologyID"] = tvDF["TopologyID"].apply(
                    lambda x: "NODATA" if type(x) != str else x)
                pass
            elif "txt" in tvFileType:
                # Assume that the user uploaded a CSV file
                tvDF = pd.read_csv(io.StringIO(tvDecoded.decode("utf-8")), sep='\t')
                tvDF.sort_values(by=["Chromosome", "Window"], inplace=True)
                validated_tvDF = tree_utils.tv_header_validation(tvDF)
                # If validated_DF returns False, raise error 
                # indicating issue with input file 
                if not validated_tvDF:
                    tvFilename = f"{tvFilename} is malformed!"
                    tvButtonColor = 'danger'
                    chromButtonColor = 'warning'
                    tvWarningLabel = "WARNING: Tree Viewer file appears to be malformed. First four headers must be ['Chromosome', 'Window', 'NewickTree', 'TopologyID']"
                    return [
                        tvDF.to_json(),
                        chromDF,
                        modalOpen,
                        tvFilename,
                        chromFilename,
                        projectName,
                        windowSize,
                        tvButtonColor,
                        chromButtonColor,
                        projectButtonColor,
                        tvWarningLabel, 
                        treeTaxaOptions,
                        treeTaxaValues,
                        gffData,
                        None,
                    ]
                if len(tvDF.columns) == 4:
                    tvDF["None"] = [pd.NA]*len(tvDF)
                tvDF["TopologyID"] = tvDF["TopologyID"].apply(
                    lambda x: "NODATA" if type(x) != str else x)
                pass
            elif "tsv" in tvFileType:
                # Assume that the user uploaded a CSV file
                tvDF = pd.read_csv(io.StringIO(tvDecoded.decode("utf-8")), sep='\t')
                tvDF.sort_values(by=["Chromosome", "Window"], inplace=True)
                validated_tvDF = tree_utils.tv_header_validation(tvDF)
                # If validated_DF returns False, raise error 
                # indicating issue with input file 
                if not validated_tvDF:
                    tvFilename = f"{tvFilename} is malformed!"
                    tvButtonColor = 'danger'
                    chromButtonColor = 'warning'
                    tvWarningLabel = "WARNING: Tree Viewer file appears to be malformed. First four headers must be ['Chromosome', 'Window', 'NewickTree', 'TopologyID']"
                    return [
                        tvDF.to_json(),
                        chromDF,
                        modalOpen,
                        tvFilename,
                        chromFilename,
                        projectName,
                        windowSize,
                        tvButtonColor,
                        chromButtonColor,
                        projectButtonColor,
                        tvWarningLabel, 
                        treeTaxaOptions,
                        treeTaxaValues,
                        gffData,
                        None,
                    ]
                if len(tvDF.columns) == 4:
                    tvDF["None"] = [pd.NA]*len(tvDF)
                tvDF["TopologyID"] = tvDF["TopologyID"].apply(
                    lambda x: "NODATA" if type(x) != str else x)
                pass
            elif "xls" in tvFileType:
                # Assume that the user uploaded an excel file
                tvDF = pd.read_excel(io.BytesIO(tvDecoded), engine="openpyxl", comment="#")
                tvDF.sort_values(by=["Chromosome", "Window"], inplace=True)
                validated_tvDF = tree_utils.tv_header_validation(tvDF)
                # If validated_DF returns False, raise error 
                # indicating issue with input file 
                if not validated_tvDF:
                    tvFilename = f"{tvFilename} is malformed!"
                    tvButtonColor = 'danger'
                    chromButtonColor = 'warning'
                    tvWarningLabel = "WARNING: Tree Viewer file appears to be malformed. First four headers must be ['Chromosome', 'Window', 'NewickTree', 'TopologyID']"
                    return [
                        tvDF.to_json(),
                        chromDF,
                        modalOpen,
                        tvFilename,
                        chromFilename,
                        projectName,
                        windowSize,
                        tvButtonColor,
                        chromButtonColor,
                        projectButtonColor,
                        tvWarningLabel, 
                        treeTaxaOptions,
                        treeTaxaValues,
                        gffData,
                        None,
                    ]
                if len(tvDF.columns) == 4:
                    tvDF["None"] = [pd.NA]*len(tvDF)
                tvDF["TopologyID"] = tvDF["TopologyID"].apply(lambda x: "NODATA" if type(x) != str else x)
                tvDF[["TopologyID"]] = tvDF[["TopologyID"]].fillna(value="NoData")
                pass
            else:
                tvFilename = f"Invalid file type - {tvFilename}"
                tvButtonColor = 'danger'
                tvWarningLabel = 'WARNING: Please provide a valid Tree Viewer file'
                return [
                    None,
                    None,
                    modalOpen,
                    tvFilename,
                    chromFilename,
                    projectName,
                    windowSize,
                    tvButtonColor,
                    chromButtonColor,
                    projectButtonColor,
                    tvWarningLabel, 
                    treeTaxaOptions,
                    treeTaxaValues,
                    None,
                    None,
                ]

            if "bed" in chromFileType:
                # Assume that the user uploaded a CSV file
                header = pd.read_csv(
                    io.StringIO(chromDecoded.decode("utf-8")),
                    sep="\t",
                    nrows=1
                )
                if [i for i in header] == ["Chromosome", "Start", "End"]:
                    chromDF = pd.read_csv(
                        io.StringIO(chromDecoded.decode("utf-8")),
                        sep="\t",
                    )
                else:
                    chromFilename = f"Invalid file - {chromFilename}"
                    chromButtonColor = 'danger'
                    tvWarningLabel = 'WARNING: Header appears malformed or missing'
                    return [
                        None,
                        None,
                        modalOpen,
                        tvFilename,
                        chromFilename,
                        projectName,
                        windowSize,
                        tvButtonColor,
                        chromButtonColor,
                        projectButtonColor,
                        tvWarningLabel, 
                        treeTaxaOptions,
                        treeTaxaValues,
                        None,
                        None,
                    ]
                pass
            else:
                chromFilename = f"Invalid file type - {chromFilename}"
                chromButtonColor = 'danger'
                tvWarningLabel = 'WARNING: Please provide a valid BED file'
                return [
                    None,
                    None,
                    modalOpen,
                    tvFilename,
                    chromFilename,
                    projectName,
                    windowSize,
                    tvButtonColor,
                    chromButtonColor,
                    projectButtonColor,
                    tvWarningLabel, 
                    treeTaxaOptions,
                    treeTaxaValues,
                    None,
                    None,
                ]
            # If Chromosome data type is int -> add chr to make categorical
            if tvDF["Chromosome"].dtype != 'object':
                tvDF["Chromosome"] = tvDF["Chromosome"].apply(lambda x: f"chr{x}")
            if chromDF["Chromosome"].dtype != 'object':
                chromDF["Chromosome"] = chromDF["Chromosome"].apply(lambda x: f"chr{x}")
            # Validate chromosome length file
            msg, validate = tree_utils.validate_chrom_lengths(chromDF, tvDF)
            if not validate:
                # Raise modal telling user there is an issue and they need to fix it
                tvButtonColor = 'warning'
                chromButtonColor = 'warning'
                tvWarningLabel = msg
                return [
                    None,
                    None,
                    modalOpen,
                    tvFilename,
                    chromFilename,
                    projectName,
                    windowSize,
                    tvButtonColor,
                    chromButtonColor,
                    projectButtonColor,
                    tvWarningLabel, 
                    treeTaxaOptions,
                    treeTaxaValues,
                    gffData,
                    None,
                ]
            # Calculate needed variables and return to hidden div
            windowSize = data_utils.get_window_size(tvDF["Window"])
            # Fix TopologyID NaN values + remove heterotachy information from NewickTrees
            tvDF[["TopologyID"]] = tvDF[["TopologyID"]].fillna(value="NoData")
            tvDF["NewickTree"] = tvDF["NewickTree"].apply(lambda x: tree_utils.remove_heterotachy_info(x))
            # Collect taxa names from first tree and build option list-dict
            init_tree = tree_utils.get_valid_init_tree(tvDF["NewickTree"])
            treeTaxaOptions = [{"label": i, "value": i} for i in tree_utils.get_taxa_from_tree(init_tree)]
            treeTaxaValues = [i['value'] for i in treeTaxaOptions]
            modalOpen = False
            return [
                tvDF.to_json(),
                chromDF.to_json(),
                modalOpen,
                tvFilename,
                chromFilename,
                projectName,
                windowSize,
                tvButtonColor,
                chromButtonColor,
                projectButtonColor,
                tvWarningLabel, 
                treeTaxaOptions,
                treeTaxaValues,
                None,
                None,
            ]
    # - Return current filenames -
    else:
        # If init new session, return no file chosen messages
        if not tvFilename:
            tvFilename = "No file chosen..."
        if not chromFilename:
            chromFilename = "No file chosen..."
        if not projectName:
            projectName = "No project.ini file chosen..."
        return [
            tvDF,
            chromDF,
            modalOpen,
            tvFilename,
            chromFilename,
            projectName,
            windowSize,
            tvButtonColor,
            chromButtonColor,
            projectButtonColor,
            tvWarningLabel, 
            treeTaxaOptions,
            treeTaxaValues,
            gffData,
            None,
        ]

# ---------------------- Dropdown Option Callbacks ---------------------
# Set topology color mapping
@app.callback(
    [Output("topology-color-chart-store", "data"),
     Output("current-color-list", "data"),
    ],
    [Input("input-data-upload", "children"),
     Input("current-color-list", "data"),
     Input("line_color", "value"),
     Input("step-up-button", "n_clicks"),
     Input("step-down-button", "n_clicks"),
    ],
)
def set_topology_colors(
    data,
    curr_colors,
    line_color,
    btn_up,
    btn_down,
):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if not data:
        return None, dash.no_update
    elif button_id == 'step-up-button':
        step_color_set = curr_colors[1:] + [curr_colors[0]]
        topo_colors = tree_utils.set_topology_colors(data, step_color_set)
        return topo_colors, step_color_set
    elif button_id == 'step-down-button':
        step_color_set = [curr_colors[-1]] + curr_colors[:-1]
        topo_colors = tree_utils.set_topology_colors(data, step_color_set)
        return topo_colors, step_color_set
    else:
        color_set = COLOR_SWATCHES[line_color]
        topo_colors = tree_utils.set_topology_colors(data, color_set)
        return topo_colors, color_set


# Set chromosome options + value
@app.callback(
    [Output("chromosome-options", "options"),
     Output("chromosome-options", "value"),
     Output("sex-chrom-dropdown", "options"),
     Output("sex-chrom-dropdown", "value"),
     Output("topo-quantile-chrom-dropdown", "options"),
     Output("topo-quantile-chrom-dropdown", "value"),
    ],
    [Input("input-data-upload", "children")]
)
def set_chromosome_options(
    tv_input_json,
):
    if not tv_input_json:
        raise PreventUpdate
    else:
        df = pd.read_json(tv_input_json)
        chromosome_options = sorted(
            [{"label": i, "value": i} for i in tree_utils.sorted_nicely(df["Chromosome"].unique())],
            key=lambda i: (i["label"], i["value"])
        )
        chrom_val = chromosome_options[0]["value"]
        sex_chrom_val = chromosome_options[-1]["value"]
        return (
            chromosome_options, chrom_val,
            chromosome_options, sex_chrom_val,
            chromosome_options, chrom_val,
        )


# Set topology dropdown options + value
@app.callback(
    [Output("topology_options", "options"),
     Output("topology_options", "value"), 
     Output("correlation-topologies", "options"),
     Output("correlation-topologies", "value"),
     Output("topo-quantile-topologies", "options"),
     Output("topo-quantile-topologies", "value"),
     Output("topology-freq-order-store", "data"),
    ],
    [Input("input-data-upload", "children"),
     Input("chromFile-data-upload", "children"),
     Input("chromosome-options", "value"),
     Input("topo-freq-order", "value"),
     Input("window-size", "children"),
     Input('current-data-range', 'data'),],
    [State("topology_options", "value"),
    ],
)
def set_topology_options(
    tv_input_json,
    chrom_json,
    chromValue,
    freq_order,
    window_size,
    dataRange,
    topo_vals,
):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if not tv_input_json:
        return [
            [{"label": "None", "value": "None", "disabled": True}], None,
            [{"label": "None", "value": "None", "disabled": True}], None,
            [{"label": "None", "value": "None", "disabled": True}], None,
            None,
        ]
    else:
        df = pd.read_json(tv_input_json)
        # If no topology values present, set according to sorting option (whole-genome vs current view)
        if not topo_vals:
            if not freq_order:
                pass
            else:
                chromosome_df = pd.read_json(chrom_json)
                chromosome_df["End"] = chromosome_df["End"].apply(lambda x: data_utils.roundUp(x, window_size))
                chromosome_df = chromosome_df[chromosome_df["Chromosome"] == chromValue]
                dataMin, dataMax = dataRange
                df = df[(df["Chromosome"] == chromValue) & (df["Window"] >= dataMin) & (df["Window"] <= dataMax)]
            # Topology order is sorted most frequent to least frequent
            value_counts = df["TopologyID"].value_counts()
            topology_options = [{"label": i, "value": i} for i in value_counts.index]
            topoOrder = list(value_counts.index)
            topoOrder.reverse()
            if topo_vals:
                return [
                    topology_options, topo_vals,
                    topology_options, [c["value"] for c in topology_options[:3]],
                    topology_options, [c["value"] for c in topology_options[:3]],
                    topoOrder,
                ]
            else:
                return [
                    topology_options, [c["value"] for c in topology_options[:3]],
                    topology_options, [c["value"] for c in topology_options[:3]],
                    topology_options, [c["value"] for c in topology_options[:3]],
                    topoOrder,
                ]
        # If topologies already loaded, return based on sorting option
        else:
            if button_id == "chromosome-options":
                return [
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                ]
            elif button_id == "current-data-range":
                return [
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                ]
            elif not freq_order:
                value_counts = df["TopologyID"].value_counts()
                topology_options = [{"label": i, "value": i} for i in value_counts.index]
                topoOrder = list(value_counts.index)
                topoOrder.reverse()
                return [
                    topology_options, [c["value"] for c in topology_options[:3]],
                    topology_options, [c["value"] for c in topology_options[:3]],
                    topology_options, [c["value"] for c in topology_options[:3]],
                    topoOrder,
                ]

            else:
                chromosome_df = pd.read_json(chrom_json)
                chromosome_df["End"] = chromosome_df["End"].apply(lambda x: data_utils.roundUp(x, window_size))
                chromosome_df = chromosome_df[chromosome_df["Chromosome"] == chromValue]
                dataMin, dataMax = dataRange
                df = df[(df["Chromosome"] == chromValue) & (df["Window"] >= dataMin) & (df["Window"] <= dataMax)]
                value_counts = df["TopologyID"].value_counts()
                topology_options = [{"label": i, "value": i} for i in value_counts.index]
                topoOrder = list(value_counts.index)
                topoOrder.reverse()
                return [
                    topology_options, [c["value"] for c in topology_options[:3]],
                    topology_options, [c["value"] for c in topology_options[:3]],
                    topology_options, [c["value"] for c in topology_options[:3]],
                    topoOrder,
                ]


# Set alt data dropdown options + value
@app.callback(
    [Output("alt_data_options", "options"),
     Output("alt_data_options", "value"),
     Output("correlation-add-data-type", "options"),
     Output("correlation-add-data-type", "value"),
     Output("data-freq-dist-dropdown", "options"),
     Output("data-freq-dist-dropdown", "value"),
     Output("topo-quantile-additional-data", "options"),
     Output("topo-quantile-additional-data", "value"),
    ],
    [Input("input-data-upload", "children")])
def set_alt_data_options(
    tv_input_json,
):
    if not tv_input_json:
        return (
            [{"label": "None", "value": "None", "disabled": True}],
            None,
            [{"label": "None", "value": "None", "disabled": True}],
            None,
            [{"label": "None", "value": "None", "disabled": True}],
            None,
            [{"label": "None", "value": "None", "disabled": True}],
            None,
        )
    else:
        df = pd.read_json(tv_input_json)
        alt_data_cols = [col for col in df.columns][4:]
        if not alt_data_cols:
            return None, None, None, None, None, None, None, None
        elif (len(alt_data_cols) == 1) and (alt_data_cols[0] == "None"):
            return (
                [{"label": "", "value": ""}], None,
                [{"label": "", "value": ""}], None,
                [{"label": "", "value": ""}], None,
                [{"label": "", "value": ""}], None,
            )
        else:
            alt_data_options = [{"label": i, "value": i} for i in alt_data_cols]
            # Filter non-numeric alt data types
            numeric_alt_data_options = [{"label": i, "value": i} for i in tree_utils.filter_numeric_dtypes(df)]
            return (
                alt_data_options, alt_data_options[0]["value"],
                alt_data_options, alt_data_options[0]["value"],
                numeric_alt_data_options, numeric_alt_data_options[0]["value"],
                numeric_alt_data_options, numeric_alt_data_options[0]["value"],
            )


# Set init chroms per whole-genome graphs
@app.callback(
    Output("whole-genome-per-graph-chrom-count", "value"),
    [Input("whole-genome-per-graph-chrom-count", "value"),
     Input("input-data-upload", "children"),
    ],
)
def set_whole_genome_chrom_count(chrom_count, tv_input_json):
    if not tv_input_json:
        raise PreventUpdate
    elif not chrom_count:
        df = pd.read_json(tv_input_json)
        chrom_count = tree_utils.get_recommended_chrom_num(len(df)) 
    return chrom_count


# ---------------------- Relayout Data Handlers ----------------------
@app.callback(
    # Outputs
    Output("init-graph-trigger", "children"),
    # Inputs
    [Input("input-data-upload", "children"),
     Input("chromFile-data-upload", "children"),
     Input("topology_options", "value"),
    ],
)
def init_graph_building(
    tv_input_json,
    chromosome_lengths,
    current_topologies,
):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if not tv_input_json:
        raise PreventUpdate
    elif not chromosome_lengths:
        raise PreventUpdate
    elif not current_topologies:
        raise PreventUpdate
    elif button_id == "init-graph-trigger":
        raise PreventUpdate
    elif button_id == "topology_options":
        raise PreventUpdate
    else:
        return True 


# ---------------------- Relayout Data Handlers ----------------------
@app.callback(
    # Outputs
    [Output("current-data-range", "data"),
     Output("range-init", "data"),
     Output("main-div-trigger", 'data'),
    ],
    # Inputs
    [Input('topologyGraph', 'relayoutData'),
     Input("chromFile-data-upload", "children"),
     Input("window-size", "children"),
     Input("chromosome-options", "value"),
     Input("range-init", "data"),
    ],
)
def update_relayout(
    relayout_data,
    chromosome_length_data,
    window_size,
    chromosome,
    range_set,
):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if not chromosome_length_data:
        raise PreventUpdate
    relayout_keys = [k for k in relayout_data.keys()]
    chromosome_df = pd.read_json(chromosome_length_data)
    chromosome_df = chromosome_df[chromosome_df["Chromosome"] == chromosome]
    chromosome_df["End"] = chromosome_df["End"].apply(lambda x: data_utils.roundUp(x, window_size))
    try:
        if button_id == "chromosome-options":
            dataRange = [0, max(chromosome_df["End"])]
            return dataRange, True, "trigger"
        elif "xaxis.range" in relayout_keys[0]:
            dataRange = [relayout_data['xaxis.range[0]'], relayout_data['xaxis.range[1]']]
            return dataRange, True, []
        elif "xaxis2.range" in relayout_keys[0]:
            dataRange = [relayout_data['xaxis2.range[0]'], relayout_data['xaxis2.range[1]']]
            return dataRange, True, []
        elif "xaxis.autorange" in relayout_keys[0]:
            dataRange = [0, max(chromosome_df["End"])]
            return dataRange, True, "trigger"
        elif "xaxis2.autorange" in relayout_keys[0]:
            dataRange = [0, max(chromosome_df["End"])]
            return dataRange, True, "trigger"
        else:
            raise KeyError
    except KeyError:
        if not range_set:
            dataRange = [0, max(chromosome_df["End"])]
            return dataRange, True, 'trigger'
        else:
            return dash.no_update, True, []
    except TypeError:
        dataRange = [0, max(chromosome_df["End"])]
        return dataRange, True, []
    except IndexError:
        dataRange = [0, max(chromosome_df["End"])]
        return dataRange, True, []

# ---------------------- Graph Building Callbacks ----------------------
@app.callback(
    # Outputs
    [Output("topologyGraph", "figure"),
     Output("additional-graphs-div", "children"),
     Output("topologyGraph", "config"),
    ],
    [Input("update-options", "n_clicks"),
     Input("current-data-range", "data"),
     Input("toggle-chrom-whole-genome", "value"),
    ],
    # Inputs
    [State("input-data-upload", "children"),
     State("chromFile-data-upload", "children"),
     State("gff-data-upload", "children"),
     State("main-graph-switches", "value"),
     State("alt-graph-switches", "value"),
     State("alt_data_options", "value"),
     State("chromosome-options", "value"),
     State("topology_options", "value"),
     State("template-option", "value"),
     State("snapshot-file-option", "value"),
     State("topology-color-chart-store", "data"),
     State("topology-freq-order-store", "data"),
     State("window-size", "children"),
     State("main-input-modal", "is_open"),
     State("alt-data-switch", "value"),
     State("axis-line-width", "value"),
     State("tv-pixel-height", "value"),
     State("tv-pixel-width", "value"),
     State("axis-gridlines", "value"),
     State("editable-graph-option", "value"),
     State("main-div-trigger", 'data'),
     State("font-family-option", "value"),
    ],
)
def home_tab_graphs(
    update_options,
    dataRange,
    view_toggle,
    #------------------
    tv_input_json,
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
    input_modal_open,
    alt_data_switch,
    axis_line_width,
    pixel_height,
    pixel_width,
    axis_gridlines,
    editable_graphs,
    main_div_trigger,
    font_family,
):
    # -- Check context --
    ctx = dash.callback_context
    button_id = ctx.triggered[0]["prop_id"].split(".")[0] if ctx.triggered else None
    # -- Check to prevent update --
    if input_modal_open:
        raise PreventUpdate
    elif not view_toggle:
        raise PreventUpdate
    elif not tv_input_json:
        raise PreventUpdate
    elif not chromosome_length_data:
        raise PreventUpdate
    elif not current_topologies:
        # If not topologies are selected, return No Data Loaded graph
        return [], None, None
    elif (button_id == 'update-options') or (button_id == 'current-data-range'):
        # -- Run callback --
        topoOrder = [e for e in topoOrder if e in current_topologies]
        whole_tv_input_df = pd.read_json(tv_input_json)
        chromosome_df = pd.read_json(chromosome_length_data)
        chromosome_df["End"] = chromosome_df["End"].apply(lambda x: data_utils.roundUp(x, window_size))
        tv_df = whole_tv_input_df[whole_tv_input_df["Chromosome"] == chromosome]
        chromosome_df = chromosome_df[chromosome_df["Chromosome"] == chromosome]
        dataMin, dataMax = dataRange
        
        # Set gridline bools
        xaxis_gridlines, yaxis_gridlines = tree_utils.get_gridline_bools(axis_gridlines)
        # Filter out data not within selected range
        tv_df = tv_df.sort_values(by=["Chromosome", "Window", "TopologyID"])
        tv_df = tv_df.reset_index(drop=True)
        # -- Set editable value
        editable_graphs = False if not editable_graphs else True
        # --- Generate main distribution ---
        topology_graphs = []
        if button_id == "main-graph-switches":
            main_div_trigger = 'trigger'
        elif button_id == "topology_options":
            main_div_trigger = 'trigger'
        elif button_id == "update-options":
            main_div_trigger = 'trigger'
        if not current_topologies:
            pass
        elif (not main_div_trigger) and (button_id != 'chromosome-options'):
            topo_dist_graph = dash.no_update
        else:
            # Only keep data for selected topologies
            tv_input_df_filtered = tv_df[tv_df["TopologyID"].isin(current_topologies)]
            # Create figure depending on toggle-toolbar, add graph to output list
            if "topo_rug_only" in main_graph_switches:
                # Build histogram + Heatmap Figures
                topo_dist_graph = tree_utils.build_rug_plot(
                    tv_input_df_filtered,
                    chromosome,
                    template,
                    current_topologies,
                    color_mapping,
                    dataRange,
                    topoOrder,
                    axis_line_width,
                    xaxis_gridlines,
                    yaxis_gridlines,
                    font_family,
                )
            elif "topo_tile_only" in main_graph_switches:
                # Build histogram + Heatmap Figures
                topo_dist_graph = tree_utils.build_tile_plot(
                    tv_input_df_filtered,
                    chromosome_df,
                    template,
                    current_topologies,
                    color_mapping,
                    dataRange,
                    window_size,
                    axis_line_width,
                    xaxis_gridlines,
                    yaxis_gridlines,
                    font_family,
                )
            else:
                # Build histogram + Heatmap Figures
                topo_dist_graph = tree_utils.build_histogram_with_rug_plot(
                    tv_input_df_filtered,
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
                    font_family,
                )

        # --- Generate alternative data graphs ---
        if not alt_data_switch:
            pass
        else:
            alt_dropdown_options = [alt_dropdown_options] if type(alt_dropdown_options) == str else alt_dropdown_options
            # Collect y-max from additional int/float data
            y_maxes = tree_utils.get_y_max_list(alt_dropdown_options, tv_df)
            # Generate graphs
            for y_max, alt_data_option in zip(y_maxes, alt_dropdown_options):
                alt_graph = tree_utils.build_alt_data_graph(
                    alt_data_option,
                    chromosome_df,
                    color_mapping,
                    tv_df,
                    window_size,
                    template,
                    dataRange,
                    y_max,
                    axis_line_width,
                    xaxis_gridlines,
                    yaxis_gridlines,
                    font_family,
                )
                topology_graphs.append(
                    html.Div(
                        dcc.Graph(
                            figure=alt_graph,
                            config=dict(
                                editable=True,
                                displaylogo=False,
                                doubleClick="reset",
                                showAxisDragHandles=True,
                                showAxisRangeEntryBoxes=True,
                                modeBarButtonsToRemove=[
                                    # "zoomIn2d",
                                    # "zoomOut2d",
                                    "autoScale2d",
                                    "lasso2d",
                                    "select2d",
                                ],
                                toImageButtonOptions=dict(
                                    format=snapshot_file_type,
                                    filename=f"{alt_data_option}",
                                ),
                            ),
                            className="alt-data-graph",
                        )
                    )
                )

        # --- Generate GFF graph ---
        if not alt_graph_switches:
            pass
        else:
            for graph in alt_graph_switches:
                if graph == "gff_switch":
                    for gff_json in gff_data:
                        try:
                            gff_df = pd.read_json(gff_json)
                        except ValueError:
                            raise PreventUpdate
                        # gff_df = gff_df[(gff_df["start"] >= dataMin) & (gff_df["end"] <= dataMax) & (gff_df["chromosome"] == chromosome)]
                        gff_df = gff_df[(gff_df["chromosome"] == chromosome)]
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
                                font_family,
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

        # --- Set configs ---
        # Set pixel default
        if (not pixel_height) or (not pixel_width):
            pixel_width = 1500
            pixel_height = 400

        tv_main_dist_config = dict(
            doubleClick="autosize",
            displaylogo=False,
            editable=editable_graphs,
            modeBarButtonsToRemove=[
                "resetScale2d",
                "pan2d",
                "lasso2d",
                "select2d",
                # "autoScale2d",
            ],
            toImageButtonOptions=dict(
                format=snapshot_file_type,
                filename="Graph_Name",
                width=int(pixel_width),
                height=int(pixel_height),
            ),
        )
        return topo_dist_graph, topology_graphs, tv_main_dist_config
    elif (view_toggle) and (button_id == "toggle-chrom-whole-genome"):
        # -- Run callback --
        topoOrder = [e for e in topoOrder if e in current_topologies]
        whole_tv_input_df = pd.read_json(tv_input_json)
        chromosome_df = pd.read_json(chromosome_length_data)
        chromosome_df["End"] = chromosome_df["End"].apply(lambda x: data_utils.roundUp(x, window_size))
        tv_df = whole_tv_input_df[whole_tv_input_df["Chromosome"] == chromosome]
        chromosome_df = chromosome_df[chromosome_df["Chromosome"] == chromosome]
        dataMin, dataMax = dataRange
        # Set gridline bools
        xaxis_gridlines, yaxis_gridlines = tree_utils.get_gridline_bools(axis_gridlines)
        # Filter out data not within selected range
        tv_df = tv_df.sort_values(by=["Chromosome", "Window", "TopologyID"])
        tv_df = tv_df.reset_index(drop=True)
        tv_input_df_filtered = tv_df[tv_df["TopologyID"].isin(current_topologies)]
        # -- Set editable value
        editable_graphs = False if not editable_graphs else True
        # --- Generate main distribution ---
        topo_dist_graph = tree_utils.build_histogram_with_rug_plot(
                tv_input_df_filtered,
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
                font_family,
            )
        # --- Generate alternative data graphs ---
        topology_graphs = []
        if not alt_data_switch:
            pass
        else:
            alt_dropdown_options = [alt_dropdown_options] if type(alt_dropdown_options) == str else alt_dropdown_options
            # Collect y-max from additional int/float data
            y_maxes = tree_utils.get_y_max_list(alt_dropdown_options, tv_df)
            # Generate graphs
            for y_max, alt_data_option in zip(y_maxes, alt_dropdown_options):
                alt_graph = tree_utils.build_alt_data_graph(
                    alt_data_option,
                    chromosome_df,
                    color_mapping,
                    tv_df,
                    window_size,
                    template,
                    dataRange,
                    y_max,
                    axis_line_width,
                    xaxis_gridlines,
                    yaxis_gridlines,
                    font_family,
                )
                topology_graphs.append(
                    html.Div(
                        dcc.Graph(
                            figure=alt_graph,
                            config=dict(
                                editable=True,
                                displaylogo=False,
                                doubleClick="reset",
                                showAxisDragHandles=True,
                                showAxisRangeEntryBoxes=True,
                                modeBarButtonsToRemove=[
                                    # "zoomIn2d",
                                    # "zoomOut2d",
                                    "autoScale2d",
                                    "lasso2d",
                                    "select2d",
                                ],
                                toImageButtonOptions=dict(
                                    format=snapshot_file_type,
                                    filename=f"{alt_data_option}",
                                ),
                            ),
                            className="alt-data-graph",
                        )
                    )
                )

        # --- Generate GFF graph ---
        if not alt_graph_switches:
            pass
        else:
            for graph in alt_graph_switches:
                if graph == "gff_switch":
                    for gff_json in gff_data:
                        try:
                            gff_df = pd.read_json(gff_json)
                        except ValueError:
                            raise PreventUpdate
                        # gff_df = gff_df[(gff_df["start"] >= dataMin) & (gff_df["end"] <= dataMax) & (gff_df["chromosome"] == chromosome)]
                        gff_df = gff_df[(gff_df["chromosome"] == chromosome)]
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
                                font_family,
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

        # --- Set configs ---
        # Set pixel default
        if (not pixel_height) or (not pixel_width):
            pixel_width = 1500
            pixel_height = 400

        tv_main_dist_config = dict(
            doubleClick="autosize",
            displaylogo=False,
            editable=editable_graphs,
            modeBarButtonsToRemove=[
                "resetScale2d",
                "pan2d",
                "lasso2d",
                "select2d",
                # "autoScale2d",
            ],
            toImageButtonOptions=dict(
                format=snapshot_file_type,
                filename="Graph_Name",
                width=int(pixel_width),
                height=int(pixel_height),
            ),
        )
        
        return topo_dist_graph, topology_graphs, tv_main_dist_config
    else:
        raise PreventUpdate


@app.callback(
    [Output("whole-genome-graphs", "children"),
     Output("whole-genome-alt-graphs", "children"),
     Output("break-init-graph-trigger", "children"),
    ],
    [
        Input("update-options", "n_clicks"),
        Input("init-graph-trigger", "children"),
        Input("toggle-chrom-whole-genome", "value"),
    ],
    [
        State("input-data-upload", "children"),
        State("chromFile-data-upload", "children"),
        State("template-option", "value"),
        State("topology-color-chart-store", "data"),
        State("whole-genome-graph-type", 'value'),
        State("topology-freq-order-store", "data"),
        State("topology_options", "value"),
        State("window-size", "children"),
        State("snapshot-file-option", "value"),
        State("axis-line-width", "value"),
        State("tv-pixel-height", "value"),
        State("tv-pixel-width", "value"),
        State("axis-gridlines", "value"),
        State("editable-graph-option", "value"),
        State("wg-squish-expand", 'value'),
        State("whole-genome-per-graph-chrom-count", "value"),
        State("font-family-option", "value"),
        State("alt-data-switch", "value"),
        State("alt_data_options", "value"),
        State("break-init-graph-trigger", "children"),
    ],
)
def whole_genome_plot(
    update_options_btn,
    init_graphs,
    view_toggle,
    tv_input_json,
    chromosome_lengths,
    template,
    color_mapping,
    chrom_plot_type,
    topoOrder,
    current_topologies,
    window_size,
    snapshot_file_type,
    axis_line_width,
    pixel_height,
    pixel_width,
    axis_gridlines,
    editable_graphs,
    wg_squish_expand,
    num_chroms_per_graph,
    font_family,
    alt_data_switch,
    alt_dropdown_options,
    break_init,
):
    
    ctx = dash.callback_context
    button_id = ctx.triggered[0]["prop_id"].split(".")[0] if ctx.triggered else None
    # -- Check to prevent update --
    if view_toggle:
        raise PreventUpdate
    elif not tv_input_json:
        raise PreventUpdate
    elif not chromosome_lengths:
        raise PreventUpdate
    elif not current_topologies:
        raise PreventUpdate
    elif (button_id == 'update-options'):
        # -- Set editable value --
        editable_graphs = False if not editable_graphs else True
        # -- Read in json data --
        df = pd.read_json(tv_input_json)
        chromosome_df = pd.read_json(chromosome_lengths)
        sorted_chromosomes = tree_utils.sorted_nicely(chromosome_df['Chromosome'].unique())
        # -- Set up graph config --
        if (not pixel_height) or (not pixel_width):
            # pixel_width, pixel_height = 1500, 1123
            config_setting = dict(
                doubleClick="reset",
                displaylogo=False,
                editable=editable_graphs,
                modeBarButtonsToRemove=[
                    "zoomin2D"
                    "autoScale2d",
                    # "resetScale2d",
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
            config_setting = dict(
                doubleClick="reset",
                displaylogo=False,
                editable=editable_graphs,
                modeBarButtonsToRemove=[
                    "zoomin2D"
                    "autoScale2d",
                    # "resetScale2d",
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
        # -- Set gridline bools --
        xaxis_gridlines, yaxis_gridlines = tree_utils.get_gridline_bools(axis_gridlines)
        df["GenomeFrequency"] = df["TopologyID"].map(df["TopologyID"].value_counts()/len(df))
        df_grouped = df.groupby(by="Chromosome")
        # -- Round up chromosome lengths for graph range --
        chromosome_df["End"] = chromosome_df["End"].apply(lambda x: data_utils.roundUp(x, window_size))
        # Remove topologies not currently chosen
        topoOrder = [e for e in topoOrder if e in current_topologies]
        # -- Set whole genome graph collection list --
        whole_genome_graphs=[]
        alt_data_graphs=[]
        # -- Group Chromosomes into groups --
        chrom_groups = tree_utils.mygrouper(num_chroms_per_graph, sorted_chromosomes)
        if chrom_plot_type == 'bar':
            topoFreqTable = tree_utils.make_topo_freq_table(df_grouped)
            # Iterate through chromosome groups + append graph
            for chromGroup in chrom_groups:
                chrom_order = [i for i in sorted_chromosomes if i in chromGroup]
                fig = tree_utils.build_whole_genome_bar_plot(
                    topoFreqTable,
                    chrom_order,
                    chromGroup,
                    template,
                    color_mapping,
                    current_topologies,
                    axis_line_width,
                    xaxis_gridlines,
                    yaxis_gridlines,
                    font_family,
                )
                whole_genome_graphs.append(
                    html.Div([
                        dcc.Graph(
                            figure=fig,
                            config=config_setting,
                            className="stats-2-graph-style",
                        )
                    ])
                )
                continue
        elif chrom_plot_type == 'pie':
            topoFreqTable = tree_utils.make_topo_freq_table(df_grouped)
            # Iterate through chromosome groups + append graph
            for chromGroup in chrom_groups:
                chrom_order = [i for i in sorted_chromosomes if i in chromGroup]
                fig = tree_utils.build_whole_genome_pie_charts(
                    topoFreqTable,
                    chrom_order,
                    template,
                    color_mapping,
                    chromGroup,
                    font_family,
                )
                whole_genome_graphs.append(
                    html.Div([
                        dcc.Graph(
                            figure=fig,
                            config=config_setting,
                            className="stats-2-graph-style",
                        )
                    ])
                )
                continue
        elif chrom_plot_type == 'rug':
            # Iterate through chromosome groups + append graph
            for chromGroup in chrom_groups:
                chrom_order = [i for i in sorted_chromosomes if i in chromGroup]
                fig = tree_utils.build_whole_genome_rug_plot(
                    df,
                    chromosome_df,
                    chrom_order,
                    chromGroup,
                    template,
                    color_mapping,
                    current_topologies,
                    topoOrder,
                    window_size,
                    axis_line_width,
                    xaxis_gridlines,
                    yaxis_gridlines,
                    wg_squish_expand,
                    font_family,
                )
                whole_genome_graphs.append(
                    html.Div([
                        dcc.Graph(
                            figure=fig,
                            config=config_setting,
                        )
                    ])
                )
                continue
        elif chrom_plot_type == 'tile':
            # Iterate through chromosome groups + append graph
            for chromGroup in chrom_groups:
                chrom_order = [i for i in sorted_chromosomes if i in chromGroup]
                fig = tree_utils.build_whole_genome_tile_plot(
                    df,
                    chromosome_df,
                    chrom_order,
                    template,
                    color_mapping,
                    current_topologies,
                    topoOrder,
                    window_size,
                    axis_line_width,
                    chromGroup,
                    xaxis_gridlines,
                    yaxis_gridlines,
                    wg_squish_expand,
                    font_family,
                )
                whole_genome_graphs.append(
                    html.Div([
                        dcc.Graph(
                            figure=fig,
                            config=config_setting,
                        )
                    ])
                )
                continue

        # -- Generate alternative data graphs --
        if not alt_data_switch:
            pass
        else:
            alt_dropdown_options = [alt_dropdown_options] if type(alt_dropdown_options) == str else alt_dropdown_options
            # Collect y-max from additional int/float data
            y_maxes = tree_utils.get_y_max_list(alt_dropdown_options, df)
            # Generate graphs
            for y_max, alt_data_option in zip(y_maxes, alt_dropdown_options):
                alt_graph = tree_utils.build_whole_genome_alt_data_graph(
                    alt_data_option,
                    chromosome_df,
                    color_mapping,
                    df,
                    window_size,
                    template,
                    y_max,
                    axis_line_width,
                    xaxis_gridlines,
                    yaxis_gridlines,
                    font_family,
                )
                alt_data_graphs.append(
                    html.Div(
                        dcc.Graph(
                            figure=alt_graph,
                            config=dict(
                                editable=True,
                                displaylogo=False,
                                doubleClick="reset",
                                modeBarButtonsToRemove=[
                                    "autoScale2d",
                                    "lasso2d",
                                    "select2d",
                                ],
                                toImageButtonOptions=dict(
                                    format=snapshot_file_type,
                                    filename=f"{alt_data_option}",
                                ),
                            ),
                            className="stats-2-graph-style",
                        ),
                        style={'borderTop': '2px solid orange', 'paddingtop': '2px'},
                    )
                )
        
        return whole_genome_graphs, alt_data_graphs, dash.no_update
    elif (init_graphs) and (not break_init):
        # -- Set editable value --
        editable_graphs = False if not editable_graphs else True
        # -- Read in json data --
        df = pd.read_json(tv_input_json)
        chromosome_df = pd.read_json(chromosome_lengths)
        sorted_chromosomes = tree_utils.sorted_nicely(chromosome_df['Chromosome'].unique())
        # -- Set up graph config --
        if (not pixel_height) or (not pixel_width):
            # pixel_width, pixel_height = 1500, 1123
            config_setting = dict(
                doubleClick="reset",
                displaylogo=False,
                editable=editable_graphs,
                modeBarButtonsToRemove=[
                    "zoomin2D"
                    "autoScale2d",
                    # "resetScale2d",
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
            config_setting = dict(
                doubleClick="reset",
                displaylogo=False,
                editable=editable_graphs,
                modeBarButtonsToRemove=[
                    "zoomin2D"
                    "autoScale2d",
                    # "resetScale2d",
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
        # -- Set gridline bools --
        xaxis_gridlines, yaxis_gridlines = tree_utils.get_gridline_bools(axis_gridlines)
        df["GenomeFrequency"] = df["TopologyID"].map(df["TopologyID"].value_counts()/len(df))
        df_grouped = df.groupby(by="Chromosome")
        # -- Round up chromosome lengths for graph range --
        chromosome_df["End"] = chromosome_df["End"].apply(lambda x: data_utils.roundUp(x, window_size))
        # Remove topologies not currently chosen
        topoOrder = [e for e in topoOrder if e in current_topologies]
        # -- Set whole genome graph collection list --
        whole_genome_graphs=[]
        alt_data_graphs=[]
        # -- Group Chromosomes into groups --
        chrom_groups = tree_utils.mygrouper(num_chroms_per_graph, sorted_chromosomes)
        for chromGroup in chrom_groups:
            chrom_order = [i for i in sorted_chromosomes if i in chromGroup]
            fig = tree_utils.build_whole_genome_rug_plot(
                df,
                chromosome_df,
                chrom_order,
                chromGroup,
                template,
                color_mapping,
                current_topologies,
                topoOrder,
                window_size,
                axis_line_width,
                xaxis_gridlines,
                yaxis_gridlines,
                wg_squish_expand,
                font_family,
            )
            whole_genome_graphs.append(
                html.Div([
                    dcc.Graph(
                        figure=fig,
                        config=config_setting,
                    )
                ])
            )
            continue
        # -- Generate alternative data graphs --
        if not alt_data_switch:
            pass
        else:
            alt_dropdown_options = [alt_dropdown_options] if type(alt_dropdown_options) == str else alt_dropdown_options
            # Collect y-max from additional int/float data
            y_maxes = tree_utils.get_y_max_list(alt_dropdown_options, df)
            # Generate graphs
            for y_max, alt_data_option in zip(y_maxes, alt_dropdown_options):
                alt_graph = tree_utils.build_whole_genome_alt_data_graph(
                    alt_data_option,
                    chromosome_df,
                    color_mapping,
                    df,
                    window_size,
                    template,
                    y_max,
                    axis_line_width,
                    xaxis_gridlines,
                    yaxis_gridlines,
                    font_family,
                )
                alt_data_graphs.append(
                    html.Div(
                        dcc.Graph(
                            figure=alt_graph,
                            config=dict(
                                editable=True,
                                displaylogo=False,
                                doubleClick="reset",
                                modeBarButtonsToRemove=[
                                    "autoScale2d",
                                    "lasso2d",
                                    "select2d",
                                ],
                                toImageButtonOptions=dict(
                                    format=snapshot_file_type,
                                    filename=f"{alt_data_option}",
                                ),
                            ),
                            className="stats-2-graph-style",
                        ),
                        style={'borderTop': '2px solid orange', 'paddingtop': '2px'},
                    )
                )
        
        return whole_genome_graphs, alt_data_graphs, True
    elif not view_toggle:
        raise PreventUpdate
    else:
        raise PreventUpdate


@app.callback(
    # Outputs
    Output("tree-graphs-div", "children"),
    # Inputs
    [Input("update-options", "n_clicks"),
     Input("input-data-upload", "children"),
     Input("chromFile-data-upload", "children"),
     Input("alt_data_options", "value"),
     Input("chromosome-options", "value"),
     Input("topology_options", "value"),
     Input("template-option", "value"),
     Input("snapshot-file-option", "value"),
     Input("topology-color-chart-store", "data"),
     Input("window-size", "children"),
     Input("main-input-modal", "is_open"),
     Input("tree-graph-switch", "value"),
     Input("tree-taxa-choices", "value"),
     Input("editable-graph-option", "value"),
     Input("toggle-chrom-whole-genome", "value"),
     Input('topologyGraph', 'clickData'),
     Input("font-family-option", "value"),
     Input("tree-shape", "value"),
    ],
)
def generate_tree_graphs(
    update_btn,
    tv_input_json,
    chromosome_length_data,
    alt_dropdown_options,
    chromosome,
    current_topologies,
    template,
    snapshot_file_type,
    color_mapping,
    window_size,
    input_modal_open,
    tree_switch,
    tree_taxa_values,
    editable_graphs,
    view_toggle,
    tv_hoverData,
    font_family,
    tree_shape,
):
    ctx = dash.callback_context
    button_id = ctx.triggered[0]["prop_id"].split(".")[0] if ctx.triggered else None
    if not tree_switch:
        return None
    elif input_modal_open:
        raise PreventUpdate
    elif not tv_input_json:
        raise PreventUpdate
    elif not chromosome_length_data:
        raise PreventUpdate
    elif not current_topologies:
        # If not topologies are selected, return No Data Loaded graph
        return None
    elif button_id != "update-options":
        raise PreventUpdate
    else:
        tv_df = pd.read_json(tv_input_json)
        chromosome_df = pd.read_json(chromosome_length_data)
        chromosome_df["End"] = chromosome_df["End"].apply(lambda x: data_utils.roundUp(x, window_size))
        topology_graphs = []
        # Put alt data option into list if only one selected.
        alt_dropdown_options = [alt_dropdown_options] if type(alt_dropdown_options) == str else alt_dropdown_options         
        # Filter out data not within selected range
        tv_df = tv_df.sort_values(by=["Chromosome", "Window", "TopologyID"])
        tv_df = tv_df.reset_index(drop=True)
        chrom_info = tv_df[(tv_df["Chromosome"] == chromosome)]
        chrom_info.reset_index(inplace=True,)
        # --- Set editable value ---
        editable_graphs = False if not editable_graphs else True
        # --- Generate tree plots ---
        tree_divs = []
        if (type(current_topologies) == str) or (type(current_topologies) == int):
            current_topologies = [current_topologies]
        if 'dark' in template:
            tree_title_className = "tree-title-dark"
        else:
            tree_title_className = "tree-title-light"

        # Create hoverData tree
        if tv_hoverData:
            skip = False
            x_pos = tree_utils.get_Treexpos(tv_hoverData, tv_df)
            poi_info = chrom_info[chrom_info["Window"] == x_pos]
            try:
                curr_tree = poi_info["NewickTree"].tolist()[0]
                curr_topo = poi_info["TopologyID"].tolist()[0]
            except IndexError:
                curr_tree = None

            try:
                pruned_tree = Tree(curr_tree)
                pruned_tree.prune(tree_taxa_values, preserve_branch_length=True)
            # No data for topology on current chromosome
            except IndexError:
                tree_divs.append(
                    dbc.Col(
                        children=[
                            html.Div(f"{curr_topo}", className=tree_title_className),
                            html.Div(
                                children=[
                                    dcc.Graph(
                                        figure=tree_utils.no_tree_data(template, "Topology not present on CHR"),
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
                        ], width=3,
                    ),
                )
                skip = True
                pass
            except NewickError:
                tree_divs.append(
                    dbc.Col(
                        children=[
                            html.Div(f"{curr_topo}", className=tree_title_className),
                            html.Div(
                                children=[
                                    dcc.Graph(
                                        figure=tree_utils.no_tree_data(template, "No Data"),
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
                        ], width=3,
                    ),
                )
                skip = True
                pass
            except TreeError:
                tree_divs.append(
                    dbc.Col(
                        children=[
                            html.Div(f"{curr_topo}", className=tree_title_className),
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
                        ], width=3,
                    ),
                )
                skip = True
                pass
            except ValueError:
                # Assumes taxa in dropdown selection
                # is not found in a particular topology/tree
                # Solution is to check list and remove taxa
                # not present in tree
                tree_taxa = pruned_tree.get_leaf_names()
                trimmed_taxa_list = [t for t in tree_taxa_values if t in tree_taxa]
                pruned_tree.prune(trimmed_taxa_list, preserve_branch_length=True)
                skip = True

            if skip:
                pass
            else:
                tree = tree_utils.DrawTree(io.StringIO(pruned_tree.write()), template, curr_topo, color_mapping, False, font_family)
                # Generate Tree figrues
                if tree_shape == 'rectangle':
                    hoverfig = tree.create_square_tree()
                elif tree_shape == 'circle':
                    hoverfig = tree.create_circular_tree()

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
                            html.Div(f"{curr_topo}: {x_pos}", className=tree_title_className),
                            html.Div(
                                children=[
                                    dcc.Graph(
                                        figure=hoverfig,
                                        style=tree_style,
                                        config=dict(
                                            doubleClick="autosize",
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
                                                "autoScale2d",
                                                "resetScale2d",
                                            ],
                                            displaylogo=False,
                                        ),
                                    )
                                ],
                                className="tree-div"
                            ),
                        ], width=3,
                    ),
                )
        else:
            pass

        for topo in current_topologies:
            # Filter data
            curr_topo_df = tv_df[tv_df["TopologyID"] == topo]

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
                                        figure=tree_utils.no_tree_data(template, "Topology not present on CHR"),
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
                        ], width=3,
                    ),
                )
                continue
            except NewickError:
                tree_divs.append(
                    dbc.Col(
                        children=[
                            html.Div(f"{topo}", className=tree_title_className),
                            html.Div(
                                children=[
                                    dcc.Graph(
                                        figure=tree_utils.no_tree_data(template, "No Data"),
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
                        ], width=3,
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
                        ], width=3,
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

            tree = tree_utils.DrawTree(io.StringIO(first_tree.write()), template, topo, color_mapping, True, font_family)

            # Generate Tree figrues
            if tree_shape == 'rectangle':
                fig = tree.create_square_tree()
            elif tree_shape == 'circle':
                fig = tree.create_circular_tree()
            # Set height dependent on number of taxa
            if len(tree_taxa_values) < 20:
                tree_height = "40vh"
                col_width = 3
            elif len(tree_taxa_values) < 80:
                tree_height = "80vh"
                col_width = 6
            elif len(tree_taxa_values) > 80:
                tree_height = "150vh"
                col_width = 6

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
                                            scale=1,
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
                    width=col_width,
                ),
            )
            continue
        topology_graphs.append(
            dbc.Row(
                children=tree_divs,
                no_gutters=True,
                style={"borderTop": "2px solid orange", "paddingTop": "2px"},
            )
        )
        return topology_graphs


@app.callback(
    [Output("normRF-fig", "figure"),
     Output("normRF-fig", "config"),
    ],
    [Input("template-option", "value"),
     Input("topology-color-chart-store", "data"),
     Input("toggle-chrom-whole-genome", "value"),
     Input("snapshot-file-option", "value"),
     Input("axis-line-width", "value"),
     Input("tv-pixel-height", "value"),
     Input("tv-pixel-width", "value"),
     Input("editable-graph-option", "value"),
     Input("topology_options", "value"),
     Input("font-family-option", "value"),
    ],
    [State("input-data-upload", "children"),
    ],
)
def NormRF(
    template,
    color_mapping,
    view_toggle,
    snapshot_file_type,
    axis_line_width,
    pixel_height, 
    pixel_width,
    editable_graphs,
    topologies,
    font_family,
    tv_input_json,
):
    def get_topo_trees(df, top_3_topos):
        trees = []
        for t in top_3_topos:
            topo_df = df[df["TopologyID"] == t]
            topo_df.reset_index(drop=True, inplace=True)
            trees.append(topo_df["NewickTree"][0])
        return trees

    def calc_rf_dist(ref, alt):
        # --- ensure tree is valid --- 
        if ref == 'NoTree':
            return 0
        elif alt == 'NoTree':
            return 0
        ref = Tree(ref)
        alt = Tree(alt)
        try:
            compare_results = ref.compare(alt)
        except TreeError:
            compare_results = ref.compare(alt, unrooted=True)
        return compare_results["norm_rf"]

    if view_toggle:
        raise PreventUpdate
    elif not tv_input_json:
        raise PreventUpdate
    elif not topologies:
        raise PreventUpdate
    else:
        # Set editable value
        editable_graphs = False if not editable_graphs else True
        # Load session data
        df = pd.read_json(tv_input_json)
        # Create pie chart
        topo_freq_df = pd.DataFrame(df["TopologyID"].value_counts()/len(df))
        topo_freq_df = topo_freq_df.reset_index()
        topo_freq_df.columns = ['TopologyID', 'Frequency']
        topo_freq_df.sort_values(by="Frequency")

        # Create RF graph
        topos = topo_freq_df["TopologyID"][0:len(topologies)]
        topo_trees = get_topo_trees(df, topos)
        rf_df = pd.DataFrame(
            {"TopologyID": topos[1:], "normRF-Distance": [calc_rf_dist(topo_trees[0], t) for t in topo_trees[1:]]}
        )
        rf_fig = tree_utils.build_rf_graph(rf_df, topos[0], template, color_mapping, axis_line_width, font_family)
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
        return rf_fig, config


@app.callback(
    Output("quantile-graphs-div", "children"),
    [Input("update-options", "n_clicks"),
     Input("input-data-upload", "children"),
     Input("chromFile-data-upload", "children"),
     Input("template-option", "value"),
     Input("topology-color-chart-store", "data"),
     Input("topology_options", "value"),
     Input("window-size", "children"),
     Input("toggle-chrom-whole-genome", "value"),
     Input("snapshot-file-option", "value"),
     Input("axis-line-width", "value"),
     Input("tv-pixel-height", "value"),
     Input("tv-pixel-width", "value"),
     Input("axis-gridlines", "value"),
     Input("editable-graph-option", "value"),
     Input("n-quantiles", "value"),
     Input("quantile-toggle", "value"),
     Input("chromosome-options", "value"),
     Input("font-family-option", "value"),
    ],
)
def quantile_graph(
    update_btn,
    tv_input_json,
    chromosome_lengths,
    template,
    color_mapping,
    current_topologies,
    window_size,
    view_toggle,
    snapshot_file_type,
    axis_line_width,
    pixel_height,
    pixel_width,
    axis_gridlines,
    editable_graphs,
    n_quantiles,
    quantile_toggle,
    chromosome,
    font_family,
):
    # ID context
    ctx = dash.callback_context
    button_id = ctx.triggered[0]["prop_id"].split(".")[0] if ctx.triggered else None

    # Run
    run = True
    if button_id == "quantile-toggle": # By pass run button if toggle is turned on/off
        run = True
    if not quantile_toggle:
        run = False
    if not tv_input_json:
        run = False
    if not chromosome_lengths:
        run = False
    # Generate graph if run = True
    if button_id == "n-quantiles":
        return dash.no_update
    elif button_id != "update-options":
        raise PreventUpdate
    elif run:
        # Set default quantile number if not provided
        if not n_quantiles:
            n_quantiles = 25
        # Set editable value
        editable_graphs = False if not editable_graphs else True
        # Load df data
        df = pd.read_json(tv_input_json)
        chromosome_df = pd.read_json(chromosome_lengths)
        # Determine if whole genome or single chromosome
        if view_toggle:
           df = df[df["Chromosome"] == chromosome]
           chromosome_df = chromosome_df[chromosome_df["Chromosome"] == chromosome]
        # Filter out topologies not currently selected
        df = df[df["TopologyID"].isin(current_topologies)]
        quantileCoordinates = tree_utils.get_quantile_coordinates(chromosome_df, n_quantiles, window_size)
        quantileFrequencies = tree_utils.calculateFrequencies(quantileCoordinates, df, chromosome_df, n_quantiles)
        xaxis_gridlines, yaxis_gridlines = tree_utils.get_gridline_bools(axis_gridlines)
        fig = tree_utils.plot_frequencies(quantileFrequencies, n_quantiles, template, color_mapping, axis_line_width, xaxis_gridlines, yaxis_gridlines)
        if (not pixel_height) or (not pixel_width):
            config = dict(
                doubleClick="autosize",
                displaylogo=False,
                editable=editable_graphs,
                modeBarButtonsToRemove=[
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
                modeBarButtonsToRemove=[
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
        graph = html.Div([
            dcc.Graph(
                figure=fig,
                config=config,
            )
        ], style={'border-top': 'solid 2px orange'})
        return [graph]
    else:
        return None


@app.callback(
    [Output("topo-quantile-graph-1", "figure"),
     Output("topo-quantile-graph-1", "config"),
     Output("topo-quantile-graph-2", "figure"),
     Output("topo-quantile-graph-2", "config"),
     Output("topo-quantile-dt-1", "data"),
     Output("topo-quantile-dt-1", "columns"),
     Output("topo-quantile-dt-2", "data"),
     Output("topo-quantile-dt-2", "columns"),
    ],
    [Input("input-data-upload", "children"),
     Input("chromFile-data-upload", "children"),
     Input("template-option", "value"),
     Input("topology-color-chart-store", "data"),
     Input("topo-quantile-topologies", "value"),
     Input("snapshot-file-option", "value"),
     Input("axis-line-width", "value"),
     Input("tv-pixel-height", "value"),
     Input("tv-pixel-width", "value"),
     Input("axis-gridlines", "value"),
     Input("editable-graph-option", "value"),
     Input("topo-n-quantiles", "value"),
     Input("chromosome-options", "value"),
     Input("topo-quantile-run-btn", "n_clicks"),
     Input("topo-quantile-additional-data", "value"),
     Input("topo-quantile-range-selection", "value"),
     Input("topo-quantile-chrom-dropdown", "value"),
     Input("sex-chrom-dropdown", "value"),
     Input("font-family-option", "value"),
    ],
)
def topology_quantile(
    tv_input_json,
    chromosome_lengths,
    template,
    color_mapping,
    current_topologies,
    snapshot_file_type,
    axis_line_width,
    pixel_height,
    pixel_width,
    axis_gridlines,
    editable_graphs,
    n_quantiles,
    chromosome,
    submit_btn,
    additional_data,
    range_selection,
    chromosome_selection,
    sex_chrom_selection,
    font_family,
):
    # ID context
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    # Run
    run = False
    if button_id == "topo-quantile-run-btn":
        run = True
    if button_id == "template-option":
        run = True
    if button_id == "topology-color-chart-store":
        run = True
    if button_id == "editable-graph-option":
        run = True
    if button_id == "snapshot-file-option":
        run = True
    if button_id == "axis-gridlines":
        run = True
    if button_id == "axis-line-width":
        run = True
    if not tv_input_json:
        run = False
    if not chromosome_lengths:
        run = False
    if not additional_data:
        run = False
    # Generate graph if run = True
    if run:
        # Set default quantile number if not provided
        if not n_quantiles:
            n_quantiles = 4
        # Set editable value
        editable_graphs = False if not editable_graphs else True
        # Load df data
        df = pd.read_json(tv_input_json)
        # Set config      
        if (not pixel_height) or (not pixel_width):
            config = dict(
                doubleClick="autosize",
                displaylogo=False,
                editable=editable_graphs,
                modeBarButtonsToRemove=[
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
                modeBarButtonsToRemove=[
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
        
        # filter data by range
        if range_selection == "per_chrom":
            quantile_df = df[df["Chromosome"] == chromosome_selection]
            freqDF = tree_utils.calculate_topo_quantile_frequencies(quantile_df, current_topologies, additional_data, n_quantiles)
            xaxis_gridlines, yaxis_gridlines = tree_utils.get_gridline_bools(axis_gridlines)
            fig = tree_utils.plot_frequencies_topo_quantile(freqDF, template, color_mapping, axis_line_width, xaxis_gridlines, yaxis_gridlines, chromosome_selection, additional_data)
            autosome_data = freqDF.to_dict("records")
            autosome_columns = [
                {'id': 'TopologyID', 'name': 'TopologyID'},
                {'id': 'Frequency', 'name': 'Frequency', 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal_or_exponent)},
                {'id': 'Quantile', 'name': 'Quantile'},
            ]
            sex_data = []
            sex_columns=[]
            return fig, config, tree_utils.init_data_graph("plotly_dark"), dict(), autosome_data, autosome_columns, sex_data, sex_columns
        elif range_selection == "whole_genome":
            freqDF = tree_utils.calculate_topo_quantile_frequencies(df, current_topologies, additional_data, n_quantiles)
            xaxis_gridlines, yaxis_gridlines = tree_utils.get_gridline_bools(axis_gridlines)
            fig = tree_utils.plot_frequencies_topo_quantile(freqDF, template, color_mapping, axis_line_width, xaxis_gridlines, yaxis_gridlines, "Whole Genome", additional_data)
            autosome_data = freqDF.to_dict("records")
            autosome_columns = [
                {'id': 'TopologyID', 'name': 'TopologyID'},
                {'id': 'Frequency', 'name': 'Frequency', 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal_or_exponent)},
                {'id': 'Quantile', 'name': 'Quantile'},
            ]
            sex_data = []
            sex_columns=[]
            return fig, config, tree_utils.init_data_graph("plotly_dark"), dict(), autosome_data, autosome_columns, sex_data, sex_columns
        elif range_selection == "a_vs_x":
            if type(sex_chrom_selection) == str:
                sex_chrom_selection = [sex_chrom_selection]
            sex_chromosome_title = "+".join(sex_chrom_selection)
            autosomes_df = df[~df["Chromosome"].isin(sex_chrom_selection)]
            sex_chrom_df = df[df["Chromosome"].isin(sex_chrom_selection)]
            autosome_freqDF = tree_utils.calculate_topo_quantile_frequencies(autosomes_df, current_topologies, additional_data, n_quantiles)
            sex_freqDF = tree_utils.calculate_topo_quantile_frequencies(sex_chrom_df, current_topologies, additional_data, n_quantiles)
            xaxis_gridlines, yaxis_gridlines = tree_utils.get_gridline_bools(axis_gridlines)
            autosome_fig = tree_utils.plot_frequencies_topo_quantile(autosome_freqDF, template, color_mapping, axis_line_width, xaxis_gridlines, yaxis_gridlines, "Autosomes", additional_data)
            sex_chrom_fig = tree_utils.plot_frequencies_topo_quantile(sex_freqDF, template, color_mapping, axis_line_width, xaxis_gridlines, yaxis_gridlines, sex_chromosome_title, additional_data)
            autosome_data = autosome_freqDF.to_dict("records")
            autosome_columns = [
                {'id': 'TopologyID', 'name': 'TopologyID'},
                {'id': 'Frequency', 'name': 'Frequency', 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal_or_exponent)},
                {'id': 'Quantile', 'name': 'Quantile'},
            ]
            sex_data = sex_freqDF.to_dict("records")
            sex_columns = [
                {'id': 'TopologyID', 'name': 'TopologyID'},
                {'id': 'Frequency', 'name': 'Frequency', 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal_or_exponent)},
                {'id': 'Quantile', 'name': 'Quantile'},
            ]
            return autosome_fig, config, sex_chrom_fig, config, autosome_data, autosome_columns, sex_data, sex_columns
    else:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

# ---------------------- DataTable Callbacks ---------------------
@app.callback(
    [Output("basic-whole-genome-stats", "data"),
     Output("basic-whole-genome-stats", "columns"),
    ],
    [Input("input-data-upload", "children"),],
)
def whole_genome_data_table_stats(
    tv_input_json,
):
    if not tv_input_json:
        raise PreventUpdate
    else:
        tv_df = pd.read_json(tv_input_json)
        # Check if addtional data types available
        if len(tv_df.columns) < 4:
            columns = [{'id': "No Additional Data", 'name': "No Additional Data"}]
            data = []
            return data, columns
        else:
            data, columns = tree_utils.whole_genome_datatable(tv_df)
    
            return data, columns


@app.callback(
    [Output("basic-stats-pie-chart", "figure"),
     Output("basic-stats-pie-chart", "config"),
     Output("basic-alt-data-stats", "columns"),
     Output("basic-alt-data-stats", "data"),
     Output("basic-stats", "columns"),
     Output("basic-stats", "data"),
     Output("stats-topology-freq-pie-graph", "figure"),
     Output("stats-topology-freq-pie-graph", "config"),
    ],
    [Input("input-data-upload", "children"),
     Input("chromFile-data-upload", "children"),
     Input("current-data-range", "data"),
     Input("chromosome-options", "value"),
     Input("template-option", "value"),
     Input("topology-color-chart-store", "data"),
     Input("snapshot-file-option", "value"),
     Input("tv-pixel-height", "value"),
     Input("tv-pixel-width", "value"),
     Input("window-size", "children"),
     Input("font-family-option", "value"),
     Input("toggle-chrom-whole-genome", "value"),
    ],
     prevent_initial_call=True,
)
def current_view_summary_stats(
    tv_input_json,
    chromosome_length_json,
    dataRange,
    chromosome,
    template,
    color_mapping,
    snapshot_file_type,
    pixel_height,
    pixel_width,
    window_size,
    font_family,
    view_toggle,
):
    if not tv_input_json:
        raise PreventUpdate
    elif not chromosome_length_json:
        raise PreventUpdate
    chromosome_df = pd.read_json(chromosome_length_json)
    whole_tv_input_df = pd.read_json(tv_input_json)
    # tv_df = whole_tv_input_df
    dataMin, dataMax = dataRange
    # Filter data to current view only
    if view_toggle:
        tv_df = whole_tv_input_df[(whole_tv_input_df["Chromosome"] == chromosome)]
        current_view_tv_input_df = tv_df[(tv_df["Window"] >= dataMin) & (tv_df["Window"] <= dataMax)]
        # number_windows_in_chromosome = len(tv_df)
    else:
        current_view_tv_input_df = whole_tv_input_df
        # number_windows_in_chromosome = len(current_view_tv_input_df)
    # Generate dataframes of data
    basic_stats_topo_freqs, basic_stats_datatable = tree_utils.basic_stats_dfs(current_view_tv_input_df)
    # Set up data table columns + data
    if view_toggle:
        current_view_title = f"Current View - {chromosome}:{int(dataMin)}-{int(dataMax)}"
    else:
        current_view_title = "Whole genome"
        chromosome = "Whole genome"

    if basic_stats_datatable.empty:
        basic_alt_data_stats_columns = []
        basic_alt_data_stats_data = []
    else:
        basic_alt_data_stats_columns = [{"name": [current_view_title, i], "id": i, 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal)} for i in basic_stats_datatable.columns]
        basic_alt_data_stats_data = basic_stats_datatable.to_dict('records')
    
    # Generate general stats data
    perc_missing_data = round(((len(whole_tv_input_df[whole_tv_input_df["TopologyID"] == "NODATA"]) / (data_utils.roundUp(chromosome_df["End"].sum(), window_size) / window_size) * 100)), 2)
    genome_length = chromosome_df["End"].sum()
    avg_chrom_len = int(chromosome_df["End"].mean())
    basic_alt_data_stats_df = pd.DataFrame(
        {
            "General Statistics": [
                "Genome length",
                # "Average chromosome length",
                "Total topologies",
                "Percent missing data (*based on total NODATA entries*)",
                "Windows in genome",
                # f"Windows in chromosome {chromosome}",
                f"Windows in current view",
            ],
            "Value": [
                genome_length,
                # avg_chrom_len,
                len(whole_tv_input_df["TopologyID"].unique()),
                f"{perc_missing_data}%",
                len(whole_tv_input_df),
                # number_windows_in_chromosome,
                len(current_view_tv_input_df)],
        }
    )
    basic_stats_columns = [{"name": "General Statistics", "id": "General Statistics"}, {"name": "Value", "id": "Value", 'type': 'numeric', 'format': Format(group=True, groups=[3])}]
    basic_stats_data = basic_alt_data_stats_df.to_dict('records')

    # Create pie chart frequencies
    if (not pixel_height) or (not pixel_width):
        pixel_width = 750
        pixel_height = 400
    pie_fig_config = dict(
        displaylogo=False,
        toImageButtonOptions=dict(
            format=snapshot_file_type,
            filename=f"basic_stats_chart_{chromosome}_{int(dataMin)}_{int(dataMax)}",
            width=int(pixel_width),
            height=int(pixel_height),
        ),
    )
    pie_fig = tree_utils.current_view_topo_freq_chart(basic_stats_topo_freqs, template, color_mapping, current_view_title)
    # Create pie chart
    topo_freq_df = pd.DataFrame(whole_tv_input_df["TopologyID"].value_counts()/len(whole_tv_input_df))
    topo_freq_df = topo_freq_df.reset_index()
    topo_freq_df.columns = ['TopologyID', 'Frequency']
    topo_freq_df.sort_values(by="Frequency")
    wg_fig_pie = tree_utils.build_topology_frequency_pie_chart(topo_freq_df, template, color_mapping, font_family)
    return (pie_fig, pie_fig_config, basic_alt_data_stats_columns, basic_alt_data_stats_data, 
            basic_stats_columns, basic_stats_data, wg_fig_pie, pie_fig_config)


@app.callback(
    [Output("correlation-stats", "data"),
     Output("correlation-stats", "columns"),
     Output("correlation-stats-heatmap", "figure"),
     Output("correlation-test-output", "data"),
     Output("correlation-test-output", "columns"),
     Output("correlation-stats", "style_data_conditional"),
     Output("mean-freq-stats", "data"),
     Output("mean-freq-stats", "columns"),
     Output("correlation-topo-alert", "is_open"),
    ],
    [Input("correlation-run-btn", "n_clicks"),
     Input("input-data-upload", "children"),
     Input("current-data-range", "data"),
     Input("chromosome-options", "value"),
    ],
    [State("correlation-range-selection", "value"),
     State("correlation-add-data-type", "value"),
     State("correlation-test-type", "value"),
     State("correlation-posthoc-type", "value"),
     State("correlation-pval-adjustment-type", "value"),
     State("template-option", "value"),
     State("correlation-topologies", "value"),
     State("alpha-input", "value"),
    ],
    prevent_initial_call=True,
)
def exploreatory_stats_testing(
    btn,
    tv_input_json,
    dataRange,
    chromosome,
    range_selection,
    additional_data_type,
    test_type,
    posthoc_type,
    pval_adjustment,
    template,
    topologies,
    alpha,
):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if not tv_input_json:
        raise PreventUpdate
    elif not additional_data_type:
        raise PreventUpdate
    elif not topologies:
        raise PreventUpdate
    elif len(topologies) < 2:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, True
    elif button_id == 'correlation-run-btn':
        tv_df = pd.read_json(tv_input_json)
        # Select data range
        if range_selection == "Current Chromosome":
            tv_df = tv_df[(tv_df["Chromosome"] == chromosome) & (tv_df["TopologyID"].isin(topologies))]
        elif range_selection == "Current View":
            dataMin, dataMax = dataRange
            tv_df = tv_df[(tv_df["Chromosome"] == chromosome) & (tv_df["Window"] >= dataMin) & (tv_df["Window"] <= dataMax) & (tv_df["TopologyID"].isin(topologies))]
        else:
            tv_df = tv_df[tv_df["TopologyID"].isin(topologies)]
        # Generate average information for addtional data per topology
        mean_freq_dt_columns = [
            {'id': 'TopologyID', 'name': 'TopologyID'},
            {'id': 'Total Windows', 'name': 'Total Windows'},
            {'id': f'Mean ({additional_data_type})', 'name': f'Mean ({additional_data_type})', 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal_or_exponent)}
        ]
        mean_freq_dt_data = tree_utils.mean_frequency_of_alt_data_per_topology(tv_df, topologies, additional_data_type)
        # Run test
        if test_type == "Kruskal-Wallis H-test":
            posthoc, data, columns, H, p = tree_utils.kruskal_wallis_H_test(tv_df, additional_data_type, posthoc_type, pval_adjustment, alpha)
            heatmap = tree_utils.stats_test_heatmap(posthoc, template)
            test_output_columns = [
                {'id': 'Kruskal-Wallis H-test', 'name': 'Kruskal-Wallis H-test'},
                {'id': 'Value', 'name': 'Value', 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal_or_exponent)}
            ]
            test_output_data = pd.DataFrame({"Kruskal-Wallis H-test": ['H', 'p-value'], "Value": [H, p]}).to_dict('records')
            cell_conditional_green = [
                {
                    'if': {
                        'filter_query': f'{{p-value}} < {alpha}',
                        'column_id': 'p-value'
                    },
                    'backgroundColor': 'green'
                }
            ]
            cell_conditional_red = [
                {
                    'if': {
                        'filter_query': f'{{p-value}} > {alpha}',
                        'column_id': 'p-value'
                    },
                    'backgroundColor': 'red'
                }
            ]
            cell_conditional = cell_conditional_red + cell_conditional_green
            
            return data, columns, heatmap, test_output_data, test_output_columns, cell_conditional, mean_freq_dt_data, mean_freq_dt_columns, False
        elif test_type == "One-Way ANOVA":
            tv_df = tv_df.dropna()
            posthoc, data, columns, F, p = tree_utils.one_way_anova(tv_df, additional_data_type, posthoc_type, pval_adjustment, alpha)
            heatmap = tree_utils.stats_test_heatmap(posthoc, template)
            test_output_columns = [
                {'id': 'One-Way ANOVA', 'name': 'One-Way ANOVA'},
                {'id': 'Value', 'name': 'Value', 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal_or_exponent)}
            ]
            test_output_data = pd.DataFrame({"One-Way ANOVA": ['F', 'p-value'], "Value": [F, p]}).to_dict('records')
            cell_conditional_green = [
                {
                    'if': {
                        'filter_query': f'{{p-value}} < {alpha}',
                        'column_id': 'p-value'
                    },
                    'backgroundColor': 'green'
                }
            ]
            cell_conditional_red = [
                {
                    'if': {
                        'filter_query': f'{{p-value}} > {alpha}',
                        'column_id': 'p-value'
                    },
                    'backgroundColor': 'red'
                }
            ]
            cell_conditional = cell_conditional_red + cell_conditional_green
            return data, columns, heatmap, test_output_data, test_output_columns, cell_conditional, mean_freq_dt_data, mean_freq_dt_columns, False
    else:
        raise PreventUpdate 


@app.callback(
    [Output("data-freq-dist-row", "children"),
     Output("freq-dist-loading", "children")
    ],
    [Input("data-freq-dist-button", "n_clicks"),
     Input("input-data-upload", "children"),
     Input("editable-graph-option", "value"),
     Input("snapshot-file-option", "value"),
     Input("template-option", "value"),
     Input("font-family-option", "value"),
    ],
    [State("data-freq-dist-dropdown", "value"),
     State("freq-dist-loading", "children"),
    ],
    prevent_initial_call=True,
)
def additional_data_frequency_distributions(
    btn,
    tv_input_json,
    editable_graphs,
    snapshot_file_type,
    template,
    font_family,
    selected_data,
    loading,
):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if not tv_input_json:
        raise PreventUpdate
    elif button_id == 'data-freq-dist-button':
        tv_df = pd.read_json(tv_input_json)
        col_lists = []
        if type(selected_data) == str:
            selected_data = [selected_data]
        elif not selected_data:
            raise PreventUpdate
        for d in selected_data:
            data = tv_df[d]
            fig = tree_utils.frequency_distribution(data, d, template)
            col_lists.append(fig)
        fig_cols = []
        for fig in col_lists:
            # Append figure to list
            config = dict(
                doubleClick="autosize",
                displaylogo=False,
                editable=editable_graphs,
                modeBarButtonsToRemove=[
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
            fig_cols.append(
                dbc.Col([
                    dcc.Graph(figure=fig, config=config),
                ], width=4)
            )
        return fig_cols, None
    else:
        raise PreventUpdate 


# ---------------------- Export Option Callbacks ---------------------
@app.callback(
    [Output("tv-current-view-download", "data"),
     Output("gff-current-view-download", "data"),],
    [Input("curr-view-export-button", "n_clicks"),
     Input("input-data-upload", "children"),
     Input("chromosome-options", "value"),
     Input("gff-data-upload", "children"),],
    [State('current-data-range', 'data'),],
    prevent_initial_call=True,
)
def current_view_summary(
    n,
    tv_input_json,
    current_chromosome,
    gff_json,
    relayout_data,
):
    if not tv_input_json:
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
        df = pd.read_json(tv_input_json)
        df = df[(df["Window"] >= dataMin) & (df["Window"] <= dataMax) & (df["Chromosome"] == current_chromosome)]
        # Gather gff information into single DF
        if not gff_json:
            gffDFs = pd.DataFrame()
        else:
            # Concat multiple GFF df's into one + output
            if len(gff_json) == 1:
                gffDFs = pd.read_json(gff_json[0])
            else:
                gffDFs = pd.DataFrame()
                for j in gff_json:
                    gff_df = pd.read_json(j)
                    gffDFs = pd.concat([gffDFs, gff_df])
            gffDFs = gffDFs[(gffDFs["start"] >= dataMin) & (gffDFs["end"] <= dataMax) & (gffDFs["chromosome"] == current_chromosome)]
            gffDFs = gffDFs.sort_values(by=['chromosome', 'start', 'end'])
        # Return according to available data
        if (len(df) == 0) and (len(gffDFs) > 1):
            return None, dcc.send_data_frame(gffDFs.to_csv, f"current_view_{current_chromosome}_{int(dataMin)}_{int(dataMax)}.gff", sep="\t", index=False)
        elif (len(df) > 1) and (len(gffDFs) == 0):
            return dcc.send_data_frame(df.to_excel, f"current_view_{current_chromosome}_{int(dataMin)}_{int(dataMax)}.xlsx", index=False), None
        else:
            return dcc.send_data_frame(df.to_excel, f"current_view_{current_chromosome}_{int(dataMin)}_{int(dataMax)}.xlsx", index=False), dcc.send_data_frame(gffDFs.to_csv, f"current_view_{current_chromosome}_{int(dataMin)}_{int(dataMax)}.gff", sep="\t", index=False, header=False)
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
    tv_input_json,
    taxa_options,
    prune_taxa_choices,
    topobinner,
):
    if not tv_input_json:
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
        df = pd.read_json(tv_input_json)
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


@app.callback(
    Output("project-ini-download", "data"),
    [Input("project-ini-download-button", "n_clicks"),],
    prevent_initial_call=True,
)
def download_project_ini_file(
    btn,
):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    
    if button_id == "project-ini-download-button":
        content = tree_utils.project_ini_template()
        return dict(content=content, filename='project.ini')
    else:
        raise PreventUpdate


@app.callback(
    Output("tree-viewer-download", "data"),
    [Input("tree-viewer-download-button", "n_clicks"),],
    prevent_initial_call=True,
)
def download_tree_viewer_template_file(
    btn,
):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    
    if button_id == "tree-viewer-download-button":
        content = tree_utils.tree_viewer_template()
        return dcc.send_data_frame(content.to_csv, filename='TreeViewerInput.tsv', sep="\t", index=False)
    else:
        raise PreventUpdate


@app.callback(
    Output("chrom-length-download", "data"),
    [Input("chrom-length-download-button", "n_clicks"),],
    prevent_initial_call=True,
)
def download_chrom_length_template_file(
    btn,
):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    
    if button_id == "chrom-length-download-button":
        content = tree_utils.chrom_len_template()
        return dcc.send_data_frame(content.to_csv, filename='ChromosomeLengths.bed', sep="\t", index=False)
    else:
        raise PreventUpdate


