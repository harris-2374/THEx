import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

from thex.app import app

layout = dbc.Container(
    children=[
        # Header
        dbc.Row(
            children=[
                dcc.Markdown("""# ***Welcome to Tree House Explorer!*** """, className="homepage-title")
            ],
            id="homepage-row",
            no_gutters=True,
            key='updates-row',
            className="homepage-title-row",
            justify='center',
            align='center',
        ),
        # Ideogram
        dbc.Row(
            dbc.Col([
                html.Div([dbc.CardImg(src=app.get_asset_url("signal_homepage_tile.png"), top=True, style={"height": "100px", "border-radius": "50px", "border": "solid 2px orange"})]),
            ], width=12),
            align='center',
            no_gutters=True,
            style={"margin-top": "5px", "padding": "10px"},
        ),
        # Body
        dbc.Row(
            children=[
                ### Tree Viewer
                dbc.Col(
                    children=[
                        # Header
                        dbc.Col(
                            dbc.Card(
                                children=[
                                    dbc.CardBody(
                                        [html.H4("Dashboards")], 
                                        style={
                                            "text-align": "center",
                                            "text-decoration": "underline",
                                            "border-radius": "3px",
                                        }
                                    ),
                                    dbc.CardFooter(
                                        [
                                            "If you are new to THEx, we reccomend starting with the Docs!"
                                        ],
                                        style={"text-align": "center"}
                                    ),
                                ],
                                color="warning",
                                outline=True,
                            ),
                            align='stretch',
                            style={"margin-bottom": "5px"},
                        ),
                        # Tree Viewer
                        dbc.Col(
                            dbc.Card(
                                children=[
                                    html.Div([dbc.CardImg(src=app.get_asset_url("tree5.svg"), top=True, className="homepage-card-image")], style={"padding": "5px"}),
                                    dbc.CardBody([dbc.Button("Tree Viewer", block=True, className="homepage-buttons", color="info", outline=True, href="/apps/tree_viewer"),]),
                                ],
                                color="warning",
                                outline=True,
                            ),
                            align='stretch',
                            style={"margin-bottom": "5px"},
                        ),
                    ],
                    align="start",
                    width=2,
                ),
                ### Middle info pane
                dbc.Col(
                    children=[
                        dbc.Jumbotron(
                            children=[
                                html.Hr(style={"background-color": "white"}),
                                dcc.Markdown("""## _- Introduction -_ """),
                                html.P(
                                    """
                                    Genomes are mosaics of evolutionary histories that reflect ancient signatures of species divergence, 
                                    as well as incomplete lineage sorting (ILS) or gene flow. In recent years, rapid advancements in massively 
                                    parallel sequencing technologies have enabled researchers to collect large volumes of genome-scale phylogenetic 
                                    data for many species. Understanding how and why phylogenetic signal varies across species genomes can yield 
                                    powerful insights into evolutionary histories and adaptive evolution. By integrating diverse data types with 
                                    local genealogies, users can differentiate genetic variation that is consistent with the species tree from that 
                                    stemming from natural selection, ILS, or gene flow. While several tools have been developed for the independent 
                                    analysis of genomic data (e.g., IGV, UCSC Genome Browser, etc.) and phylogenetic tree visualization (Iroki, FigTree, etc.), 
                                    a tool that can simultaneously analyze phylogenetic signal variation with other chromosomal and gene-based annotations 
                                    has yet to be developed for the field of phylogenomics. Tree House Explorer (THEx) is a novel genome browser that allows 
                                    users to integrate phylogenomic data and genomic annotations into a single interactive platform for combined analysis. 
                                    THEx allows users to visualize genome-wide variation in evolutionary histories as well as genetic divergence on a 
                                    chromosome-by-chromosome basis, with continuous sliding window comparisons to gene annotations, recombination rates, 
                                    GC-content, and other user-specified, highly customizable feature annotations. THEx provides a new resource for 
                                    interactive data visualization in phylogenomics and a novel approach to analyze and interpret the diverse evolutionary 
                                    histories woven throughout genomes.
                                    """
                                ),
                                html.Hr(style={"background-color": "white"}),
                                dcc.Markdown("""## _- Dashboards -_ """),
                                dcc.Markdown(
                                    """
                                    Tree House Explorer offers a growing collection of dashboards that provide a new way to visualize a multitude of
                                    phylogenomic and genomic data concurrently. Each dashboard provides highly interactive graphs that allow for smooth surfing
                                    of data at a chromosome and genome-wide level. In addition to being interactive, each dashboard provides a growing set of 
                                    graph themes, colors, and other attributes that brings your figures to as close to publication ready as possible.  

                                    Tree Viewer:  
                                    * Windowed based approach to visualizing phylogenetic signal across the genome with the added benefit of being able to
                                    visualize additional data types like recombination rate, GC-content, and gene annotations concurrently.  

                                    Signal Tracer:  
                                    * Visualizer for window-based calculations from multi-alignment sequences. Signal Tracer provides an interactive way to 
                                    dissect and extract local information with ease. 
                                    """
                                ),
                                html.Hr(style={"background-color": "white"}),
                            ],
                            className="homepage-jumbotron",
                            key="homepage-jumbotron",
                        ),
                    ],
                    width=8,
                    align='start',
                ),
                ### Signal Tracer
                dbc.Col(
                    children=[
                        # Header
                        dbc.Col(
                            dbc.Card(
                                children=[
                                    dbc.CardBody(
                                        [html.H4("Dashboards")], 
                                        style={
                                            "text-align": "center",
                                            "text-decoration": "underline",
                                            "border-radius": "3px",
                                        }
                                    ),
                                    dbc.CardFooter(
                                        [
                                            "Check out THEx's GitHub page for example data sets and updates!"
                                        ],
                                        style={"text-align": "center"}
                                    ),
                                ],
                                color="warning",
                                outline=True,
                            ),
                            align='stretch',
                            style={"margin-bottom": "5px"},
                        ),
                        # Signal Tracer
                        dbc.Col(
                            dbc.Card(
                                children=[
                                    html.Div([dbc.CardImg(src=app.get_asset_url("p-distance-logo-V7-71pt.svg"), top=True, className="homepage-card-image")], style={"padding": "5px"}),
                                    dbc.CardBody(
                                        children=[
                                            dbc.Button("Signal Tracer", block=True, className="homepage-buttons", color="info", outline=True, href="/apps/signal_tracer"),
                                        ],
                                    ),
                                ],
                                color="warning",
                                outline=True,
                            ),
                            align='stretch',
                            style={"margin-bottom": "5px"},
                        ),
                    ],
                    align="start",
                    width=2,
                ),
            ],
            justify="center",
            align="center",
            id="homepage-row",
            no_gutters=True,
            key='updates-row',
        ),
        # Footer
        dbc.Row(
            dbc.NavbarSimple(
                [dbc.Button(
                    children=[
                        html.Img(src=app.get_asset_url("docs4.svg"), height="50px"),
                        html.Div("Docs")
                    ],
                    className="logo-png",
                    href="/apps/documentation",
                    outline=True,
                ),
                dbc.Button(
                    children=[
                        html.Img(src=app.get_asset_url("GitHub-Mark-64px.png"), height="50px"),
                        html.Div("GitHub")
                    ],
                    className="logo-png",
                    href="https://github.com/harris-2374/THEx",
                    target="_blank",
                    outline=True,
                ),],
                fixed='bottom',
                fluid=True,
                color='dark',
            ),
        ),
    ],
    fluid=True,
    key='home_container',
)


