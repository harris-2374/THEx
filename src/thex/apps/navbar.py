import dash_bootstrap_components as dbc
import dash_html_components as html

def Navbar(current_program=None):
    if not current_program:
        brand_title = "Tree House Explorer"
    else:
        brand_title = f"Tree House Explorer: {current_program}"
    navbar = dbc.NavbarSimple(
        children=[
            html.Div(
                children=[],
                style={'background-color': 'black'}
            ),
            html.Div(
                children=[
                    dbc.DropdownMenu(
                        nav=True,
                        in_navbar=True,
                        label="Program Toggle",
                        color="warning",
                        right=True,
                        children=[
                        dbc.DropdownMenuItem(
                            "Tree Viewer",
                            href="/apps/tree_viewer"
                        ),
                        dbc.DropdownMenuItem(
                            "p-Distance Tracer",
                            href="/apps/p_distance_tracer",
                        ),
                        dbc.DropdownMenuItem(
                            "Data Preparation",
                            href="/apps/data_prep"
                        ),
                        ],
                        style={
                            "color":'black',
                            "font-size": "20px",
                            "border": "2px black solid",
                            "border-radius": "5px",
                        },
                    ),
                ],
            ),
        ],
        brand=brand_title,
        brand_href="/",
        brand_style={
            "color":'black',
            "font-size": "35px",
            "border":"2px black solid",
            "border-radius": "5px",
            "padding-left": "10px",
            "padding-right": "10px",
            "padding-top": "0px",
            "padding-bottom": "0px",
        },
        sticky="top",
        fluid=True,
        color="warning"
    )
    return navbar
    