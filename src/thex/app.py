import dash
import dash_bootstrap_components as dbc
from flask import cli
cli.show_server_banner = lambda *_: None


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(
    __name__,
    external_stylesheets = [dbc.themes.DARKLY],
    # external_stylesheets = [dbc.themes.SUPERHERO],
    suppress_callback_exceptions = True,
    serve_locally = True,
)

app.title = "THEx"
server = app.server

