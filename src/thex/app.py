import dash
import dash_bootstrap_components as dbc
from flask import cli
cli.show_server_banner = lambda *_: None

app = dash.Dash(
    __name__,
    external_stylesheets = [dbc.themes.DARKLY],
    suppress_callback_exceptions = True,
    serve_locally = True,
)

app.title = "THEx"
server = app.server

