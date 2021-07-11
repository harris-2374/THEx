"""
Author: Andrew Harris
Python Version: Python3.8.3
"""
import os

from gevent import monkey
monkey.patch_all()

import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from gevent.pywsgi import WSGIServer
from thex.app import app, server
from thex.apps import homepage
from thex.apps.tree_viewer import tv_layout
from thex.apps.signal_tracer import signalTracer_layout
from thex.apps.docs import docs_layout

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])

@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/':
        return homepage.layout
    elif pathname == '/apps/signal_tracer':
        return signalTracer_layout()
    elif pathname == '/apps/tree_viewer':
        return tv_layout()
    elif pathname == '/apps/documentation':
        return docs_layout()
    else:
        return '404 - Page not found'


def main():
    app.run_server(debug=True)

    # print("Tree House Explorer running on http://127.0.0.1:8050/")
    # http_server = WSGIServer(('127.0.0.1', 8050), application=server, log=None)
    # http_server.serve_forever()
