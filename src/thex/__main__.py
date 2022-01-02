"""
Author: Andrew Harris
Python Version: Python3.8.3
"""
import argparse
from pathlib import Path

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
    parser = argparse.ArgumentParser(description='Tree House Explorer Genome Browser')
    parser.add_argument(
        '--host',
        type=str,
        action='store',
        default='127.0.0.1',
        help="Host address",
    )
    parser.add_argument(
        '--port',
        type=int,
        action='store',
        default=8050,
        help='Port number',
    )
    parser.add_argument(
        '--dev',
        action='store_true',
        default=False,
        help='Run in development mode',
    )
    args = parser.parse_args()
    if args.dev:
        app.run_server(debug=True, port=args.port, host=args.host)
    else:
        print(f"Tree House Explorer running on http://{args.host}:{args.port}/")
        http_server = WSGIServer((args.host, args.port), application=server, log=None)
        http_server.serve_forever()
