import math
import itertools

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from Bio import Phylo
from ete3 import Tree
from plotly.subplots import make_subplots

from thex.apps.utils import data_utils

# -------------------------------------------------------------------------------------
# --------------------------------------- Classes -------------------------------------
class DrawTree():
    def __init__(self, newicktree, template, topology, color_map):
        self.newicktree = Phylo.read(newicktree, "newick")
        self.template = template
        self.topology = topology
        self.color_map = color_map

    def create_square_tree(self):

        def get_x_coordinates(tree):
            """Associates to each clade an x-coord.
            returns dict {clade: x-coord}
            """
            xcoords = tree.depths(unit_branch_lengths=True)
            # xcoords = tree.depths(unit_branch_lengths=False)

            # first_node = [k for k in xcoords][0]
            # xcoords[first_node] = -0.01

            # tree.depth() maps tree clades to depths (by branch length).
            # returns a dict {clade: depth} where clade runs over all Clade instances of the tree, and depth is the distance from root to clade

            #  If there are no branch lengths, assign unit branch lengths
            if not max(xcoords.values()):
                xcoords = tree.depths(unit_branch_lengths=True)
            return xcoords

        def get_y_coordinates(tree, dist=1):
            """
            returns  dict {clade: y-coord}
            The y-coordinates are  (float) multiple of integers (i*dist below)
            dist depends on the number of tree leafs
            """
            maxheight = tree.count_terminals()  # Counts the number of tree leafs.
            # Rows are defined by the tips/leafs
            ycoords = dict(
                (leaf, maxheight - i * dist)
                for i, leaf in enumerate(reversed(tree.get_terminals()))
            )
            def calc_row(clade):
                for subclade in clade:
                    if subclade not in ycoords:
                        calc_row(subclade)
                ycoords[clade] = (ycoords[clade.clades[0]] +
                                  ycoords[clade.clades[-1]]) / 2

            if tree.root.clades:
                calc_row(tree.root)
            return ycoords

        def get_clade_lines(
            orientation="horizontal",
            y_curr=0,
            x_start=0,
            x_curr=0,
            y_bot=0,
            y_top=0,
            line_color="white",
            line_width=2,
            root_clade = False
        ):
            """define a shape of type 'line', for branch
            """
            branch_line = dict(
                type="line", layer="below", line=dict(color=line_color, width=line_width)
            )

            if root_clade:
                branch_line.update(x0=-0.01, y0=y_curr, x1=-0.01, y1=y_curr)
                return branch_line
            elif orientation == "horizontal":
                branch_line.update(x0=x_start, y0=y_curr, x1=x_curr, y1=y_curr)
            elif orientation == "vertical":
                branch_line.update(x0=x_curr, y0=y_bot, x1=x_curr, y1=y_top)
            else:
                raise ValueError("Line type can be 'horizontal' or 'vertical'")
            return branch_line

        def draw_clade(
            clade,
            x_start,
            line_shapes,
            line_color="white",
            line_width=2,
            x_coords=0,
            y_coords=0,
            init_clade=False,
        ):
            """Recursively draw the tree branches, down from the given clade"""

            x_curr = x_coords[clade]
            y_curr = y_coords[clade]

            # Draw a horizontal line from start to here
            if init_clade:
                branch_line = get_clade_lines(
                    orientation="horizontal",
                    y_curr=y_curr,
                    x_start=x_start,
                    x_curr=x_curr,
                    line_color=line_color,
                    line_width=line_width,
                    root_clade=True,
                )
            else:
                branch_line = get_clade_lines(
                    orientation="horizontal",
                    y_curr=y_curr,
                    x_start=x_start,
                    x_curr=x_curr,
                    line_color=line_color,
                    line_width=line_width,
                    root_clade=False,
                )

            line_shapes.append(branch_line)

            if clade.clades:
                # Draw a vertical line connecting all children
                y_top = y_coords[clade.clades[0]]
                y_bot = y_coords[clade.clades[-1]]

                line_shapes.append(
                    get_clade_lines(
                        orientation="vertical",
                        x_curr=x_curr,
                        y_bot=y_bot,
                        y_top=y_top,
                        line_color=line_color,
                        line_width=line_width,
                    )
                )

                # Draw descendants
                for child in clade:
                    draw_clade(child, x_curr, line_shapes,
                               x_coords=x_coords, y_coords=y_coords, 
                               line_color=line_color)


        if 'dark' in self.template:
            text_color = 'white'
        else:
            text_color = 'black'
        line_color = self.color_map[self.topology]
        tree = self.newicktree
        tree.ladderize()

        x_coords = get_x_coordinates(tree)
        y_coords = get_y_coordinates(tree)
        line_shapes = []
        
        draw_clade(
            tree.root,
            0,
            line_shapes,
            line_color=line_color,
            line_width=2,
            x_coords=x_coords,
            y_coords=y_coords,
            init_clade=True,
        )
        my_tree_clades = x_coords.keys()
        X = []
        Y = []
        text = []

        for cl in my_tree_clades:
            X.append(x_coords[cl])
            Y.append(y_coords[cl])
            # Add confidence values if internal node
            if not cl.name:
                if not cl.name:
                    text.append(" ")
                else:
                    text.append(cl.name)
            else:
                text.append(cl.name)

        axis = dict(
            showline=False,
            visible=False,
            zeroline=False,
            showgrid=False,
            showticklabels=False,
            title="",  # y title
        )

        label_legend = ["Tree_1"]
        nodes = []

        for elt in label_legend:
            node = dict(
                type="scatter",
                x=X,
                y=Y,
                mode="markers+text",
                marker=dict(color=text_color, size=5),
                text=text,  # vignet information of each node
                textposition='middle right',
                textfont=dict(color=text_color, size=15),
                showlegend=False,
                name=elt,
            )
            nodes.append(node)
        # Set graph x-range
        if max(x_coords.values()) < 0.1:
            x_range = [-0.1, (max(x_coords.values())+(max(x_coords.values())*1.25))]
            show_xaxis = True
        elif max(x_coords.values()) < 0.5:
            x_range = [-0.1, 0.5]
            show_xaxis = True
        elif max(x_coords.values()) < 1:
            x_range = [-0.1, 1]
            show_xaxis = True
        elif max(x_coords.values()) == 1:
            x_range = [-0.1, max(x_coords.values())+2]
            show_xaxis = False
        else:
            x_range = [-0.1, max(x_coords.values())+2]
            show_xaxis = False

        layout = dict(
            autosize=True,
            template=self.template,
            dragmode="pan",
            margin=dict(t=20, b=10, r=20, l=10),
            showlegend=False,
            xaxis=dict(
                showline=True,
                zeroline=False,
                visible=show_xaxis,
                showgrid=False,
                showticklabels=True,
                range=x_range,
            ),
            yaxis=axis,
            hovermode="closest",
            shapes=line_shapes,
        )

        fig = go.Figure(data=nodes, layout=layout)
        return fig

    def create_angular_tree(self):

        def get_x_coordinates(tree):
            """Associates to each clade an x-coord.
            returns dict {clade: x-coord}
            """
            # xcoords = tree.depths(unit_branch_lengths=True)
            # print("===========================")
            # nodes = [n for n in tree.find_clades()]
            # nodes = tree.get_terminals() + tree.get_nonterminals()
            # print(tree.root.clades)
            # root_xcoord = {tree.root.clades[1]:0}
            terminal_nodes = tree.get_terminals()
            internal_nodes = tree.get_nonterminals()
            terminal_xcoords = dict((leaf, i) for i, leaf in enumerate(terminal_nodes))
            internal_xcoords = dict(
                (leaf, i+0.5) for leaf, i in zip(internal_nodes, range(1, len(internal_nodes)))
            )

            xcoords = {**terminal_xcoords, **internal_xcoords}


            # print(xcoords)
            # print("===========================")
            # tree.depth() maps tree clades to depths (by branch length).
            # returns a dict {clade: depth} where clade runs over all Clade instances of the tree, and depth
            # is the distance from root to clade

            #  If there are no branch lengths, assign unit branch lengths
            if not max(xcoords.values()):
                xcoords = tree.depths(unit_branch_lengths=True)

            return xcoords

        def get_y_coordinates(tree, dist=1):
            """
            returns  dict {clade: y-coord}
            The y-coordinates are  (float) multiple of integers (i*dist below)
            dist depends on the number of tree leafs
            """
            maxheight = tree.count_terminals()  # Counts the number of tree leafs.
            # Rows are defined by the tips/leafs
            # root_ycoord = {tree.root:maxheight}
            terminal_nodes = tree.get_terminals()
            internal_nodes = tree.get_nonterminals()
            terminal_ycoords = dict((leaf, 1) for _, leaf in enumerate(terminal_nodes))
            internal_ycoords = dict(
                (leaf, i) for leaf, i in zip(internal_nodes, reversed(range(1, len(internal_nodes))))
            )
            ycoords = {**terminal_ycoords, **internal_ycoords}
            # print("***********************")
            # print(ycoords)
            # print("***********************")

            def calc_row(clade):
                for subclade in clade:
                    if subclade not in ycoords:
                        calc_row(subclade)
                ycoords[clade] = (ycoords[clade.clades[0]] +
                                  ycoords[clade.clades[-1]]) / 2

            if tree.root.clades:
                calc_row(tree.root)
            return ycoords

        def get_clade_lines(
            orientation="horizontal",
            y_curr=0,
            last_y_curr=0,
            x_start=0,
            x_curr=0,
            y_bot=0,
            y_top=0,
            line_color="rgb(25,25,25)",
            line_width=0.5,
            init_flag=False,
        ):
            """define a shape of type 'line', for branch
            """
            branch_line = dict(
                type="line", layer="below", line=dict(color=line_color, width=line_width)
            )
            if orientation == "horizontal":
                if init_flag:
                    branch_line.update(x0=x_start, y0=y_curr,
                                       x1=x_curr, y1=y_curr)
                else:
                    branch_line.update(
                        x0=x_start, y0=last_y_curr, x1=x_curr, y1=y_curr)
            elif orientation == "vertical":
                branch_line.update(x0=x_curr, y0=y_bot, x1=x_curr, y1=y_top)
            else:
                raise ValueError("Line type can be 'horizontal' or 'vertical'")
            return branch_line

        def draw_clade(
            clade,
            x_start,
            line_shapes,
            line_color="rgb(15,15,15)",
            line_width=1,
            x_coords=0,
            y_coords=0,
            last_clade_y_coord=0,
            init_flag=True
        ):
            """Recursively draw the tree branches, down from the given clade"""
            x_curr = x_coords[clade]
            y_curr = y_coords[clade]

            # Draw a horizontal line from start to here
            branch_line = get_clade_lines(
                orientation="horizontal",
                y_curr=y_curr,
                last_y_curr=last_clade_y_coord,
                x_start=x_start,
                x_curr=x_curr,
                line_color=line_color,
                line_width=line_width,
                init_flag=init_flag,
            )

            line_shapes.append(branch_line)

            if clade.clades:
                # Draw descendants
                for child in clade:
                    draw_clade(child, x_curr, line_shapes, x_coords=x_coords,
                               y_coords=y_coords, last_clade_y_coord=y_coords[clade], 
                               init_flag=False, line_color=line_color)
        if 'dark' in self.template:
            text_color = 'white'
        else:
            text_color = 'black'
        line_color = self.color_map[self.topology]
        # Load in Tree object and ladderize
        tree = self.newicktree
        tree.ladderize()

        # Get coordinates + put into dictionary
        # dict(keys=clade_names, values=)
        x_coords = get_x_coordinates(tree)
        y_coords = get_y_coordinates(tree)
        # y_coords = get_x_coordinates(tree)
        # x_coords = get_y_coordinates(tree)
        line_shapes = []

        draw_clade(
            tree.root,
            0,
            line_shapes,
            line_color=line_color,
            line_width=2,
            x_coords=x_coords,
            y_coords=y_coords,
        )
        #
        my_tree_clades = x_coords.keys()
        X = []
        Y = []
        text = []

        for cl in my_tree_clades:
            X.append(x_coords[cl])
            Y.append(y_coords[cl])
            # Add confidence values if internal node
            if not cl.name:
                text.append(cl.confidence)
            else:
                text.append(cl.name)

        axis = dict(
            showline=False,
            zeroline=False,
            showgrid=False,
            visible=False,
            showticklabels=False,
        )

        label_legend = ["Tree_1"]
        nodes = []

        for elt in label_legend:
            node = dict(
                type="scatter",
                x=X,
                y=Y,
                mode="markers+text",
                marker=dict(color=text_color, size=5),
                text=text,  # vignet information of each node
                textposition='right',
                textfont=dict(color=text_color, size=25),
                showlegend=False,
                name=elt,
            )
            nodes.append(node)

        layout = dict(
            template=self.template,
            dragmode="select",
            autosize=True,
            showlegend=True,
            xaxis=dict(
                showline=True,
                zeroline=False,
                visible=False,
                showgrid=False,
                showticklabels=True,
                range=[-1, (max(x_coords.values())+2)]
            ),
            yaxis=axis,
            hovermode="closest",
            shapes=line_shapes,
            legend={"x": 0, "y": 1},
            font=dict(family="Open Sans"),
        )

        fig = dict(data=nodes, layout=layout)
        return fig


class RFDistance():

    def __init__(self, t1, t2):
        self.t1 = Tree(t1)
        self.t2 = Tree(t2)
        self.compare = self.t1.compare(self.t2)

    def NormRF(self):
        return self.compare['norm_rf']

    def RF(self):
        return self.compare['rf']

    def MaxRF(self):
        return self.compare['max_rf']


# -------------------------------------------------------------------------------------
# ------------------------------ Alt Data Graph Functions -----------------------------

def make_alt_data_str_figure(
    alt_data_to_graph,
    chromosome_df,
    color_mapping,
    topology_df,
    window_size,
    template,
    dataRange,
    axis_line_width,
    xaxis_gridlines,
    yaxis_gridlines,
):
    # sort dataframe
    topology_df.sort_values(by=["Window"], inplace=True)
    topology_df.fillna("NULL", inplace=True)

    # Build graph
    fig = px.histogram(
        topology_df,
        x="Window",
        y=[1]*len(topology_df),
        color=alt_data_to_graph,
        color_discrete_sequence=list(color_mapping.values()),
        nbins=int(chromosome_df["End"].max()/window_size),
    )
    # Update layout
    fig.update_layout(
        template=template,
        margin=dict(
            l=60,
            r=50,
            b=40,
            t=40,
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0
        ),
        title={
            'text': str(alt_data_to_graph),
            'y':0.95,
            'x':0.501,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        hovermode="x unified",
    )
    fig.update_xaxes(
        title="Position",
        range=dataRange,
        showline=True,
        showgrid=xaxis_gridlines,
        linewidth=axis_line_width,
    )
    fig.update_yaxes(
        title="y-axis",
        range=[0, 1],
        nticks=1,
        showline=True,
        showgrid=yaxis_gridlines,
        linewidth=axis_line_width,
    )
    return fig


def make_alt_data_int_figure(
    alt_data_to_graph,
    color_mapping,
    topology_df,
    template,
    dataRange,
    y_max,
    axis_line_width,
    xaxis_gridlines,
    yaxis_gridlines,
):
    # sort dataframe
    topology_df = topology_df.sort_values(by=["Window"])
    y_range = [0, (y_max*1.1)]
    # Build graph
    fig = px.line(
        topology_df,
        x="Window",
        y=alt_data_to_graph,
        color_discrete_sequence=list(color_mapping.values()),
    )
    # Update layout
    fig.update_layout(
        template=template,
        margin=dict(
            l=60,
            r=50,
            b=40,
            t=40,
        ),
        title={
            'text': str(alt_data_to_graph),
            'y':0.95,
            'x':0.501,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        hovermode="x unified",
    )
    fig.update_xaxes(
        title="Position",
        range=dataRange,
        showline=True,
        showgrid=xaxis_gridlines,
        linewidth=axis_line_width,
    )
    fig.update_yaxes(
        fixedrange=True,
        linewidth=axis_line_width,
        range=y_range,
        showgrid=yaxis_gridlines,
        showline=True,
        title="Edit me",
    )
    return fig


# ----------------------------------------------------------------------------------------
# ----------------------------- Distribution Graph Functions -----------------------------
def build_histogram_with_rug_plot(
    topology_df,
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
):
    # --- Set up topology data ---
    # Extract current topology data
    if (type(current_topologies) == str) or (type(current_topologies) == int):
        wanted_rows = topology_df[topology_df["TopologyID"] == current_topologies]
    elif type(current_topologies) == list:
        wanted_rows = topology_df[topology_df["TopologyID"].isin(current_topologies)]

    # Add in psuedodata for missing current_topologies (fixes issue where topology is dropped from legend)
    if len(wanted_rows['TopologyID'].unique()) < len(current_topologies):
        missing_topologies = [t for t in current_topologies if t not in wanted_rows['TopologyID'].unique()]
        for mt in missing_topologies:
            missing_row_data = [chromosome, 0, 'NA', mt] + ['NULL']*(len(wanted_rows.columns)-4)
            missing_row = pd.DataFrame(data={i:j for i,j in zip(wanted_rows.columns, missing_row_data)}, index=[0])
            wanted_rows = pd.concat([wanted_rows, missing_row])


    # Group data by topology ID
    grouped_topology_df = wanted_rows.groupby(by='TopologyID')

    # Set row heights based on number of current_topologies being shown
    if len(current_topologies) <= 6:
        subplot_row_heights = [1, 1]
    elif len(current_topologies) <= 8:
        subplot_row_heights = [4, 2]
    else:
        subplot_row_heights = [8, 2]
    # Build figure
    fig = make_subplots(rows=2, cols=1, row_heights=subplot_row_heights, vertical_spacing=0.05, shared_xaxes=True)
    for topology, data in grouped_topology_df:
        # fig.add_trace(
        #     go.Scattergl(
        #         x=data['Window'],
        #         y=data['TopologyID'],
        #         name=topology,
        #         legendgroup=topology,
        #         mode='markers',
        #         marker_symbol='line-ns-open',
        #         marker_size=int(50/len(grouped_topology_df)),
        #         marker_line_width=1,
        #         marker_color=[color_mapping[topology]]*len(data),
        #     ),
        #     row=1, col=1,
        # )
        fig.add_trace(
            go.Scatter(
                x=data['Window'],
                y=data['TopologyID'],
                name=topology,
                legendgroup=topology,
                mode='markers',
                marker_symbol='line-ns-open',
                marker_size=int(50/len(grouped_topology_df)),
                marker_line_width=1,
                marker_color=[color_mapping[topology]]*len(data),
            ),
            row=1, col=1,
        )
        fig.add_trace(
            go.Histogram(
                x=data['Window'],
                y=[1]*len(data),
                histfunc="min",
                nbinsx=int(chromosome_df["End"].max() / window_size),
                name=topology,
                legendgroup=topology,
                showlegend=False,
                marker_color=color_mapping[topology],
            ),
            row=2, col=1
        )
    # Update layout + axes
    fig.update_layout(
        template=template,
        legend_title_text='Topology',
        margin=dict(
            l=60,
            r=50,
            b=40,
            t=40,
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            traceorder='normal',
        ),
        hovermode="x unified",
    )
    fig.update_xaxes(
        rangemode="tozero",
        range=dataRange,
        linewidth=axis_line_width,
        showgrid=xaxis_gridlines,
        row=1,
        col=1
    )
    fig.update_xaxes(
        rangemode="tozero",
        range=dataRange,
        linewidth=axis_line_width,
        title='Position',
        showgrid=xaxis_gridlines,
        row=2,
        col=1,
    )
    fig.update_yaxes(
        rangemode="tozero",
        categoryarray=topoOrder,
        linewidth=axis_line_width,
        showgrid=yaxis_gridlines,
        showticklabels=False,
        fixedrange=True,
        ticklen=0,
        title="",
        type='category',
        row=1,
        col=1,
    )
    fig.update_yaxes(
        rangemode="tozero",
        fixedrange=True,
        linewidth=axis_line_width,
        nticks=1,
        showgrid=yaxis_gridlines,
        showticklabels=False,
        ticklen=0,
        title="",
        row=2,
        col=1,
    )
    return fig


def build_rug_plot(
    topology_df,
    chromosome,
    template,
    current_topologies,
    color_mapping,
    dataRange,
    topoOrder,
    axis_line_width,
    xaxis_gridlines,
    yaxis_gridlines,
):
    # --- Group wanted data ---
    if (type(current_topologies) == str) or (type(current_topologies) == int):
        wanted_rows = topology_df[topology_df["TopologyID"] == current_topologies]
    elif type(current_topologies) == list:
        wanted_rows = topology_df[topology_df["TopologyID"].isin(current_topologies)]
    # Add in psuedodata for missing current_topologies (fixes issue where topology is dropped from legend)
    if len(wanted_rows['TopologyID'].unique()) < len(current_topologies):
        missing_topologies = [t for t in current_topologies if t not in wanted_rows['TopologyID'].unique()]
        for mt in missing_topologies:
            missing_row_data = [chromosome, 0, 'NA', mt] + ['NULL']*(len(wanted_rows.columns)-4)
            missing_row = pd.DataFrame(data={i:j for i,j in zip(wanted_rows.columns, missing_row_data)}, index=[0])
            wanted_rows = pd.concat([wanted_rows, missing_row])
    else:
        pass
    # --- Group data by topology ID
    grouped_topology_df = wanted_rows.groupby(by='TopologyID')
    # --- Build figure ---
    fig = go.Figure()
    for topology, data in grouped_topology_df:
        fig.add_trace(go.Scatter(
            x=data['Window'],
            y=data['TopologyID'],
            name=topology,
            legendgroup=topology,
            mode='markers',
            marker_symbol='line-ns-open',
            marker_size=int(100/len(grouped_topology_df)),
            marker_line_width=1,
            marker_color=[color_mapping[topology]]*len(data),
        ))
    # Update figure layout + axes
    fig.update_layout(
        template=template,
        legend_title_text='Topology',
        xaxis_title_text='Position',
        margin=dict(
            l=60,
            r=60,
            b=40,
            t=40,
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            traceorder='normal',
        ),
        hovermode="x unified",
    )
    fig.update_xaxes(
        fixedrange=True,
        rangemode="tozero",
        range=dataRange,
        linewidth=axis_line_width,
        showgrid=xaxis_gridlines,
        showline=True,
    )
    fig.update_yaxes(
        fixedrange=True,
        title="",
        showline=True,
        showgrid=yaxis_gridlines,
        linewidth=axis_line_width,
        showticklabels=False,
        type='category',
        categoryarray=topoOrder,
    )
    fig.for_each_annotation(lambda a: a.update(text=""))
    return fig


def build_tile_plot(
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
):
    # Extract current topology data
    if (type(current_topologies) == str) or (type(current_topologies) == int):
        wanted_rows = topology_df_filtered[topology_df_filtered["TopologyID"] == current_topologies]
    elif type(current_topologies) == list:
        wanted_rows = topology_df_filtered[topology_df_filtered["TopologyID"].isin(current_topologies)]

    # fig = px.histogram(
    #     wanted_rows,
    #     x="Window",
    #     y=[1]*len(wanted_rows),
    #     color="TopologyID",
    #     color_discrete_map=color_mapping,
    #     nbins=int(chromosome_df["End"].max()/window_size)
    # )
    grouped_topology_df = wanted_rows.groupby(by='TopologyID')

    # Build figure
    fig = go.Figure()
    for topology, data in grouped_topology_df:
        fig.add_trace(
            go.Scatter(
                x=data['Window'],
                y=[1]*len(data),
                name=topology,
                legendgroup=topology,
                mode='markers',
                marker_symbol='line-ns-open',
                marker_size=225,
                # marker_line_width=2,
                marker_color=[color_mapping[topology]]*len(data),
                # showlegend = False
            ),
        )
    
    # Update layout + axes
    fig.update_layout(
        template=template,
        legend_title_text='Topology',
        margin=dict(
            l=60,
            r=50,
            b=40,
            t=40,
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            traceorder='normal',
        ),
        hovermode="x unified",
    )
    fig.update_xaxes(
        fixedrange=True,
        linewidth=axis_line_width,
        rangemode="tozero",
        range=dataRange,
        showgrid=xaxis_gridlines,
    )
    fig.update_yaxes(
        fixedrange=True,
        linewidth=axis_line_width,
        # range=[0, 1],
        showline=False,
        showgrid=yaxis_gridlines,
        showticklabels=False,
        ticklen=0,
        title="",
    )

    return fig


def build_alt_data_graph(
    alt_data_to_graph,
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
):
    # Check input type and graph accordingly
    try:
        input_type = type(topology_df[alt_data_to_graph].dropna().to_list()[0])
    except IndexError:
        return no_data_graph(template)
    if input_type == str:
        alt_data_graph_data = make_alt_data_str_figure(
            alt_data_to_graph,
            chromosome_df,
            color_mapping,
            topology_df,
            window_size,
            template,
            dataRange,
            axis_line_width,
            xaxis_gridlines,
            yaxis_gridlines,
        )
    else:
        alt_data_graph_data = make_alt_data_int_figure(
            alt_data_to_graph,
            color_mapping,
            topology_df,
            template,
            dataRange,
            y_max,
            axis_line_width,
            xaxis_gridlines,
            yaxis_gridlines,
        )
    return alt_data_graph_data


def build_RF_heatmap(
    rfdist_df,
    template,
    color_mapping,
    axis_line_width,
):
    # Reset index and use to get proper window placement
    rfdist_df.reset_index(inplace=True)
    # Create figure + return figure
    fig = px.bar(
        rfdist_df,
        x='index',
        y='NormRF',
        color='TopologyID',
        color_discrete_map=color_mapping,
        text='NormRF',
        range_y=[0, 1.2],
    )

    fig.update_traces(
        textposition='outside', 
        texttemplate='%{text:.2f}',
    )
    fig.update_xaxes(
        dtick=1,
        title='Window',
        # NOTE: Replace index labels with Window values
        tickmode = 'array',
        tickvals = rfdist_df['index'],
        ticktext = [str(s) for s in rfdist_df['Window']],
        linewidth=axis_line_width,
        )
    fig.update_layout(
        title={
            'text': "Normalized RF-Distance",
            'y':0.9,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        template=template,
    )
    return fig


def build_gff_figure(
    data,
    dataRange,
    template,
    axis_line_width,
    xaxis_gridlines,
    yaxis_gridlines,
):
    regionStart, regionEnd = dataRange
    # Show gene names if showing less than 1Mb of data
    # if abs(regionEnd - regionStart) <= 10000000:
    if abs(regionEnd - regionStart) <= 10000000:
        show_gene_names = True
    else:
        show_gene_names = False
    # Separate 
    # group data by feature and gene name
    attr_group = data.groupby(by=['feature', 'attribute', 'strand'])
    positive_text_pos = "top center"
    negative_text_pos = "top center"
    fig = go.Figure()
    for fg, gene_data in attr_group:
        feature, gene, strand = fg
        feature_strand = f"{feature} ({strand})"
        x_values = sorted(gene_data['start'].to_list() + gene_data['end'].to_list())
        # Set color, y-values, and arrow direction
        if strand == '+':
            colorValue = 'red'
            y_values = [feature_strand]*len(x_values)
            markerSymbol = ['square']*(len(x_values)-1) + ['triangle-right']
            text_pos = positive_text_pos
            text_val = [gene] + ['']*(len(x_values)-1)
            if positive_text_pos == "top center":
                positive_text_pos = "bottom center"
            elif positive_text_pos == "bottom center":
                positive_text_pos = "top center"
        else:
            colorValue = '#009BFF'
            y_values = [feature_strand]*len(x_values)
            markerSymbol = ['triangle-left'] + ['square']*(len(x_values)-1)
            text_pos = negative_text_pos
            text_val = ['']*(len(x_values)-1) + [gene]
            if negative_text_pos == "top center":
                negative_text_pos = "bottom center"
            elif negative_text_pos == "bottom center":
                negative_text_pos = "top center"
        if show_gene_names:
            fig.add_trace(go.Scatter(
                x=x_values,
                y=y_values,
                name=gene,
                # legendgroup=gene,
                mode='markers+lines+text',
                marker_symbol=markerSymbol,
                marker_size=8,
                marker_color=colorValue,
                text=text_val,
                textposition=text_pos,
                textfont=dict(
                    size=10,
                ),
                # hoverinfo=['all'],
                hovertemplate=None,
            ))
        else:
            fig.add_trace(go.Scatter(
                x=x_values,
                y=y_values,
                name=gene,
                # legendgroup=gene,
                mode='markers+lines',
                marker_symbol=markerSymbol,
                marker_size=8,
                marker_color=colorValue,
                # hoverinfo=['all'],
                hovertemplate=None,
            ))
    fig.update_layout(
        hovermode="x unified",
        # hovermode="x",
        showlegend=False,
        template=template,
        title='',
        margin=dict(
            l=20,
            r=50,
            b=20,
            t=40,
        ),
    )
    fig.update_xaxes(
        range=dataRange,
        title='Position',
        matches="x",
        rangemode="tozero",
        linewidth=axis_line_width,
        showgrid=xaxis_gridlines,
    )
    fig.update_yaxes(
        automargin=True,
        categoryorder='array',
        fixedrange=True,
        showticklabels=True,
        showgrid=yaxis_gridlines,
        title='',
        linewidth=axis_line_width,
    )
    return fig


# ---------------------------------------------------------------------------------
# -------------------------------- Stats Tab Graphs -------------------------------

def build_topology_frequency_pie_chart(df, template, color_mapping):
    fig = px.pie(
        df,
        values='Frequency',
        names='TopologyID',
        template=template,
        color_discrete_map=color_mapping,
        hover_data=['TopologyID'],
        title='Whole Genome Frequencies'
    )
    fig.update_traces(textposition='inside')
    fig.update_layout(
        uniformtext_minsize=12,
        uniformtext_mode='hide',
        showlegend=False,
        margin=dict(l=20, r=20, t=40, b=10),
        title_x=0.5,
    )
    return fig


def build_rf_graph(df, ref_topo, template, color_mapping, axis_line_width):
    fig = px.bar(
        df, x="TopologyID", y="normRF-Distance",
        color="TopologyID", color_discrete_map=color_mapping,
        text='normRF-Distance')
    fig.update_traces(texttemplate='%{text:.2f}', textposition='inside')
    fig.update_layout(
        title=f"Normalized RF-Distance from {ref_topo}",
        title_x=0.5,
        template=template,
    )
    fig.update_xaxes(linewidth=axis_line_width)
    fig.update_yaxes(linewidth=axis_line_width, range=[0, 1])
    return fig


def build_stats_chromosome_freq_rug_plot(
    df,
    chrom_df,
    chromGroup,
    template,
    color_mapping,
    currTopologies,
    topoOrder,
    window_size,
    axis_line_width,
    xaxis_gridlines,
    yaxis_gridlines,
):
    df = df[df['TopologyID'].isin(currTopologies)]
    df = df[df['Chromosome'].isin(chromGroup)]

    grouped_topology_df = df.groupby(by='TopologyID')
    num_chroms = len(df['Chromosome'].unique())
    chrom_row_dict = {chrom:i for chrom, i in zip(sorted(df['Chromosome'].unique()), range(1, len(df['Chromosome'].unique())+1, 1))}
    row_heights = [1]*num_chroms
    chrom_shapes = []
    # --- Build figure ---
    # If chromosome name longer than 5 characters, use subplot titles 
    # instead of row ittles
    if df.Chromosome.map(len).max() > 5:
        fig = make_subplots(
            rows=num_chroms,
            # row_heights=row_heights,
            subplot_titles=df['Chromosome'].unique(),
            shared_xaxes=True,
            vertical_spacing=0.01,
            cols=1,
        )
    else:
        fig = make_subplots(
            rows=num_chroms,
            # row_heights=row_heights,
            row_titles=[c for c in df['Chromosome'].unique()],
            shared_xaxes=True,
            vertical_spacing=0.01,
            cols=1,
        )
    for topo, data in grouped_topology_df:
        add_legend = True
        for chrom in chrom_row_dict.keys():
            chrom_data = data[data["Chromosome"] == chrom]
            chrom_length_data = chrom_df[chrom_df['Chromosome'] == chrom]
            chrom_length = chrom_length_data['End'].max()
            if len(chrom_data) == 0:
                fig.add_trace(
                    go.Scatter(
                        x=[0],
                        y=[topo],
                        name=topo,
                        legendgroup=topo,
                        mode='markers',
                        marker_symbol='line-ns-open',
                        marker_color=[color_mapping[topo]]*len(chrom_data),
                        showlegend = False,
                    ),
                    row=chrom_row_dict[chrom], col=1,
                )
            elif add_legend:
                fig.add_trace(
                    go.Scatter(
                        x=chrom_data['Window'],
                        y=chrom_data['TopologyID'],
                        name=topo,
                        legendgroup=topo,
                        mode='markers',
                        # marker_size=int(25/len(grouped_topology_df)),
                        marker_symbol='line-ns-open',
                        marker_color=[color_mapping[topo]]*len(chrom_data),
                    ),
                    row=chrom_row_dict[chrom], col=1,
                )
                chrom_shapes.append(dict(type="line", xref="x", yref="y", x0=chrom_length, x1=chrom_length, y0=-1, y1=len(currTopologies), line_width=2))
                add_legend = False
            else:
                fig.add_trace(
                    go.Scatter(
                        x=chrom_data['Window'],
                        y=chrom_data['TopologyID'],
                        name=topo,
                        legendgroup=topo,
                        mode='markers',
                        # marker_size=int(25/len(grouped_topology_df)),
                        marker_symbol='line-ns-open',
                        marker_color=[color_mapping[topo]]*len(chrom_data),
                        showlegend = False,
                    ),
                    row=chrom_row_dict[chrom], col=1,
                )
                chrom_ref = chrom_row_dict[chrom]
                chrom_shapes.append(dict(type="rect", xref=f"x{chrom_ref}", yref=f"y{chrom_ref}", x0=chrom_length, x1=chrom_length, y0=-1, y1=len(currTopologies), line_width=2))

    # Update layout + axes
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.update_xaxes(
        rangemode="tozero",
        range=[0, (chrom_df['End'].max()+(2*window_size))],
        fixedrange=True,
        linewidth=axis_line_width,
        ticklen=0,
        matches="x",
        showgrid=xaxis_gridlines,
    )
    fig.update_yaxes(
        fixedrange=True,
        title="",
        showgrid=yaxis_gridlines,
        showticklabels=False,
        linewidth=axis_line_width,
        categoryarray=topoOrder,
        # matches='y'
    )
    fig.update_layout(
        template=template,
        legend_title_text='Topology',
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            traceorder='normal',
            itemsizing='constant',
        ),
        # hovermode="x unified",
        height=125*num_chroms,
        shapes=chrom_shapes,
        title_x=0.5,
    )
    # Rotate chromosome names to 0-degrees
    for annotation in fig['layout']['annotations']: 
        annotation['textangle']=0
        annotation['align']="center"
    return fig


def build_stats_chromosome_freq_bar_plot(
    df,
    template,
    color_mapping,
    currTopologies,
    axis_line_width,
    chromGroup,
    xaxis_gridlines,
    yaxis_gridlines,
):
    # Filter df to chromosomes in group
    df = df[df['Chromosome'].isin(chromGroup)]
    df = df[df['TopologyID'].isin(currTopologies)]
    number_of_chrom_rows = len(df["Chromosome"].unique()) // 3
    fig = px.bar(
        df,
        x='TopologyID',
        y='Frequency',
        facet_col='Chromosome',
        facet_col_wrap=3,
        facet_row_spacing=0.05,
        color='TopologyID',
        template=template,
        color_discrete_map=color_mapping,
        text='Frequency',
        height=int(500*number_of_chrom_rows),
    )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.update_traces(texttemplate='%{text:.2}', textposition='outside')
    # Remove y-axis labels
    for axis in fig.layout:
        if type(fig.layout[axis]) == go.layout.YAxis:
            fig.layout[axis].title.text = ''
    fig.update_layout(
        uniformtext_minsize=12, 
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            traceorder='normal',
        ),
        margin=dict(l=10, r=10, t=10, b=10),
        title="",
        annotations = list(fig.layout.annotations) + 
        [go.layout.Annotation(
                x=-0.07,
                y=0.5,
                font=dict(
                    size=12,
                    # color='white',
                ),
                showarrow=False,
                text="Frequency",
                textangle=-90,
                xref="paper",
                yref="paper"
            )
        ],
        title_x=0.5,
    )
    fig.update_xaxes(
        title="",
        linewidth=axis_line_width,
        showgrid=xaxis_gridlines,
    )
    fig.update_yaxes(
        range=[0, 1.1],
        matches='y',
        linewidth=axis_line_width,
        showgrid=yaxis_gridlines,
    )
    return fig


def build_stats_chromosome_freq_pie_plot(
    df,
    template,
    color_mapping,
    chromGroup,
):
    # Filter df to chromosomes in group
    df = df[df['Chromosome'].isin(chromGroup)]
    number_of_chrom_rows = (len(df["Chromosome"].unique()) // 3)+(math.ceil(len(df["Chromosome"].unique()) % 3))
    specs = [[{'type':'domain'}, {'type':'domain'}, {'type':'domain'}] for _ in range(number_of_chrom_rows)]
    fig = make_subplots(
        rows=number_of_chrom_rows,
        cols=3,
        specs=specs,
        vertical_spacing=0.03,
        horizontal_spacing=0.001,
        subplot_titles=df["Chromosome"].unique(),
        # row_heights=[2]*number_of_chrom_rows,
        column_widths=[2]*3,
    )
    col_pos = 1
    row_num = 1
    for c in df['Chromosome'].unique():
        chrom_df = df[df["Chromosome"] == c]
        fig.add_trace(go.Pie(labels=chrom_df["TopologyID"], values=chrom_df['Frequency']), row=row_num, col=col_pos)
        if col_pos == 3:
            col_pos = 1
            row_num += 1
        else:
            col_pos += 1
    fig.update_traces(textposition='inside')
    fig.update_layout(
        uniformtext_minsize=12, 
        # uniformtext_mode='hide',
        
        showlegend=True,
        # margin=dict(l=20, r=20, t=40, b=10),
        template=template,
        height=int(200*number_of_chrom_rows),
    )
    return fig


def build_stats_chromosome_freq_tile_plot(
    df,
    chrom_df,
    template,
    color_mapping,
    currTopologies,
    topoOrder,
    window_size,
    axis_line_width,
    chromGroup,
    xaxis_gridlines,
    yaxis_gridlines,
):
    """
    Max chromosomes per graph if # current_topologies <= 3: 20
    Max chromosomes per graph if # current_topologies > 3: 20/2

    Returns: List of figures to display
    """
    df = df[df['TopologyID'].isin(currTopologies)]
    df = df[df['Chromosome'].isin(chromGroup)]
    grouped_topology_df = df.groupby(by='TopologyID')
    num_chroms = len(df['Chromosome'].unique())
    chrom_row_dict = {chrom:i for chrom, i in zip(sorted(df['Chromosome'].unique()), range(1, len(df['Chromosome'].unique())+1, 1))}
    row_heights = [1]*num_chroms
    chrom_shapes = []
    # --- Build figure ---
    # If longest chromosome name longer 
    # than 5 characters, use subplot titles 
    # instead of row ittles
    if df.Chromosome.map(len).max() > 5:
        fig = make_subplots(
            rows=num_chroms,
            # row_heights=row_heights,
            subplot_titles=df['Chromosome'].unique(),
            vertical_spacing=0.01,
            cols=1,
            shared_xaxes=True,
        )
    else:
        fig = make_subplots(
            rows=num_chroms,
            # row_heights=row_heights,
            cols=1,
            shared_xaxes=True,
            vertical_spacing=0.001,
            row_titles=[c for c in df['Chromosome'].unique()]
        )
    for topo, data in grouped_topology_df:
        add_legend = True
        for chrom in chrom_row_dict.keys():
            chrom_data = data[data["Chromosome"] == chrom]
            chrom_length_data = chrom_df[chrom_df['Chromosome'] == chrom]
            chrom_length = chrom_length_data['End'].max()
            if add_legend:
                fig.add_trace(
                    go.Scatter(
                        x=chrom_data['Window'],
                        y=[1]*len(chrom_data),
                        name=topo,
                        legendgroup=topo,
                        mode='markers',
                        marker_symbol='line-ns-open',
                        marker_size=50,
                        marker_line_width=1,
                        marker_color=[color_mapping[topo]]*len(chrom_data),
                        # showlegend = False
                    ),
                    row=chrom_row_dict[chrom], col=1,
                )
                chrom_shapes.append(dict(type="line", xref="x", yref="y", x0=chrom_length, x1=chrom_length, y0=-1, y1=len(currTopologies), line_width=2))
                add_legend = False
            else:
                fig.add_trace(
                    go.Scatter(
                        x=chrom_data['Window'],
                        y=[1]*len(chrom_data),
                        name=topo,
                        legendgroup=topo,
                        mode='markers',
                        marker_symbol='line-ns-open',
                        marker_size=50,
                        marker_line_width=1,
                        marker_color=[color_mapping[topo]]*len(chrom_data),
                        showlegend = False
                    ),
                    row=chrom_row_dict[chrom], col=1,
                )
                chrom_ref = chrom_row_dict[chrom]
                chrom_shapes.append(dict(type="rect", xref=f"x{chrom_ref}", yref=f"y{chrom_ref}", x0=chrom_length, x1=chrom_length, y0=-1, y1=len(currTopologies), line_width=2))
    # Set graph height
    if num_chroms < 5:
        graph_height = 100*num_chroms
    else:
        graph_height = 75*num_chroms

    # Update layout + axes
    fig.update_layout(
        template=template,
        legend_title_text='Topology',
        margin=dict(
            l=60,
            r=50,
            b=40,
            t=40,
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            traceorder='normal',
            itemsizing='constant',
        ),
        hovermode="x unified",
        height=graph_height,
        shapes=chrom_shapes,
        title_x=0.5,
    )
    fig.update_xaxes(
        linewidth=axis_line_width,
        fixedrange=True,
        rangemode="tozero",
        range=[0, chrom_df['End'].max()],
        ticklen=0,
        showgrid=xaxis_gridlines,
    )
    fig.update_yaxes(
        # categoryarray=topoOrder,
        fixedrange=True,
        linewidth=axis_line_width,
        showgrid=yaxis_gridlines,
        showticklabels=False,
        title="",
        ticklen=0,
    )
    # Rotate chromosome names to 0-degrees
    for annotation in fig['layout']['annotations']: 
        annotation['textangle']=0
        annotation['align']="center"
    return fig

# ---------------------------------------------------------------------------------
# ------------------------- Graph Customization Functions -------------------------

def set_topology_colors(data, color):
    df = pd.read_json(data) 
    # Set colors to current_topologies
    sorted_topologies = df.assign(freq=df.groupby('TopologyID')['TopologyID'].transform('count')).sort_values(by=['freq','TopologyID'],ascending=[False,True]).loc[:,['TopologyID']]
    unique_topos = sorted_topologies["TopologyID"].unique()
    color_list = (color * ((len(unique_topos) // len(color))))+ color[:len(unique_topos) % len(color)]
    output_dict = dict()
    for s, c in zip(unique_topos, color_list):
        output_dict[s] = c
    return output_dict


def get_RFxpos(hoverdata, df):
    hoverdata = hoverdata['points'][0]
    if ('customdata' in hoverdata.keys()) or ('marker.color' in hoverdata.keys()):
        return int(hoverdata['x'])  # Window
    else:
        return df.loc[hoverdata['binNumber']]['Window']

# ---------------------------------------------------------------------------------
# ------------------------- Init + Empty Graph Functions --------------------------

def no_data_graph(template):
    """This function returns a blank figure with a "NO DATA" watermark"""
    fig = go.Figure()
    fig.update_layout(
        template=template,
        title='',
        annotations=[
            dict(
                name="draft watermark",
                text="NO DATA",
                textangle=0,
                opacity=0.5,
                font=dict(color="white", size=50),
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
            )
        ],
    )
    fig.update_xaxes(showgrid=False, range=[0.2, 1], zeroline=False, visible=False)
    fig.update_yaxes(showgrid=False, range=[0.2, 1], zeroline=False, visible=False)
    return fig


def init_data_graph(template):
    """
    This function returns a blank figure with a "NO DATA" watermark.
    """
    fig = go.Figure()
    fig.update_layout(
        template=template,
        annotations=[
            dict(
                name="draft watermark",
                text="No Data Loaded",
                textangle=0,
                opacity=0.9,
                font=dict(color="white", size=50),
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
            )
        ],
    )
    fig.update_xaxes(range=[0.2, 1], showgrid=False, visible=False, zeroline=False)
    fig.update_yaxes(range=[0.2, 1], showgrid=False, visible=False, zeroline=False)
    return fig


def loading_data_graph(template):
    """
    This function returns a blank figure with a "NO DATA" watermark.
    """
    fig = go.Figure()
    fig.update_layout(
        template=template,
        annotations=[
            dict(
                name="draft watermark",
                text="GATHERING DATA...",
                textangle=0,
                opacity=0.9,
                font=dict(color="white", size=100),
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
            )
        ],
    )
    fig.update_xaxes(showgrid=False, range=[0.2, 1], zeroline=False, visible=False)
    fig.update_yaxes(showgrid=False, range=[0.2, 1], zeroline=False, visible=False)
    return fig


def init_RF_graph(template):
    """
    This function returns a blank figure with a "NO DATA" watermark.
    """
    fig = go.Figure()
    fig.update_layout(
        template=template,
        annotations=[
            dict(
                name="draft watermark",
                text="Hover Over Data to Activate",
                textangle=0,
                opacity=0.9,
                font=dict(color="white", size=100),
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
            )
        ],
    )
    fig.update_xaxes(showgrid=False, range=[0.2, 1], zeroline=False, visible=False)
    fig.update_yaxes(showgrid=False, range=[0.2, 1], zeroline=False, visible=False)
    return fig


def no_tree_data(template, msg):
    """
    This function returns a blank figure with a "NO DATA" watermark.
    """
    fig = go.Figure()
    fig.update_layout(
        template=template,
        annotations=[
            dict(
                name="draft watermark",
                text=msg,
                textangle=0,
                opacity=0.9,
                font=dict(size=25),
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
            )
        ],
    )
    fig.update_xaxes(showgrid=False, range=[0.2, 1], zeroline=False, visible=False)
    fig.update_yaxes(showgrid=False, range=[0.2, 1], zeroline=False, visible=False)
    return fig


def zoom_in_gff(template):
    """
    This function returns a blank figure with a "NO DATA" watermark.
    """
    fig = go.Figure()
    fig.update_layout(
        template=template,
        annotations=[
            dict(
                name="draft watermark",
                text="Zoom in to view",
                textangle=0,
                opacity=0.9,
                font=dict(color="white", size=100),
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                showarrow=False,
            )
        ],
    )
    fig.update_xaxes(showgrid=False, range=[0.2, 1], zeroline=False, visible=False)
    fig.update_yaxes(showgrid=False, range=[0.2, 1], zeroline=False, visible=False)
    return fig


# ---------------------------------------------------------------------------------
# --------------------------- Input File Verification -----------------------------
def validate_chrom_lengths(chromDF, tvDF):
    """Ensure all chromosomes in chromDF are present in tvDF.
       Chromosome length file can contain for chromosomes than TV file,
       but not the other way around.
       Return True if all are found, False if not."""
    chrom_names = chromDF['Chromosome'].unique()
    tv_chrom_names = tvDF['Chromosome'].unique()
    missing_chromosomes = []
    valid = True
    issue_files = []
    # Check chromosome length file against TV file
    for c in chrom_names:
        if c not in tv_chrom_names:
            missing_chromosomes.append(c)
            valid = False
            issue_files.append("Chromosome Length File")
            continue
        else:
            continue
    # Check TV file against chromosome length file
    for c in tv_chrom_names:
        if c not in chrom_names:
            missing_chromosomes.append(c)
            valid = False
            issue_files.append("Tree Viewer File")
            continue
        else:
            continue
    try:
        if not valid:
            missing_chroms = "-".join(missing_chromosomes)
            if len(issue_files) > 1:
                missing_files = " & ".join(list(set(issue_files)))
            else:
                missing_files = issue_files[0]
            msg = f"ERROR: Chromosome(s) {missing_chroms} is missing from {missing_files}, please validate consistency of chromosomes between files"
            return msg, False
        else:
            return None, True
    except UnboundLocalError:
        return None, True


def get_taxa_from_tree(tree):
    """Collect leaf names from tree"""
    tree = Tree(tree)
    taxa = []
    for leaf in tree.iter_leaves():
        taxa.append(leaf.name)
    return sorted(taxa)


def validate_gff_gtf_filename(f):
    """Ensure file extension is gff or gtf"""
    if "gtf" in f.lower():
        return True
    elif "gff" in f.lower():
        return True
    else:
        return False


def get_y_max_list(alt_dropdown_options, topology_df):
    """Generate list of max y-values for additinal data"""
    y_maxes = []
    for i in alt_dropdown_options:
        try:
            data_type = type(topology_df[i][0])
        except KeyError:
            data_type = str
        if data_type == str:
            y_maxes.append(1)
        else:
            y_maxes.append(topology_df[i].max())
    return y_maxes


def make_basic_info_df(df, topo_freq_df, chrom_lengths, window_size):
    genome_length = chrom_lengths["End"].sum()
    avg_chrom_len = int(chrom_lengths["End"].mean())
    missing_data_coverage = round((100-(len(df) / (data_utils.roundUp(chrom_lengths["End"].sum(), window_size) / window_size) * 100)), 2)
    number_of_topologies = len(df["TopologyID"].unique())
    avg_topo_freq = round(topo_freq_df['Frequency'].mean(), 5)
    median_topo_freq = round(topo_freq_df['Frequency'].median(), 5)
    df = pd.DataFrame(
        {
            "": [
                "Genome Length",
                "Average Chromosome Length",
                "Total Missing Data",
                "Number of Topologies",
                "Average Topology Frequency",
                "Median Topology Frequency",
            ],
            " ": [
                f"{genome_length:,}-bp",
                f"{avg_chrom_len:,}-bp",
                f"{missing_data_coverage}%",
                f"{number_of_topologies:,}",
                f"{avg_topo_freq}",
                f"{median_topo_freq}",
            ]
        }
    )
    return df


def validate_tree_viewer_input(df):
    """Return False when required headers are not present/correct"""
    def fix_column_names(columns):
        """ Fix column names """
        if columns[:4] == ["Chromosome", "Window", "NewickTree", "TopologyID"]:
            return columns
        else:
            return ["Chromosome", "Window", "NewickTree", "TopologyID"] + columns[4:]

    def check_newick(df):
        """Check if string contains basic newick characters"""
        if "(" not in df["NewickTree"][0]:
            return False
        elif ")" not in df["NewickTree"][0]:
            return False
        elif ";" not in df["NewickTree"][0]:
            return False
        else:
            return True

    def check_window(df):
        """Return False if row type is not int"""
        if type(df["Window"][0]) == np.int32:
            return True
        elif type(df["Window"][0]) == np.int64:
            return True
        else:
            return False

    # Fix required headers if needed
    cols = fix_column_names(list(df.columns))
    df.columns = cols
    # Check reuqired column types
    newick_check = check_newick(df)
    window_check = check_window(df)

    if not newick_check:
        return False
    elif not window_check:
        return False
    else:
        return df


def tv_header_validation(df):
    """Return False if first four required column headers are not valid"""
    required_cols = list(df.columns[:4])
    try:
        assert required_cols == ["Chromosome", "Window", "NewickTree", "TopologyID"]
        return True
    except AssertionError:
        return False


# ---------------------------------------------------------------------------------
# --------------------------- Tree Prune Export Tools -----------------------------
def prune_tree(x, prune_taxa_choices):
    tree = Tree(x)
    try:
        tree.prune(prune_taxa_choices, preserve_branch_length=True)
    except ValueError:
        # Assumes taxa in dropdown selection
        # is not found in a particular topology/tree
        # Solution is to check list and remove taxa
        # not present in tree
        tree_taxa = tree.get_leaf_names()
        trimmed_taxa_list = [t for t in prune_taxa_choices if t in tree_taxa]
        tree.prune(trimmed_taxa_list, preserve_branch_length=True)
    return tree.write()


def remove_heterotachy_info(l):
    """Remove any information in bracketsete3 
       does not support this format of newick"""
    if ("[" not in l) and ("]" not in l):
        return l
    open_brackets = [i for i, x in enumerate(l) if x == "["]
    close_brackets = [i for i, x in enumerate(l) if x == "]"]
    final_string = f'{l[:open_brackets[0]]}'
    for ob, cb in zip(open_brackets[1:], close_brackets[:-1]):
        final_string += l[cb+1:ob]
    final_string += l[close_brackets[-1]+1:]
    return final_string


def tv_topobinner(df):
    """Bin tree topologies that have RF-distance of 0"""
    trees = df['NewickTree']
    topologies = dict()
    topoCount = 1
    
    for n, t in enumerate(trees):
        if len(topologies.keys()) == 0:
            topologies[n] = {'count': 1, 'idx': [n]}
            continue
        else:
            # Iterate through topology list
            # add new topology if no rf == 0
            # increase count if rf == 0 with topology 
            new_topology = True
            for idx in topologies.keys():
                t1 = Tree(remove_heterotachy_info(t))
                t2 = Tree(remove_heterotachy_info(df.at[idx, 'NewickTree']))
                comparison = t1.compare(t2)
                rf = comparison['rf']
                if rf == 0:
                    topologies[idx]['count'] += 1
                    topologies[idx]['idx'].append(n)
                    new_topology = False
                    break
                else:
                    continue
            if new_topology:
                topologies[n] = {'count': 1, 'idx': [n]}
                continue
            else:
                continue
    # Sort topologies dictionary by 'count'
    topologies = {k: v for k, v in sorted(topologies.items(), key=lambda item: item[1]['count'], reverse=True)}

    # Update DataFrame TopologyID column with results
    for topo in topologies.keys():
        idx = topologies[topo]['idx']
        topoName = f'topo{topoCount}'
        for i in idx:
            df.at[i, 'TopologyID'] = topoName
            continue
        topoCount += 1
    return df


def mygrouper(n, iterable):
    args = [iter(iterable)] * n
    return ([e for e in t if e != None] for t in itertools.zip_longest(*args))


def make_topo_freq_table(df_grouped):
    dataTableDF = pd.DataFrame(columns=["Chromosome", "TopologyID", 'Frequency'], index=range(len(df_grouped)))
    idx = 0
    for chrom, data in df_grouped:
        chromFreqs = data["TopologyID"].value_counts()/len(data)
        freqTopoOrder = [i for i in chromFreqs.index]
        freqs = [f for f in chromFreqs]
        for t, f in zip(freqTopoOrder, freqs):
            dataTableDF.at[idx, 'Chromosome'] = chrom
            dataTableDF.at[idx, 'TopologyID'] = t
            dataTableDF.at[idx, 'Frequency'] = round(f, 3)
            idx += 1
            continue
    return dataTableDF


def get_gridline_bools(axis_gridlines):
    """If gridlines ON, return True else False"""
    if 'xaxis' in axis_gridlines:
        xaxis_gridlines = True
    else:
        xaxis_gridlines = False
    if 'yaxis' in axis_gridlines:
        yaxis_gridlines = True
    else:
        yaxis_gridlines = False
    return xaxis_gridlines, yaxis_gridlines

