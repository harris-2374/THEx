import math
import itertools
import re

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scipy.stats as ss
import scikit_posthocs as sp
from dash_table.Format import Format, Scheme
from Bio import Phylo
from ete3 import Tree
from plotly.subplots import make_subplots

# -------------------------------------------------------------------------------------
# --------------------------------------- Classes -------------------------------------
class DrawTree():
    def __init__(self, newicktree, template, topology, color_map, branch_len, font_family):
        self.newicktree = Phylo.read(newicktree, "newick")
        self.template = template
        self.topology = topology
        self.color_map = color_map
        self.branch_len = branch_len
        self.font_family = font_family

    def create_square_tree(self):

        def get_x_coordinates(tree):
            """Associates to each clade an x-coord.
            returns dict {clade: x-coord}
            """
            if self.branch_len:
                xcoords = tree.depths(unit_branch_lengths=True)
            else:
                max_depth = max(tree.depths().values())
                xcoords = {k: abs(v-max_depth) for k, v in tree.depths().items()}
                # xcoords = 

            # tree.depth() maps tree clades to depths (by branch length).
            # returns a dict {clade: depth} where clade runs over all Clade instances of the tree, and depth is the distance from root to clade

            #  If there are no branch lengths, assign unit branch lengths
            if not max(xcoords.values()):
                xcoords = tree.depths(unit_branch_lengths=True)
            return xcoords

        def get_y_coordinates(tree, dist=1.3):
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
                # This is intermediate placement of internal nodes
                ycoords[clade] = (ycoords[clade.clades[0]] + ycoords[clade.clades[-1]]) / 2

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
                textfont=dict(color=text_color, size=12),
                showlegend=False,
                name=elt,
            )
            nodes.append(node)
        # Set graph x-range
        if self.branch_len:
            layout = dict(
                autosize=True,
                showlegend=False,
                template=self.template,
                dragmode="pan",
                margin=dict(t=0, b=0, r=0, l=0),
                xaxis=dict(
                    showline=True,
                    zeroline=False,
                    visible=False,
                    showgrid=False,
                    showticklabels=True,
                    range=[-1, max(x_coords.values())+2]
                ),
                yaxis=axis,
                hovermode="closest",
                shapes=line_shapes,
                font=dict(family=self.font_family,size=14),
            )
        else:
            layout = dict(
                autosize=True,
                showlegend=False,
                template=self.template,
                dragmode="pan",
                margin=dict(t=0, b=0, r=0, l=0),
                xaxis=dict(
                    showline=True,
                    zeroline=False,
                    visible=True,
                    showgrid=False,
                    showticklabels=True,
                    autorange='reversed',
                ),
                yaxis=axis,
                hovermode="closest",
                shapes=line_shapes,
                font=dict(family=self.font_family,size=14),
            )
        fig = go.Figure(data=nodes, layout=layout)
        return fig

    def create_circular_tree(self):
        def get_circular_tree_data(tree, order='level', dist=1, start_angle=0, end_angle=360, start_leaf='first'):
            """Define  data needed to get the Plotly plot of a circular tree
               Source code found at: https://chart-studio.plotly.com/~empet/14834.embed
            """
            # tree:  an instance of Bio.Phylo.Newick.Tree or Bio.Phylo.PhyloXML.Phylogeny
            # order: tree  traversal method to associate polar coordinates to its nodes
            # dist:  the vertical distance between two consecutive leafs in the associated rectangular tree layout
            # start_angle:  angle in degrees representing the angle of the first leaf mapped to a circle
            # end_angle: angle in degrees representing the angle of the last leaf
            # the list of leafs mapped in anticlockwise direction onto circles can be tree.get_terminals() 
            # or its reversed version tree.get_terminals()[::-1]. 
            # start leaf: is a keyword with two possible values"
            # 'first': to map  the leafs in the list tree.get_terminals() onto a circle,
            #         in the counter-clockwise direction
            # 'last': to map  the leafs in the  list, tree.get_terminals()[::-1] 
            
            start_angle *= np.pi/180 # conversion to radians
            end_angle *= np.pi/180
            
            def get_radius(tree):
                """
                Associates to  each clade root its radius, equal to the distance from that clade to the tree root
                returns dict {clade: node_radius}
                """
                if self.branch_len:
                    node_radius = tree.depths(unit_branch_lengths=True)
                else:
                    node_radius = tree.depths()
                
                #  If the tree did not record  the branch lengths  assign  the unit branch length
                #  (ex: the case of a newick tree "(A, (B, C), (D, E))")
                if not np.count_nonzero(node_radius.values()):
                    node_radius = tree.depths(unit_branch_lengths=True)
                return node_radius
        
            
            def get_vertical_position(tree):
                """
                returns a dict {clade: ycoord}, where y-coord is the cartesian y-coordinate 
                of a  clade root in a rectangular phylogram
                
                """
                n_leafs = tree.count_terminals() # Counts the number of tree leafs.
                
                # Assign y-coordinates to the tree leafs
                if start_leaf == 'first':
                    node_ycoord = dict((leaf, k) for k, leaf in enumerate(tree.get_terminals()))
                elif start_leaf == 'last':
                    node_ycoord = dict((leaf, k) for k, leaf in enumerate(reversed(tree.get_terminals())))
                else:
                    raise ValueError("start leaf can be only 'first' or 'last'")
                    
                def assign_ycoord(clade):#compute the y-coord for the root of this clade
                    for subclade in clade:
                        if subclade not in node_ycoord: # if the subclade root hasn't a y-coord yet
                            assign_ycoord(subclade)
                    node_ycoord[clade] = 0.5 * (node_ycoord[clade.clades[0]] + node_ycoord[clade.clades[-1]])

                if tree.root.clades:
                    assign_ycoord(tree.root)
                return node_ycoord

            node_radius = get_radius(tree)
            node_ycoord = get_vertical_position(tree)
            y_vals = node_ycoord.values()
            ymin, ymax = min(y_vals), max(y_vals)
            ymin -= dist # this dist subtraction is necessary to avoid coincidence of the  first and last leaf angle
                        # when the interval  [ymin, ymax] is mapped onto [0, 2pi],
                        
            def ycoord2theta(y):
                # maps an y in the interval [ymin-dist, ymax] to the interval [radian(start_angle), radian(end_angle)]
                
                return start_angle + (end_angle - start_angle) * (y-ymin) / float(ymax-ymin)


            def get_points_on_lines(linetype='radial', x_left=0, x_right=0, y_right=0,  y_bot=0, y_top=0):
                """
                - define the points that generate a radial branch and the circular arcs, perpendicular to that branch
                
                - a circular arc (angular linetype) is defined by 10 points on the segment of ends
                (x_bot, y_bot), (x_top, y_top) in the rectangular layout,
                mapped by the polar transformation into 10 points that are spline interpolated
                - returns for each linetype the lists X, Y, containing the x-coords, resp y-coords of the
                line representative points
                """
            
                if linetype == 'radial':
                    theta = ycoord2theta(y_right) 
                    X = [x_left*np.cos(theta), x_right*np.cos(theta), None]
                    Y = [x_left*np.sin(theta), x_right*np.sin(theta), None]
                
                elif linetype == 'angular':
                    theta_b = ycoord2theta(y_bot)
                    theta_t = ycoord2theta(y_top)
                    t = np.linspace(0,1, 10)# 10 points that span the circular arc 
                    theta = (1-t) * theta_b + t * theta_t
                    X = list(x_right * np.cos(theta)) + [None]
                    Y = list(x_right * np.sin(theta)) + [None]
                
                else:
                    raise ValueError("linetype can be only 'radial' or 'angular'")
            
                return X,Y   
                

            def get_line_lists(clade,  x_left,  xlines, ylines, xarc, yarc):
                """Recursively compute the lists of points that span the tree branches"""
                
                # xlines, ylines  - the lists of x-coords, resp y-coords of radial edge ends
                # xarc, yarc - the lists of points generating arc segments for tree branches
                
                x_right = node_radius[clade]
                y_right = node_ycoord[clade]
        
                X,Y = get_points_on_lines(linetype='radial', x_left=x_left, x_right=x_right, y_right=y_right)
        
                xlines.extend(X)
                ylines.extend(Y)
        
                if clade.clades:
                
                    y_top = node_ycoord[clade.clades[0]]
                    y_bot = node_ycoord[clade.clades[-1]]
            
                    X,Y = get_points_on_lines(linetype='angular',  x_right=x_right, y_bot=y_bot, y_top=y_top)
                    xarc.extend(X)
                    yarc.extend(Y)
            
                    # get and append the lists of points representing the  branches of the descedants
                    for child in clade:
                        get_line_lists(child, x_right, xlines, ylines, xarc, yarc)

            xlines = []
            ylines = []
            xarc = []
            yarc = []
            get_line_lists(tree.root,  0, xlines, ylines, xarc, yarc)  
            xnodes = []
            ynodes = []

            for clade in tree.find_clades(order='preorder'): #it was 'level'
                theta = ycoord2theta(node_ycoord[clade])
                xnodes.append(node_radius[clade]*np.cos(theta))
                ynodes.append(node_radius[clade]*np.sin(theta))
                
            return xnodes, ynodes,  xlines, ylines, xarc, yarc

        if 'dark' in self.template:
            text_color = 'white'
        else:
            text_color = 'black'
        line_color = self.color_map[self.topology]
        
        tree = self.newicktree
        tree.ladderize()

        traverse_order = 'preorder'

        all_clades=list(tree.find_clades(order=traverse_order))
        for k in range(len((all_clades))):
            all_clades[k].id=k

        xnodes, ynodes,  xlines, ylines, xarc, yarc = get_circular_tree_data(tree, order=traverse_order, start_leaf='last')

        tooltip=[]
        clade_names=[]
        color=[]
        
        for clade in tree.find_clades(order=traverse_order):
            if self.branch_len:
                branch_length = 1
            else:
                branch_length = clade.branch_length
            if clade.name and clade.confidence and clade.branch_length:
                tooltip.append(f"name: {clade.name}<br>branch-length: {branch_length}\
                            <br>confidence: {int(clade.confidence)}")
                color.append[clade.confidence.value]
                clade_names.append(clade.name)
            elif clade.name is None and clade.branch_length is not None and clade.confidence is not None: 
                color.append(clade.confidence)
                clade_names.append(clade.name)
                tooltip.append(f"branch-length: {branch_length}\
                            <br>confidence: {int(clade.confidence)}")
            elif clade.name and clade.branch_length and clade.confidence is None:
                tooltip.append(f"name: {clade.name}<br>branch-length: {branch_length}")
                color.append(-1)
                clade_names.append(clade.name)
            else: 
                tooltip.append('')
                color.append(-1)
                clade_names.append(clade.name)

        trace_nodes=dict(type='scatter',
            x=xnodes,
            y= ynodes, 
            mode='markers+text',
            marker=dict(color=text_color, size=8),
            text=clade_names,
            textposition='top center',
            textfont=dict(color=text_color, size=12),
            hoverinfo='text',
            hovertemplate=tooltip,
        )

        trace_radial_lines=dict(type='scatter',
            x=xlines,
            y=ylines, 
            mode='lines',
            line=dict(color=line_color, width=1),
            hoverinfo='none',
        )

        trace_arcs=dict(type='scatter',
            x=xarc,
            y=yarc,
            mode='lines',
            line=dict(color=line_color, width=1, shape='spline'),
            hoverinfo='none',
        )
        layout=dict(
            font=dict(family=self.font_family,size=14),
            autosize=True,
            showlegend=False,
            template=self.template,
            xaxis=dict(visible=False),
            yaxis=dict(visible=False), 
            hovermode='closest',
            margin=dict(t=20, b=10, r=20, l=10, pad=20),
        )
        fig = go.Figure(data=[trace_radial_lines, trace_arcs, trace_nodes], layout=layout)
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
                range=[0, (max(x_coords.values())+2)]
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
    font_family,
    whole_genome,
):
    # sort dataframe
    topology_df.sort_values(by=["Window"], inplace=True)
    topology_df.fillna("NULL", inplace=True)

    # Build graph
    if whole_genome:
        fig = px.histogram(
            topology_df,
            x="Window",
            y=[1]*len(topology_df),
            category_orders={"Chromosome": chromosome_df['Chromosome']},
            color=alt_data_to_graph,
            color_discrete_sequence=list(color_mapping.values()),
            nbins=int(chromosome_df["End"].max()/window_size),
            facet_row="Chromosome",
        )
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
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
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top',
            },
            hovermode="x unified",
            font=dict(family=font_family,),
            height=100*len(topology_df["Chromosome"].unique())
        )
    else:
        fig = px.histogram(
            topology_df,
            x="Window",
            y=[1]*len(topology_df),
            color=alt_data_to_graph,
            color_discrete_sequence=list(color_mapping.values()),
            nbins=int(chromosome_df["End"].max()/window_size),
        )
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
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top',
            },
            hovermode="x unified",
            font=dict(family=font_family,),
        )
    
    if dataRange:
        fig.update_xaxes(
            title="Position",
            range=dataRange,
            showline=True,
            showgrid=xaxis_gridlines,
            linewidth=axis_line_width,
        )
    else:
        fig.update_xaxes(
            title="Position",
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
    chromosome_df,
    template,
    dataRange,
    y_max,
    axis_line_width,
    xaxis_gridlines,
    yaxis_gridlines,
    font_family,
    whole_genome,
):
    # sort dataframe
    topology_df = topology_df.sort_values(by=["Window"])
    y_range = [0, (y_max*1.1)]
    # Build graph
    if whole_genome:
        fig = px.line(
            topology_df,
            x="Window",
            y=alt_data_to_graph,
            category_orders={"Chromosome": chromosome_df['Chromosome']},
            color_discrete_sequence=list(color_mapping.values()),
            facet_row="Chromosome",
        )
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
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
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top',
            },
            hovermode="x unified",
            font=dict(family=font_family,),
            height=100*len(topology_df["Chromosome"].unique()),
        )
    else:
        fig = px.line(
            topology_df,
            x="Window",
            y=alt_data_to_graph,
            color_discrete_sequence=list(color_mapping.values()),
        )
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
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top',
            },
            hovermode="x unified",
            font=dict(family=font_family,),
        )

    # Update X-axis
    if dataRange:
        fig.update_xaxes(
            title="Position",
            range=dataRange,
            showline=True,
            showgrid=xaxis_gridlines,
            linewidth=axis_line_width,
            matches=None,
        )
    else:
        fig.update_xaxes(
            title="Position",
            range=[0, chromosome_df['End'].max()],
            showline=True,
            showgrid=xaxis_gridlines,
            linewidth=axis_line_width,
            matches=None,
        )
    if y_max < 0.1:
        fig.update_yaxes(
            # fixedrange=True,
            linewidth=axis_line_width,
            range=y_range,
            showgrid=yaxis_gridlines,
            showline=True,
            title="Value",
            showexponent = 'all',
            exponentformat = 'e',
        )
    else:
        fig.update_yaxes(
            # fixedrange=True,
            linewidth=axis_line_width,
            range=y_range,
            showgrid=yaxis_gridlines,
            showline=True,
            title="Value",
        )
    return fig


# ----------------------------------------------------------------------------------------
# -------------------------- Single Chromosome Graph Functions ---------------------------
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
    font_family,
):
    # --- Set up topology data ---
    # Extract current topology data
    if (type(current_topologies) == str) or (type(current_topologies) == int):
        wanted_rows = topology_df[topology_df["TopologyID"] == current_topologies]
    elif type(current_topologies) == list:
        wanted_rows = topology_df[topology_df["TopologyID"].isin(current_topologies)]

    # Group data by topology ID
    grouped_topology_df = wanted_rows.sort_values(['TopologyID'], ascending=False).groupby(by='TopologyID')

    # Set row heights based on number of current_topologies being shown
    if len(current_topologies) <= 6:
        subplot_row_heights = [1, 1]
    elif len(current_topologies) <= 8:
        subplot_row_heights = [4, 2]
    else:
        subplot_row_heights = [8, 2]
    # Build figure
    fig = make_subplots(rows=2, cols=1, vertical_spacing=0.05, shared_xaxes=True)
    for topology, data in grouped_topology_df:
        fig.add_trace(
            go.Scatter(
                x=data['Window'],
                y=data['TopologyID'],
                name=topology,
                legendgroup=topology,
                mode='markers',
                marker_symbol='line-ns-open',
                marker_line_width=1,
                marker_color=[color_mapping[topology]]*len(data),
            ),
            row=1, col=1,
        )
        fig.add_trace(
            go.Bar(
                x=data['Window'],
                y=[1]*len(data),
                name=topology,
                legendgroup=topology,
                showlegend=False,
                marker_color=color_mapping[topology],
                marker_line_width=0,
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
            itemsizing='constant'
        ),
        hovermode="x unified",
        font=dict(family=font_family,),
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
    font_family,
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
            # marker_size=int(100/len(grouped_topology_df)),
            marker_line_width=1,
            marker_color=[color_mapping[topology]]*len(data),
            showlegend=True,
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
        font=dict(family=font_family,),
    )
    fig.update_xaxes(
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
        ticklen=0,
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
    font_family,
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
        font=dict(family=font_family,),
    )
    fig.update_xaxes(
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
    font_family,
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
            font_family,
            False,
        )
    else:
        alt_data_graph_data = make_alt_data_int_figure(
            alt_data_to_graph,
            color_mapping,
            topology_df,
            chromosome_df,
            template,
            dataRange,
            y_max,
            axis_line_width,
            xaxis_gridlines,
            yaxis_gridlines,
            font_family,
            False,
        )
    return alt_data_graph_data


def build_gff_figure(
    gff_df,
    min_range,
    max_range,
    template="plotly_dark",
    axis_line_width=1,
    xaxis_gridlines=False,
    yaxis_gridlines=False,
    font_family="Arial",
):
    """
    Supported features: gene, mRNA, CDS, exon
    """
    def gff_groups(gff_df):
        group_list=[]
        prev=None
        idx=0
        if len(gff_df.feature.unique()) == 1:
            [group_list.append(i) for i in range(len(gff_df))]
        else:
            for i in gff_df.itertuples():
                if not prev:
                    prev = i.feature
                    group_list.append(idx)
                elif prev != i.feature:
                    idx += 1
                    prev = i.feature
                    group_list.append(idx)
                else:
                    group_list.append(idx)
        return group_list


    def category_order_sorting(feature_list):
        built_in = ['gene', 'mRNA', 'exon', 'CDS']
        for i in built_in:
            try:
                old_index = feature_list.index(i)
                feature_list.insert(built_in.index(i), feature_list.pop(old_index))
            except ValueError:
                continue
        reversed(feature_list)
        return feature_list


    def create_marker_symbol_list(x_values, feature, strand):
        if strand == '+':
            if feature == 'gene':
                marker_symbol = ['square']*(len(x_values)-1) + ['triangle-right']
            elif (feature == 'exon') or (feature == 'CDS'):
                marker_symbol = []
                markers = {0: "triangle-right", 1: "triangle-left"}
                for n, _ in enumerate(x_values, 1):
                    marker_symbol.append(markers[0]) if n % 2 != 0 else marker_symbol.append(markers[1])
            else:
                marker_symbol = ['square']*(len(x_values)-1) + ['triangle-right']
            pass
        else:
            if feature == 'gene':
                marker_symbol = ['triangle-left'] + ['square']*(len(x_values)-1)
            elif (feature == 'exon') or (feature == 'CDS'):
                marker_symbol = []
                markers = {0: "triangle-right", 1: "triangle-left"}
                for n, _ in enumerate(x_values, 1):
                    marker_symbol.append(markers[0]) if n % 2 != 0 else marker_symbol.append(markers[1])
            else:
                marker_symbol = ['triangle-left'] + ['square']*(len(x_values)-1)
            pass
        return marker_symbol


    def get_attributes(data):
        attribute_list=[]
        genes=[]
        for attr in data.attribute:
            for i in attr.split(";"):
                # Split depending on GFF or GTF
                if '"' in i:
                    info = i.replace(" ", "").split('"')
                else:
                    info = i.replace(" ", "").split('=')
                
                if len(info) < 2:
                    continue
                else:
                    attribute_list.append(f"{info[0]} = {info[1]}<br>")
                    if (info[0] == 'gene') or (info[0] == 'ID'):
                        genes.append(info[1])
                    continue
        return attribute_list, genes
    

    def chunker(seq, size):
        return (seq[pos:pos + size] for pos in range(0, len(seq), size))

    features_graphed = list()
    features_order=list()
    fig = go.Figure()
    gff_df.insert(0, 'groups', gff_groups(gff_df))
    for _, data in gff_df.groupby(by=['groups']):
        feature = data.feature.unique()[0]
        strand = data.strand.unique()[0]
        attribute_list, genes = get_attributes(data)
        x_values = sorted(data['start'].to_list() + data['end'].to_list())
        # Set legend show if feature in list already
        features_graphed.append(feature) if feature not in features_graphed else None
        features_order.append(feature) if feature not in features_graphed else None
        y_values = [feature]*len(x_values)
        marker_symbol = create_marker_symbol_list(x_values, feature, strand)
        color_value = 'red' if strand == '+' else '#009BFF'
        fig.add_trace(go.Scatter(
            x=x_values if feature.lower() not in ['exon', 'cds'] else [min(x_values), max(x_values)],
            y=y_values if feature.lower() not in ['exon', 'cds'] else [feature, feature],
            name=feature,
            legendgroup=feature,
            mode='markers+lines' if feature.lower() not in ['exon', 'cds'] else 'lines',
            marker_symbol=marker_symbol,
            marker_size=8,
            marker_color=color_value,
            text="".join(attribute_list) if feature.lower() not in ['exon', 'cds'] else None,
            hovertemplate= None if feature.lower() in ['exon', 'cds'] else "",
            # text="".join(attribute_list),
            textfont=dict(size=10),
            showlegend=False,
        ))
        if feature.lower() in ['exon', 'cds']:
            for group, attr in zip(chunker(x_values, 2), chunker(attribute_list, 2)):
                fig.add_trace(
                    go.Scatter(
                        x=[group[0], group[1]], 
                        y=[feature, feature], 
                        fill="toself",
                        fillcolor=color_value,
                        marker_symbol='line-ns',
                        mode='markers+lines',
                        name=feature,
                        legendgroup=feature,
                        line_color=color_value,
                        line_width=10,
                        customdata=["".join(attr), "".join(attr)],
                        hovertemplate='%{customdata}',
                        textfont=dict(size=10),
                        showlegend=False,
                    )
                )
        # ----------------------------------------------
        # TODO: Resolve how to include gene names without
        # exponentially increasing load time. 
        # import statistics
        # if (genes) and (feature == 'gene'):
        #     fig.add_annotation(
        #         x=statistics.median(x_values),
        #         y=y_values[0],
        #         text=genes[0],
        #         showarrow=False,
        #         yshift=-15,
        #         # visible=False,
        #     )
    fig.update_layout(
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            traceorder='normal',
        ),
        margin=dict(t=40, b=10),
        template=template,
        title='',
        title_x=0.5,
        height=100*len(features_graphed),
        font=dict(family=font_family),
        hovermode='closest',
    )
    fig.update_xaxes(
        range=[min_range, max_range],
        title='Position',
        matches="x",
        rangemode="tozero",
        linewidth=axis_line_width,
        showgrid=xaxis_gridlines,
    )
    fig.update_yaxes(
        fixedrange=True,
        showgrid=yaxis_gridlines,
        title='',
        linewidth=axis_line_width,
        categoryorder='array',
        categoryarray=category_order_sorting(list(gff_df.feature.unique())),
    )
    return fig


# ----------------------------------------------------------------------------------------
# ------------------------------- Quantile Graph Functions -------------------------------
def get_quantile_coordinates(
    chromLengths,
    QUANTILES,
    WINDOWSIZE,
):
    quantileCoordinates = pd.DataFrame(columns=chromLengths["Chromosome"], index=range(1, QUANTILES+1))
    for row in chromLengths.itertuples(index=False):
        chrom, _, end = row
        chunkSize = end // QUANTILES
        for i in range(QUANTILES):
            q = i + 1
            if q == 1:
                quantileCoordinates.at[q, chrom] = [0, chunkSize]
            else:
                quantileCoordinates.at[q, chrom] = [chunkSize*(q-1) + WINDOWSIZE, chunkSize*q]
    return quantileCoordinates


def calculateFrequencies(
    quantileCoordinates,
    input_df,
    chromLengths,
    QUANTILES,
):
    quantileFrequencies = pd.DataFrame(columns=chromLengths["Chromosome"], index=range(1, QUANTILES+1))
    topos = input_df["TopologyID"].unique()
    for chrom in quantileCoordinates.columns:
        for q, quantile in enumerate(quantileCoordinates[chrom], 1):
            quantileData = input_df[(input_df['Window'] >= quantile[0]) & (input_df['Window'] <= quantile[1]) & (input_df['Chromosome'] == chrom)]
            topoQD = quantileData['TopologyID'].value_counts().to_dict()
            # Add missing topologies as count=0
            for i in topos:
                if i not in topoQD.keys():
                    topoQD[i] = 0
            quantileFrequencies.at[q, chrom] = topoQD
            continue
    return quantileFrequencies


def plot_frequencies(
    quantileFrequencies,
    n_quantiles,
    template,
    color_mapping,
    axis_line_width,
    xaxis_gridlines,
    yaxis_gridlines,
):
    def reorganizeDF(df):
        new_df = pd.DataFrame(columns=['Chr', 'Quantile', 'TopologyID', 'Frequency'])
        nidx = 0
        for c in df.columns:
            for idx in df.index:
                chromTotal = sum([v for v in df.at[idx, c].values()])
                for topo, freq in zip(df.at[idx, c].keys(), df.at[idx, c].values()):
                    new_df.at[nidx, 'TopologyID'] = topo
                    new_df.at[nidx, 'Chr'] = c
                    new_df.at[nidx, 'Quantile'] = idx
                    try:
                        new_df.at[nidx, 'Frequency'] = int(freq)/chromTotal
                    except ZeroDivisionError:
                        new_df.at[nidx, 'Frequency'] = 0.0
                    nidx += 1
        return new_df
    # Organize DataFrame
    organizedDF= reorganizeDF(quantileFrequencies)
    # Create line graph
    fig = px.line(
        organizedDF,
        x='Quantile',
        y='Frequency',
        color='TopologyID',
        facet_col='Chr',
        facet_col_wrap=1,
        facet_row_spacing=0.01,
        color_discrete_map=color_mapping,
    )
    fig.update_traces(texttemplate='%{text:.3}', textposition='top center')
    if len(organizedDF["Chr"].unique()) == 1:
        fig.update_layout(
            uniformtext_minsize=12,
            template=template,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="left",
                x=0,
                traceorder='normal',
            ),
            height=300,
        )
    else:
        fig.update_layout(
            uniformtext_minsize=12,
            template=template,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="left",
                x=0,
                traceorder='normal',
            ),
            height=100*len(organizedDF["Chr"].unique()),
        )
    fig.update_xaxes(
        range=[1, n_quantiles],
        rangemode="tozero",
        linewidth=axis_line_width,
        showgrid=xaxis_gridlines,
    )
    fig.update_yaxes(
        range=[0, 1],
        fixedrange=True,
        showgrid=yaxis_gridlines,
        linewidth=axis_line_width,
    )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    return fig


def calculate_topo_quantile_frequencies(df, current_topologies, additional_data, n_quantiles):
    final_df = pd.DataFrame(columns=["TopologyID", "Frequency", "Quantile"])
    for topology in current_topologies:
        topo_df = pd.DataFrame(columns=["TopologyID", "Frequency", "Quantile"])
        tidx = 0
        df = df.sort_values(by=additional_data)
        df = df.assign(Quantile = pd.qcut(df[additional_data].rank(method='first'), q=n_quantiles, labels=False))
        df['Quantile'] = df['Quantile'].apply(lambda x: x+1)
        df_group = df.groupby(by="Quantile")
        for rank, data in df_group:
            counts = data["TopologyID"].value_counts()
            for t, f in zip(counts.index, counts):
                if t == topology:
                    topo_df.at[tidx, "TopologyID"] = t
                    topo_df.at[tidx, "Frequency"] = (f/len(data))*100
                    topo_df.at[tidx, "Quantile"] = rank
                    tidx += 1
                    break
                else:
                    continue
            # -- Concat dfs -- 
        final_df = pd.concat([final_df, topo_df])
    return final_df


def plot_frequencies_topo_quantile(
    final_df,
    template,
    color_mapping,
    axis_line_width,
    xaxis_gridlines,
    yaxis_gridlines,
    graph_title,
    additional_data
):
    fig = px.line(
        final_df,
        x="Quantile", y="Frequency",
        color="TopologyID",
        color_discrete_map=color_mapping,
        markers=True,
    )
    fig.update_layout(
        template=template,
        title=graph_title,
        title_x=0.5,
        margin=dict(
            t=80
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            # itemsizing='constant'
        ),
    )
    fig.update_xaxes(
        title=f"{additional_data} Quantiles",
        linewidth=axis_line_width,
        showgrid=xaxis_gridlines,
        tick0=0,
        dtick=1,
    )
    fig.update_yaxes(
        rangemode="tozero",
        linewidth=axis_line_width,
        showgrid=yaxis_gridlines,
        title='% Windows Observed',        
    )
    return fig

# ---------------------------------------------------------------------------------
# -------------------------------- Whole Genome Graph Functions -------------------------------
def build_whole_genome_alt_data_graph(
    alt_data_to_graph,
    chromosome_df,
    color_mapping,
    topology_df,
    window_size,
    template,
    y_max,
    axis_line_width,
    xaxis_gridlines,
    yaxis_gridlines,
    font_family,
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
            None,
            axis_line_width,
            xaxis_gridlines,
            yaxis_gridlines,
            font_family,
            True,
        )
    else:
        alt_data_graph_data = make_alt_data_int_figure(
            alt_data_to_graph,
            color_mapping,
            topology_df,
            chromosome_df,
            template,
            None,
            y_max,
            axis_line_width,
            xaxis_gridlines,
            yaxis_gridlines,
            font_family,
            True,
        )
    return alt_data_graph_data


def build_topology_frequency_pie_chart(
    df,
    template,
    color_mapping,
    font_family,
):
    """Returns pie graph for whole genome topology frequencies"""
    fig = px.pie(
        df,
        values='Frequency',
        names='TopologyID',
        color="TopologyID",
        color_discrete_map=color_mapping,
        template=template,
        title='Whole Genome Topology Frequencies',
    )
    fig.update_traces(textposition='inside')
    fig.update_layout(
        margin=dict(l=120, r=20, t=40, b=10),
        uniformtext_minsize=12,
        uniformtext_mode='hide',
        legend=dict(itemclick=False, itemdoubleclick=False),
        title_x=0.5,
        font=dict(family=font_family,),
    )
    return fig


def build_rf_graph(
    df,
    ref_topo,
    template,
    color_mapping,
    axis_line_width,
    font_family,
):
    fig = px.bar(
        df, x="TopologyID", y="normRF-Distance",
        color="TopologyID", color_discrete_map=color_mapping,
        text='normRF-Distance')
    fig.update_traces(texttemplate='%{text:.2f}', textposition='inside')
    fig.update_layout(
        title=f"Normalized RF-Distance from {ref_topo}",
        title_x=0.5,
        template=template,
        font=dict(family=font_family,),
    )
    fig.update_xaxes(linewidth=axis_line_width)
    fig.update_yaxes(linewidth=axis_line_width, range=[0, 1])
    return fig


def build_whole_genome_rug_plot(
    df,
    chrom_df,
    chrom_order,
    chromGroup,
    template,
    color_mapping,
    currTopologies,
    topoOrder,
    window_size,
    axis_line_width,
    xaxis_gridlines,
    yaxis_gridlines,
    wg_squish_expand,
    font_family,
):
    grouped_topology_df = df.groupby(by='TopologyID')
    num_chroms = len(df['Chromosome'].unique())
    chrom_row_dict = {chrom:i for chrom, i in zip(chrom_order, range(1, len(df['Chromosome'].unique())+1, 1))}
    chrom_shapes = []
    row_height = [1]*num_chroms
    # --- Build figure ---
    if df.Chromosome.map(len).max() > 5:
        fig = make_subplots(
            rows=num_chroms,
            subplot_titles=[c for c in chrom_row_dict.keys()],
            shared_xaxes=True,
            cols=1,
            row_heights=row_height,
        )
    else:
        fig = make_subplots(
            rows=num_chroms,
            row_titles=[c for c in chrom_row_dict.keys()],
            shared_xaxes=True,
            cols=1,
            row_heights=row_height,
        )
    for topology, data in grouped_topology_df:
        add_legend = True
        for chrom in chrom_row_dict.keys():
            chrom_data = data[data["Chromosome"] == chrom]
            chrom_length_data = chrom_df[chrom_df['Chromosome'] == chrom]
            chrom_length = chrom_length_data['End'].max()
            if len(chrom_data) == 0:
                fig.add_trace(
                    go.Scatter(
                        x=[0],
                        y=[topology],
                        name=topology,
                        legendgroup=topology,
                        mode='markers',
                        marker_symbol='line-ns-open',
                        marker_color=[color_mapping[topology]]*len(chrom_data),
                        showlegend = False,
                    ),
                    row=chrom_row_dict[chrom], col=1,
                )
            elif add_legend:
                fig.add_trace(
                    go.Scatter(
                        x=chrom_data['Window'],
                        y=chrom_data['TopologyID'],
                        name=topology,
                        legendgroup=topology,
                        mode='markers',
                        # marker_size=int(25/len(grouped_topology_df)),
                        marker_symbol='line-ns-open',
                        marker_color=[color_mapping[topology]]*len(chrom_data),
                    ),
                    # go.Box(
                    #     x=chrom_data['Window'],
                    #     y=chrom_data['TopologyID'],
                    #     boxpoints='all',
                    #     jitter=0,
                    #     legendgroup=topology,
                    #     marker_symbol='line-ns-open',
                    #     marker_color=color_mapping[topology],
                    #     name=topology,
                    # ),
                    row=chrom_row_dict[chrom], col=1,
                )
                chrom_shapes.append(dict(type="line", xref="x", yref="y", x0=chrom_length, x1=chrom_length, y0=-1, y1=len(currTopologies), line_width=2))
                add_legend = False
            else:
                fig.add_trace(
                    go.Scatter(
                        x=chrom_data['Window'],
                        y=chrom_data['TopologyID'],
                        name=topology,
                        legendgroup=topology,
                        mode='markers',
                        # marker_size=int(25/len(grouped_topology_df)),
                        marker_symbol='line-ns-open',
                        marker_color=[color_mapping[topology]]*len(chrom_data),
                        showlegend = False,
                    ),
                    # go.Box(
                    #     x=chrom_data['Window'],
                    #     y=chrom_data['TopologyID'],
                    #     boxpoints='all',
                    #     jitter=0,
                    #     marker_symbol='line-ns-open',
                    #     marker_color=color_mapping[topology],
                    #     legendgroup=topology,
                    #     showlegend = False,
                    #     name=topology,
                    # ),
                    row=chrom_row_dict[chrom], col=1,
                )
                chrom_ref = chrom_row_dict[chrom]
                chrom_shapes.append(dict(type="rect", xref=f"x{chrom_ref}", yref=f"y{chrom_ref}", x0=chrom_length, x1=chrom_length, y0=-1, y1=len(currTopologies), line_width=2))

    # Update layout + axes
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    if wg_squish_expand == 'expand':
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
            margin=dict(
                t=20,
                b=30,
            ),
            height=100+(125*num_chroms),
            shapes=chrom_shapes,
            title_x=0.5,
            font=dict(family=font_family,),
        )
    elif wg_squish_expand == 'squish':
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
            margin=dict(
                t=20,
                b=30,
            ),
            height=100+(100*num_chroms),
            shapes=chrom_shapes,
            title_x=0.5,
            font=dict(family=font_family,),
        )
    else:
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
            # height=25*num_chroms,
            height=40+(20*num_chroms),
            shapes=chrom_shapes,
            title_x=0.5,
            margin=dict(
                t=20,
                b=30,
            ),
            font=dict(family=font_family,),
        )
    fig.update_xaxes(
        rangemode="tozero",
        range=[0, chrom_df['End'].max()],
        # fixedrange=True,
        linewidth=axis_line_width,
        matches=None,
        showgrid=xaxis_gridlines,
    )
    fig.update_yaxes(
        fixedrange=True,
        range=[-1, len(topoOrder)],
        title="",
        showgrid=yaxis_gridlines,
        showticklabels=False,
        ticklen=0,
        linewidth=axis_line_width,
        categoryarray=topoOrder,
    )
    # Rotate chromosome names to 0-degrees
    for annotation in fig['layout']['annotations']: 
        annotation['textangle']=0
        annotation['align']="center"
    return fig


def build_whole_genome_tile_plot(
    df,
    chrom_df,
    chrom_order,
    template,
    color_mapping,
    currTopologies,
    topoOrder,
    window_size,
    axis_line_width,
    chromGroup,
    xaxis_gridlines,
    yaxis_gridlines,
    wg_squish_expand,
    font_family,
):
    """
    Max chromosomes per graph if # current_topologies <= 3: 20
    Max chromosomes per graph if # current_topologies > 3: 20/2

    Returns: List of figures to display
    """
    grouped_topology_df = df.groupby(by='TopologyID')
    num_chroms = len(chrom_order)
    chrom_row_dict = {chrom:i for chrom, i in zip(chrom_order, range(1, len(chrom_order)+1, 1))}
    chrom_shapes = []
    # --- Build figure ---
    # If longest chromosome name longer 
    # than 5 characters, use subplot titles 
    # instead of row titles
    if df.Chromosome.map(len).max() > 5:
        fig = make_subplots(
            rows=num_chroms,
            cols=1,
            shared_xaxes=True,
            subplot_titles=chrom_row_dict.keys(),
            vertical_spacing=0.03,
        )
    else:
        fig = make_subplots(
            rows=num_chroms,
            cols=1,
            shared_xaxes=True,
            row_titles=[c for c in chrom_row_dict.keys()],
            vertical_spacing=0.001,
        )
    for topology, data in grouped_topology_df:
        add_legend = True
        for chrom in chrom_row_dict.keys():
            chrom_data = data[data["Chromosome"] == chrom]
            chrom_length_data = chrom_df[chrom_df['Chromosome'] == chrom]
            chrom_length = chrom_length_data['End'].max()
            if add_legend:
                fig.add_trace(
                    go.Histogram(
                        x=chrom_data['Window'],
                        y=[1]*len(chrom_data),
                        nbinsx=int(chrom_length/window_size),
                        name=topology,
                        legendgroup=topology,
                        marker_line_width=0,
                        marker_color=color_mapping[topology],
                    ),
                    row=chrom_row_dict[chrom], col=1,
                )
                chrom_shapes.append(dict(type="line", xref="x", yref="y", x0=chrom_length, x1=chrom_length, y0=0, y1=1, line_width=2))
                add_legend = False
            else:
                fig.add_trace(
                    go.Histogram(
                        x=chrom_data['Window'],
                        y=[1]*len(chrom_data),
                        nbinsx=int(chrom_length/window_size),
                        name=topology,
                        legendgroup=topology,
                        marker_line_width=0,
                        marker_color=color_mapping[topology],
                        showlegend = False
                    ),
                    row=chrom_row_dict[chrom], col=1,
                )
                chrom_ref = chrom_row_dict[chrom]
                chrom_shapes.append(dict(type="rect", xref=f"x{chrom_ref}", yref=f"y{chrom_ref}", x0=chrom_length, x1=chrom_length, y0=0, y1=1, line_width=2))

    # Update layout + axes
    if wg_squish_expand == 'expand':
        if num_chroms < 5:
                fig.update_layout(
                    barmode="relative",
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
                    height=125*num_chroms,
                    shapes=chrom_shapes,
                    title_x=0.5,
                    font=dict(family=font_family,),
                )
        else:
            fig.update_layout(
                barmode="relative",
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
                height=100*num_chroms,
                shapes=chrom_shapes,
                title_x=0.5,
                font=dict(family=font_family,),
            )
    elif wg_squish_expand == 'squish':
        if num_chroms < 5:
            fig.update_layout(
                barmode="relative",
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
                # height=75*num_chroms,
                height=200+(75*num_chroms),
                shapes=chrom_shapes,
                title_x=0.5,
                font=dict(family=font_family,),
            )
        else:
            fig.update_layout(
                barmode="relative",
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
                height=50*num_chroms,
                shapes=chrom_shapes,
                title_x=0.5,
                font=dict(family=font_family,),
            )
    else:
        fig.update_layout(
            barmode="relative",
            template=template,
            legend_title_text='Topology',
            margin=dict(
                t=20,
                b=30,
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
            height=40+(20*num_chroms),
            shapes=chrom_shapes,
            title_x=0.5,
            font=dict(family=font_family,),
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
        range=[0, 1],
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


def build_whole_genome_bar_plot(
    df,
    chrom_order,
    chromGroup,
    template,
    color_mapping,
    currTopologies,
    axis_line_width,
    xaxis_gridlines,
    yaxis_gridlines,
    font_family,
):
    # Filter df to chromosomes in group
    df = df[df['Chromosome'].isin(chromGroup)]
    df = df[df['TopologyID'].isin(currTopologies)]
    number_of_chrom_rows = len(df["Chromosome"].unique()) // 3
    fig = px.bar(
        df,
        x='TopologyID',
        y='Frequency',
        category_orders={"Chromosome": chrom_order},
        facet_col='Chromosome',
        facet_col_wrap=3,
        facet_row_spacing=0.025,
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
        margin=dict(l=40, r=40, t=40, b=40),
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
        font=dict(family=font_family,),
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


def build_whole_genome_pie_charts(
    df,
    chrom_order,
    template,
    color_mapping,
    chromGroup,
    font_family,
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
        horizontal_spacing=0.0001,
        subplot_titles=chrom_order,
        column_widths=[2]*3,
    )
    col_pos = 1
    row_num = 1
    for c in chrom_order:
        chrom_df = df[df["Chromosome"] == c]
        colors = [color_mapping[i] for i in chrom_df["TopologyID"]]
        fig.add_trace(go.Pie(
            labels=chrom_df["TopologyID"],
            values=chrom_df['Frequency'],
            marker_colors=colors,
            ),
            row=row_num,
            col=col_pos,  
        )
        if col_pos == 3:
            col_pos = 1
            row_num += 1
        else:
            col_pos += 1
    fig.update_traces(textposition='inside')
    fig.update_layout(
        uniformtext_minsize=12, 
        showlegend=True,
        template=template,
        height=int(200*number_of_chrom_rows),
        font=dict(family=font_family,),
    )
    return fig


def get_recommended_chrom_num(df_len):
    """Divide len of DF by 500k"""
    return int(500000 // df_len)

# ---------------------------------------------------------------------------------
# --------------------------- Stats DataFrame Generators --------------------------
def _get_valid_cols(topology_df):
        valid_cols = list()
        for i in topology_df.columns[4:]:
            data = topology_df[i].unique()
            flag = None
            for j in data:
                if type(j) == str:
                    flag = False
                    break
                else:
                    flag = True
            if flag:
                valid_cols.append(i)
            else:
                continue
        return valid_cols


def basic_stats_dfs(topology_df):
    """Generate dataframes of basic statistics

    :param topology_df: Current View Tree Viewer input file dataframe
    :type topology_df: Object
    """
    # Calculate current view topologies
    topo_freq_df = pd.DataFrame(topology_df["TopologyID"].value_counts()/len(topology_df))
    if len(topo_freq_df) > 25:  # If more than 25 topologies loaded, just show top 25
        topo_freq_df = topo_freq_df.head(25)
        remainder_freq = 1.0 - sum(topo_freq_df['TopologyID'])
        topo_freq_df.at["Other", "TopologyID"] = remainder_freq
    topo_names = [i for i in topo_freq_df.index]
    topo_freqs = [round(i, 4) for i in topo_freq_df["TopologyID"]]
    # Calculate median + average of additional data
    if len(topology_df.columns) > 4:
        valid_cols = _get_valid_cols(topology_df)
        additional_dt_names = [i for i in valid_cols]
        additional_dt_avg = [topology_df[i].mean() for i in valid_cols]
        additional_dt_std = [topology_df[i].std() for i in valid_cols]
        topo_freq_df = pd.DataFrame(
            {
                "TopologyID": topo_names,
                "Frequency": topo_freqs,
            }
        )
        additional_data_df = pd.DataFrame(
            {
                "Additional Data": additional_dt_names,
                "Average": additional_dt_avg,
                "Std Dev": additional_dt_std,
            }
        )
        return topo_freq_df, additional_data_df
    else:  # No additional data types present in file
        topo_freq_df = pd.DataFrame(
            {
                "TopologyID": topo_names,
                "Frequency": topo_freqs,
            }
        )
        return topo_freq_df, pd.DataFrame()


def current_view_topo_freq_chart(basic_stats_topo_freqs, template, color_mapping, current_view_title):
    """Return pie chart figure object for local topology frequencies

    :param basic_stats_topo_freqs: Dataframe of topology frequencies
    :type basic_stats_topo_freqs: DataFrame
    :return: Plotly express pie chart
    :rtype: Figure object
    """
    if "Other" in basic_stats_topo_freqs["TopologyID"].to_list():
        fig = px.bar(
            basic_stats_topo_freqs,
            x='TopologyID',
            y="Frequency",
            color="TopologyID",
            color_discrete_map=color_mapping,
            text="Frequency",
        )
        fig.update_layout(
            title=current_view_title,
            title_x=0.5,
            template=template,
            uniformtext_minsize=12, 
            uniformtext_mode='hide',
        )
        fig.update_traces(textposition='outside')
        return fig
    else:
        fig = px.pie(
            basic_stats_topo_freqs,
            values="Frequency",
            names="TopologyID",
            color="TopologyID",
            color_discrete_map=color_mapping,
            template=template,
            title="Current View Topology Frequencies",
        )
        fig.update_layout(
            legend=dict(itemclick=False, itemdoubleclick=False),
            margin=dict(l=120, r=20, t=40, b=10),
            uniformtext_minsize=12, 
            uniformtext_mode='hide',
            title_x=0.5,
        )
        fig.update_traces(textposition='inside')
        
        return fig


def whole_genome_datatable(tv_df):
    valid_cols = _get_valid_cols(tv_df[4:])
    for i in tv_df.columns.to_list()[4:]:
        if i in valid_cols:
            continue
        else:
            tv_df.drop(labels=i, axis=1, inplace=True)
    df_group = tv_df.groupby(by="TopologyID")
    out_df = pd.DataFrame(columns=["TopologyID", "Additional Data", "Num. Windows", "Average", "Std Dev"])
    idx = 0
    for topology, data in df_group:
        additional_datatypes = [i for i in data.columns[4:]]
        for datatype in additional_datatypes:
            dt_data = data[datatype]
            mean = dt_data.mean()
            stdev = dt_data.std()
            out_df.at[idx, "TopologyID"] = topology
            out_df.at[idx, "Additional Data"] = datatype
            out_df.at[idx, "Num. Windows"] = len(dt_data)
            out_df.at[idx, "Average"] = mean
            out_df.at[idx, "Std Dev"] = stdev
            idx += 1
            continue
    columns = [{'id': c, 'name': ["Per-Topology Whole Genome Comparison", c], 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal)} for c in out_df.columns]
    data = out_df.to_dict('records')
    return data, columns

# --- post-hoc tests ---
def mann_whitney_posthoc(tv_df, additional_data_type, pval_adjustment):
    return sp.posthoc_mannwhitney(tv_df, val_col=additional_data_type, group_col='TopologyID', p_adjust=pval_adjustment)


def dunns_test_posthoc(tv_df, additional_data_type, pval_adjustment):
    return sp.posthoc_dunn(tv_df, val_col=additional_data_type, group_col='TopologyID', p_adjust=pval_adjustment)


def tukeyHSD_posthoc(tv_df, additional_data_type, pval_adjustment, alpha):
    return sp.posthoc_tukey_hsd(tv_df[additional_data_type], tv_df["TopologyID"], alpha=alpha)

# --- Significance tests ---
def kruskal_wallis_H_test(tv_df, additional_data_type, posthoc_type, pval_adjustment, alpha):
    """Return dataframe with Kruskal-Wallis H test information for each topology
    """
    d = [tv_df.loc[ids, additional_data_type].values for ids in tv_df.groupby('TopologyID').groups.values()]
    H, p = ss.kruskal(*d, nan_policy='omit')
    if posthoc_type == "Mann-Whitney rank test":
        posthoc = mann_whitney_posthoc(tv_df, additional_data_type, pval_adjustment)
        posthoc_df = pd.DataFrame(columns=[posthoc_type, "p-value"])
        idx = 0
        for c1 in posthoc.columns:
            for c2, pval in zip(posthoc.index, posthoc[c1]):
                if c1 == c2:  # Remove self-self comparisons
                    continue
                posthoc_df.at[idx, posthoc_type] = f"{c1} vs {c2}"
                posthoc_df.at[idx, "p-value"] = float(pval)
                idx += 1
        data = posthoc_df.to_dict('records')
        columns = [
            {'id': posthoc_type, 'name': posthoc_type},
            {'id': 'p-value', 'name': 'p-value', 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal_or_exponent)},
        ]
    elif posthoc_type == "Dunn's test":
        posthoc = dunns_test_posthoc(tv_df, additional_data_type, pval_adjustment)
        posthoc_df = pd.DataFrame(columns=[posthoc_type, "p-value"])
        idx = 0
        for c1 in posthoc.columns:
            for c2, pval in zip(posthoc.index, posthoc[c1]):
                if c1 == c2:  # Remove self-self comparisons
                    continue
                posthoc_df.at[idx, posthoc_type] = f"{c1} vs {c2}"
                posthoc_df.at[idx, "p-value"] = float(pval)
                idx += 1
        data = posthoc_df.to_dict('records')
        columns = [
            {'id': posthoc_type, 'name': posthoc_type},
            {'id': 'p-value', 'name': 'p-value', 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal_or_exponent)},
        ]
    elif posthoc_type == "TukeyHSD":
        posthoc = tukeyHSD_posthoc(tv_df, additional_data_type, pval_adjustment, alpha)
        posthoc_df = pd.DataFrame(columns=[posthoc_type, "p-value"])
        idx = 0
        for c1 in posthoc.columns:
            for c2, pval in zip(posthoc.index, posthoc[c1]):
                if c1 == c2:  # Remove self-self comparisons
                    continue
                posthoc_df.at[idx, posthoc_type] = f"{c1} vs {c2}"
                posthoc_df.at[idx, "p-value"] = float(pval)
                idx += 1
        data = posthoc_df.to_dict('records')
        columns = [
            {'id': posthoc_type, 'name': posthoc_type},
            {'id': 'p-value', 'name': 'p-value', 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal_or_exponent)},
        ]
    else:
        pass
    
    return posthoc, data, columns, H, p


def one_way_anova(tv_df, additional_data_type, posthoc_type, pval_adjustment, alpha):
    d = [tv_df.loc[ids, additional_data_type].values for ids in tv_df.groupby('TopologyID').groups.values()]
    F, p = ss.f_oneway(*d)
    if posthoc_type == "Mann-Whitney rank test":
        posthoc = mann_whitney_posthoc(tv_df, additional_data_type, pval_adjustment)
        posthoc_df = pd.DataFrame(columns=[posthoc_type, "p-value"])
        idx = 0
        for c1 in posthoc.columns:
            for c2, pval in zip(posthoc.index, posthoc[c1]):
                posthoc_df.at[idx, posthoc_type] = f"{c1} vs {c2}"
                posthoc_df.at[idx, "p-value"] = float(pval)
                idx += 1
        data = posthoc_df.to_dict('records')
        columns = [
            {'id': posthoc_type, 'name': posthoc_type},
            {'id': 'p-value', 'name': 'p-value', 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal_or_exponent)},
        ]
    elif posthoc_type == "Dunn's test":
        posthoc = dunns_test_posthoc(tv_df, additional_data_type, pval_adjustment)
        posthoc_df = pd.DataFrame(columns=[posthoc_type, "p-value"])
        idx = 0
        for c1 in posthoc.columns:
            for c2, pval in zip(posthoc.index, posthoc[c1]):
                posthoc_df.at[idx, posthoc_type] = f"{c1} vs {c2}"
                posthoc_df.at[idx, "p-value"] = float(pval)
                idx += 1
        data = posthoc_df.to_dict('records')
        columns = [
            {'id': posthoc_type, 'name': posthoc_type},
            {'id': 'p-value', 'name': 'p-value', 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal_or_exponent)},
        ]
    elif posthoc_type == "TukeyHSD":
        posthoc = tukeyHSD_posthoc(tv_df, additional_data_type, pval_adjustment, alpha)
        posthoc_df = pd.DataFrame(columns=[posthoc_type, "p-value"])
        idx = 0
        for c1 in posthoc.columns:
            for c2, pval in zip(posthoc.index, posthoc[c1]):
                posthoc_df.at[idx, posthoc_type] = f"{c1} vs {c2}"
                posthoc_df.at[idx, "p-value"] = float(pval)
                idx += 1
        data = posthoc_df.to_dict('records')
        columns = [
            {'id': posthoc_type, 'name': posthoc_type},
            {'id': 'p-value', 'name': 'p-value', 'type': 'numeric', 'format': Format(precision=4, scheme=Scheme.decimal_or_exponent)},
        ]
    else:
        pass
    
    return posthoc, data, columns, F, p


def stats_test_heatmap(posthoc, template):
    fig = go.Figure(data=go.Heatmap(
        z=posthoc.values,
        x=posthoc.columns,
        y=posthoc.index,
        zmin=0,
        zmax=1,
        colorscale='Viridis',
        colorbar=dict(title='p-value'),
        hovertemplate = 'p-value: %{z}<extra></extra>',
    ))
    fig.update_layout(
        template=template,
        coloraxis_colorbar=dict(title="log(p-value)"),
        margin=dict(
            t=60,
        ),
    )
    return fig


def frequency_distribution(data, name, template):
    """Return frequency density distribution"""
    fig = px.histogram(data, x=name, histnorm='density')
    fig.update_layout(template=template, margin=dict(t=20, pad=30))
    return fig


def mean_frequency_of_alt_data_per_topology(tv_df, topologies, additional_data_type):
    out_df = pd.DataFrame(columns=["TopologyID", "Total Windows", f"Mean ({additional_data_type})"])
    idx = 1
    for i in topologies:
        topo_df = tv_df[tv_df["TopologyID"] == i]
        additional_data_mean = topo_df[f"{additional_data_type}"].mean()
        out_df.at[idx, "TopologyID"] = i
        out_df.at[idx, "Total Windows"] = len(topo_df)
        out_df.at[idx, f"Mean ({additional_data_type})"] = additional_data_mean
        idx += 1
        continue
    return out_df.to_dict('records')

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
        return int(hoverdata['x'])
    else:
        return df.loc[hoverdata['binNumber']]['Window']


def get_Treexpos(hoverdata, df):
    hoverdata = hoverdata['points'][0]
    if ('customdata' in hoverdata.keys()) or ('marker.color' in hoverdata.keys()):
        return int(hoverdata['x'])
    else:
        return int(hoverdata['x'])

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
    This function returns a blank figure with a "NO DATA LOADED" watermark.
    """
    fig = go.Figure()
    fig.update_layout(
        template=template,
        # annotations=[
        #     dict(
        #         name="draft watermark",
        #         # text="NO DATA LOADED",
        #         textangle=0,
        #         opacity=0.9,
        #         font=dict(color="white", size=50),
        #         xref="paper",
        #         yref="paper",
        #         x=0.5,
        #         y=0.5,
        #         showarrow=False,
        #     )
        # ],
    )
    fig.update_xaxes(range=[0.2, 1], showgrid=False, visible=False, zeroline=False)
    fig.update_yaxes(range=[0.2, 1], showgrid=False, visible=False, zeroline=False)
    return fig


def init_stats_graph(template):
    """
    This function returns a blank figure with a "NO DATA" watermark.
    """
    fig = go.Figure()
    fig.update_layout(
        template=template,
        annotations=[
            dict(
                name="draft watermark",
                text="NO DATA",
                textangle=0,
                opacity=0.9,
                font=dict(color="white", size=35),
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
        height=300,
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
        height=300,
        template=template,
        annotations=[
            dict(
                name="draft watermark",
                text="Zoom in to view - max 5Mb",
                textangle=0,
                opacity=0.9,
                font=dict(color="white", size=25),
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
    # for c in chrom_names:
    #     if c not in tv_chrom_names:
    #         missing_chromosomes.append(c)
    #         valid = False
    #         issue_files.append("Chromosome Length File")
    #         continue
    #     else:
    #         continue
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
            missing_chroms = ", ".join(missing_chromosomes)
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
    if tree == "NoTree":
        return "NoTree"
    tree = Tree(tree)
    taxa = []
    for leaf in tree.iter_leaves():
        taxa.append(leaf.name)
    return sorted(taxa)


def get_valid_init_tree(trees):
    """Returns first NewickTree entry that is not NoTree"""
    for i in range(len(trees)):
        if trees[i] == "NoTree":
            continue
        else:
            return trees[i]


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


def valid_gff_gtf(i):
    """Ensure gff/gtf file has 9 columns of data

    :param i: gff/gtf file
    :type i: str
    :return: True or False
    :rtype: bool
    """
    gffdf = pd.read_csv(i, sep="\t", comment="#", nrows=1)
    try:
        assert len(gffdf.loc[0]) == 9
        return True
    except AssertionError:
        return False

# ---------------------------------------------------------------------------------
# --------------------------- Tree Prune Export Tools -----------------------------
def prune_tree(x, prune_taxa_choices):
    if x == "NoTree":
        return "NoTree"
    else:
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
    """Remove any information in brackets - ete3 
       does not support this format of newick"""
    # --- Ensure tree is NaN value, if so return NoTree ---
    if type(l) == float:
        return "NoTree"
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
        if t == "NoTree":
            continue
        elif len(topologies.keys()) == 0:
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
    for topology in topologies.keys():
        idx = topologies[topology]['idx']
        topoName = f'topology{topoCount}'
        for i in idx:
            df.at[i, 'TopologyID'] = topoName
            continue
        topoCount += 1
    return df


def mygrouper(n, iterable):
    args = [iter(iterable)] * n
    return list([e for e in t if e != None] for t in itertools.zip_longest(*args))


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


# ---------------------------------------------------------------------------------
# ----------------------------- Template Generaters -------------------------------
def project_ini_template():
    content = """[MAIN]\nProjectDir = /path/to/Project\nTreeViewerFile = /path/to/TreeViewerInput.xlsx\nChromLengths = /path/to/ChromosomeLengths.bed\n\n[ADDITIONAL]\n# Load multiple gff/gtf files by listing them with ";" separating the files\nGFF_GTF = None"""
    return content


def tree_viewer_template():
    content = pd.DataFrame(columns=["Chromosome", "Window", "NewickTree", "TopologyID"])
    return content


def chrom_len_template():
    content = pd.DataFrame({"Chromosome": ["chr1", "chr2", "chr3"], "Start": [0, 0, 0], "End": [1000000, 1500000, 2000000]})
    return content
# ---------------------------------------------------------------------------------
# ------------------------------- Misc. Functions ---------------------------------
def divide_input_into_cpu_size_chunks(l, n):
    """Divides chromosomes into sets of size n, where n
    is the number of cores available to use"""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def filter_numeric_dtypes(df):
    filtered_names = []
    for name, data_type in zip(df.dtypes.index[4:], df.dtypes[4:]):
        if str(data_type) == 'object':
            continue
        else:
            filtered_names.append(name)
    return filtered_names


def sorted_nicely(l): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
