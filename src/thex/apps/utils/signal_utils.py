from pathlib import Path

import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html

import plotly
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# -------------------- Graphing Functions --------------------
def single_chromosome_graph_line(
    df,
    chromosome,
    chosen_template,
    marker_width,
    colors,
    font_size,
    xaxis_gridlines,
    yaxis_gridlines,
    font_family,
    samples,
):
    """ Filter out current chromosome and set x- and y-max"""
    curr_chrom_data = df[df["Chromosome"] == chromosome]
    y_max = float(curr_chrom_data["Value"].max())
    fig = px.line(
        curr_chrom_data,
        x='Window',
        y='Value',
        category_orders={"Sample": samples},
        color='Sample',
        color_discrete_sequence=colors,
        height=500,
    )
    fig.update_layout(
        font=dict(
            size=font_size,
            family=font_family,
        ),
        legend=dict(
            itemsizing='trace',
            orientation="h",
            xanchor="left",
            x=0,
            y=1.02,
            yanchor="bottom",
        ),
        showlegend=True,
        template=chosen_template,
        title_x=0.5,
    )
    fig.update_xaxes(
        title="Position",
        rangemode='tozero',
        showgrid=xaxis_gridlines,
    )
    fig.update_yaxes(
        title="Value",
        range=[0, y_max],
        # fixedrange=True,
        showgrid=yaxis_gridlines,
    )
    fig.update_traces(
        line=dict(width=float(marker_width)),
    )
    return fig

def single_chromosome_graph_scatter(
    df,
    chromosome,
    chosen_template,
    marker_width,
    colors,
    font_size,
    xaxis_gridlines,
    yaxis_gridlines,
    font_family,
    samples,
):
    """ Filter out current chromosome and set x- and y-max"""
    curr_chrom_data = df[df["Chromosome"] == chromosome]
    y_max = float(curr_chrom_data["Value"].max())
    fig = px.scatter(
        curr_chrom_data,
        x='Window',
        y='Value',
        category_orders={"Sample": samples},
        color='Sample',
        color_discrete_sequence=colors,
        height=500,
    )
    fig.update_layout(
        font=dict(
            size=font_size,
            family=font_family,
        ),
        legend=dict(
            itemsizing='trace',
            orientation="h",
            xanchor="left",
            x=0,
            y=1.02,
            yanchor="bottom",
        ),
        showlegend=True,
        template=chosen_template,
        title_x=0.5,
    )
    fig.update_xaxes(
        title="Position",
        rangemode='tozero',
        showgrid=xaxis_gridlines,
    )
    fig.update_yaxes(
        title="Value",
        range=[0, y_max],
        # fixedrange=True,
        showgrid=yaxis_gridlines,
    )
    fig.update_traces(
        marker=dict(size=float(marker_width)),
    )
    return fig

def whole_genome_line(
    df,
    chromosomes,
    samples,
    colors,
    marker_width,
    template,
    font_size,
    y_max,
    x_max,
    xaxis_gridlines,
    yaxis_gridlines,
    font_family,
):
    fig = make_subplots(
        rows=len(chromosomes),
        cols=1,
        x_title="Position",
        y_title="value",
        row_titles=chromosomes,
        row_heights=[2]*len(chromosomes),
    )
    for n, sample in enumerate(samples):
        legend_flag = True
        for row, current_chromosome in enumerate(chromosomes, start=1):
            filt = (df['Chromosome'] == current_chromosome) & (df["Sample"] == sample)
            sample_chromosome_data = df[filt]
            # Make figure
            fig.add_trace(
                go.Scatter(
                    x=sample_chromosome_data['Window'],
                    y=sample_chromosome_data['Value'],
                    mode='lines',
                    legendgroup=str(sample),
                    name=sample,
                    line=dict(
                        color=colors[n],
                        width=float(marker_width)
                    ),
                    showlegend=legend_flag,
                ),
                row=row,
                col=1
            )
            legend_flag = False
            continue
    # --- Update Figure ---
    fig.update_layout(
        font=dict(size=font_size, family=font_family),
        height=125*len(chromosomes),
        hovermode="x unified",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            itemsizing='trace',
            title="",
        ),
        margin=dict(
            l=60,
            r=50,
            b=60,
            t=10,
        ),
        template=template,
        title_x=0.5,
        font_family="Arial",
    )
    fig.update_xaxes(
        fixedrange=True,
        range=[0, x_max],
        showgrid=xaxis_gridlines,
    )
    fig.update_yaxes(
        range=[0.0, y_max],
        fixedrange=True,
        showgrid=yaxis_gridlines,
    )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    # Rotate chromosome names to 0-degrees
    for annotation in fig['layout']['annotations']:
        if annotation['text'] == "value":
            continue
        annotation['textangle']=0
        annotation['align']="center"
    return fig


def whole_genome_scatter(
    df,
    chromosomes,
    samples,
    colors,
    marker_width,
    template,
    font_size,
    y_max,
    x_max,
    xaxis_gridlines,
    yaxis_gridlines,
    font_family,
):
    # fig = make_subplots(
    #     rows=len(chromosomes),
    #     cols=1,
    #     x_title="Position",
    #     y_title="Edit Me!",
    #     row_titles=chromosomes,
    #     row_heights=[2]*len(chromosomes),
    # )
    # for n, sample in enumerate(samples):
    #     legend_flag = True
    #     for row, current_chromosome in enumerate(chromosomes, start=1):
    #         filt = (df['Chromosome'] == current_chromosome) & (df["Sample"] == sample)
    #         sample_chromosome_data = df[filt]
    #         # Make figure
    #         fig.add_trace(
    #             go.Scatter(
    #                 x=sample_chromosome_data['Window'],
    #                 y=sample_chromosome_data['Value'],
    #                 mode='markers',
    #                 legendgroup=str(sample),
    #                 name=sample,
    #                 line=dict(
    #                     color=colors[n],
    #                     width=float(marker_width)
    #                 ),
    #                 showlegend=legend_flag,
    #             ),
    #             row=row,
    #             col=1
    #         )
    #         legend_flag = False
    #         continue

    fig = px.scatter(
        df,
        x='Window',
        y='Value',
        category_orders={"Sample": samples},
        color='Sample',
        color_discrete_sequence=colors,
        # height=500,
        facet_row="Chromosome",
    )
    # --- Update Figure ---
    fig.update_layout(
        font=dict(size=font_size, family=font_family),
        height=125*len(chromosomes),
        hovermode="x unified",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            itemsizing='trace',
            title="",
        ),
        margin=dict(
            l=60,
            r=50,
            b=60,
            t=10,
        ),
        template=template,
        title_x=0.5,
        font_family=font_family,
    )
    fig.update_xaxes(
        fixedrange=True,
        range=[0, x_max],
        showgrid=xaxis_gridlines,
    )
    fig.update_yaxes(
        range=[0.0, y_max],
        fixedrange=True,
        showgrid=yaxis_gridlines,
        title='',
    )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.update_traces(marker=dict(size=float(marker_width)))
    # Rotate chromosome names to 0-degrees
    for annotation in fig['layout']['annotations']:
        if annotation['text'] == "value":
            continue
        annotation['textangle']=0
        annotation['align']="center"
    return fig


# -------------------- File Validation --------------------

def validate_signal_tracer_headers(df):
    """Validate that headers are correct"""
    expected_headers = ["Chromosome", "Window", "Sample", "Value"]
    try:
        assert list(df.columns) == expected_headers
        return True
    except AssertionError:
        return False


def validate_signal_tracer_values(xlsx_df):
    """Return False if value column data are not int or float"""
    try:
        assert xlsx_df['Value'].dtype != "object"
        return True
    except AssertionError:
        return False


def validate_file_type(filename):
    """Return False if file type is not valid """
    valid_filetypes = ['.tsv', '.csv', '.xlsx', '.txt']
    filetype = Path(filename).suffix
    if filetype not in valid_filetypes:
        return False
    else:
        return True

