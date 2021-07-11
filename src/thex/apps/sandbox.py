# from pathlib import Path
# from io import StringIO
# import json

# import dash
# import dash_bootstrap_components as dbc
# import dash_core_components as dcc
# import dash_html_components as html
# from dash.dependencies import Input, Output, State
# from dash.exceptions import PreventUpdate
# import dash_cytoscape
# import dash_daq as daq
# import dash_bio as dashbio

# import plotly.express as px

# from Bio import Phylo

# import pandas as pd

# from apps import navbar, tree_utils, data_utils
# from app import app


# # annotations = pd.read_json('https://eweitz.github.io/ideogram/data/annotations/SRR562646.json')
# # print(annotations)

# LOCAL_ORGS = [{'label': i.stem, 'value': i.name} for i in Path("apps/app_data/tree_viewer/bands/").iterdir() if i.is_file()]


# def gff2DashAnnots(
#     gff_file,
# ):
#     """
#     This function will take in a GFF file and parse it out into JSON format for input to the ideogram.
#     This function will only pull gene information as of now. 
#     It pulls...
#         1. Gene Name
#         2. Chromosome
#         3. Start Position
#         4. Stop Position
#     """

#     def _get_gene_name(attribute):
#         attribute_split = attribute.split(";")
#         gene_name = [attr for attr in attribute_split if "gene=" in attr]
#         return gene_name[0].strip("gene=")


#     gff_cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

#     gff_path = Path(gff_file)

#     annotation_collect_list = []

#     """Read in GFF file through context manager to skip inital metadata lines"""
#     with open(gff_path, 'r') as fh:
#         while "#" in fh.readline():
#             pass
#         read_gff = pd.read_csv(fh, sep='\t', header=None)
#         read_gff.columns = gff_cols
#         read_gff = read_gff[read_gff['feature'] == "gene"]
#         read_gff.drop(labels=["source", "feature", "score", "strand", "frame"], axis=1, inplace=True)

#         for chrom, start, stop, attribute in zip(read_gff['seqname'], read_gff['start'], read_gff['end'], read_gff['attribute']):
#             # # May need to remove, this will skip any gene with "LOC" in it.
#             if "LOC" in str(_get_gene_name(attribute)):
#                 continue
#             annotation_dict=dict(name=str(_get_gene_name(attribute)), chr=str(chrom), start=int(start), stop=int(stop))
#             annotation_collect_list.append(annotation_dict)
#             continue
#     return annotation_collect_list


# TEST_GFF = "C:\\Users\\andrewharris.2374\\Desktop\\GCF_000181335.3_Felis_catus_9.0_genomic.gff"

# test_out = gff2DashAnnots(TEST_GFF)
# # print(test_out[0])

# body = dbc.Container(
#     children=[
#         dbc.Row(
#             dbc.Col(
#                 id="col1",
#                 children=[
#                     html.P("HELLO"),
#                     # dcc.Dropdown(
#                     #     id='local_orgs',
#                     #     options=LOCAL_ORGS,
#                     #     value=LOCAL_ORGS[0]['value']
#                     # ),
#                     dashbio.Ideogram(
#                         id='ideogram-annotations',
#                         organism="felis-catus",
#                         # dataDir=None,
#                         # localOrganism={
#                         #     "chrBands": [
#                         #         "A1 p 1 0 3713 0 89880064",
#                         #         "A1 q 1 3714 10000.0 89880064 242100913",
#                         #         "A2 p 1 0 2523 0 61080734",
#                         #         "A2 q 1 2524 7082.655941904689 61080734 171471747",
#                         #         "A3 p 1 0 2142 0 51856977",
#                         #         "A3 q 1 2143 5914.9882264177995 51856977 143202405",
#                         #         "B1 p 1 0 1623 0 39294480",
#                         #         "B1 q 1 1624 8600.25211883443 39294480 208212889",
#                         #         "B2 p 1 0 1246 0 30165856",
#                         #         "B2 q 1 1247 6414.789439476422 30165856 155302638",
#                         #         "B3 p 1 0 1225 0 29663499",
#                         #         "B3 q 1 1226 6185.5119480693575 29663499 149751809",
#                         #         "B4 p 1 0 1851 0 44808485",
#                         #         "B4 q 1 1852 5969.770754230902 44808485 144528695",
#                         #         "C1 p 1 0 4498 0 108907037",
#                         #         "C1 q 1 4499 9202.366865919254 108907037 222790142",
#                         #         "C2 p 1 0 3148 0 76225444",
#                         #         "C2 q 1 3149 6658.097567769189 76225444 161193150",
#                         #         "D1 p 1 0 1371 0 33184684",
#                         #         "D1 q 1 1372 4859.462384596625 33184684 117648028",
#                         #         "D2 p 1 0 860 0 20812855",
#                         #         "D2 q 1 861 3725.16810789557 20812855 90186660",
#                         #         "D3 p 1 0 1313 0 31780724",
#                         #         "D3 q 1 1314 4001.8108481895733 31780724 96884206",
#                         #         "D4 p 1 0 1343 0 32525998",
#                         #         "D4 q 1 1344 3986.835522590532 32525998 96521652",
#                         #         "E1 p 1 0 1068 0 25862966",
#                         #         "E1 q 1 1069 2622.6538435235066 25862966 63494689",
#                         #         "E2 p 1 0 1030 0 24929418",
#                         #         "E2 q 1 1031 2657.5816754561347 24929418 64340295",
#                         #         "E3 p 1 0 753 0 18233305",
#                         #         "E3 q 1 754 1844.2013888646509 18233305 44648284",
#                         #         "F1 p 1 0 41 0 1000001",
#                         #         "F1 q 1 42 2960.0980067349024 1000001 71664243",
#                         #         "F2 p 1 0 41 0 1000001",
#                         #         "F2 q 1 42 3542.012912607232 1000001 85752456",
#                         #         "X p 1 0 2128 0 51521357",
#                         #         "X q 1 2129 5392.669006580739 51521357 130557009"
#                         #     ]
#                         # },
#                         container="col1",
#                         # organism="danio-rerio",
#                         # annotationsLayout="track",
#                         annotationHeight=5,
#                         filterable=True,
#                         chromosomes=[1, 2],
#                         # annotations=test_out,
#                         # annotationsPath='https://eweitz.github.io/ideogram/data/annotations/SRR562646.json',
#                         style={"color": "black", "background": "white"},
#                         showBandLabels=True,
#                         # orientation='horizontal',
#                         # fullChromosomeLabels=True,
#                     )  
#                 ],
#                 width=12,
#             )
#         ),
#     ],
#     fluid=True
# )



# nav = navbar.Navbar()
# layout = html.Div([nav, body])


# # @app.callback(
# #     dash.dependencies.Output('ideogram-annotations', 'localOrganism'),
# #     [dash.dependencies.Input('local_orgs', 'value')]
# # )
# # def update_ideogram_rotated(value):
# #     print(value)
# #     full_path = Path("apps/app_data/tree_viewer/bands/") / value
# #     with open(full_path) as fh:
# #         json_out = json.load(fh)
# #     print((json_out))
# #     return json_out

# @app.callback(
#     dash.dependencies.Output('ideogram-annotations', 'annotations'),
#     [dash.dependencies.Input('ideogram-annotations', 'chromosomes')]
# )
# def update_ideogram_rotated(value):
#     return test_out





