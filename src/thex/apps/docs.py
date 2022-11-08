import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt
import pandas as pd
from dash.dependencies import Input, Output

from thex.app import app

docs_sidebar = html.Div(
    [
        html.Div("Contents", className="tv-docs-title"),
        html.Hr(style={"background-color": "black"}),
        dbc.Nav(
            [
                dbc.Button("Tree Viewer", id="docs-treeviewer-button", className="tv-docs-buttons", color="info", outline=True),
                dbc.Button("Signal Tracer", id="docs-pdist-button", className="tv-docs-buttons", color="info", outline=True),
                dbc.Button("THExBuilder", id="docs-thexb-button", className="tv-docs-buttons", color="info", outline=True),
            ],
            # vertical=True,
            pills=True,
        ),
    ],
    style={'width': "100%"}
)

docs_tv_df = pd.DataFrame({
    "Chromosome": ["Chr1", "Chr2", "Chr3", "Chr4"], "Window": ["1000", "2000", "3000", "4000"], 
    "NewickTree": ["(A,(B,C));", "(B,(A,C));", "(C,(A,B));", "(A,(B,C));"], 
    "TopologyID": ["Tree1", "Tree2", "Tree3", "Tree1"],
    "ALT_DATA1": [0.5, 0.25, 0.5, 0.35],
    "ALT_DATA2": ["LOW", "HIGH", "LOW", "HIGH"],
    }
)

docs_tv_chrom_df = pd.DataFrame({"Chromosome": ["Chr1", "Chr2", "Chr3", "Chr4"], "Start": [0, 0, 0, 0], "End": [1000, 2000, 3000, 4000]})

docs_pdist_df = pd.DataFrame(
    {"Chromosome": ["Chr1", "Chr1", "Chr2", "Chr2"],
     "Window": [100, 100, 200, 200],
     "Sample": ["Tiger1", "Tiger2", "Tiger1", "Tiger2"],
     "Value": [0.03, 0.05, 0.02, 0.05]}
)

docs_treeviewer_contents = html.Div(
    children=[
        dbc.Col([
            # html.Div(children=[
            #     html.Video(
            #         controls = True,
            #         id = 'movie_player',
            #         src = "https://www.w3schools.com/html/mov_bbb.mp4",
            #         autoPlay=False,
            #         style={"width": "100%", "height": "90vh"}
            #     ),
            # ]),
            html.Br(),
            # Purpose
            html.H3(
                children=["Purpose"],
                className="tv-docs-content-title"
            ),
            dcc.Markdown(
                children=[
                """
                Tree Viewer is a dashboard designed for window-based approaches to visualizing phylogenetic 
                signal across a reference genome alongside additional data types like recombination rate, GC-content, 
                and gene annotations. Through simultaneous visualization of phylogenetic signal and additional data types, 
                users are offered an all-in-one experience for identifying and understanding the implications of phylogenetic 
                signal variation and its underlying genomic context. Through Tree Viewer, users are able to make more impactful 
                findings and gain a better understanding of their data all in one place. Download the example data sets from 
                GitHub and refer to Li et al. 2019 for an example clade that illustrates Tree Viewer (2). 
                """
                ],
                className="tv-docs-content-text"
            ),
            html.Hr(style={"background-color": "black"}),
            # Input File Structure
            html.H3(
                children=[
                """
                Input File Structures
                """
                ],
                className="tv-docs-content-title"
            ),
            dcc.Markdown(
                """
                There are two required input files to run Tree Viewer: 1. Tree Viewer input file and 2. chromosome length BED file. 
                These two files provide the required information to visualize phylogenetic signal across the genome in proper scale 
                to the chromosome length. Additional window-based data types can be added to the Tree Viewer input file providing 
                the ability to visualize a multitude of data types concurrently. These additional data types are added as a new column 
                in the Tree Viewer file keeping your data succinct and organized.
                """,
                className="tv-docs-content-text",
            ),
            html.Br(),
            dcc.Markdown(
                children=[
                """
                ### _- Tree Viewer Input File -_
                The input file for Tree Viewer is designed to be simple to make and even easier to incorporate new window-based data. 
                Tree Viewer takes a tab or comma delimited file where the first four columns are Chromosome, Window (i.e., 100,000 - which 
                covers bases 1-100,000 for 100kb windows), NewickTree, and TopologyID. The first four columns are required and must have 
                the appropriate headers in the order given in the example input below. Tree Viewer accepts four different file extensions 
                (.csv, .tsv, .txt, .xlsx) for the input file. Note there are column and per-cell limitations to [Excel files](https://support.microsoft.com/en-us/office/excel-specifications-and-limits-1672b34d-7043-467e-8e27-269d656771c3) (.xlsx) files, so 
                large datasets may be better off in a flat file format like .csv, .tsv, or .txt. Ensure the headers of your file match the headers shown in the example below.

                """
                ],
                className="tv-docs-content-text"
            ),
            html.Div([
                dt.DataTable(
                    columns=[{"name": i, "id": i} for i in ["Chromosome", "Window", "NewickTree", "TopologyID"]],
                    data=docs_tv_df.to_dict('records'),
                    style_header={
                        'fontWeight': 'bold',
                        'backgroundColor': 'rgb(230, 230, 230)',
                        'textAlign': 'center'
                    },
                    style_cell={
                        'color': 'black',
                        'textAlign': 'center'
                    },
                    style_data_conditional=[                
                        {
                            "if": {"state": "selected"},
                            "backgroundColor": "inherit !important",
                            "border": "inherit !important",
                        }
                    ],
                ),
            ],className="tv-docs-content-text"),
            html.Br(),
            # Additional Data
            dcc.Markdown(
                children=[
                """
                ### _- Additional Data -_
                Additional numeric and categorical data (i.e., GC-content, recombination rate low/high regions, etc.) can be easily incorporated 
                into your Tree Viewer file by adding an additional column on the right-most side of the Tree Viewer input file. When you load your 
                Tree Viewer file into a new session, the additional data features will load as options in the “Additional Data" dropdown in the 
                “Single Chromosome” tab of the toolbar.
                """
                ],
                className="tv-docs-content-text"
            ),
            html.Div([
                dt.DataTable(
                    columns=[{"name": i, "id": i} for i in docs_tv_df.columns],
                    data=docs_tv_df.to_dict('records'),
                    style_header={
                        'fontWeight': 'bold',
                        'backgroundColor': 'rgb(230, 230, 230)',
                        'textAlign': 'center'
                    },
                    style_cell={
                        'color': 'black',
                        'textAlign': 'center'
                    },
                    style_data_conditional=[                
                        {
                            "if": {"state": "selected"},
                            "backgroundColor": "inherit !important",
                            "border": "inherit !important",
                        }
                    ],
                ),
            ],className="tv-docs-content-text"),
            html.Br(),
            # Chrom length File
            dcc.Markdown(
                children=[
                """
                ### _- Chromosome Length BED File -_
                The second file required to run Tree Viewer is a BED file (.bed) of chromosome lengths. Ensure that the chromosome length bed 
                file contains all chromosomes that are found in the Tree Viewer main input file, otherwise you will be prompted with an error 
                and asked to update your BED file with the missing data. It is also important that you create your file with headers as shown below.
                """
                ],
                className="tv-docs-content-text"
            ),
            html.Div([
                dt.DataTable(
                    columns=[{"name": i, "id": i} for i in ["Chromosome", "Start", "End"]],
                    data=docs_tv_chrom_df.to_dict('records'),
                    style_header={
                        'fontWeight': 'bold',
                        'backgroundColor': 'rgb(230, 230, 230)',
                        'textAlign': 'center'
                    },
                    style_cell={
                        'color': 'black',
                        'textAlign': 'center'
                    },
                    style_data_conditional=[                
                        {
                            "if": {"state": "selected"},
                            "backgroundColor": "inherit !important",
                            "border": "inherit !important",
                        }
                    ],
                ),
            ],className="tv-docs-content-text"),
            html.Br(),
            # Chrom length File
            dcc.Markdown(
                children=[
                """
                ### _- Standardized File Types -_
                Tree Viewer currently only offers the ability to load GFFv3/GTF gene annotation files. Provided one of these files, users can 
                investigate underlying genes, coding regions, and other annotations concurrently with phylogenetic signal and all other user 
                provided window-based data types. More standard data types will be added in future updates.
                """
                ],
                className="tv-docs-content-text"
            ),
            html.Hr(style={"background-color": "black"}),
            html.H3(
                children=["Usage"],
                className="tv-docs-content-title"
            ),
            dcc.Markdown(
                children=[
                """
                Upon opening Tree Viewer, a pop-up will appear prompting you to select two input files. The two required input files are a 
                Tree Viewer input file and a chromosome length BED file (formats are described in input file structure section). Clicking each 
                button will bring up a file system browser allowing you to navigate to the location of your input files. When a file is selected, 
                the upload button will turn green if the file type is valid, or red if it is invalid. After selecting both a Tree Viewer input file 
                as well as a chromosome length BED file, you click submit. Tree Viewer then compares the input files and returns an error message 
                if one or both files are incorrectly formatted. If both files pass the validation checks, they will automatically be uploaded into 
                the session and the toolbar dropdowns will be populated with the data provided in the input files.  

                The dashboard is separated into two sections, the “Single Chromosome” tab and “Whole Genome” tab. Within the “Single Chromosome” tab, 
                the main distribution graph is the master graph for this view, meaning all zooming and panning performed on the master graph will update 
                the x-axis range for all other additional data graphs being shown. This allows for smooth, simultaneous investigation of phylogenetic 
                signal alongside additional data types like recombination rates, GC-content, and gene annotations without the need to align the data yourself.
                “Whole Genome” view allows for a large-scale overview of the phylogenetic signal, making it easy to pick out interesting chromosomal regions 
                to investigate in further detail in the “Single Chromosome” view.  

                “Single Chromosome” view provides the most interactive place to explore your phylogenetic signal. This is where you will add/remove topologies, 
                load additional data types into graphs, and visualize tree topologies. The dropdowns for the topologies, additional data, and tree taxa allow the 
                user to customize which data are being shown at a given time. You may plot one topology or all topologies, with the caveat that as one increases 
                the number of topologies, this complicates visualization and analysis. By default, Tree Viewer orders the topologies by their genome-wide frequencies. 
                However, it may be useful to look at them by chromosome frequency, and there is a switch next to the Topologies header that allows you to change between 
                these alternatives. Turning the Additional Data switch to “On” will load all additional data types that are selected within its respective dropdown, 
                allowing you to load multiple data types at once. Lastly, turning on the Trees switch will load simple representations of the tree topologies selected 
                in the Topologies dropdown. Currently, Tree Viewer plots the trees without branch lengths since branch lengths can change from window-to-window. 
                By selecting/removing taxa names from the Tree’s dropdown, you can prune the tree topology to zoom in on a specific clade or set of taxa. Note that this 
                function does not re-bin topologies, it simply visually prunes the trees. If you are wanting to prune and re-bin topologies, see the File Pruning function 
                in the Tree Viewer export options section.  

                On the right most side of the single chromosome toolbar, there are a few extra options. Under Graph Switches you will find toggles for the main topology graph 
                as well as a switch for GFF/GTF gene annotation files. The main distribution options allow you to independently preview the rug plot and tiled histogram chromosome 
                displays. These previews can be exported using the snapshot button in the top-right of each graph. To upload a GFF/GTF file, go to “File - > GFF/GTF” and upload the 
                file in the same manner as the main input files. Once the file is loaded into the session, the GFF/GTF switch will activate allowing you to load the GFF/GTF graph.
                """
                ],
                className="tv-docs-content-text"
            ),
            html.Hr(style={"background-color": "black"}),
            # Graph Customization + Export Options
            html.H3(
                children=["Graph Customization + Export Options"],
                className="tv-docs-content-title"
            ),
            html.Br(),
            dcc.Markdown(
                children=[
                """
                ### _- Graph Customization -_
                Graph Customization allow you to customize the look of the topology or data distribution graphs. Dropdowns provide lists of themes, 
                color palettes, gridline toggles, and several output formats to choose from. Making changes to the graph customizations will update 
                all graphs currently loaded on the browser. Through these customizations, one can create publication-quality graphs that can be dropped 
                directly into your manuscript. We have allowed users the choice to specify the output dimensions in pixels, however leaving the inputs 
                blank will revert each graph to default dimensions that are optimized for a word document or .pdf file.
                """
                ],
                className="tv-docs-content-text"
            ),
            html.Br(),
            dcc.Markdown(
                children=[
                """
                ### _- Export Options -_
                Tree Viewer offers two export options: current view and file pruning. The current view option allows you to extract local information for 
                a chromosomal region. For example, you are in single chromosome view and looking at Chromosome 1 and you zoom into range 1.5Mb-2Mb. If you 
                click the current view export option, a download prompt will appear allowing you to choose a place to download a new Tree Viewer file that 
                only contains the information for the current range. If you also have a GFF/GTF gene annotation file loaded into the session, it will also 
                extract the information from the given range.  

                The file pruning export option allows you to select a subset of taxa and create a new Tree Viewer input file with a pruned set of Newick trees. 
                This process will also clear the TopologyID column, enabling you to re-bin the new tree topologies. There is an option to run the binning process 
                before downloading the new file, but please note that re-binning large trees across large genomes (e.g, >2 Gbp) may require considerable run time 
                and it will lock your session until it completes. In such circumstances we recommend running the --topobinner command in THExBuilder.    

                ##### ***_Walk-Through: Start New Session + New Load Input Files_***  
                1.	Click “File + New Session” – A pop-up will appear titled “New Session Upload”.  
                2.	Click “Tree Viewer Input” button – navigate to and select Tree Viewer input file.  
                3.	Click “Chromosome Lengths” button – navigate to and select Chromosome Length BED file.  
                4.	Assuming both input buttons have turned green, click “Submit”.  

                ##### ***_Walk-Through: Load GFF/GTF Files_***  
                1.	Click “File + GFF/GTF” – A pop-up will appear titled “Alternative Input Files”.  
                2.	Click “GFF/GTF File” button – navigate to and select GFF/GTF input file.  
                3.	Once the button turns green, click “Submit” – you should see the GFF/GTF switch now active.  

                ##### ***_Walk-Through: Zooming + Panning_***  
                1.	Zoom in by clicking and dragging across a region of interest.  
                2.	You can zoom in or out in steps by using the plus and minus buttons in the toolbar at the top right corner of the graph.  
                3.	Reset the x-axis range by double clicking the graph.  
                4.	Pan by clicking and dragging the x-axis of the graph  
            
                """
                ],
                className="tv-docs-content-text"
            ),
            html.Br(),
            html.Hr(style={"background-color": "black"}),
            # References
            html.H3(
                children=["References"],
                className="tv-docs-content-title"
            ),
            dcc.Markdown(
                children=[
                """
                1.	Edelman, N. B. et al. Genomic architecture and introgression shape a butterfly radiation. Science (80-. ). 366, 594–599 (2019).

                2.	Li, G., Figueiró, H. V, Eizirik, E., Murphy, W. J. & Yoder, A. Recombination-Aware Phylogenomics Reveals the Structured Genomic Landscape of Hybridizing Cat Species. Mol. Biol. Evol. 36, 2111–2126 (2019).

                3.	Nelson, T. C. et al. Ancient and recent introgression shape the evolutionary history of pollinator adaptation and speciation in a model monkeyflower radiation (Mimulus section Erythranthe). PLOS Genet. 17, e1009095 (2021).

                4.	Small, S. T. et al. Radiation with reticulation marks the origin of a major malaria vector. Proc. Natl. Acad. Sci. 202018142 (2020) doi:10.1073/pnas.2018142117.

                5.	trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Salvador Capella-Gutierrez; Jose M. Silla-Martinez; Toni Gabaldon. Bioinformatics 2009 25: 1972-1973.

                6.	Foley et al. Zoonomia Phylogenetic Analysis (unpublished)

                7.	L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015) IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies.. Mol. Biol. Evol., 32:268-274. 

                8.	S.M. Crotty, B.Q. Minh, N.G. Bean, B.R. Holland, J. Tuke, L.S. Jermiin, A. von Haeseler (2019) GHOST: Recovering historical signal from heterotachously-evolved sequence alignments. Syst. Biol., in press. https://doi.org/10.1093/sysbio/syz051

                9.	Robinson, D. F. & Foulds, L. R. Comparison of phylogenetic trees. Math. Biosci. 53, 131–147 (1981).

                10.	ETE 3: Reconstruction, analysis and visualization of phylogenomic data. Jaime Huerta-Cepas, Francois Serra and Peer Bork. Mol Biol Evol 2016; doi: 10.1093/molbev/msw046

                11.	Bredemeyer, K. R. et al. Ultracontinuous Single Haplotype Genome Assemblies for the Domestic Cat (Felis catus) and Asian Leopard Cat (Prionailurus bengalensis). J. Hered. 112, 165–173 (2021).

                """
                ],
                className="tv-docs-content-text"
            ),
            # html.Hr(style={"background-color": "black"}),
            ## Important comments + Performance
            # html.H3(
            #     children=["Notes"],
            #     className="tv-docs-content-title"
            # ),
            # dcc.Markdown(
            #     children=[
            #     """
            #     #### _- Important Comments -_
            #     Tree Viewer was designed with the mindset that users have chromosome level assemblies, however, it is still possible to use Tree Viewer
            #     with scaffold level assemblies. The only stipulation is that Tree Viewer will treat separate scaffolds as independent sequences, so there
            #     is currently no way to link separate scaffolds together into psuedo-chromosomes within the dashboard.  

            #     #### _- Performance -_
            #     Currently, Tree Viewer reliably works with a ~2.5-3.0Gb genome at 50-100kb with little-to-no latency. However, as the size of the input grows,
            #     it takes Tree Viewer longer to load, parse, and store data, increasing the latency. There are updates in the works that should resolve this issue, but the
            #     best current workaround is to split the input file by chromosome (or scaffold) and load the data one chromosome at a time. This does
            #     limit the users ability to visualize the data in the whole genome view, but it will provide you the best user experience for the time being.
            #     You can easily parse your input file by using the "--parse_treeviewer" argument in THExBuilder.  

            #     """
            #     ],
            #     className="tv-docs-content-text"
            # ),
            # html.Hr(style={"background-color": "black"}),
        ], width={"size": 8, "offset": 2}, style={"background-color": "white"}),
    ],
)

docs_pdist_contents = html.Div(
    children=[
        dbc.Col([
            html.Br(),
            # html.Hr(style={"background-color": "black"}),
            html.H3(
                children=["Purpose"],
                className="tv-docs-content-title"
            ),
            dcc.Markdown(
                children=[
                """
                Signal Tracer (ST) is a simple dashboard that visualizes window-based values calculated from multiple-sequence alignments. 
                It provides an easy way to investigate variation in variables like genetic distance, branch length, or divergence times, at 
                single chromosome and whole genome levels. Signal Tracer was originally developed to validate the phasing of F1 hybrid long-reads 
                using the trio-binning approach for single haplotype genome assembly references. By visualizing genetic distance of reference 
                and non-reference species, we were able to visually validate that there were no regions that indicate potential improper phasing. 
                If you would like to learn more about this process, refer to our [paper](https://academic.oup.com/jhered/article/112/2/165/6030926).  

                """
                ],
                className="tv-docs-content-text"
            ),
            html.Hr(style={"background-color": "black"}),
            html.H3(
                children=["Input File Structure"],
                className="tv-docs-content-title"
            ),
            dcc.Markdown(
                children=[
                """
                The input file for Signal Tracer is a tab or comma delimited file consisting of four columns: _Chromosome_, _Window_, _Sample_, and _Value_. 
                The Value column can only contain numerical values, not categorical values. Ensure the headers of your input file match the headers 
                listed in italics above, an error message will appear if they do not match exactly. 
                """
                ],
                className="tv-docs-content-text"
            ),
            html.Div([
                dt.DataTable(
                    columns=[{"name": i, "id": i} for i in docs_pdist_df.columns],
                    data=docs_pdist_df.to_dict('records'),
                    style_header={
                        'fontWeight': 'bold',
                        'backgroundColor': 'rgb(230, 230, 230)',
                        'textAlign': 'center'
                    },
                    style_cell={
                        'color': 'black',
                        'textAlign': 'center'
                    },
                    style_data_conditional=[                
                        {
                            "if": {"state": "selected"},
                            "backgroundColor": "inherit !important",
                            "border": "inherit !important",
                        }
                    ],
                ),
            ],className="tv-docs-content-text"),
            html.Hr(style={"background-color": "black"}),
            html.H3(
                children=["Usage"],
                className="tv-docs-content-title"
            ),
            dcc.Markdown(
                children=[
                """
                ### _- New Session -_
                Upon opening the Signal Tracer dashboard, a pop-up will appear promting you to choose a file. 
                Selecting a file will automatically load the file into the session and load the first chromosome 
                into the plot. The dashboard will start in single chromosome view, but can be easily changed to 
                whole-genome view by using the "View" dropdown and choosing "Whole Genome". 

                If you are already loaded into a session and wish to switch to another file, all you need to do is 
                go to "File -> New Session" and the input pop-up will appear again. Select the new file and it will 
                load into the session replacing the old file.

                ### _ - Graph Interactability -_
                The Signal Tracer graph is highly interactive allowing easy zooming, panning, and manipulation of traces. 
                You can zoom by click+dragging over a region and releasing, or remove a trace by clicking it's respective 
                name in the legend. Whole genome viewing allows for zooming of individual chromosomes without affecting others, 
                making it easy to highlight an important regions of one chromosome, while maintaining perspective of the other 
                chromosomes. 

                ### _- p-Distance Calculator -_
                This dashboard was originally designed to visualize raw genetic distances variation across the genome, so THExBuilder 
                comes will a “--pdistance” command to generate input files for Signal Tracer from multi-alignment fasta files. Visit the 
                Signal Tracer Pipeline in the THExBuilder documentation for more information.

                ##### ***_Walk-Through: Load Input Files_***
                1.	Click “File + New Session” – A pop-up will appear titled “New Session Upload”.
                2.	Click “Signal Tracer File” button – navigate to and select Tree Viewer input file.
                3.	Assuming the input button turned green, click “Submit”.
                
                ##### ***_Walk-Through: Load GFF/GTF Files_***
                1.	Click “File + GFF/GTF” – A pop-up will appear titled “Alternative Input Files”.
                2.	Click “GFF/GTF File” button – navigate to and select GFF/GTF input file.
                3.	Once the button turns green, click “Submit” – you should see the GFF/GTF switch now active.
                
                ##### ***_Walk-Through: Zooming + Panning_***
                1.	Zoom in by clicking and dragging across a region, then release.
                2.	You can zoom in or out in steps by using the plus and minus buttons in the top right corner of the graph.
                3.	Reset the view by double clicking the graph.
                4.	Pan by clicking and dragging the x-axis of the graph


                """
                ],
                className="tv-docs-content-text"
            ),
            html.Hr(style={"background-color": "black"}),
            # Output Options
            html.H3(
                children=["Graph Customization + Export Options"],
                className="tv-docs-content-title"
            ),
            dcc.Markdown(
                children=[
                """
                ### _- Graph Customization -_
                Graph Customization allow you to customize the look of the chromosome graphs in an easy-to-use way. 
                Dropdowns provide lists of themes, color palettes, gridline toggles, and several output formats to 
                choose from. Making changes to the graph customizations will update all graphs currently loaded on 
                the browser. Through these customizations, one can create publication-quality graphs that can be 
                dropped directly into your in-progress manuscript. Although Signal Tracer allows users to specify 
                the output dimensions in pixels, leaving the inputs blank will revert each graph to default dimensions 
                that are optimized for a word document or .pdf file.
                
                ### _- Export Options -_
                Signal Tracer offers a single export option that allows you to extract local information from your input 
                file. Similar to Tree Viewer, zooming into a region of a chromosome in single chromosome view and selecting 
                “File -> Current View” will bring up a download prompt where you can specify the file name and download location.  

                ##### ***_Walk-Through: Current View Export_***
                1.	Click “File + Current View”.
                2.	A file-system prompt will appear and ask you to select a location and file name. Once you have selected a location and file name, click “Save”.  
                    _- If a GFF/GTF file is loaded into the session, then a second file-system prompt will appear to download a subset of the GFF/GTF file too_
                """
                ],
                className="tv-docs-content-text"
            ),
            html.Br(),
            html.Hr(style={"background-color": "black"}),
            # References
            html.H3(
                children=["References"],
                className="tv-docs-content-title"
            ),
            dcc.Markdown(
                children=[
                """
                1.	Edelman, N. B. et al. Genomic architecture and introgression shape a butterfly radiation. Science (80-. ). 366, 594–599 (2019).

                2.	Li, G., Figueiró, H. V, Eizirik, E., Murphy, W. J. & Yoder, A. Recombination-Aware Phylogenomics Reveals the Structured Genomic Landscape of Hybridizing Cat Species. Mol. Biol. Evol. 36, 2111–2126 (2019).

                3.	Nelson, T. C. et al. Ancient and recent introgression shape the evolutionary history of pollinator adaptation and speciation in a model monkeyflower radiation (Mimulus section Erythranthe). PLOS Genet. 17, e1009095 (2021).

                4.	Small, S. T. et al. Radiation with reticulation marks the origin of a major malaria vector. Proc. Natl. Acad. Sci. 202018142 (2020) doi:10.1073/pnas.2018142117.

                5.	trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Salvador Capella-Gutierrez; Jose M. Silla-Martinez; Toni Gabaldon. Bioinformatics 2009 25: 1972-1973.

                6.	Foley et al. Zoonomia Phylogenetic Analysis (unpublished)

                7.	L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015) IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies.. Mol. Biol. Evol., 32:268-274. 

                8.	S.M. Crotty, B.Q. Minh, N.G. Bean, B.R. Holland, J. Tuke, L.S. Jermiin, A. von Haeseler (2019) GHOST: Recovering historical signal from heterotachously-evolved sequence alignments. Syst. Biol., in press. https://doi.org/10.1093/sysbio/syz051

                9.	Robinson, D. F. & Foulds, L. R. Comparison of phylogenetic trees. Math. Biosci. 53, 131–147 (1981).

                10.	ETE 3: Reconstruction, analysis and visualization of phylogenomic data. Jaime Huerta-Cepas, Francois Serra and Peer Bork. Mol Biol Evol 2016; doi: 10.1093/molbev/msw046

                11.	Bredemeyer, K. R. et al. Ultracontinuous Single Haplotype Genome Assemblies for the Domestic Cat (Felis catus) and Asian Leopard Cat (Prionailurus bengalensis). J. Hered. 112, 165–173 (2021).

                """
                ],
                className="tv-docs-content-text"
            ),
            # html.Hr(style={"background-color": "black"}),
        ], width={"size": 8, "offset": 2}, style={"background-color": "white"}),
    ],
)

docs_thexb_contents = html.Div(
    children=[
        dbc.Col([
            html.Br(),
            # html.Hr(style={"background-color": "black"}),
            # Purpose
            html.H3(
                children=["Purpose"],
                className="tv-docs-content-title"
            ),
            dcc.Markdown(
                children=[
                """
                THExBuilder is a suite of pipelines and tools designed to generate, modify, and analyze THEx input files. 
                It provides a customized solution for generating Tree Viewer input files. It also contains a single command 
                pipeline for calculating raw p-distance values from multiple-sequence fasta files in non-overlapping sliding 
                windows for Signal Tracer. THExBuilder was designed to be flexible and scalable, taking advantage of 
                multiprocessing and partitioned workflows. For example, the Tree Viewer Pipeline is designed to be run either 
                with a configuration file or by command line arguments, as well as being able to run the pipeline step-by-step 
                or from start-to-finish with a single command. This provides the user the freedom and flexibility to tailor their 
                runs and track the progress through pipelines with ease.
	            
                NOTE: Due to limitations with a few dependencies in the Tree Viewer Pipeline, THExBuilder is only supported for 
                Linux and OSX systems. The Windows installation of THEx only comes with the browser.

                """
                ],
                className="tv-docs-content-text"
            ),
            html.Hr(style={"background-color": "black"}),
            # Usage
            html.H3(
                children=["Usage"],
                className="tv-docs-content-title"
            ),
            dcc.Markdown(
                children=[
                """
                In addition to being a novel interactive genome browser, THEx comes with a built-in command line suite 
                that contains tools and pipelines that help simplify the process of building and manipulating input files 
                for the dashboards. There are two main pipelines, one for Tree Viewer and one for Signal Tracer. To use THExBuilder, 
                activate the conda environment that THEx is installed within and enter the command _thexb_ to begin. This command will 
                bring up the help section, providing an overview of all the options and available tools.  

                THExBuilder's output is contained in a single THExBuilderOutput directory. By default, the output directory is written 
                to the current working directory, so if you do not provide an output location be sure to use the same working directory 
                as you move through the pipelines. Within the output directory you will find the intermediate output files for each stage 
                of the different pipelines as well as log files to help track what has been performed, to what files, in what order, and 
                with what input parameters. 

                """
                ],
                className="tv-docs-content-text"
            ),
            html.Hr(style={"background-color": "black"}),
            # Input Files
            html.H3(
                children=["Input Files"],
                className="tv-docs-content-title"
            ),
            dcc.Markdown(
                children=[
                """
                All of the pipelines for THExBuilder start with multiple-sequence alignment fasta files, organized within a 
                single directory. Each fasta file ideally represents the multiple-sequence alignment of a single chromosome, 
                but scaffold level assemblies can be used too with the caveat that each sequence will be treated as independent 
                sequences and there is no way to link scaffolds into pseudo-chromosomes within the pipeline or application. 

                """
                ],
                className="tv-docs-content-text"
            ),
            html.Hr(style={"background-color": "black"}),
            # Pipelines
            html.H3(
                children=["Pipelines"],
                className="tv-docs-content-title"
            ),
            html.Br(),
            # Tree Viewer Pipeline
            dcc.Markdown(
                children=[
                """
                ### _- Tree Viewer Pipeline -_
                The Tree Viewer pipeline is designed to take multi-alignment fasta files, run them through several filtration 
                steps, and then through IQ-Tree to generate per-window maximum likelihood phylogenies. These phylogenies along 
                with their chromosome-window coordinates are put together in a Tree Viewer file where the topologies can then 
                be binned and labeled for viewing in Tree Viewer. The examples below show how to run each stage with and without 
                a configuration file, providing the user flexibility to run the pipeline with all arguments conveniently located 
                in a single file. It is recommended to use a configuration file as this helps improve reproducibility and also 
                reduces the length of the commands, but the option is left to the user. Since the configuration file has to be 
                built with a specific format, you can generate a blank configuration file running “thexb --tv_config_template”.

                The pipeline is designed to be flexible, providing the user complete control over how the pipeline 
                is run. Although not recommended, one can run the entire Tree Viewer pipeline start-to-finish by providing a 
                configuration file and calling the “--tv_all” argument. This will run through each of the pipeline steps described 
                below, but will not stop after each step is completed. Although the option to run the pipeline start-to-finish is 
                available, it is recommended to run each step individually to inspect and validate the intermediate results.  
                

                ```py
                $ thexb --tv_all -c config.ini
                
                ```
                ### _Stages:_
                                
                1.	_Fasta Windowing_   
                
                    The fasta windowing stage will take a directory containing multi-alignment fastas (ideally whole chromosomes) 
                    and will return a directory of subdirectories, labeled by the original sequence names (i.e., chromosome names), 
                    where each subdirectory contains the windowed fasta files. Note that this step can return thousands of files when 
                    using a small window size (i.e., 5kb for a 2.5Gb genome), so be mindful when opening these directories in a file 
                    browser.

                    _With configuration file_: 
                    ```py
                    $ thexb --minifastas -c config.ini
                    ```

                    _Without configuration file_: 
                    ```py
                    $ thexb --minifastas --window_size 10kb -i dir_of_multi_alignment_files/ 
                    ```
                    
                    _Descriptions_:
                    ```
                    -i = "Pathway to directory containing per-chromosome multi-alignment fasta files"
                    --window_size = "Size of window as 5bp/kb/mb/gb"
                    ```

                2. _Trimal Gap Trimming_  

                    The Trimal stage of the pipeline uses the gap-threshold feature that removes regions of a window where a taxon 
                    has missing data that exceeds the provided threshold for allowed gaps. A second threshold is provided to remove resulting 
                    sequences that do not meet a minimum sequence length.

                    _With configuration file_:
                    ```py
                    $ thexb --trimal -c config.ini
                    ```
                    _Without configuration file_:    
                    ```py
                    $ thexb --trimal --trimal_gap_threshold 0.9 --trimal_min_seq_len 1kb -i THExBuilderOutput/windowed_fastas
                    ```
                    _Descriptions_:
                    ```py
                    -i = "Input file path - Cannot be used when running Tree Viewer pipeline with configuration file"
                    --trimal_gap_threshold = "Percentage of alignment with no gaps (i.e., removes sequence with gaps in 10 percent or more of the sequences)" (default: 0.1)  
                    --trimal_min_seq_len = "Minimum sequence length to retain window in study. If above, window is excluded" (default: 1kb)  
                    ```
                        
                3.	_Pairwise Distance Filter_  

                    The pairwise distance and coverage filter walks through each window of the genome and conducts a 
                    sliding sub-window analysis with a step size of n-bp to mask local regions of extreme divergence 
                    that may be caused by misalignment or undetected paralogy. It also masks entire windows where a 
                    single taxon’s base coverage (valid bases include A, T, G, and C) is below the provided threshold. 
                    This step was integrated from Foley et al. 2021 and updated to run on Python 3.

                    _With configuration file_: 
                    ```py 
                    $ thexb --pw_filter --pw_subwindow_size 100bp --pw_step 10bp --pw_min_seq_len 1000bp --pw_min_pdist 0.15 --pw_zscore_cutoff 2 --pw_seq_coverage 0.9 -i THExBuilderOutput/trimal_filtered_windows
                    ```
                    
                    _Without configuration file_:
                    ```py
                    $ thexb --pw_filter -c config.ini
                    ```

                    _Descriptions_:  
                    ```py
                    -i = "Pathway to directory containing the output from the Trimal stage" (i.e., THExBuilderOutput/trimal_filtered_windows/)
                    --pw_subwindow_size = "Sub-window size" (default: 100bp)
                    --pw_step 10bp = "Sub-window step size" (default: 10bp)
                    --pw_min_seq_len = "Minimum length of sequence in a given window to be retained" (default: 1000bp)
                    --pw_max_pdist = "Maximum p-distance for a given taxon to be considered non-extreme. Any sequence above this value in a given sub-window will be masked. Note this value is dataset and taxon dependent." (default: 0.15)
                    --pw_zscore_cutoff = "Maximum number of standard deviations away from the mean (calculated as z-score) allowed before entire window is masked." (default: 2)
                    --pw_seq_coverage = "Minimum valid, non-missing sequence coverage. Individual sequences below this value are masked." (default: 0.9)
                    ```
                        
                4.	_IQ-Tree Maximum Likelihood Tree Construction_  

                    IQ-Tree is used to infer maximum likelihood phylogenies providing the basis for the distribution 
                    of phylogenetic signal visualized in Tree Viewer. The filtered fasta files generated by the pairwise 
                    filtration step is passed to IQ-Tree and the resulting Newick files are collected and organized into the 
                    basis of the Tree Viewer input file. IQ-Tree creates several output files per-window, but we are only 
                    interested in the trees so the other data can be saved or discarded.

                    There are two approaches built into the Tree Viewer pipeline for the IQ-Tree analysis. The first is to have 
                    THExBuilder call IQ-Tree for each filtered window by using the “--iqtree” command, a model of sequence evolution, 
                    and number of bootstrap replications. The second approach is to take the filtered fasta files from the pairwise 
                    filtration step and run IQ-Tree externally on a different local machine, server, or cluster. If you choose to run 
                    IQ-Tree externally, you will need to organize the resulting “.treefile” Newick trees into a single directory and 
                    pass the directory to the “-iqtree_external” command rather than the “--iqtree” command. You can also optionally 
                    create sub-directories per chromosome and organize the “.treefile” by chromosome and pass it to the “--iqtree_external” 
                    command, but it is not required.

                    _With configuration file_: 
                    ```py
                    $ thexb --iqtree -c config.ini
                    ```

                    _Without configuration file_: 
                    ```py
                    $ thexb –iqtree --iqtree_model “GTR*H4” --iqtree_bootstrap 1000 -i dir_of_windowed_fastas
                    ```

                    _Descriptions_:
                    ```py
                    -i = "Pathway to directory containing the output from the Pairwise Distance + Coverage Filter step" 
                    --iqtree_model = "Any valid model of nucleotide sequence evolution supported by IQ-Tree" (default: GTR*H4)
                    --iqtree_bootstrap = "Number of bootstrap replicates to run" (default: 1000)
                    ```
                    

                5.	_Topobinner_  

                    Topobinner is a tool used to organize and label equal tree topologies based on RF-distance. Trees are treated 
                    as equal topologies when their RF-distance is equal to 0 (9-10). This is what allows you to visualize discordant 
                    topologies across the genome and identify regions of interest. The trees are binned and then labeled by their 
                    whole genome frequency, meaning the most frequent topology in the genome will be labeled Tree1, then Tree 2 for 
                    the second most frequent topology, and so on. Once this step is completed, an updated Tree Viewer file with binned 
                    topologies will be generated and placed in the THExBuilderOutput directory ready for visualization in Tree Viewer. 
                    It should also be noted that users can also pass a Tree Viewer input file produced by the File Pruning export option 
                    (see Tree Viewer export options section for more details) in Tree Viewer to bin new trees that have been pruned into 
                    a subset from a larger Tree Viewer input file. No special steps need to be taken to run this, simply provide the 
                    Tree Viewer input file as you would normally.
                    
                    _With configuration file_:
                    ```py
                    $ thexb --topobinner -c config.ini
                    ```

                    _Without configuration file_: 
                    ```py
                    $ thexb --topobinner -i TreeViewer_file_with_blank_TopologyID_col.xlsx
                    ```
                    
                    _Descriptions_:
                    ```py
                    -i = "Tree Viewer input file from IQ-Tree step or File Pruning export option with TopologyID column blank"
                    ```

                    _Topobinner Alternative_
                    An alternative to using topobinner is [PhyBin](https://github.com/rrnewton/PhyBin). PhyBin is not hosted on conda, 
                    so you would need download it from GitHub and set it up externally. 
                    
                    It is important that you use PhyBin's "--bin" option and not it’s clustering algorithm. When you run PhyBin, it will 
                    generate binned results in a new output directory named either phybin_output or whatever you provide as an output location. 
                    Passing the output directory generated by PhyBin and a Tree Viewer file with the TopologyID column blank to THExBuilder's 
                    "--phybin_external" command will produce a complete Tree Viewer input file just like topobinner. 

                    _Example Command_:
                    ```py
                    $ thexb --phybin_external -i pathway_to_phybin_output_directory/ --tv_file_name ./THExBuilderOutput/TreeViewer_input_file.xlsx
                    ```
                
                """
                ],
                className="tv-docs-content-text"
            ),
            html.Br(),
            html.Hr(style={"background-color": "black"}), 
            dcc.Markdown(
                children=[
                """
                ### _- Signal Tracer Pipeline -_
                The Signal Tracer pipeline is a single script that calculates raw, uncorrected p-distance between a 
                reference and other taxon of one or more multi-alignment fasta files. There are two ways to run the pipeline, 
                you can provide a single multi-alignment fasta file or a directory of multi-alignment fasta files. The directory 
                of files is typically the input directory used for the entire Tree Viewer pipeline that is ideally a multi-alignment 
                fasta file for each chromosome. Since some chromosomes can be quite large, this process may take several hours to 
                run, if not a day or two depending on the system configuration being used (i.e., processor, amount of memory, etc.). 
                The only two additional arguments that are needed to run the pipeline are the window size and missing data threshold. 
                The missing data threshold is a value between 0.0 and 1.0 that indicates the maximum amount of missing data 
                (gaps or masked sequence) allowed for a single taxon, and if that threshold is surpassed, then the values for all 
                taxa in the window are set to NULL. The output from this pipeline is formatted to be directly uploaded into the Signal 
                Tracer dashboard for quick visualization.
                
                _Example Command_:
                ```py
                $ thexb --pdistance --pdist_threshold 0.75 --window_size 100kb
                ```
                _Descriptions_:  
                ```py
                --pdist_threshold = "Maximum frequency of missing data allowed in a single taxon before the entire window is masked." (default: 0.75)  
                --window_size = "Window size (bp/kb/mb)" (default: 100kb)  
                ```
                """
                ],
                className="tv-docs-content-text"
            ),
            html.Br(),
            html.Hr(style={"background-color": "black"}),
            # References
            html.H3(
                children=["References"],
                className="tv-docs-content-title"
            ),
            dcc.Markdown(
                children=[
                """
                1.	Edelman, N. B. et al. Genomic architecture and introgression shape a butterfly radiation. Science (80-. ). 366, 594–599 (2019).

                2.	Li, G., Figueiró, H. V, Eizirik, E., Murphy, W. J. & Yoder, A. Recombination-Aware Phylogenomics Reveals the Structured Genomic Landscape of Hybridizing Cat Species. Mol. Biol. Evol. 36, 2111–2126 (2019).

                3.	Nelson, T. C. et al. Ancient and recent introgression shape the evolutionary history of pollinator adaptation and speciation in a model monkeyflower radiation (Mimulus section Erythranthe). PLOS Genet. 17, e1009095 (2021).

                4.	Small, S. T. et al. Radiation with reticulation marks the origin of a major malaria vector. Proc. Natl. Acad. Sci. 202018142 (2020) doi:10.1073/pnas.2018142117.

                5.	trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Salvador Capella-Gutierrez; Jose M. Silla-Martinez; Toni Gabaldon. Bioinformatics 2009 25: 1972-1973.

                6.	Foley et al. Zoonomia Phylogenetic Analysis (unpublished)

                7.	L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015) IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies.. Mol. Biol. Evol., 32:268-274. 

                8.	S.M. Crotty, B.Q. Minh, N.G. Bean, B.R. Holland, J. Tuke, L.S. Jermiin, A. von Haeseler (2019) GHOST: Recovering historical signal from heterotachously-evolved sequence alignments. Syst. Biol., in press. https://doi.org/10.1093/sysbio/syz051

                9.	Robinson, D. F. & Foulds, L. R. Comparison of phylogenetic trees. Math. Biosci. 53, 131–147 (1981).

                10.	ETE 3: Reconstruction, analysis and visualization of phylogenomic data. Jaime Huerta-Cepas, Francois Serra and Peer Bork. Mol Biol Evol 2016; doi: 10.1093/molbev/msw046

                11.	Bredemeyer, K. R. et al. Ultracontinuous Single Haplotype Genome Assemblies for the Domestic Cat (Felis catus) and Asian Leopard Cat (Prionailurus bengalensis). J. Hered. 112, 165–173 (2021).

                """
                ],
                className="tv-docs-content-text"
            ),
        ], width={"size": 8, "offset": 2}, style={"background-color": "white"}),
    ],
)

docs_content = html.Div(id="docs-content", className="docs-content-div")

layout = dbc.Container(
    children=[
        # Nav
        dbc.Row(
            children=[
                dbc.Col(
                    children=[
                        html.Div([
                            dbc.Button(
                                children=[
                                    html.Img(src=app.get_asset_url("docs4.svg"), height='30px'),
                                    html.Div("Docs")
                                ],
                                className="logo-png",
                                href="/",
                                outline=True,
                            ),
                        ]),
                    ],
                    width=2,
                    className="homepage-title-col"
                ),
                dbc.Col(
                    children=[
                        html.H2("Contents | ", style={'vertical-align': 'middle', 'paddingTop': '4px'}),
                        dbc.Button("Tree Viewer", id="docs-treeviewer-button", className="tv-docs-buttons", color="primary"),
                        dbc.Button("Signal Tracer", id="docs-pdist-button", className="tv-docs-buttons", color="primary"),
                        dbc.Button("THExBuilder", id="docs-thexb-button", className="tv-docs-buttons", color="primary"),
                    ],
                    width=8,
                    className="homepage-title-col"
                ),
                dbc.Col(className="homepage-title-col", width=2)
            ],
        ),
        # Body
        dbc.Row(
            children=[
                # Sidebar
                # dbc.Col(
                #     children=[docs_sidebar],
                #     width=12,
                #     style={"border": "orange solid 2px", "background-color": "#575151", "border-radius": "5px"},
                # ),
                # Content
                dbc.Col(
                    children=[docs_content],
                    width=12,
                    style={"border": "orange solid 2px", "border-radius": "5px"},
                ),
            ],
            style={"padding": "2px"}
        ),
    ],
    fluid=True,
)

def docs_layout():
    return layout

@app.callback(
    [Output("docs-content", "children"),
     Output("docs-treeviewer-button", "active"),
     Output("docs-pdist-button", "active"),
     Output("docs-thexb-button", "active"),],
    [Input("docs-treeviewer-button", "n_clicks"),
     Input("docs-pdist-button", "n_clicks"),
     Input("docs-thexb-button", "n_clicks"),],)
def render_page_content(btn1, btn2, btn3):
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if button_id == "docs-treeviewer-button":
        page_content = html.P(
            children=[docs_treeviewer_contents],
        )
        return page_content, True, False, False
    elif button_id == "docs-pdist-button":
        page_content = html.P(
            children=[docs_pdist_contents],
        )
        return page_content, False, True, False
    elif button_id == "docs-thexb-button":
        page_content = html.P(
            children=[docs_thexb_contents],
        )
        return page_content, False, False, True
    else:
        page_content = html.P(
            children=[docs_treeviewer_contents],
        )
        return page_content, True, False, False