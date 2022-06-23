
from math import log10

from dash import Dash, html, dcc, callback, Input, Output, State
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
import plotly.express as px
import plotly.graph_objects as go
# from plotly.subplots import make_subplots

from chr6_project.analysis.analyze import analyze
from chr6_project.analysis import network
from chr6_project.analysis.phenotype_homogeneity import make_phenotype_homogeneity_table

from chr6_project.dashboard.dashboard_general import head_layout

# @callback(
#     Output("login_modal", "is_open"),
#     Output("username_input", "value"),
#     Output("password_input", "value"),
#     Input("open_login", "n_clicks"),
#     State("login_modal", "is_open"),
#     )
# def toggle_modal(n1, is_open):
#     if n1 == 0:
#         raise PreventUpdate
#     return not is_open, "", ""
#
#
# @callback(
#     Output("login_div", "style"),
#     Output("logout_div", "style"),
#     Input("user_credentials", "resources"),
#     Input("logout_button", "n_clicks"),
#     )
# def toggle_login_buttons(creds, logout):
#     if creds == {} and logout == 0:
#         raise PreventUpdate
#     if creds != {}:
#         return {'display': 'none'}, {'display': 'block'}
#     return {'display': 'block'}, {'display': 'none'}
#
#
# @callback(
#     Output("user_credentials", "resources"),
#     Output("open_login", "n_clicks"),
#     Output("username_input", "invalid"),
#     Output("password_input", "invalid"),
#     Input("login_button", "n_clicks"),
#     Input("password_input", "n_submit"),
#     Input("logout_button", "n_clicks"),
#     State("username_input", "value"),
#     State("password_input", "value"),
#     State("open_login", "n_clicks")
#     )
# def update_credentials(login_clicks, submits, logouts, user, password, click_count):
#     if login_clicks == 0 and submits == 0 and click_count == 0 and logouts == 0:
#         raise PreventUpdate
#     if ctx.triggered_id == "logout_button":
#         return {}, 0, False, False
#     credentials = {"hostname": "c-head", "username": user, "password": password}
#     if not test_connection(credentials):
#         return {}, 0, True, True
#     click_count = click_count + 1
#     return credentials, click_count, False, False
#
#
# @callback(
#     Output("historical_coverage", "resources"),
#     Output("historical_insert_len", "resources"),
#     Input("user_credentials", "resources"),
#     )
# def fetch_historical_info(login):
#     if login == {}:
#         return {}, {}
#     agg = "/home/tmedina/QC_Project/aggregate/"
#     hist_coverage = stream_json_table(login, agg + "Clinical_2022_WGS_Germline_aggregate_coverage.binned.json.gz")
#     hist_insert = stream_tsv_table(login, agg + "Clinical_2022_DNA_aggregate_insert_len.tsv.gz")
#     hist_insert = hist_insert.to_json()
#     return hist_coverage, hist_insert
#
#
# @callback(
#     Output(f"navbar-collapse", "is_open"),
#     Input(f"navbar-toggler", "n_clicks"),
#     State(f"navbar-collapse", "is_open"),
#     )
# def toggle_navbar_collapse(n, is_open):
#     if n:
#         return not is_open
#     return is_open
#
#
# @callback(
#     Output("current_page_layout", "children"),
#     Input("url", "pathname")
#     )
# def display_page(pathname):
#     if pathname == "/case-explorer":
#         return case_explorer_layout
#     elif pathname == '/historical-trends':
#         return historical_trends_layout
#     else:
#         return case_explorer_layout


comparison, _, ontology = analyze(
    genotypes="/home/tyler/Documents/Chr6_docs/PatientData/2022-Feb-21/c6_array_2022-02-21_18_28_06.csv",
    phenotypes="/home/tyler/Documents/Chr6_docs/PatientData/2022-Feb-21/c6_questionnaire_2022-02-21_10_18_09.csv",
    patient_hpo="/home/tyler/Documents/Chr6_docs/PatientData/2022-Feb-21/c6_research_patients_2022-02-21_10_20_53.csv",
    geneset_gtf="/home/tyler/Documents/Chr6_docs/GeneSets/hg19.ensGene.chr6.gtf.gz",
    drop_list_file="/home/tyler/Documents/Chr6_docs/PatientData/drop_list.txt",
    expand_hpos=False
    )

with open("/home/tyler/Documents/Chr6_docs/Phenotype_Homogeneity/selected_phenotypes_hpos.txt") as infile:
    selected_hpos = infile.readlines()
selected_hpos = [ontology[x.strip().split("\t")[1]] for x in selected_hpos]

nx_graph = network.make_nx_graph_from_comparisons(comparison)


class DefaultSlider(dcc.Slider):
    def __init__(self, slider_id):
        super().__init__(
            id=slider_id, min=0, max=100, value=0,
            marks={str(sim): str(sim) for sim in range(0, 101, 10)},
            tooltip={"placement": "bottom"}
            )


app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "ClinGen QC Dashboard"

app.layout = html.Div([
    head_layout,
    html.Button(id="option_pane"),
    dbc.Offcanvas(id="off_canvas", is_open=False, children=[
        html.H3("Similarity Settings"),
        html.Hr(),
        html.Label("Size Similarity"),
        DefaultSlider("length-slider"),
        html.Br(),
        html.Label("Overlap"),
        DefaultSlider("overlap-slider"),
        html.Br(),
        html.Label("Gene Similarity"),
        DefaultSlider("gene-slider"),
        html.Br(),
        html.Label("HI Gene Similarity"),
        DefaultSlider("hi-slider"),
        ]),
    # Panes
    html.Div([
        # Left Pane
        html.Div([
            dcc.Graph(id="graph-with-slider"),
            html.Div(id="subsize-5-count"),
            html.Div(id="link-5-count"),
            html.Hr(),

            ], style={"padding": 10, "flex": 1}),
        # Right Pane
        html.Div([
            dcc.Graph(id="homogeneity-graph"),
            dcc.Graph(id="homogeneity-heatmap"),
            html.Label("Prevalence Level Threshold"),
            dcc.Slider(id="homogeneity-slider",
                       min=0, max=100, value=20,
                       marks={str(sim): str(sim) for sim in range(0, 101, 10)},
                       tooltip={"placement": "bottom"}),
            html.Br(),
            html.Label("Group Size Threshold"),
            dcc.Slider(id="group-size-slider",
                       min=0, max=100, value=5,
                       marks={str(sim): str(sim) for sim in range(0, 101, 10)},
                       tooltip={"placement": "bottom"}),
            ], style={"padding": 10, "flex": 1}),
        ], style={"display": "flex", "flex-direction": "row"}),
    ])


@app.callback(Output("homogeneity-graph", "figure"),
              Output("homogeneity-heatmap", "figure"),
              Input("homogeneity-slider", "value"),
              Input("group-size-slider", "value"),
              Input("length-slider", "value"),
              Input("overlap-slider", "value"),
              Input("gene-slider", "value"),
              Input("hi-slider", "value"))
def update_homogeneity_figures(homogeneity_threshold, group_size_threshold,
                               length_threshold, overlap_threshold,
                               gene_threshold, hi_threshold):
    all_homogens, upper_homogens, lower_homogens = comparison.test_all_homogeneities(
        selected_hpos, length_threshold/100, overlap_threshold/100,
        gene_threshold/100, hi_threshold/100, hpo_similarity=0,
        group_size_threshold=group_size_threshold
        )
    table = make_phenotype_homogeneity_table(all_homogens, selected_hpos, homogeneity_threshold / 100, 2)
    upper_homogens = [100*homogen.calculate_homogeneity(homogeneity_threshold/100, 2)
                      for homogen in upper_homogens.values()]
    lower_homogens = [100*homogen.calculate_homogeneity(homogeneity_threshold/100, 2)
                      for homogen in lower_homogens.values()]
    histogram = go.Figure()
    histogram.add_trace(
        trace=go.Histogram(
            x=upper_homogens,
            name=f"Group size >= {group_size_threshold}",
            xbins=dict(start=0,
                       end=105,
                       size=5),
            autobinx=False,
            hovertemplate="Homogeneity Score: %{x}%<br>Count: %{y}<extra></extra>"
            )
        )
    histogram.add_trace(
        trace=go.Histogram(
            x=lower_homogens,
            name=f"Group size < {group_size_threshold}",
            xbins=dict(start=0,
                       end=105,
                       size=5),
            autobinx=False,
            hovertemplate="Homogeneity Score: %{x}%<br>Count: %{y}<extra></extra>"
            )
        )
    histogram.update_layout(
        transition_duration=500,
        title_text="Distribution of Homogeneity Scores",
        xaxis_title_text="Score",
        yaxis_title_text="Count",
        barmode="stack"
        )
    heatmap = px.imshow(table, aspect="auto")
    return histogram, heatmap


@app.callback(Output("graph-with-slider", "figure"),
              Output("subsize-5-count", "children"),
              Output("link-5-count", "children"),
              Input("length-slider", "value"),
              Input("overlap-slider", "value"),
              Input("gene-slider", "value"),
              Input("hi-slider", "value"))
def update_figure(length_threshold, overlap_threshold, gene_threshold, hi_threshold):
    filtered_graph = network.filter_graph_edges(nx_graph, length_threshold/100,
                                                overlap_threshold/100, gene_threshold/100,
                                                hi_threshold/100)
    subsizes = network.get_subnet_sizes(filtered_graph)
    subsize_gt_5 = network.count_subnets_over_size_n(filtered_graph, 5)
    subsize_gt_5 = f"Clusters with at least 5 individuals: {subsize_gt_5}"

    links = sorted(list(network.get_node_degrees(filtered_graph).values()))
    links_gt_5 = network.count_nodes_over_n_degree(filtered_graph, 5)
    links_gt_5 = f"Individuals with at least 5 links: {links_gt_5}"

    # fig1 = make_subplots(rows=1, cols=2)
    fig1 = go.Figure()
    fig1.add_trace(
        # row=1,
        # col=1,
        trace=go.Histogram(
            x=subsizes,
            xbins=dict(start=0,
                       end=500,
                       size=1),
            autobinx=False,
            name="Subnet Sizes",
            hovertemplate="Subnet Size: %{x}<br>Count: %{y}<extra></extra>"
            )
        )
    fig1.add_trace(
        # row=1,
        # col=2,
        trace=go.Histogram(
            x=links,
            xbins=dict(start=0,
                       end=500,
                       size=1),
            autobinx=False,
            name="Node Degrees",
            hovertemplate="Links: %{x}<br>Count: %{y}<extra></extra>"
            )
        )
    fig1.update_layout(
        transition_duration=500,
        title_text="Distribution of Subnet Sizes and Node Degrees",
        xaxis_title_text="Size",
        yaxis_title_text="Count",
        legend_orientation="h"
        )
    fig1.update_yaxes(type="log", dtick=log10(2))

    return fig1, subsize_gt_5, links_gt_5


@app.callback(
    Output("off_canvas", "is_open"),
    Input("option_pane", "n_clicks"),
    State("off_canvas", "is_open"),
    )
def toggle_offcanvas(n1, is_open):
    if n1:
        return not is_open
    return is_open


if __name__ == "__main__":
    app.run_server(debug=True, use_reloader=False)
