
from math import log10
import sys

from dash import Dash, dash_table, html, dcc, Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from numpy import linspace
import plotly.express as px
import plotly.graph_objects as go

from chr6_project.analysis.analyze import analyze_online
from chr6_project.analysis.hpo import get_default_termset_yaml_path
from chr6_project.analysis import network
# from chr6_project.analysis.phenotype_homogeneity import phenotype_homo_test
from chr6_project.analysis.plotting import plot_precision_stats, ph_histogram
from chr6_project.analysis.phenotype_prediction import PredictionDatabase

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

# c6_dir = "/home/tyler/Documents/Chr6_docs/"
# comparison, _, ontology, termset = analyze(
#     genotypes=f"{c6_dir}/PatientData/2022-Oct-25/c6_array_2022-10-25_18_58_17.csv",
#     phenotypes=f"{c6_dir}/PatientData/2022-Oct-25/c6_questionnaire_2022-10-25_19_00_59.csv",
#     patient_hpo=f"{c6_dir}/PatientData/2022-Oct-25/c6_research_patients_2022-10-25_19_49_06.csv",
#     geneset_gtf=f"{c6_dir}/GeneSets/hg19.ensGene.chr6.gtf.gz",
#     drop_list_file=f"{c6_dir}/PatientData/drop_list.txt",
#     hpo_termset_yaml=f"{c6_dir}/Phenotype_Homogeneity/selected_phenotypes_hpos.yaml",
#     expand_hpos=False
#     )

comparison, ontology, termset = analyze_online(
    username=sys.argv[1], password=sys.argv[2], drop_list_file=sys.argv[3],
    hpo_termset_yaml=get_default_termset_yaml_path()
    )

patient_ids = sorted(comparison.patient_db.list_ids())
prediction_db = PredictionDatabase(dict())
nx_graph = network.make_nx_graph_from_comparisons(comparison)


class DefaultSlider(dcc.Slider):
    def __init__(self, slider_id, minimum=0, maximum=1, value=0.0, show_dec_as_perc=True,
                 nticks=11, step=0.01):
        if show_dec_as_perc:
            marks = {i: str(int(i*100)) for i in linspace(minimum, maximum, nticks)}
        else:
            marks = {i: str(int(i)) for i in linspace(minimum, maximum, nticks)}
        super().__init__(
            id=slider_id, min=minimum, max=maximum, value=value,
            marks=marks,
            step=step,
            tooltip={"placement": "bottom"}
            )


app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "The Chromosome 6 Project"

app.layout = dbc.Container(fluid=True, children=[
    head_layout,
    dbc.Offcanvas(id="off_canvas", is_open=False, placement="end", children=[
        dcc.Tabs(children=[
            dcc.Tab(label="Genotype Similarity", children=[
                html.Label("Size Similarity"),
                DefaultSlider("length-slider"),
                html.Br(),
                html.Label("Overlap"),
                DefaultSlider("loci-slider"),
                html.Br(),
                html.Label("Gene Similarity"),
                DefaultSlider("gene-slider"),
                html.Br(),
                html.Label("HI Gene Similarity"),
                DefaultSlider("hi-gene-slider", value=0.75),
                html.Br(),
                dbc.Switch("dom-gene-switch", label="Group by Dominant-Effect Genes",
                           value=True),
                html.Br(),
                html.Label("Group Size Threshold"),
                DefaultSlider("group-size-slider", 1, 100, 5, False, step=1)
                ]),
            dcc.Tab(label="Phenotype Homogeneity", children=[
                html.Label("Relative Prevalence Threshold"),
                DefaultSlider("rel-prevalence-slider", value=0.2),
                html.Br(),
                html.Label("Absolute Prevalence Threshold"),
                DefaultSlider("abs-prevalence-slider", 1, 10, 2, False, 10, 1),
                # html.Label("Group Size Threshold"),
                # DefaultSlider("group-size-slider"),
                ]),
            dcc.Tab(label="Predictions", children=[
                html.Br(),
                html.Label(),
                dbc.Switch(id="adjusted-freq-switch", label="Adjust Population Freq",
                           value=True)
                ])
            ]),
        # html.H3("Similarity Settings"),
        # html.Hr(),
        # html.Label("Size Similarity"),
        # DefaultSlider("length-slider"),
        # html.Br(),
        # html.Label("Overlap"),
        # DefaultSlider("overlap-slider"),
        # html.Br(),
        # html.Label("Gene Similarity"),
        # DefaultSlider("gene-slider"),
        # html.Br(),
        # html.Label("HI Gene Similarity"),
        # DefaultSlider("hi-slider"),
        # html.Br(),
        # html.H3("Phenotype Homogeneity Settings"),
        # html.Hr(),
        # html.Label("Prevalence Level Threshold"),
        # DefaultSlider("homogeneity-slider"),
        # html.Br(),
        # html.Label("Group Size Threshold"),
        # DefaultSlider("group-size-slider"),
        ]),
    # Panes
    dcc.Tabs(children=[
        dcc.Tab(label="Network", children=[
            dcc.Graph(id="connectivity-bargraph")
            ]),
        dcc.Tab(label="Homogeneity", children=[
            dcc.Graph(id="homogeneity-graph"),
            dcc.Graph(id="homogeneity-heatmap")
            ]),
        dcc.Tab(label="Predictions", children=[
            html.Br(),
            dcc.Tabs(children=[
                dcc.Tab(label="Summary", children=[
                    dcc.Graph(id="precision-stats-graph"),
                    ]),
                dcc.Tab(label="Patients", children=[
                    html.Br(),
                    dcc.Dropdown(id="prediction-selector", options=patient_ids),
                    html.Br(),
                    html.Div(id="prediction-table", children=[dash_table.DataTable()]),
                    ])
                ]),
            ])
        ]),
    # dbc.Row(children=[
    #     dbc.Col(width=6, children=[
    #         dcc.Graph(id="connectivity-bargraph")
    #         ]),
    #     dbc.Col(width=6, children=[
    #         dcc.Graph(id="homogeneity-graph")
    #         ])
    #     ]),
    # dbc.Row(children=[
    #     dbc.Col(children=[
    #         dcc.Graph(id="homogeneity-heatmap")
    #         ])
    #     ]),
    # html.Div([
    #     # Left Pane
    #     html.Div([
    #         dcc.Graph(id="connectivity-bargraph"),
    #         ], style={"padding": 10, "flex": 1}),
    #     # Right Pane
    #     html.Div([
    #         dcc.Graph(id="homogeneity-graph"),
    #         dcc.Graph(id="homogeneity-heatmap"),
    #         ], style={"padding": 10, "flex": 1}),
    #     ], style={"display": "flex", "flex-direction": "row"}),
    ])


@app.callback(Output("connectivity-bargraph", "figure"),
              Input("length-slider", "value"),
              Input("loci-slider", "value"),
              Input("gene-slider", "value"),
              Input("hi-gene-slider", "value"),
              Input("dom-gene-switch", "value"),
              Input("group-size-slider", "value"))
def update_connectivity_plot(length_similarity, loci_similarity, gene_similarity,
                             hi_gene_similarity, dom_gene_match, group_size_thresholds):
    filtered_graph = network.filter_graph_edges(nx_graph,
                                                length_similarity,
                                                loci_similarity,
                                                gene_similarity,
                                                hi_gene_similarity,
                                                dom_gene_match,
                                                hpo_similarity=0)
    subsizes = network.get_subnet_sizes(filtered_graph)
    subsize_gt_5 = network.count_subnets_over_size_n(filtered_graph,
                                                     group_size_thresholds)
    subsize_gt_5 = (f"Clusters with at least {group_size_thresholds} individuals: "
                    f"{subsize_gt_5:>5}")

    links = sorted(list(network.get_node_degrees(filtered_graph).values()))
    links_gt_5 = network.count_nodes_over_n_degree(filtered_graph,
                                                   group_size_thresholds)
    links_gt_5 = (f"Individuals with at least {group_size_thresholds} links:    "
                  f"{links_gt_5:>5}")

    params = dict(xbins=dict(start=0, end=500, size=1),
                  autobinx=False)

    fig1 = go.Figure()
    fig1.add_trace(
        trace=go.Histogram(
            x=subsizes,
            name="Subnet Sizes",
            hovertemplate="Subnet Size: %{x}<br>Count: %{y}<extra></extra>",
            **params
            )
        )
    fig1.add_trace(
        trace=go.Histogram(
            x=links,
            name="Node Degrees",
            hovertemplate="Links: %{x}<br>Count: %{y}<extra></extra>",
            **params
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
    fig1.add_annotation(text=subsize_gt_5 + "<br>" + links_gt_5, showarrow=False,
                        align="right", xref="paper", yref="paper", x=0.95, y=0.95)

    return fig1


@app.callback(Output("homogeneity-graph", "figure"),
              Output("homogeneity-heatmap", "figure"),
              Input("length-slider", "value"),
              Input("loci-slider", "value"),
              Input("gene-slider", "value"),
              Input("hi-gene-slider", "value"),
              Input("dom-gene-switch", "value"),
              Input("group-size-slider", "value"),
              Input("rel-prevalence-slider", "value"),
              Input("abs-prevalence-slider", "value"))
def update_homogeneity_figures(length_similarity, loci_similarity, gene_similarity,
                               hi_gene_similarity, dom_gene_match, group_size_threshold,
                               rel_threshold, abs_threshold):
    ph_database = comparison.compare_all_patient_pheno_prevalences(
        list(termset), length_similarity, loci_similarity, gene_similarity,
        hi_gene_similarity, dom_gene_match, hpo_similarity=0
        )
    histogram = ph_histogram(existing_ph_database=ph_database,
                             rel_threshold=rel_threshold,
                             abs_threshold=abs_threshold,
                             min_size=group_size_threshold+1)

    table = ph_database.make_phenotype_homogeneity_table(rel_threshold, abs_threshold,
                                                         min_size=group_size_threshold+1)
    heatmap = px.imshow(table, aspect="auto", height=1000)
    return histogram, heatmap


@app.callback(
    Output("precision-stats-graph", "figure"),
    Output("prediction-selector", "options"),
    Input("length-slider", "value"),
    Input("loci-slider", "value"),
    Input("gene-slider", "value"),
    Input("hi-gene-slider", "value"),
    Input("dom-gene-switch", "value"),
    Input("group-size-slider", "value"),
    Input("rel-prevalence-slider", "value"),
    Input("abs-prevalence-slider", "value"),
    Input("adjusted-freq-switch", "value")
    )
def update_predictions(length_similarity, loci_similarity, gene_similarity,
                       hi_gene_similarity, dom_gene_match, group_size_threshold,
                       rel_threshold, abs_threshold, use_adjusted_frequency):
    prediction_db = comparison.test_all_phenotype_predictions(
        length_similarity, loci_similarity, gene_similarity,
        hi_gene_similarity, dom_gene_match,
        )
    precision_plot = plot_precision_stats(prediction_db, comparison.patient_db,
                                          list(termset), rel_threshold, abs_threshold,
                                          use_adjusted_frequency, group_size_threshold)
    prediction_ids = sorted(prediction_db.predictions.keys())
    return precision_plot, prediction_ids


@app.callback(
    Output("prediction-table", "children"),
    Input("prediction-selector", "value"),
    State("length-slider", "value"),
    State("loci-slider", "value"),
    State("gene-slider", "value"),
    State("hi-gene-slider", "value"),
    State("dom-gene-switch", "value"),
    )
def update_prediction_table(patient_id, length_similarity, loci_similarity,
                            gene_similarity, hi_gene_similarity, dom_gene_match):
    if prediction_db is None or patient_id is None:
        raise PreventUpdate
    prediction = comparison.test_phenotype_prediction(
        patient_id, length_similarity, loci_similarity, gene_similarity,
        hi_gene_similarity, dom_gene_match, 0
        )
    table = prediction.convert_patient_predictions_to_df()
    table = dash_table.DataTable(table.to_dict("records"),
                                 sort_action="native",
                                 page_action="native", page_size=20)
    return table


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
