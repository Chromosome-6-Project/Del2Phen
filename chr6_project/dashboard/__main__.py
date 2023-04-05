
import sys

from dash import Dash, dash_table, html, dcc, Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from numpy import linspace
import plotly.express as px

from chr6_project.analysis import network, plotting
from chr6_project.analysis.analyze import analyze_online
from chr6_project.analysis.hpo import get_default_termset_yaml_path
from chr6_project.analysis.phenotype_prediction import PredictionDatabase

from chr6_project.dashboard.dashboard_general import head_layout


comparison, geneset, ontology, termset = analyze_online(
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
            dcc.Graph(id="connectivity-bargraph"),
            dcc.Graph(id="size-vs-hi"),
            ]),
        dcc.Tab(label="Homogeneity", children=[
            dcc.Graph(id="homogeneity-histogram"),
            dcc.Graph(id="homogeneity-revcum"),
            dcc.Graph(id="homogeneity-heatmap")
            ]),
        dcc.Tab(label="Predictions", children=[
            html.Br(),
            dcc.Tabs(children=[
                dcc.Tab(label="Summary", children=[
                    dcc.Graph(id="precision-stats-graph"),
                    dcc.Graph(id="min-precision-graph")
                    ]),
                dcc.Tab(label="Patients", children=[
                    html.Br(),
                    html.Label("Select patient:"),
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
              Output("size-vs-hi", "figure"),
              Input("length-slider", "value"),
              Input("loci-slider", "value"),
              Input("gene-slider", "value"),
              Input("hi-gene-slider", "value"),
              Input("dom-gene-switch", "value"),
              Input("group-size-slider", "value"))
def update_connectivity_plot(length_similarity, loci_similarity, gene_similarity,
                             hi_gene_similarity, dom_gene_match, group_size_threshold):
    connectivity = plotting.connectivity_histogram(
        comparison,
        length_similarity,
        loci_similarity,
        gene_similarity,
        hi_gene_similarity,
        dom_gene_match,
        hpo_similarity=0,
        min_size=group_size_threshold
        )
    size_vs_hi = plotting.plot_min_degree_count_vs_hi_score(
        comparison=comparison,
        hi_scores=[x/100 for x in range(0, 101, 10)],
        min_degrees=[group_size_threshold],
        dom_gene_match=dom_gene_match
        )
    return connectivity, size_vs_hi


# @app.callback(Output("size-vs-hi", "figure"),
#               Input("dom-gene-switch", "value"),
#               Input("group-size-slider", "value"))
# def update_size_vs_hi(dom_gene_match, group_size_threshold):
#
#     return size_vs_hi


@app.callback(Output("homogeneity-histogram", "figure"),
              Output("homogeneity-revcum", "figure"),
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
    histogram = plotting.ph_histogram(existing_ph_database=ph_database,
                                      rel_threshold=rel_threshold,
                                      abs_threshold=abs_threshold,
                                      min_size=group_size_threshold+1)

    hi_scores = [i.round(2)
                 for i in linspace(hi_gene_similarity-0.15, hi_gene_similarity+0.15, 7)
                 if i.round(2) > 0]
    revcum = plotting.patients_per_ph(comparison, termset, hi_scores, dom_gene_match,
                                      rel_threshold, abs_threshold, group_size_threshold+1)

    table = ph_database.make_phenotype_homogeneity_table(rel_threshold, abs_threshold,
                                                         min_size=group_size_threshold+1)
    heatmap = px.imshow(table, aspect="auto", height=1000)
    return histogram, revcum, heatmap


@app.callback(
    Output("precision-stats-graph", "figure"),
    Output("min-precision-graph", "figure"),
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
    precision_plot = plotting.plot_precision_stats(prediction_db, comparison.patient_db,
                                                   list(termset), rel_threshold, abs_threshold,
                                                   use_adjusted_frequency, group_size_threshold,
                                                   geneset)
    hi_scores = [i.round(2)
                 for i in linspace(hi_gene_similarity - 0.15, hi_gene_similarity + 0.15, 7)
                 if i.round(2) > 0]
    min_precision_plot = plotting.minimum_precision(comparison, list(termset),
                                                    rel_threshold, abs_threshold,
                                                    use_adjusted_frequency,
                                                    group_size_threshold, hi_scores,
                                                    dom_gene_match)
    prediction_ids = sorted(prediction_db.predictions.keys())
    return precision_plot, min_precision_plot, prediction_ids


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
