import pandas as pd
import plotly.graph_objs as go
import dash
import dash_core_components as dcc
import dash_html_components as html
import re

# Program to compare before and after molecule states for the Sizing Selection Step.
# Comparison based on log files provided by SizingStep.py test.

app = dash.Dash()


def provide_info(info):
    return info

def generate_graphs(input_log, output_log, info):

    # Input data frame - adding seq length and tail length columns
    input_df = pd.read_csv(input_log)
    input_df['seq_length'] = input_df.apply(lambda row: len(row['sequence']), axis=1)

    # Output data frame - filtering out retained molecules and adding seq length and tail length columns
    output_df = pd.read_csv(output_log)
    retained_df = output_df[output_df["note"] != "removed"]
    retained_df['seq_length'] = retained_df.apply(lambda row: len(row['sequence']), axis=1)

    print("Total Sequences Input: " + str(len(input_df)))
    print("Total Sequences Output: " + str(len(retained_df)))

    # Histogram comparing sequence lengths before and after step
    data_seq_lengths = [
        go.Histogram(
            x=retained_df['seq_length'],
            opacity=0.75,
            name="Retained Molecules",
            xbins = dict(start=0, end=1000, size=20)
        ),
        go.Histogram(
            x=input_df['seq_length'],
            opacity=0.75,
            name="Input Molecules",
            xbins = dict(start=0, end=1000, size=20)
        )
    ]

    figure_seq_lengths = {
        'data': data_seq_lengths,
        'layout': go.Layout(
                barmode='overlay',
                title='Sequence Lengths',
                hovermode='closest')
        }


    app.layout = html.Div([
        html.H1("Sizing Selection Step Data", style={"text-align": "center", "marginTop": 5, "marginBottom": 5}),
        html.P("Molecules: " +  info['molecule_info']),
        html.P("Model: " + info['model_info']),
        html.P("Parameters: " + info['param_info']),
        dcc.Graph(id='seq_lengths', figure=figure_seq_lengths),
    ])


if __name__ == '__main__':
    info = {
        "molecule_info": "100K, Size 25-1000K uniform distribution",
        "model_info": "Retention - Normal/Line Chimera",
        "param_info": "mean - 250, std dev - 50"
    }
    generate_graphs("../../data/sizing_step_input_data.log","../../data/sizing_step_output_data.log", info)
    app.run_server()