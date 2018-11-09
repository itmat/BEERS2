import pandas as pd
import plotly.graph_objs as go
import dash
import dash_core_components as dcc
import dash_html_components as html
import re


app = dash.Dash()


def provide_info(info):
    return info

def generate_graphs(input_log, output_log, info):

    # Input data frame - adding seq length and tail length columns
    input_df = pd.read_csv(input_log)
    input_df['seq_length'] = input_df.apply(lambda row: len(row['sequence']), axis=1)

    # Output data frame - filtering out retained molecules and adding seq length and tail length columns
    output_df = pd.read_csv(output_log)
    output_df['seq_length'] = output_df.apply(lambda row: len(row['sequence']), axis=1)

    print("Total Sequences Input: " + str(len(input_df)))
    print("Total Sequences Output: " + str(len(output_df)))

    # Histogram comparing sequence lengths before and after step
    data_seq_lengths = [
        go.Histogram(
            x=output_df['seq_length'],
            opacity=0.75,
            name="Post PCR Amplification Molecules",
            xbins = dict(start=0, end=1000, size=5)
        ),
        go.Histogram(
            x=input_df['seq_length'],
            opacity=0.75,
            name="Pre PCR Amplification Molecules",
            xbins = dict(start=0, end=1000, size=5)
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
        html.H1("PCR Amplication Step Data", style={"text-align": "center", "marginTop": 5, "marginBottom": 5}),
        html.P("Molecules: " +  info['molecule_info']),
        dcc.Graph(id='seq_lengths', figure=figure_seq_lengths),
    ])


if __name__ == '__main__':
    info = {
        "molecule_info": "original 10K - actual transcripts post sizing",
    }
    generate_graphs("../../data/tests/pcr_amplification_step_input_data.log","../../data/tests/pcr_amplification_step_output_data.log", info)
    app.run_server()