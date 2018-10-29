import pandas as pd
import plotly.graph_objs as go
import dash
import dash_core_components as dcc
import dash_html_components as html
import re

app = dash.Dash()


def provide_info(info):
    return info

def generate_graphs(info):

    # Original sample data frame - adding seq length column
    original_df = pd.read_csv("../../data/original_sample.log")
    original_df['seq_length'] = original_df.apply(lambda row: len(row['sequence']), axis=1)

    # Post polyA data
    polyA_df = pd.read_csv("../../data/polyA_selection_step.log")
    polyA_retained_df = polyA_df[polyA_df["note"] != "removed"]
    polyA_retained_df['seq_length'] = polyA_retained_df.apply(lambda row: len(row['sequence']), axis=1)

    # Post Fragmentation data
    frag_df = pd.read_csv("../../data/fragment_step.log")
    frag_df['seq_length'] = frag_df.apply(lambda row: len(row['sequence']), axis=1)

    # Post Sizing data
    sizing_df = pd.read_csv("../../data/sizing_step.log")
    sizing_retained_df = sizing_df[sizing_df["note"] != "removed"]
    sizing_retained_df['seq_length'] = sizing_retained_df.apply(lambda row: len(row['sequence']), axis=1)

    # Post PCR Amplification data
    pcr_amp_df = pd.read_csv("../../data/pcr_amplification_step.log")
    pcr_amp_df['seq_length'] = pcr_amp_df.apply(lambda row: len(row['sequence']), axis=1)

    print(f"Total Sequences Originally: {len(original_df)}")
    print(f"Total Sequences Post Poly A: {len(polyA_retained_df)}")
    print(f"Total Sequences Post fragment: {len(frag_df)}")
    print(f"Total Sequences Post sizing: {len(sizing_retained_df)}")
    print(f"Total Sequences Post PCR Amplification: {len(pcr_amp_df)}")


    # Histogram comparing sequence lengths before and after step
    data_seq_lengths = [
        go.Histogram(
            x=original_df['seq_length'],
            opacity=0.75,
            name="Starting Material",
            xbins = dict(start=0, end=6500, size=20)
        ),
        go.Histogram(
            x=polyA_retained_df['seq_length'],
            opacity=0.75,
            name="Post PolyA Selection Step",
            xbins = dict(start=0, end=6500, size=20)
        ),
        go.Histogram(
            x=frag_df['seq_length'],
            opacity=0.75,
            name="Post Fragment Step",
            xbins=dict(start=0, end=6500, size=20)
        ),
        go.Histogram(
            x=sizing_retained_df['seq_length'],
            opacity=0.75,
            name="Post Sizing Step",
            xbins=dict(start=0, end=6500, size=20)
        ),
        go.Histogram(
            x=pcr_amp_df['seq_length'],
            opacity=0.75,
            name="Post PCR Amplification Step",
            xbins=dict(start=0, end=6500, size=20)
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
        html.H1("Pipeline Sequence Length Data", style={"text-align": "center", "marginTop": 5, "marginBottom": 5}),
        html.P("Starter Molecules: " +  info['molecule_info']),
        dcc.Graph(id='seq_lengths', figure=figure_seq_lengths),
    ])


if __name__ == '__main__':
    info = {
        "molecule_info": "50K, Size 5k, Tail: 200, normal distribution",
    }
    generate_graphs(info)
    app.run_server()