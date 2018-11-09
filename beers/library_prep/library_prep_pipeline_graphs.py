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
    original_df = pd.read_csv("../../data/library_prep/output/original_sample.log", quotechar="'")
    original_df['seq_length'] = original_df.apply(lambda row: len(row['sequence']), axis=1)
    print(f"Original:\n{original_df.head(1)}\n")

    # Post polyA data
    polyA_df = pd.read_csv("../../data/library_prep/output/polya_selection_step.log", quotechar="'")
    polyA_retained_df = polyA_df[polyA_df["note"] != "removed"]
    polyA_retained_df['seq_length'] = polyA_retained_df.apply(lambda row: len(row['sequence']), axis=1)
    print(f"Poly A:\n{polyA_df.head(1)}\n")

    # Post Fragmentation data
    frag_df = pd.read_csv("../../data/library_prep/output/fragment_step.log", quotechar="'")
    frag_df['seq_length'] = frag_df.apply(lambda row: len(row['sequence']), axis=1)
    print(f"Fragmented:\n{frag_df.head(1)}\n")

    # Post First Strand Prime data
    first_strand_prime_df = pd.read_csv("../../data/library_prep/output/first_strand_prime_step.log", quotechar="'")
    first_strand_prime_df['seq_length'] = first_strand_prime_df.apply(lambda row: len(row['sequence']), axis=1)
    print(f"First Strand Prime:\n{first_strand_prime_df.head(1)}\n")

    # Post First Strand Synthesis data
    first_strand_syn_df = pd.read_csv("../../data/library_prep/output/first_strand_synthesis_step.log", quotechar="'")
    first_strand_syn_df['seq_length'] = first_strand_syn_df.apply(lambda row: len(row['sequence']), axis=1)
    print(f"First Strand Synthesis:\n{first_strand_syn_df.head(1)}\n")

    # Post Second Strand Prime data
    second_strand_prime_df = pd.read_csv("../../data/library_prep/output/second_strand_prime_step.log", quotechar="'")
    second_strand_prime_df['seq_length'] = second_strand_prime_df.apply(lambda row: len(row['sequence']), axis=1)
    print(f"Second Strand Prime:\n{second_strand_prime_df.head(1)}\n")

    # Post Second Strand Synthesis data
    second_strand_syn_df = pd.read_csv("../../data/library_prep/output/second_strand_synthesis_step.log", quotechar="'")
    second_strand_syn_df['seq_length'] = second_strand_syn_df.apply(lambda row: len(row['sequence']), axis=1)
    print(f"Second Strand Synthesis:\n{second_strand_syn_df.head(1)}\n")

    # Post Sizing data
    sizing_df = pd.read_csv("../../data/library_prep/output/sizing_step.log", quotechar="'")
    sizing_retained_df = sizing_df[sizing_df["note"] != "removed"]
    sizing_retained_df['seq_length'] = sizing_retained_df.apply(lambda row: len(row['sequence']), axis=1)
    print(f"Sizing:\n{sizing_df.head(1)}\n")

    # Post PCR Amplification data
    pcr_amp_df = pd.read_csv("../../data/library_prep/output/pcr_amplification_step.log", quotechar="'")
    pcr_amp_df['seq_length'] = pcr_amp_df.apply(lambda row: len(row['sequence']), axis=1)
    print(f"PCR Amplification:\n{pcr_amp_df.head(1)}\n")

    print(f"Total Sequences Originally: {len(original_df)}")
    print(f"Total Sequences Post Poly A: {len(polyA_retained_df)}")
    print(f"Total Sequences Post fragment: {len(frag_df)}")
    print(f"Total Sequences Post first strand prime: {len(first_strand_prime_df)}")
    print(f"Total Sequences Post first strand synthesis: {len(first_strand_syn_df)}")
    print(f"Total Sequences Post second strand prime: {len(second_strand_prime_df)}")
    print(f"Total Sequences Post second strand synthesis: {len(second_strand_syn_df)}")
    print(f"Total Sequences Post sizing: {len(sizing_retained_df)}")
    print(f"Total Sequences Post PCR Amplification: {len(pcr_amp_df)}")


    # Histogram comparing sequence lengths before and after step
    polya_seq_lengths = [
        go.Histogram(
            x=original_df['seq_length'],
            opacity=0.75,
            name="Pre PolyA Selection Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
        go.Histogram(
            x=polyA_retained_df['seq_length'],
            opacity=0.75,
            name="Post PolyA Selection Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
    ]

    frag_seq_lengths = [
        go.Histogram(
            x=polyA_retained_df['seq_length'],
            opacity=0.75,
            name="Pre Fragmentation Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
        go.Histogram(
            x=frag_df['seq_length'],
            opacity=0.75,
            name="Post Fragmentation Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
    ]

    first_strand_prime_seq_lengths = [
        go.Histogram(
            x=frag_df['seq_length'],
            opacity=0.75,
            name="Pre First Strand Prime Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
        go.Histogram(
            x=first_strand_prime_df['seq_length'],
            opacity=0.75,
            name="Post First Stand Prime Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
    ]

    first_strand_syn_seq_lengths = [
        go.Histogram(
            x=first_strand_prime_df['seq_length'],
            opacity=0.75,
            name="Pre First Strand Synthesis Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
        go.Histogram(
            x=first_strand_syn_df['seq_length'],
            opacity=0.75,
            name="Post First Stand Synthesis Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
    ]

    second_strand_prime_seq_lengths = [
        go.Histogram(
            x=first_strand_syn_df['seq_length'],
            opacity=0.75,
            name="Pre Second Strand Prime Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
        go.Histogram(
            x=second_strand_prime_df['seq_length'],
            opacity=0.75,
            name="Post Second Stand Prime Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
    ]

    second_strand_syn_seq_lengths = [
        go.Histogram(
            x=second_strand_prime_df['seq_length'],
            opacity=0.75,
            name="Pre Second Strand Synthesis Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
        go.Histogram(
            x=second_strand_syn_df['seq_length'],
            opacity=0.75,
            name="Post Second Stand Synthesis Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
    ]

    sizing_seq_lengths = [
        go.Histogram(
            x=second_strand_syn_df['seq_length'],
            opacity=0.75,
            name="Pre Sizing Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
        go.Histogram(
            x=sizing_retained_df['seq_length'],
            opacity=0.75,
            name="Post Sizing Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
    ]

    pcr_amp_seq_lengths = [
        go.Histogram(
            x=sizing_retained_df['seq_length'],
            opacity=0.75,
            name="Pre PCR Amplification Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
        go.Histogram(
            x=pcr_amp_df['seq_length'],
            opacity=0.75,
            name="Post PCR Amplication Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
    ]

    data_seq_lengths = [
        go.Histogram(
            x=original_df['seq_length'],
            opacity=0.75,
            name="Starting Material",
            xbins = dict(start=0, end=4000, size=20)
        ),
        go.Histogram(
            x=polyA_retained_df['seq_length'],
            opacity=0.75,
            name="Post PolyA Selection Step",
            xbins = dict(start=0, end=4000, size=20)
        ),
        go.Histogram(
            x=frag_df['seq_length'],
            opacity=0.75,
            name="Post Fragment Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
        #go.Histogram(
        #    x=first_strand_prime_df['seq_length'],
        #    opacity=0.75,
        #    name="Post First Strand Prime Step",
        #    xbins=dict(start=0, end=4000, size=20)
        #),
        go.Histogram(
            x=sizing_retained_df['seq_length'],
            opacity=0.75,
            name="Post Sizing Step",
            xbins=dict(start=0, end=4000, size=20)
        ),
        go.Histogram(
            x=pcr_amp_df['seq_length'],
            opacity=0.75,
            name="Post PCR Amplification Step",
            xbins=dict(start=0, end=4000, size=20)
        )
    ]

    figure_polya_lengths = {
        'data': polya_seq_lengths,
        'layout': go.Layout(
                height=350,
                barmode='overlay',
                title='Pre-Post Poly A Sequence Lengths',
                hovermode='closest')
        }

    figure_frag_lengths = {
        'data': frag_seq_lengths,
        'layout': go.Layout(
            height=350,
            barmode='overlay',
            title='Pre-Post Fragmentation Sequence Lengths',
            hovermode='closest')
    }

    figure_first_strand_prime_lengths = {
        'data': first_strand_prime_seq_lengths,
        'layout': go.Layout(
            height=350,
            barmode='overlay',
            title='Pre-Post First Strand Prime Sequence Lengths',
            hovermode='closest')
    }

    figure_first_strand_syn_lengths = {
        'data': first_strand_syn_seq_lengths,
        'layout': go.Layout(
            height=350,
            barmode='overlay',
            title='Pre-Post First Strand Synthesis Sequence Lengths',
            hovermode='closest')
    }

    figure_second_strand_prime_lengths = {
        'data': second_strand_prime_seq_lengths,
        'layout': go.Layout(
            height=350,
            barmode='overlay',
            title='Pre-Post Second Strand Prime Sequence Lengths',
            hovermode='closest')
    }

    figure_second_strand_syn_lengths = {
        'data': second_strand_syn_seq_lengths,
        'layout': go.Layout(
            height=350,
            barmode='overlay',
            title='Pre-Post Second Strand Synthesis Sequence Lengths',
            hovermode='closest')
    }

    figure_sizing_lengths = {
        'data': sizing_seq_lengths,
        'layout': go.Layout(
            height=350,
            barmode='overlay',
            title='Pre-Post Sizing Sequence Lengths',
            hovermode='closest')
    }

    figure_pcr_amp_lengths = {
        'data': pcr_amp_seq_lengths,
        'layout': go.Layout(
            height=350,
            barmode='overlay',
            title='Pre-Post PCR Amplification Sequence Lengths',
            hovermode='closest')
    }

    figure_seq_lengths = {
        'data': data_seq_lengths,
        'layout': go.Layout(
                height=500,
                barmode='overlay',
                title='Sequence Lengths',
                hovermode='closest')
        }


    app.layout = html.Div([
        html.H1("Pipeline Sequence Length Data", style={"text-align": "center", "marginTop": 5, "marginBottom": 5}),
        html.P("Starter Molecules: " +  info['molecule_info']),
        html.Div([
            html.Div([
                dcc.Graph(id='polya_lengths', figure=figure_polya_lengths)
            ], style= {'width': '45%', 'display': 'inline-block'}),
            html.Div([
                dcc.Graph(id='frag_lengths', figure=figure_frag_lengths)
            ], style= {'width': '45%', 'display': 'inline-block'})
        ]),
        html.Div([
            html.Div([
                dcc.Graph(id='first_strand_prime_lengths', figure=figure_first_strand_prime_lengths)
            ], style={'width': '45%', 'display': 'inline-block'}),
            html.Div([
                dcc.Graph(id='first_strand_syn_lengths', figure=figure_first_strand_syn_lengths)
            ], style={'width': '45%', 'display': 'inline-block'})
        ]),
        html.Div([
            html.Div([
                dcc.Graph(id='second_strand_prime_lengths', figure=figure_second_strand_prime_lengths)
            ], style={'width': '45%', 'display': 'inline-block'}),
            html.Div([
                dcc.Graph(id='second_strand_syn_lengths', figure=figure_second_strand_syn_lengths),
            ], style={'width': '45%', 'display': 'inline-block'})
        ]),
        html.Div([
            html.Div([
                dcc.Graph(id='sizing_lengths', figure=figure_sizing_lengths)
            ], style={'width': '45%', 'display': 'inline-block'}),
            html.Div([
                dcc.Graph(id='pcr_amp_lengths', figure=figure_pcr_amp_lengths)
            ], style={'width': '45%', 'display': 'inline-block'})
        ]),
        dcc.Graph(id='seq_lengths', figure=figure_seq_lengths),
    ])


if __name__ == '__main__':
    info = {
        "molecule_info": "10K actual transcripts",
    }
    generate_graphs(info)
    app.run_server()