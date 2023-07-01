from collections import Counter
from dataclasses import dataclass, field

import numpy as np
import pandas as pd 
from Bio import SeqIO
import plotly.graph_objects as go
from plotly.subplots import make_subplots



NT_COLOR = { 
    "A": "#e30808",
    "T": "#01b201",
    "C": "#0871fd",
    "G": "#fda400",
    "N": "lightgray",
    "#": "purple",
    "*": "black",
    "R" : "cyan",
    "Y" : "cyan",
    "S" : "cyan",
    "W" : "cyan",
    "K" : "cyan",
    "M" : "cyan",
    "B" : "cyan",
    "D" : "cyan",
    "H" : "cyan",
    "V" : "cyan"}


@dataclass
class Report:
    '''Visualization for each sample'''
    reference: str
    sample_id: str 
    
    path: str = field(init=False)
    reference_file: str = field(init=False)
    coverage_file: str = field(init=False)
    consensus: pd.DataFrame = field(init=False)
    distrib: pd.DataFrame = field(init=False)
    reference_sequence: str = field(init=False)
    output: str = field(init=False)
    

    def __post_init__(self) -> None:
        self.path = f"/home/Nanopore/input/output/{self.sample_id}"
        self.coverage_file = f"{self.path}/coverage"
        self.output = f"{self.path}/report.html"
        self.reference_file = f"{self.path}/{self.reference}"
        self.construct_data()
        self.generate_report()
        

    def analyze_seq(self, sequence:str) -> dict[str, int]:
        bps = {k:0 for k in NT_COLOR.keys()}
        cnt = Counter(sequence.upper())
        bps_percent = {k: v/len(sequence) for k, v in cnt.items()}
        bps.update(bps_percent)
        return bps


    def construct_data(self) -> None:
        self.consensus = pd.read_csv(f"{self.path}/consensus.pileup", sep=r"\t", header=None, engine="python")
        self.consensus.columns = [
            "reference_name", 
            "reference_position", 
            "base_at_pos", 
            "map_count", 
            "consensus_bp", 
            "consensus_confidence", 
            "sequences", 
            "quality"]

        df_coverage = pd.read_csv(self.coverage_file, sep=r"\t", header=None, engine="python")
        self.consensus["coverage"] = df_coverage.iloc[:,2].values

        mean_phred = [np.mean([ord(x)-33 for x in qual]) for qual in self.consensus["quality"].values]
        self.consensus["mean_phred"] = mean_phred
        
        for record in SeqIO.parse(self.reference_file, "fasta"):
            self.reference_sequence = str(record.seq).upper()

        dist = [self.analyze_seq(cons) for cons in self.consensus["sequences"].values]
        self.distrib = pd.DataFrame(dist)


    def generate_report(self) -> None:
        heights = [.2, .2, .26, .07, .07, .2]
        heights=[.26, .07, .07, .2, .2, .2]
        colscale = "Emrld"
        titles = [
            "Sequencing result",
            "Consensus sequence",
            "Reference sequence", 
            "Consensus confidence",
            "Mean Phred score", 
            "Coverage"]
        
        ref_name = self.reference_file.split('/')[-1]

        fig = make_subplots(
            rows=6, 
            cols=1, 
            vertical_spacing=0.05, 
            shared_xaxes=True,
            subplot_titles=titles,
            row_heights=heights)

        for base, col in NT_COLOR.items(): 
            fig.add_trace(
                go.Scatter(
                    name=base,
                    y=self.distrib[base],
                    marker_line_width=0,
                    marker_color=col,
                    mode="lines",
                ), 1,1)

        fig.add_trace(
            go.Bar(
                y=[1]*self.consensus.shape[0],
                text=self.consensus["consensus_bp"].values,
                textposition="inside",
                marker=dict(
                    line_width = 0,
                    color=[NT_COLOR[x.upper()] for x in self.consensus["consensus_bp"].values]),
            ), 2, 1)

        fig.add_trace(
            go.Bar(
                y=[1]*self.consensus.shape[0],
                text=list(self.reference_sequence),
                textposition="inside",
                marker=dict(
                    line_width=0,
                    color=[NT_COLOR[x] for x in self.reference_sequence])
            ), 3, 1)
        
        fig.add_trace(
            go.Scatter(
                y=self.consensus["consensus_confidence"].values,
                line=dict(color="#008000"), #074050
                mode='lines',
            ), 4, 1)
        
        fig.add_trace(
            go.Bar(
                y=self.consensus["mean_phred"].values,
                marker=dict(
                    line_width=0,
                    color=self.consensus["mean_phred"].values,
                    colorscale=colscale,
                    cmin=10,
                    cmax=20)
            ), 5, 1)

        fig.add_trace(
            go.Bar(
                y=self.consensus["coverage"].values,
                marker=dict(
                    line_width = 0,
                    color=self.consensus["coverage"].values,
                    colorscale=colscale,
                    cmin=0,
                    cmax=200),
            ), 6, 1)

        fig.update_annotations(font_size=14)

        fig.update_layout(
            title_text=f"ONT Sequencing Report<br>Barcode: <b>{self.sample_id}</b>, Reference: <b>{ref_name}</b><br>",
            title_font_size=20,
            title_x=0.5,
            margin=dict(l=40, r=40, b=20, t=120),
            plot_bgcolor = "white",
            xaxis6_rangeslider_visible=True, 
            xaxis6_rangeslider_thickness=.05,
            yaxis2_showticklabels=False,
            yaxis3_showticklabels=False,
            yaxis4_gridcolor = "rgba(0,0,0,.08)",
            yaxis5_gridcolor = "rgba(0,0,0,.08)",
            yaxis6_gridcolor = "rgba(0,0,0,.08)",
            yaxis1_tickformat="2%",
            bargap=0.0,
            bargroupgap=0.0,
            showlegend=False,
            yaxis=dict(
                gridcolor = "rgba(0,0,0,.08)",
                showline = True, 
                linewidth = 1),
            xaxis=dict(
                showline=True,
                tickmode='array',
                tickvals=list(range(self.consensus.shape[0])),
                ticktext=list(self.consensus["consensus_bp"].values))
            )

        fig.write_html(self.output)
    


if __name__ == "__main__":
    pass