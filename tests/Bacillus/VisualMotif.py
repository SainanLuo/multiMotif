import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Set non-interactive backend to prevent X11 window creation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

def plot_motifs_to_single_chart(file_path, output_file, display_both_directions=False):
    df = pd.read_csv(file_path, sep='\t', header=None, dtype={0: str, 1: int, 2: str, 3: int, 4: int, 5: str}, usecols=[0, 1, 2, 3, 4, 5], skiprows=1, names=["Sequence ID", "Length", "Motif", "Start", "End", "Strand"])

    grouped_df = df.groupby("Sequence ID")

    num_sequences = len(grouped_df)
    height_interval = 0.5
    max_length = df["Length"].max()
    line_length = 15
    multiplier = line_length / max_length

    fig, ax = plt.subplots(figsize=(line_length, num_sequences * 0.5))
    ax.axis('off')

    color_list = ['#008000','#FF0000', '#0000FF', '#FF00FF', '#00FFFF', '#FFFF00', '#800000', '#00FF00', '#000080', '#808000', '#800080', '#008080', '#808080', '#C0C0C0', '#FFA500']
    color_map = {}

    for idx, (seq_id, seq_group) in enumerate(grouped_df):
        y_position = idx * 2
        seq_length = seq_group.iloc[0]["Length"]
        line = seq_length * multiplier
        ax.plot([0, line], [y_position, y_position], color='black')
        ax.text(-0.5, y_position, seq_id, ha='right', va='center', fontsize=10)

        for _, row in seq_group.iterrows():
            start = row["Start"]
            end = row["End"]
            direction = row["Strand"]
            motif = row["Motif"]
            arrow_length = (end-start)*line_length/max_length*0.2
            shift_amount = arrow_length

            if not display_both_directions and direction == '-':
                continue

            motif_length = end - start

            if motif not in color_map:
                color_map[motif] = color_list[len(color_map) % 15]

            color = color_map[motif]

            if direction == '+':
                rectangle = plt.Rectangle((start / seq_length * line, y_position), (end - start) / seq_length * line, height_interval, facecolor=color, edgecolor='none')
            else:  # direction == '-'
                rectangle = plt.Rectangle(((start - shift_amount) / seq_length * line, y_position), (end - start) / seq_length * line, height_interval, facecolor=color, edgecolor='none')

            ax.add_patch(rectangle)

            if direction == '+':
                arrow_tail_x = (start + motif_length - shift_amount) / seq_length * line
                arrow_head_x = arrow_tail_x + arrow_length
            else:  # direction == '-'
                arrow_tail_x = (start + shift_amount) / seq_length * line
                arrow_head_x = arrow_tail_x - arrow_length

            arrow_polygon = np.array([[arrow_head_x, y_position + height_interval / 2], [arrow_tail_x, y_position], [arrow_tail_x, y_position + height_interval]])
            arrow = plt.Polygon(arrow_polygon, closed=True, edgecolor=color, facecolor=color)
            ax.add_patch(arrow)

    ax.set_xlim(0, line_length)  # 设置 x 轴的范围为 0 到 15，表示图中横线的长度
    ax.set_ylim(-1, num_sequences * 2)
    ax.set_xticks(np.arange(0, line_length + 1))  # 设置 x 轴的刻度为整数，从 0 到 15
    ax.set_xlabel('Position (cm)')
    ax.set_title('Motif Search Results')
    ax.set_yticks([])

    legend_elements = [Patch(facecolor=color, edgecolor='black', label=motif) for motif, color in color_map.items()]
    ax.legend(handles=legend_elements, title='Motif', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    
    #Save as PDF 
    output_file_pdf = f'{output_file}.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight',format='pdf')

    #Save as PNG
    output_file_png = f'{output_file}.png'
    plt.savefig(output_file_png, bbox_inches='tight',format='png')


def main():
    parser = argparse.ArgumentParser(description='Generate motif location chart from a table file.')
    parser.add_argument('-t', '--table', dest='table_file', required=True, help='Input table file.')
    parser.add_argument('-r', '--display_both_directions', action='store_true', help='Display motifs from both + and - strands.')
    parser.add_argument('-o', '--output_prefix', dest='output_prefix', required=True, help='Output file prefix.')
    args = parser.parse_args()

    plot_motifs_to_single_chart(args.table_file, args.output_prefix, args.display_both_directions)

if __name__ == '__main__':
    main()

