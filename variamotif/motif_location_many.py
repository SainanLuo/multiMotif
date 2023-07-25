import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

def plot_motifs_to_single_chart(file_path, display_both_directions=False):
    df = pd.read_csv(file_path, sep='\t', header=None, dtype={0: str, 1: str, 2: int, 3: int, 4: str}, usecols=[0, 1, 2, 3, 4], skiprows=1)

    # 将结果按照Sequence ID进行分组
    grouped_df = df.groupby(0)  # 使用列索引（从0开始）进行分组

    num_sequences = len(grouped_df)
    height_interval = 0.5  # 方块高度
    line_length = 15  # 线的长度
    arrow_length = 0.15  # 箭头长度
    shift_amount = 0.15  # 移动的距离

    # 通过减小高度来设置较低的纵横比
    fig, ax = plt.subplots(figsize=(line_length, num_sequences * 0.5))

    # 删除图表的坐标轴和边框线。
    ax.axis('off')

    # 自定义颜色列表
    color_list = ['#008000','#FF0000', '#0000FF', '#FF00FF', '#00FFFF', '#FFFF00', '#800000', '#00FF00', '#000080', '#808000', '#800080', '#008080', '#808080', '#C0C0C0', '#FFA500']

    # 使用字典映射颜色到motif
    color_map = {}

    # 绘制每个Sequence ID的线和箭头
    for idx, (seq_id, seq_group) in enumerate(grouped_df):
        y_position = idx * 2
        ax.plot([0, line_length], [y_position, y_position], color='black')
        ax.text(-0.5, y_position, seq_id, ha='right', va='center', fontsize=10)  # 左边显示Sequence ID

        for _, row in seq_group.iterrows():
            start = row[2]
            end = row[3]
            direction = row[4]
            motif = row[1]

            if not display_both_directions and direction == '-':
                continue  # Skip if only displaying one direction and the current row is for the negative strand

            motif_length = end - start

            # 获取或生成motif对应的颜色
            if motif not in color_map:
                color_map[motif] = color_list[len(color_map) % 15]

            color = color_map[motif]

            # 绘制矩形方块
            if direction == '+':
                rectangle = plt.Rectangle((start / 400 * line_length, y_position), (end - start) / 400 * line_length, height_interval, facecolor=color, edgecolor='none')
            else:  # direction == '-'
                rectangle = plt.Rectangle(((start - shift_amount) / 400 * line_length, y_position), (end - start) / 400 * line_length, height_interval, facecolor=color, edgecolor='none')

            ax.add_patch(rectangle)

            # 绘制箭头
            if direction == '+':
                arrow_tail_x = (start + motif_length - shift_amount) / 400 * line_length
                arrow_head_x = arrow_tail_x + arrow_length
            else:  # direction == '-'
                arrow_tail_x = (start + shift_amount) / 400 * line_length
                arrow_head_x = arrow_tail_x - arrow_length

            arrow_polygon = np.array([[arrow_head_x, y_position + height_interval / 2], [arrow_tail_x, y_position], [arrow_tail_x, y_position + height_interval]])
            arrow = plt.Polygon(arrow_polygon, closed=True, edgecolor=color, facecolor=color)
            ax.add_patch(arrow)

    ax.set_xlim(0, line_length)
    ax.set_ylim(-1, num_sequences * 2)
    ax.set_xticks(range(0, line_length + 1, 1))
    ax.set_xlabel('Position (cm)')
    ax.set_title('Motif Search Results')
    ax.set_yticks([])  # 隐藏y轴刻度

    # 生成图例
    legend_elements = [Patch(facecolor=color, edgecolor='black', label=motif) for motif, color in color_map.items()]
    ax.legend(handles=legend_elements, title='Motif', bbox_to_anchor=(1.05, 1), loc='upper left')

    output_file = 'motif_location.pdf'

    plt.tight_layout()
    plt.savefig(output_file)
    #plt.show()

def main():
    parser = argparse.ArgumentParser(description='Generate motif location chart from a table file.')
    parser.add_argument('-t', '--table', dest='table_file', required=True, help='Input table file.')
    parser.add_argument('-d', '--display_both_directions', action='store_true', help='Display motifs from both + and - strands.')
    args = parser.parse_args()

    plot_motifs_to_single_chart(args.table_file, args.display_both_directions)

if __name__ == '__main__':
    main()

