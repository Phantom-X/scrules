# -*-coding:utf-8-*-
"""
@Author: Phantom
@Time:2023/9/7 16:54
@Email: 2909981736@qq.com
"""
import pandas as pd
import matplotlib.pyplot as plt


def draw_marker_gene_expression(marker_genes_Expression_dic, savepath):
    marker_genes_Expression_df = pd.DataFrame(marker_genes_Expression_dic)

    df_sorted = marker_genes_Expression_df.sort_values(by='Expression', ascending=True)

    num_bars = len(df_sorted)
    if num_bars > 300:
        df_sorted = df_sorted[-100:]
        num_bars = 100
    fig_width = max(8, num_bars * 0.05)
    fig_height = max(6, num_bars * 0.15)
    plt.figure(figsize=(fig_width, fig_height))

    bar_width = 0.8
    num_bars = len(df_sorted)

    plt.barh(range(num_bars), df_sorted['Expression'], height=bar_width, align='center', color='blue', alpha=0.8)

    plt.yticks(range(num_bars), df_sorted['marker_genes'])
    plt.xlabel('Expression')
    plt.ylabel('marker_gene')
    plt.title(' Expression of Marker Genes (Sorted)')

    x_positions = [0.2, 0.4, 0.6, 0.8, 1.0]
    for x in x_positions:
        plt.axvline(x=x, color='gray', linestyle='dashed')
    plt.tight_layout()
    plt.savefig(savepath)
    plt.show()
