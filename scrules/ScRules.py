# -*-coding:utf-8-*-
"""
@Author: Phantom
@Time:2023/9/4 7:59
@Email: 2909981736@qq.com
"""
import os
import fim
import igraph as ig
import pandas as pd
from time import time
import matplotlib.pyplot as plt
import requests
from .utils.DataLoader import load_sc_transactions_data
from .utils.eval_calculate import calculate_eval_vectorization
from .utils.RegNet import RegNet
from .utils.draw import draw_marker_gene_expression
from collections import Counter


class ScRules:
    def __init__(self, datafile):
        """
        Initialization method
        :param datafile: a file that looks like this:
                item1, item2, item3,...
                item2,item6,item4,item5,...
                item3, item4, item5,...
                ...
        """
        self.data = load_sc_transactions_data(datafile)

    def fpgrowth(self, supp=-4, conf=50, algo='c', other_evaluation=False):
        """
        Calculate association rules using the fp-growth algorithm
        :param supp:minimum support of an item set         (default: -4)
                    (positive: percentage, negative: absolute number)
        :param conf:minimum confidence of an assoc. rule   (default: 50%)
        :param algo:algorithm variant to use               (default: c)
                s     simple     simple  tree nodes (only link and parent)
                c     complex    complex tree nodes (children and siblings)
                d     single     top-down processing on a single prefix tree
                t     topdown    top-down processing of the prefix trees
                Variant d does not support closed/maximal item set mining.
        :param other_evaluation: conviction,leverage,zhangs_metric,chi_square   (default: False)
        :return:rules_df:data frame of association rules
        """
        t1 = time()
        report = 'xyclQ'
        target = 'r'
        zmin = 2
        zmax = 2
        fpgrowth_results = fim.fpgrowth(self.data, target=target, report=report, supp=supp, zmin=zmin, zmax=zmax,
                                        conf=conf, algo=algo)
        rules_df = {'antecedents': [], 'consequents': [], 'antecedents_support': [], 'consequents_support': [],
                    'confidence': [], 'lift': [], 'total_transactions': []}
        for fpr in fpgrowth_results:
            rules_df['antecedents'].append(fpr[1][0])
            rules_df['consequents'].append(fpr[0])
            rules_df['antecedents_support'].append(fpr[2])
            rules_df['consequents_support'].append(fpr[3])
            rules_df['confidence'].append(fpr[4])
            rules_df['lift'].append(fpr[5])
            rules_df['total_transactions'].append(fpr[6])

        rules_df = pd.DataFrame(rules_df)

        if other_evaluation:
            rules_df = calculate_eval_vectorization(rules_df)

        print('The association rule calculation is completed and takes:', time() - t1, 's')

        return rules_df

    def rules_match(self, rules_df, matched_network='reactome'):
        """
        The method of matching with the existing human genome regulatory
        network can be matched with three existing regulatory networks,
        namely: consensus, reactome, regnewworks.
        :param rules_df:The data frame obtained by the fpgrowth method
        :param matched_network:Optional matching network, the value can
                be: consensus, reactome(default), regnewworks
        :return:Matched association rule data frame
        """
        t1 = time()
        try:
            package_path = os.path.dirname(os.path.abspath(__file__))
            matched_network_path = os.path.join(package_path, "data", f"{matched_network}.tsv")

            if not os.path.exists(matched_network_path):
                url = f'http://eptohubgene.litjxxy.com/downloads/{matched_network}.tsv'
                header = {
                    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/116.0.0.0 Safari/537.36 Edg/116.0.1938.69'}
                resp = requests.get(url, headers=header)
                if resp.status_code == 200:
                    with open(matched_network_path, 'wb') as file:
                        file.write(resp.content)
                        print(f'{matched_network}.tsv downloads successfully')
                else:
                    print(f'{matched_network}.tsv downloads failed')
                    return

            matched_network_df = pd.read_csv(matched_network_path, delimiter='\t', header=None,
                                             names=['antecedents', 'consequents'])
        except Exception as e:
            print(str(e))
            print("Please enter one of ‘consensus’, ‘reactome’, ‘regnewworks’ for the matched_network parameter")
            print("Please turn off proxy")
            return

        matched_rules_df = pd.merge(rules_df, matched_network_df, on=['antecedents', 'consequents'], how='inner')

        print('Matching of known regulatory relationships is completed, which takes:', time() - t1, 's')

        return matched_rules_df

    def not_merched_rules(self, rules_df, matched_rules_df):
        """
        Mining unknown association rules that do not appear in regulatory networks
        :param rules_df:data frame of association rules
        :param matched_rules_df:Matched association rule data frame
        :return:No matching association rule data frame
        """

        t1 = time()

        merge_df = pd.merge(rules_df, matched_rules_df, on=['antecedents', 'consequents'], how='outer',
                            indicator=True)

        not_merched_rules_df = merge_df.query("_merge != 'both'")
        not_merched_rules_df = not_merched_rules_df.drop(columns=['_merge'])
        to_be_drop = list()
        to_be_rename = dict()
        for col_name in not_merched_rules_df.columns:
            if col_name[-2:] == "_y":
                to_be_drop.append(col_name)
            elif col_name[-2:] == "_x":
                to_be_rename[col_name] = col_name[:-2]

        not_merched_rules_df = not_merched_rules_df.drop(columns=to_be_drop)
        not_merched_rules_df = not_merched_rules_df.rename(columns=to_be_rename)

        print('Unknown regulatory relationship matching is completed, which takes time:', time() - t1, 's')

        return not_merched_rules_df

    def mining_regulatory_networks_by_tree_digging(self, rules_df, root_genes, digging_depth=3, with_weights=False,
                                                   weights_type=None):
        """
        Use A in the (A, B) rule as the root node, B as the child node, dig down,
        the next layer uses B as the root node, and the node (B, C) associated with B as the child node,
        continue down Mining to obtain the regulatory network with B as the starting point.
        In addition, the starting point node can be multiple and stored in a list.
        :param rules_df:data frame of association rules
        :param root_genes:genes as root node, a list
        :param digging_depth:The depth to dig down, the default is 3
        :param with_weights:Whether to set weight
        :param weights_type:Type of weight
        :return:The data frame of the control network can also be used as the edge of the directed graph.
        """
        t1 = time()
        if with_weights:
            if weights_type is None:
                rules_tuple_list = list(rules_df[['antecedents', 'consequents', 'lift']].to_records(index=False))
                rules_tuple_list = [tuple(record) for record in rules_tuple_list]
            else:
                rules_tuple_list = list(rules_df[['antecedents', 'consequents', weights_type]].to_records(index=False))
                rules_tuple_list = [tuple(record) for record in rules_tuple_list]
        else:
            rules_tuple_list = list(rules_df[['antecedents', 'consequents']].to_records(index=False))
            rules_tuple_list = [tuple(record) for record in rules_tuple_list]

        root_genes = set(root_genes)
        regulatory_networks = set()

        for i in range(digging_depth):
            next_root_genes = set()
            for root_gene in root_genes:
                for rule_tuple in rules_tuple_list:
                    if rule_tuple[0] == root_gene:
                        regulatory_networks.add(rule_tuple)
                        next_root_genes.add(rule_tuple[1])
            root_genes = next_root_genes.copy()

        if with_weights:
            if weights_type is None:
                regulatory_networks_df = pd.DataFrame(list(regulatory_networks), columns=['Source', 'Target', 'lift'])
            else:
                regulatory_networks_df = pd.DataFrame(list(regulatory_networks),
                                                      columns=['Source', 'Target', weights_type])
        else:
            regulatory_networks_df = pd.DataFrame(list(regulatory_networks), columns=['Source', 'Target'])

        print('The regulatory network mining through tree mining is completed and takes time:', time() - t1)
        return regulatory_networks_df

    def mining_regulatory_networks_by_graph_traversal(self, rules_df, start_gene):
        """
        Starting from a certain gene and not specifying the depth, we look for the deepest regulatory network.
        :param rules_df:data frame of association rules
        :param start_gene:A string of gene names
        :return:A two-dimensional list, each row is a control relationship
        """
        t1 = time()
        regnet = RegNet()
        for rule_df in rules_df:
            regnet.add_edge(rule_df[0], rule_df[1])
        deepest_regulatory_networks_by_start_gene = regnet.find_regulatory_network(start_gene)
        print(deepest_regulatory_networks_by_start_gene)

        print(
            f'Mining the deepest regulatory network starting from the {start_gene} gene '
            f'through graph traversal is completed, time-consuming:',
            time() - t1)

        return deepest_regulatory_networks_by_start_gene

    def save_rules(self, rules_df, savepath):
        """
        Save the association rule data frame to the specified directory
        :param rules_df:data frame of association rules, As long as it is a data frame
        :param savepath: The path to save the csv must end with the csv file name
        :return:
        """
        rules_df.to_csv(savepath, index=False)

        print('Saved successfully')

    def draw_regulatory_networks(self, regulatory_networks_df, savepath, directed=True, graph_layout='kk',
                                 figsize=(30, 30), bbox=(3000, 3000)):
        """
        Draw a regulatory network diagram
        :param regulatory_networks_df:A data frame of association rules or regulatory networks.
                The first two columns are the previous and subsequent items respectively.
        :param savepath:save route, png
        :param directed:Is there any direction
        :param graph_layout:drawing layout
        :param figsize:Image size
        :param bbox:Specify the drawing area bounding box
        :return:none
        """
        print("Warning: The drawing effect is not good. It is recommended to export the result generated by "
              "mining_regulatory_networks_by_tree_digging as csv and send it to Cytoscape for the visualization"
              " of the directed graph: https://cytoscape.org/")
        graph = ig.Graph(directed=directed)
        vertices = list(set(regulatory_networks_df.iloc[:, :2].values.flatten()))
        edges = list(zip(regulatory_networks_df.iloc[:, 0], regulatory_networks_df.iloc[:, 1]))
        graph.add_vertices(vertices)
        graph.vs['label'] = vertices
        graph.add_edges(edges)
        fig, ax = plt.subplots(figsize=figsize)
        layout = graph.layout(graph_layout)
        ig.plot(graph, layout=layout, bbox=bbox, target=ax)  # matplotlib backend
        plt.savefig(savepath)
        print("Picture saved successfully")
        plt.show()

    def compute_marker_gene_ratio(self, marker_genes, savepath):
        """
        Draw a bar graph of the expression ratio of the marker gene in the cell data set
        :param marker_genes: Cell type marker genes,a list
        :param savepath:The path to save the picture
        :return:none
        """
        cell_count = len(self.data)
        counter = Counter()
        for row in self.data:
            counter.update(row)
        marker_genes_Expression_dic = {'marker_genes': [], 'Expression': []}
        for marker_gene in marker_genes:
            Expression = counter[marker_gene]
            marker_genes_Expression_dic['marker_genes'].append(marker_gene)
            marker_genes_Expression_dic['Expression'].append(round(Expression/cell_count, 3))

        draw_marker_gene_expression(marker_genes_Expression_dic, savepath)
        print("Picture saved successfully")
