# -*-coding:utf-8-*-
"""
@Author: Phantom
@Time:2023/9/5 16:56
@Email: 2909981736@qq.com
"""

class RegNet:
    def __init__(self):
        self.adjacency_list = {}
        self.result = []

    def add_edge(self, gene1, gene2):
        if gene1 in self.adjacency_list:
            self.adjacency_list[gene1].append(gene2)
        else:
            self.adjacency_list[gene1] = [gene2]

    def find_regulatory_network(self, seed_gene):
        if seed_gene not in self.adjacency_list:
            return []
        else:
            visited = []
            self._dfs(seed_gene, visited)
            return self.result

    def _dfs(self, gene, visited):
        visited.append(gene)
        if gene in self.adjacency_list:
            for next_gene in self.adjacency_list[gene]:
                if next_gene not in visited:
                    self._dfs(next_gene, visited)
                else:
                    visited.append(next_gene)
                    self.result.append(visited.copy())
                    visited.pop()
            visited.pop()
        else:
            self.result.append(visited.copy())
            visited.pop()
