#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 09:25:54 2024

@author: galeanojav
"""

#from Bio import Phylo
import numpy as np
import networkx as nx
import math
from collections import Counter


def connectance(B):
    num_edges = B.number_of_edges()
    num_nodes_top = len([n for n, d in B.nodes(data=True) if d['bipartite'] == 0])
    num_nodes_bottom = len([n for n, d in B.nodes(data=True) if d['bipartite'] == 1])
    num_possible_edges = num_nodes_top * num_nodes_bottom
    return num_edges / num_possible_edges

def nodf(B):
    # Convertir la red bipartita a una matriz de adyacencia binaria
    top_nodes = [n for n, d in B.nodes(data=True) if d['bipartite'] == 0]
    bottom_nodes = [n for n, d in B.nodes(data=True) if d['bipartite'] == 1]
    adjacency_matrix = nx.bipartite.biadjacency_matrix(B, row_order=top_nodes, column_order=bottom_nodes).toarray()
    
    def nodf_pairs(matrix):
        n_rows, n_cols = matrix.shape
        nestedness = 0
        
        # Calcular NODF para filas
        for i in range(n_rows):
            for j in range(i + 1, n_rows):
                min_deg = min(matrix[i].sum(), matrix[j].sum())
                if min_deg > 0:
                    overlap = np.dot(matrix[i], matrix[j])
                    nestedness += overlap / min_deg
        
        # Calcular NODF para columnas
        for i in range(n_cols):
            for j in range(i + 1, n_cols):
                min_deg = min(matrix[:, i].sum(), matrix[:, j].sum())
                if min_deg > 0:
                    overlap = np.dot(matrix[:, i], matrix[:, j])
                    nestedness += overlap / min_deg
        
        return nestedness
    
    # Normalizar el NODF por el número de pares de filas y columnas comparados
    n_rows, n_cols = adjacency_matrix.shape
    max_pairs = (n_rows * (n_rows - 1) / 2) + (n_cols * (n_cols - 1) / 2)
    nodf_value = (2 * nodf_pairs(adjacency_matrix)) / max_pairs
    
    return nodf_value



#def binary_nestedness(B):
#    # Convertir la red bipartita a una matriz de adyacencia
#    top_nodes = [n for n, d in B.nodes(data=True) if d['bipartite'] == 0]
#    bottom_nodes = [n for n, d in B.nodes(data=True) if d['bipartite'] == 1]
#    adjacency_matrix = nx.bipartite.biadjacency_matrix(B, row_order=top_nodes, column_order=bottom_nodes).toarray()
    
    # Calcular nestedness
#    nestedness = Phylo.BaseTree.TreeMixin.calculate_birnbaum_nestedness(adjacency_matrix)
#    return nestedness

def assortativity(B):
    return nx.degree_assortativity_coefficient(B)

# def eigenvector_centralization(B):
#     eigenvector_centrality = nx.eigenvector_centrality_numpy(B)
#     max_centrality = max(eigenvector_centrality.values())
#     centralization = sum(max_centrality - val for val in eigenvector_centrality.values()) / (len(B.nodes()) - 1)
#     return centralization


def eigenvector_centralization(B, max_iter=1000, tol=1e-06):
    try:
        eigenvector_centrality = nx.eigenvector_centrality(B, max_iter=max_iter, tol=tol)
        max_centrality = max(eigenvector_centrality.values())
        centralization = sum(max_centrality - val for val in eigenvector_centrality.values()) / (len(B.nodes()) - 1)
        return centralization
    except nx.PowerIterationFailedConvergence:
        print("Eigenvector centrality did not converge.")
        return None

# def interaction_evenness(B):
#     edges = B.edges()
#     interaction_counts = Counter(edges)
#     total_interactions = sum(interaction_counts.values())
#     proportions = [count / total_interactions for count in interaction_counts.values()]
#     shannon_diversity = -sum(p * math.log(p) for p in proportions if p > 0)
#     evenness = shannon_diversity / math.log(len(interaction_counts))
#     return evenness



def interaction_evenness(B):
    # Contar las interacciones (aristas) y sus pesos
    interaction_counts = Counter()
    for u, v, data in B.edges(data=True):
        weight = data.get('weight', 1)  # Si no hay peso, usar 1 como valor predeterminado
        interaction_counts[(u, v)] += weight
    
    # Calcular el total de interacciones
    total_interactions = sum(interaction_counts.values())
    
    # Calcular las proporciones de cada interacción
    proportions = [count / total_interactions for count in interaction_counts.values()]
    
    # Calcular la diversidad de Shannon
    shannon_diversity = -sum(p * math.log(p) for p in proportions if p > 0)
    
    # Calcular la evenness (equitatividad)
    evenness = shannon_diversity / math.log(len(proportions)) if len(proportions) > 1 else 0
    
    return evenness

#print("Connectance:", connectance(B))
#print("Nestedness:", binary_nestedness(B))
#print("Assortativity:", assortativity(B))
#print("Eigenvector Centralization:", eigenvector_centralization(B))
#print("Interaction Evenness:", interaction_evenness(B))

