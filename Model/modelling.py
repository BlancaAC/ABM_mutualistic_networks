#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 12:04:42 2024

@author: galeanojav
"""

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import networkx as nx


def initial_pollinators_random(dist_pol, n_spe_pol, n_pols, xmin, xmax, ymin, ymax, 
                               classify=True, random_distribution=False):
    """
    Genera la distribución inicial de polinizadores con opciones personalizables.
    
    Parámetros:
    - dist_pol: Distribución de abundancia de polinizadores.
    - n_spe_pol: Número de especies de polinizadores.
    - n_pols: Número total de polinizadores.
    - xmin, xmax, ymin, ymax: Límites espaciales para generar coordenadas.
    - classify: Booleano. Si es True, divide los polinizadores en generalistas y especialistas.
    - random_distribution: Booleano. Si es True, genera una distribución aleatoria de probabilidades para np.random.choice.
    
    Retorna:
    - generalistas: Lista de especies clasificadas como generalistas (si classify=True).
    - df_final: DataFrame con información de los polinizadores.
    """
    
    generalistas = []  # Inicializar lista de generalistas para el caso en que classify=False
    
    if classify:
        # Step 1: Calcular etiquetas de generalista y especialista
        polinizadores = dist_pol.values.reshape(-1, 1)
        n_clusters = 2
        
        # Calcular centros a partir de los datos de distribución de frecuencia
        c1 = np.mean(dist_pol.values[:-2])
        c2 = np.mean(dist_pol.values[2:])
        centroides_iniciales = np.array([[c1], [c2]])
        
        # Aplicar K-Means con ubicaciones iniciales personalizadas
        kmeans = KMeans(n_clusters=n_clusters, init=centroides_iniciales, random_state=42)
        kmeans.fit(polinizadores)

        # Obtener etiquetas de cluster y centros
        etiquetas = kmeans.labels_
        
        # Crear DataFrame intermedio
        df_intermedio = pd.DataFrame(data={'Abundancia': polinizadores.flatten(), 'Etiqueta': etiquetas}, index=dist_pol.index)
        df_intermedio['Etiqueta'] = df_intermedio['Etiqueta'].map({0: 'Especialista', 1: 'Generalista'})

        # Generar lista de generalistas
        generalistas = df_intermedio.loc[df_intermedio['Etiqueta'] == 'Generalista'].index.tolist()

    # Paso 2: Generar distribución de probabilidades para np.random.choice
    if random_distribution:
        probabilities = np.random.dirichlet(np.ones(len(dist_pol)), size=1).flatten()
    else:
        probabilities = dist_pol

    # Paso 3: Seleccionar especies de polinizadores
    pollinators = np.random.choice(dist_pol.index, n_pols, p=probabilities)
    
    # Paso 4: Generar atributos de los polinizadores
    indiceList = np.arange(1000, 1000 + n_pols)
    xList = np.round(xmin + np.random.rand(n_pols) * (xmax - xmin), decimals=3)
    yList = np.round(ymin + np.random.rand(n_pols) * (ymax - ymin), decimals=3)
    
    # Generar radios basados en generalistas/especialistas si classify=True
    if classify:
        radList = np.where(np.isin(pollinators, generalistas), 15, 5)
    else:
        radList = np.random.gamma(10, 2, size=n_pols)  # Radio aleatorio como alternativa

    # Paso 5: Construir DataFrame final
    df_final = pd.DataFrame({
        'Pol_id': indiceList,
        'Specie': pollinators,
        'x': xList,
        'y': yList,
        'Radius': radList
    }, index=indiceList)
    
    # Añadir columna 'Tipo' si classify=True
    if classify:
        df_final['Tipo'] = np.where(df_final['Specie'].isin(generalistas), 'Generalista', 'Especialista')
    else:
        df_final['Tipo'] = 'Sin clasificar'

    return generalistas, df_final




def initial_pollinators(dist_pol,n_spe_pol, n_pols, xmin, xmax, ymin, ymax):
    
    # Step 1: Calculate generalist and specialist labels tags
    polinizadores = dist_pol.values.reshape(-1, 1)
    n_clusters = 2
    
    # Calcular centros a partir de los datos de distribución de frecuencia
    c1 = np.mean(dist_pol.values[:-2])
    c2 = np.mean(dist_pol.values[2:])
    centroides_iniciales = np.array([[c1], [c2]])
    
    # Aplicar K-Means con ubicaciones iniciales personalizadas
    kmeans = KMeans(n_clusters=n_clusters, init=centroides_iniciales, random_state=42)
    kmeans.fit(polinizadores)

    # Obtener etiquetas de cluster y centros
    etiquetas = kmeans.labels_
    
    # Paso 2: Crear DataFrame intermedio
    df_intermedio = pd.DataFrame(data={'Abundancia': polinizadores.flatten(), 'Etiqueta': etiquetas}, index=dist_pol.index)
    df_intermedio['Etiqueta'] = df_intermedio['Etiqueta'].map({0: 'Especialista', 1: 'Generalista'})

    # Paso 3: Generar distribución
    generalistas = df_intermedio.loc[df_intermedio['Etiqueta'] == 'Generalista'].index.tolist()
    
    #pollinators = np.random.choice(dist_pol.index, n_pols, p=dist_pol)
    pollinators = np.random.choice(dist_pol.index, n_pols, p=dist_pol)


    indiceList = np.arange(1000, 1000 + n_pols)
    xList = np.round(xmin + np.random.rand(n_pols) * (xmax - xmin), decimals=3)
    yList = np.round(ymin + np.random.rand(n_pols) * (ymax - ymin), decimals=3)
    #radList = np.where(np.isin(pollinators, generalistas), pd.Series(np.random.gamma(15, 2, size=n_pols)),
                       #pd.Series(np.random.gamma(4, 2, size=n_pols)))
    
    radList = np.where(np.isin(pollinators, generalistas), 15,5)

    
    # Paso 5: Construir DataFrame final
    df_final = pd.DataFrame({
        'Pol_id': indiceList,
        'Specie': pollinators,
        'x': xList,
        'y': yList,
        'Radius': radList
    }, index=indiceList)
    
     # Añadir columna 'Tipo' indicando si la especie es generalista o especialista
    df_final['Tipo'] = np.where(df_final['Specie'].isin(generalistas), 'Generalista', 'Especialista')
    
    return generalistas,df_final

# Network Functions

def initial_network(pollinators, plants):
    """Initial Graph is a complete digraph between the list of pollinators and the plants.
    The network is bipartite and with weight =0"""
    B = nx.DiGraph()
    
    B.add_nodes_from(pollinators, bipartite=0)
    B.add_nodes_from(plants, bipartite=1)
    B.add_weighted_edges_from((u, v, 0) for u in pollinators for v in plants)
    
    return B

def remove_zero(B):
    """Remove edges with zero weight from the graph."""
    
    edge_list = [(u, v) for (u, v, w) in B.edges(data=True) if w['weight'] == 0]
    B.remove_edges_from(edge_list)


def degree_dist(df):
    """Calculate the degree distribution of each set of nodes 
    in a bipartite network.
    Input .- Pandas DataFrame 
    Output.- Tuple of 2 Series representing degree distributions for every set of nodes"""
    
    # Calculate degree distribution
    plant_degree = df.astype(bool).sum(axis=1)
    pol_degree = df.astype(bool).sum(axis=0)
    
    return pol_degree, plant_degree

## Evolution programs

def plant_pol(NewAgent, env):
    """Find plants within a certain radius.
    The agents in the returned list are all plants."""
    
    x, y, radio = NewAgent.x, NewAgent.y, NewAgent.radioAccion
    neighbors = [nb for nb in env.plant_list if (x - nb.x)**2 + (y - nb.y)**2 < radio**2]
    
    return neighbors


def update(envpol,evenp,B,xmin,xmax,ymin,ymax):
    indAgent = np.random.choice(len(envpol.pol_list))
    NewAgent = envpol.pol_list[indAgent]

    NewAgent.random_xy_pol(xmin,xmax,ymin,ymax)
    neigh=plant_pol(NewAgent,evenp)
    
    
    if neigh:
        #B[NewAgent.id][neigh[0].id]['weight'] += 1
        # randomly select one neighbor
        #selected_neigh = np.random.choice(neigh)
        #B[NewAgent.id][selected_neigh.id]['weight'] += 1
        
        # calculate the distances to each neighbor
        distances = [(n, ((NewAgent.x - n.x) ** 2 + (NewAgent.y - n.y) ** 2) ** 0.5) for n in neigh]
        # select the closest neighbor
        selected_neigh = min(distances, key=lambda x: x[1])[0]
        B[NewAgent.id][selected_neigh.id]['weight'] += 1


def update_totalinks(tlink, envpol,evenp, B,xmin,xmax,ymin,ymax):
    """This function is not tested yet"""
    total_links = sum(nx.get_edge_attributes(B,'weight').values())
    
    while total_links < tlink:
        update(envpol,evenp, B,xmin,xmax,ymin,ymax)
        total_links += 1  # increase the total links after each update
        # print(total_links)





