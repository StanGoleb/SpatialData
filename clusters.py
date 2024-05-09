import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt

import networkx as nx

from sklearn.neighbors import radius_neighbors_graph
from sklearn.cluster import AgglomerativeClustering

def clusterize(input_path, chosen_cells, cluster_size):

    if1_cell_types = ["other", "CD15+Tumor", "CD15-Tumor", "Tcell", "Bcell", "BnTcell", "Neutrophil", "Macrophage", "DC"]

    all_df = pd.read_csv(input_path)
    chosen_df = all_df[all_df['cell_type'].isin(chosen_cells)]

    # koordynaty
    all_xy = all_df[['nucleus.x', 'nucleus.y']].to_numpy()
    chosen_xy = chosen_df[['nucleus.x', 'nucleus.y']]

    # sąsiedzi, bliżej niż 30
    all_graph = radius_neighbors_graph(all_xy, 30, mode='distance')
    chosen_graph = radius_neighbors_graph(chosen_xy, 30, mode='distance')

    # grafy spójne
    all_G = nx.convert_matrix.from_scipy_sparse_array(all_graph)
    all_connected_components = list(nx.connected_components(all_G))
    
    chosen_G = nx.convert_matrix.from_scipy_sparse_array(chosen_graph)
    chosen_connected_components = list(nx.connected_components(chosen_G))

    # słownik indeksów chosen
    chosen_to_all = {}
    i = 0
    for row_ind, row in chosen_df.iterrows():
        chosen_to_all[i] = row_ind
        i+= 1
    
    # wybieramy najlepsze grafy spójne
    chosen_connected_components = [component for component in chosen_connected_components if len(component) > 20]

    expanded_components_vectors = []
    expanded_components_ids = []

    for i, component in enumerate(chosen_connected_components):
        
        # rozszerzamy graf spójny o sąsiadów każdej z komórek
        component_w_neighbors = set()
        for chosen in component:
            component_w_neighbors.add(chosen)
            # znajdujemy komórkę z grafu w liście sąsiadów wszystkich typów komórek
            chosen_real_ind = chosen_to_all[chosen]
            chosen_neighbors = [n for n in all_graph[chosen_real_ind].indices]
            for neighbor in chosen_neighbors:
                component_w_neighbors.add(neighbor)
        
        component_w_neighbors_df = all_df[all_df['cell.ID'].isin(component_w_neighbors)]
        expanded_components_ids += [component_w_neighbors]

        # tworzymy wektor %-owego uczestnictwa typów komórek w rozszerzonym grafie spójnym
        cell_type_percent_dict = {cell_type: 0 for cell_type in if1_cell_types}
        for cell_ind, cell in component_w_neighbors_df.iterrows():
            cell_type_percent_dict[cell["cell_type"]] += 1
        
        values_sum = sum(cell_type_percent_dict.values())
        for key in cell_type_percent_dict.keys():
            cell_type_percent_dict[key] = cell_type_percent_dict[key]/max(values_sum, 1)

        cell_type_vector = list(cell_type_percent_dict.values())
        cell_type_keys = list(cell_type_percent_dict.keys()).sort()
        expanded_components_vectors += [cell_type_vector]

    hierarchical_cluster = AgglomerativeClustering(n_clusters=cluster_size, metric='euclidean', linkage='complete')

    if len(expanded_components_vectors) > 1:

        component_clusters = hierarchical_cluster.fit_predict(expanded_components_vectors)
        component_cluster_ids = [i+1 for i in component_clusters]
        
        plot_data = {cluster_id + 1 : [[], []] for cluster_id in range(cluster_size)}
        plot_data[0] = [[], []]
        for cell_id, (x, y) in enumerate(all_xy):
            for component_ind, expanded_component in enumerate(expanded_components_ids):
                if cell_id in expanded_component:

                    plot_data[component_cluster_ids[component_ind]][0].append(x)
                    plot_data[component_cluster_ids[component_ind]][1].append(y)
                    break
                else:
                    plot_data[0][0].append(x)
                    plot_data[0][1].append(y)
    
    else:
        plot_data = {0: [[all_xy[:][0]],[all_xy[:][1]]]}

                
    return plot_data, expanded_components_vectors, component_clusters


# overview danych - wielu pacjentów, jaki jest rozkłąd % zawartości b-cell i t-cell wśród pacjentów. 
# ogólna funkcja sklrean radius neighbors graph która zbuduje graf na tych komórkach dla dowolnego fenotypu (na przykład tylko CK i coś innego)
# zbudować graf sąsiedztwa dla wszystkich komórek, obok graf dla komórek typu B (np. CD20)
#   w dataframe wybieram tylko B-celle (CD20+).
#   najpierw graf sąsiedztwa z parametrem jak blisko mają być na wszystkim i potem drugi zbudowany tylko na punktach z CD20+
    # elementy CD20+ blisko siebie są najprawdopodobniej TLS
#   przejść przez komórki CD20+ funckją connected components. szukamy sąsiadów całego spójnego grafu CD20+

# taki graf spójny razem z sasiadami będzie miał wektor [20% tcelli, 50%bcelli,...]. takie wektory chcemy sklastrować np. na średnio zawsze 60% bcelli. potem zaznaczamy je w obrębie raka


# na tym zrobić histogram  z zaznaczonymi TLS

#github
# raport