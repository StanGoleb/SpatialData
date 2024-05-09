import streamlit as st
import matplotlib.pyplot as plt
import os
import numpy as np
from itertools import groupby

from clusters import clusterize
from colors import color_types


def visualize(input_path, dict_path, selected_types, cluster_num, max_neighbor_distance, min_graph_size):
 
    # Graf z pokolorowanymi typami komórek na próbce z biopsji
    xs, ys, colors, legend_handles, color_dict = color_types(input_path, dict_path) #patrz colors.py
    fig1, ax1 = plt.subplots(figsize=(14, 10))
    ax1.scatter(xs, ys, c=colors, alpha=0.5)
    ax1.set_title('Scatter plot')
    ax1.set_xlabel('nucleus.x')
    ax1.set_ylabel('nucleus.y')
    ax1.legend(handles=legend_handles)
    st.pyplot(fig1)

    plot_data, components_vectors, component_clusters = clusterize(selected_file, selected_types, cluster_num, max_neighbor_distance, min_graph_size) #patrz clusters.py
    
    # Graf z zaznaczonymi na próbce z biopsji grafami spójnymi z otoczką. Pokolorowane według klastra.
    fig2, ax2 = plt.subplots(figsize=(14, 10))
    for cluster_ind, cluster in sorted(plot_data.items()):        
        ax2.scatter(cluster[0], cluster[1], label=f"Cluster {cluster_ind}")
        
    ax2.set_xlabel("nucleus.x")
    ax2.set_ylabel("nucleus.y")
    ax2.legend()
    st.pyplot(fig2)

    # Wykres słupkowy %-udziału typów komórek w danym grafie spójnym (rozszerzonym o otoczkę).
    fig3, ax3 = plt.subplots(figsize=(14, 10))

    hist_data = np.array(components_vectors)
    x = np.arange(len(hist_data))
    for i in range(hist_data.shape[1]):
        ax3.bar(x, hist_data[:, i], bottom=np.sum(hist_data[:, :i], axis=1), color=list(color_dict.values())[i], label=list(color_dict.keys())[i])

    ax3.set_xlabel("Grafy spójne z otoczką")
    ax3.set_ylabel("Udział procentowy komórek")
    ax3.legend(loc='upper right')
    st.pyplot(fig3)

 
    # Wykres posortowanych słupków.
    fig4, ax4 = plt.subplots(figsize=(14, 10))

    # Grupowanie klastrów - słupki w jednym klastrze są stawiane obok siebie.
    original_indices = list(range(len(component_clusters)))
    sorted_data = [[vector, original_index] for _,vector, original_index in sorted(zip(component_clusters,components_vectors, original_indices))]
    sorted_vectors, sorted_indices = [i[0] for i in sorted_data], [i[1] for i in sorted_data]
    hist_data = np.array(sorted_vectors)

    x = np.arange(len(hist_data))
    for i in range(hist_data.shape[1]):
        ax4.bar(x, hist_data[:, i], bottom=np.sum(hist_data[:, :i], axis=1), color=list(color_dict.values())[i], label=list(color_dict.keys())[i])

    for list_ind, sorted_index in enumerate(sorted_indices):
        if component_clusters[sorted_index] != component_clusters[sorted_indices[list_ind-1]] and list_ind != 0:
            ax4.axvline(x=list_ind-0.5, color='black')

    # Dane słupki zachowują swoje indeksy pomimo przesortowania.
    ax4.set_xticks(x)
    ax4.set_xticklabels(sorted_indices, rotation=45)
    st.pyplot(fig4)



# ========POCZĄTEK PROGRAMU========

current_dir = os.getcwd()

# Ścieżka do katalogu plików z biopsją.
if_data_path = os.path.join(current_dir, 'if_data')

# Pole wyboru opcji IF1-3.
selected_option = st.selectbox("Choose an option", ["IF1", "IF2", "IF3"])
dict_path = os.path.join(current_dir, f'{selected_option}_phen_to_cell_mapping.csv')

# Pole wyboru nazwy pliku z biopsją.
files = [f for f in os.listdir(if_data_path) if os.path.isfile(os.path.join(if_data_path, f)) and f.endswith('.csv') and selected_option in f]
files.sort()
selected_file = st.selectbox("Choose a CSV file", files)
selected_file = os.path.join(if_data_path, selected_file)

# Pola wyboru typów komórek, które będą tworzyć rdzeń grafów spójnych.
cell_types = ["other", "CD15+Tumor", "CD15-Tumor", "Tcell", "Bcell", "BnTcell", "Neutrophil", "Macrophage", "DC"]
checkbox_status = {cell_type: st.checkbox(cell_type) for cell_type in cell_types}
selected_cell_types = [cell_type for cell_type, status in checkbox_status.items() if status]

# Pole do wybrania ilości klastrów.
cluster_num = st.number_input("Wprowadź liczbę klastrów", 1, 100)

# Pole do wybrania minimalnej wielkości początkowego grafu spójnego.
min_graph_size = st.number_input("Wprowadź minimalną wielkość początkowych grafów spójnych wybranych komórek", 1, 100, 20)

# Pole do wybrania maksymalnej odległości pomiędzy sąsiadami w grafach.
max_neighbor_distance = st.number_input("Wprowadź maksymalna odległość pomiędzy wybranymi komórkami", 1, 100, 30)

run_button = st.button("Start")

# Po naciśnięciu przycisku odpalić wizualizację.
if run_button:
    if selected_file is not None:
        visualize(selected_file, dict_path, selected_cell_types, cluster_num, max_neighbor_distance, min_graph_size)
