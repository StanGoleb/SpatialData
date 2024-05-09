import streamlit as st
import matplotlib.pyplot as plt
import os
import numpy as np
from itertools import groupby

from clusters import clusterize
from colors import color_types


def visualize(input_path, dict_path, cluster_num, selected_types):
 
    # Pobierz i wykreśl graf z pokolorywanymi typami komórek
    xs, ys, colors, legend_handles, color_dict = color_types(input_path, dict_path)
    
    fig1, ax1 = plt.subplots(figsize=(14, 10))
    ax1.scatter(xs, ys, c=colors, alpha=0.5)
    ax1.set_title('Scatter plot')
    ax1.set_xlabel('nucleus.x')
    ax1.set_ylabel('nucleus.y')
    ax1.legend(handles=legend_handles)
    st.pyplot(fig1)

    plot_data, components_vectors, component_clusters = clusterize(selected_file, selected_types, cluster_num)
    fig2, ax2 = plt.subplots(figsize=(14, 10))
    for cluster_ind, cluster in sorted(plot_data.items()):        
        ax2.scatter(cluster[0], cluster[1], label=f"Cluster {cluster_ind}")
        
    ax2.set_xlabel("nucleus.x")
    ax2.set_ylabel("nucleus.y")
    ax2.legend()
    st.pyplot(fig2)

    fig3, ax3 = plt.subplots(figsize=(14, 10))
    # print(components_vectors)
    hist_data = np.array(components_vectors)
    x = np.arange(len(hist_data))

    # Plot the bars
    for i in range(hist_data.shape[1]):
        ax3.bar(x, hist_data[:, i], bottom=np.sum(hist_data[:, :i], axis=1), color=list(color_dict.values())[i], label=list(color_dict.keys())[i])


    ax3.set_xlabel("Grafy spójne z otoczką")
    ax3.set_ylabel("Udział procentowy komórek")
    ax3.legend(loc='upper right')

    st.pyplot(fig3)

    fig4, ax4 = plt.subplots(figsize=(14, 10))

    # Create a list of original indices
    original_indices = list(range(len(component_clusters)))

    # Sort the data and the original indices based on component_clusters
    sorted_data = [[vector, original_index] for _,vector, original_index in sorted(zip(component_clusters,components_vectors, original_indices))]

    sorted_vectors, sorted_indices = [i[0] for i in sorted_data], [i[1] for i in sorted_data]
    hist_data = np.array(sorted_vectors)

    x = np.arange(len(hist_data))

    # Plot the bars
    for i in range(hist_data.shape[1]):
        ax4.bar(x, hist_data[:, i], bottom=np.sum(hist_data[:, :i], axis=1), color=list(color_dict.values())[i], label=list(color_dict.keys())[i])

    for list_ind, sorted_index in enumerate(sorted_indices):
        if component_clusters[sorted_index] != component_clusters[sorted_indices[list_ind-1]] and list_ind != 0:
            ax4.axvline(x=list_ind-0.5, color='black')

    # Set the x-ticks to the sorted indices
    ax4.set_xticks(x)
    ax4.set_xticklabels(sorted_indices)
    
    st.pyplot(fig4)



# ========POCZĄTEK PROGRAMU========

current_dir = os.getcwd()

# Define the relative path to the 'if_data_modified' folder
if_data_path = os.path.join(current_dir, 'if_data_modified')

# Get a list of all CSV files in the 'if_data_modified' directory that contain the selected option
selected_option = st.selectbox("Choose an option", ["IF1", "IF2", "IF3"])
files = [f for f in os.listdir(if_data_path) if os.path.isfile(os.path.join(if_data_path, f)) and f.endswith('.csv') and selected_option in f]
files.sort()

# Create a dropdown list for the user to select a file
selected_file = st.selectbox("Choose a CSV file", files)
selected_file = os.path.join(if_data_path, selected_file)

# Define the relative path to the '_phen_to_cell_mapping.csv' file
dict_path = os.path.join(current_dir, f'{selected_option}_phen_to_cell_mapping.csv')


# Define the cell types
cell_types = ["other", "CD15+Tumor", "CD15-Tumor", "Tcell", "Bcell", "BnTcell", "Neutrophil", "Macrophage", "DC"]

# Create a dictionary to store the checkbox status of each cell type
checkbox_status = {cell_type: st.checkbox(cell_type) for cell_type in cell_types}

# Get the list of selected cell types
selected_cell_types = [cell_type for cell_type, status in checkbox_status.items() if status]

# pole do wybrania ilości klastrów
cluster_num = st.number_input("Wprowadź liczbę klastrów", 1, 100)

run_button = st.button("Run")

# Run the program when the button is pressed
if run_button:
    if selected_file is not None:
        visualize(selected_file, dict_path, cluster_num, selected_cell_types)


# zbudować macierz sąsiedztwa przy użyciu funckji  sklearn.neighbors.radius_neighbors_graph
# napisać funkcję pomocniczą która powiązuje komórkę i jej typ z typami komórek jej sąsiadów
# histogram
# dwa grafy - jeden rdzeń TLS, drugi do odpytywania pozostałych typów
