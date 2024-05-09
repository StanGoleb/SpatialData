import pandas as pd
import json
import matplotlib.patches as mpatches

def color_types(input_path, dict_path):
    df = pd.read_csv(input_path)
    dict_file = pd.read_csv(dict_path)

    # Stworzenie odpowiednio przesortowanego słownika do znajdowania typu komórki po jej markerach (tylko dla IF1).
    sorted_dict = {}
    for key_ind, key in enumerate(dict_file['phenotype']):
        old = key.split('C')[1:]
        new_key = 'C'+'C'.join([old[1], old[5], old[4], old[0], old[3], old[2]])
        sorted_dict[new_key] = dict_file.at[key_ind, 'celltype']
    with open("IF1_dict.json", "w") as jfile:
        json.dump(sorted_dict, jfile)    
    xs = [x for x in df['nucleus.x']]
    ys = [y for y in df['nucleus.y']]
    labels = [sorted_dict[phenotype] for phenotype in df['phenotype']]

    def rgb_to_hex(rgb):
        return '#' + '%02x%02x%02x' % rgb

    IF1_cell_mapping = {"other": rgb_to_hex((190, 190, 190)), 
                        "CD15+Tumor": rgb_to_hex((73, 176, 248)),
                        "CD15-Tumor": rgb_to_hex((138, 79, 45)),
                        "Tcell": rgb_to_hex((235, 74, 148)),
                        "Bcell": rgb_to_hex((204, 49, 31)),
                        "BnTcell": rgb_to_hex((236, 95, 42)),
                        "Neutrophil": rgb_to_hex((0, 40, 245)),
                        "Macrophage": rgb_to_hex((97, 209, 62)),
                        "DC": rgb_to_hex((49, 113, 30))}

    # Zczytanie odpowiedinch kolorów i nazw typów dla komórek
    colors = [IF1_cell_mapping[celltype] for celltype in labels]
    legend_handles = [mpatches.Patch(color=color, label=label) for label, color in IF1_cell_mapping.items()]
    
    return xs, ys, colors, legend_handles, IF1_cell_mapping
