import os
import pandas as pd
import matplotlib.pyplot as plt
import csv
from collections import Counter


#==============================================================================
# Program		:	additional functions for detectCOs
# Author		:	Mohamad Yassine [mohamad.yassine@inrae.fr]
# Date			:	25.02.2024
# Last update	:	26.04.2024
# Version		:	2.0
#==============================================================================



def detect_delimiter(file_path, num_lines=10):
    """
    Détecte le séparateur utilisé dans un fichier CSV donné en lisant les premières 'num_lines' lignes.
    La fonction vérifie les séparateurs courants tels que les virgules (,), les points-virgules (;), les tabulations (\t) 
    et les espaces (' ').
    
    :param file_path: Chemin d'accès au fichier CSV.
    :param num_lines: Nombre de lignes à lire pour détecter le séparateur.
    :return: Le séparateur détecté.
    """
    delimiters = [',', ';', '\t', '|', ' ']
    with open(file_path, 'r', newline='') as file:
        # Lire les premières 'num_lines' lignes
        lines = [file.readline() for _ in range(num_lines)]
        
        # Compter les occurrences des séparateurs potentiels
        delimiter_counts = {delim: 0 for delim in delimiters}
        for line in lines:
            counter = Counter(line)
            for delim in delimiters:
                delimiter_counts[delim] += counter[delim]
        
        # Détecter le séparateur ayant le plus grand nombre d'occurrences
        detected_delimiter = max(delimiter_counts, key=delimiter_counts.get)
        
    return detected_delimiter


"""
def assign_window_colors(data):
    data = data.copy()
    # Création d'une colonne 'color' pour stocker les couleurs
    data['color'] = ''

    # Parcours des chromosomes uniques dans les données
    chromosomes = data['chr_window'].str.split('_', expand=True)[0].unique()
    for chromosome in chromosomes:
        # Sélection des données pour le chromosome actuel
        chromosome_data = data[data['chr_window'].str.startswith(chromosome)].copy()

        # Assignation des couleurs en fonction du génotype
        for _, row in chromosome_data.iterrows():
            window_index = row.name
            genotype = row['genotype']

            # Assignation de la couleur en fonction du génotype
            if genotype == 'Col':
                color = 'blue'
            elif genotype == 'Ct':
                color = 'red'
            elif genotype == 'Col/Ct': 
                color = 'green'
            else:
                color = 'gray'

            # Mise à jour de la couleur de la fenêtre
            data.loc[window_index, 'color'] = color

    return data
"""

def assign_window_colors(data, genotype_colors):
    """
    Attribue des couleurs aux fenêtres en fonction du génotype, en utilisant un dictionnaire de correspondance entre génotypes et couleurs.
    Gère les valeurs NaN en les convertissant en 'NA'.

    :param data: DataFrame contenant les données, y compris une colonne 'genotype' pour laquelle les couleurs seront attribuées.
    :param genotype_colors: Dictionnaire où les clés sont des génotypes et les valeurs sont les couleurs correspondantes.
    """
    data = data.copy()
    # Convertir toutes les valeurs de 'genotype' en chaînes de caractères, et gérer NaN explicitement
    data['genotype'] = data['genotype'].apply(lambda x: 'NA' if pd.isnull(x) else str(x).strip())
    # Attribuer des couleurs en fonction du dictionnaire 'genotype_colors', avec 'gray' comme couleur par défaut
    data['color'] = data['genotype'].apply(lambda x: genotype_colors.get(x, 'gray'))
    return data



def find_files_in_directories(directory, filename, max_dirs=None):
    """
    Parcourt les sous-répertoires de 'directory' à la recherche de 'filename'.
    Retourne un dictionnaire avec les noms des répertoires où le fichier a été trouvé et le chemin complet du fichier.
    Le paramètre 'max_dirs' permet de limiter le nombre de répertoires à prendre en compte.
    
    :param directory: Répertoire parent où chercher les sous-répertoires.
    :param filename: Nom du fichier à chercher dans les sous-répertoires.
    :param max_dirs: Nombre maximal de répertoires à considérer (None pour illimité).
    :return: Dictionnaire des chemins de fichiers avec les noms des répertoires comme clés.
    """
    file_paths = {}
    dir_count = 0
    # Parcourir les répertoires et sous-répertoires
    for subdir, dirs, files in os.walk(directory):
        # Ignorer le répertoire racine
        if subdir == directory:
            continue
        # Vérifier si le fichier désiré est dans le répertoire actuel
        if filename in files:
            # Ajouter le chemin complet du fichier à la liste
            file_paths[os.path.basename(subdir)] = os.path.join(subdir, filename)
            dir_count += 1
            # Si le nombre maximal de répertoires est atteint, arrêter la recherche
            if max_dirs is not None and dir_count >= max_dirs:
                break
    return file_paths

def visualsize_genotope(file_paths, outdir, genotype_colors, snps_per_window, ratio_min_EMS):
    """
    Génère des visualisations de la taille génomique pour chaque chromosome 
    à partir de fichiers spécifiés et les enregistre dans un répertoire de sortie.
    """
    chromosome_data_all_files = {}

    for directory_name, file_path in file_paths.items():
        delimiter = detect_delimiter(file_path)
        print(f"The detected delimiter in {directory_name} is: '{delimiter}'")
        data = pd.read_csv(file_path, sep=delimiter)
        print(f"Data from {directory_name}:\n{data.head()}")

        chromosomes = data['chr_window'].str.split('_', expand=True)[0].unique()
        for chromosome in chromosomes:
            chromosome_data = data[data['chr_window'].str.startswith(chromosome)]
            chromosome_data_with_colors = assign_window_colors(chromosome_data, genotype_colors)
            chromosome_data_with_colors = chromosome_data_with_colors.reset_index(drop=True)

            if chromosome not in chromosome_data_all_files:
                chromosome_data_all_files[chromosome] = []
            chromosome_data_all_files[chromosome].append((directory_name, chromosome_data_with_colors))

    for chromosome, data_list in chromosome_data_all_files.items():
        if len(data_list) <= 50:
            plt.figure(figsize=(20, 10))
            fontsize = 12
        elif len(data_list) <= 80:
            plt.figure(figsize=(20, 15))
            fontsize = 14
        else:
            plt.figure(figsize=(60, 30))
            fontsize = 18

        y_positions = [y - (len(data_list) - 1) / 2.0 for y in range(len(data_list))]
        for i, (directory_name, chromosome_data_with_colors) in enumerate(data_list):
            y_pos = y_positions[i]
            plt.scatter(
                range(len(chromosome_data_with_colors)),
                [y_pos] * len(chromosome_data_with_colors),
                c=chromosome_data_with_colors['color'],
                marker='s',
                s=500,
                label=f'{directory_name}'
            )
        plt.xlabel('Index de la fenêtre', fontsize=fontsize)
        n = len(chromosome_data_with_colors)
        step = max(1, n // 10)
        plt.xticks(range(0, n, step), chromosome_data_with_colors['chr_window'].iloc[::step], rotation=90, fontsize=fontsize)
        plt.ylim(-(len(data_list) / 2.0), (len(data_list) / 2.0))
        plt.yticks(y_positions, [name for name, _ in data_list], fontsize=fontsize)
        plt.tight_layout()
        plt.subplots_adjust(right=0.75)
        
        # Création dynamique des éléments de la légende basée sur 'genotype_colors'
        legend_elements = [plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=color, markersize=10, label=genotype) for genotype, color in genotype_colors.items()]
        plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=20)
        

        plt.title(f'Chromosome {chromosome} | Windows of {snps_per_window} snps with ratio {ratio_min_EMS}', fontsize=fontsize)
        plt.savefig(f'{outdir}/{chromosome}_{snps_per_window}_snps_ratio_{ratio_min_EMS}.png', dpi=100, bbox_inches='tight')
        #plt.show()

