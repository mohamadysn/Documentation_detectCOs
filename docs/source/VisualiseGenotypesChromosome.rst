VisualiseGenotypesChromosome.py
===============================


This Python script includes several utility functions for handling and visualizing genomic data. The script supports operations like delimiter detection in CSV files, color assignment based on genotype, searching for specific files in directories, and creating visual representations of genomic data.


detect_delimiter
----------------

.. autofunction:: module_name.detect_delimiter

.. code-block:: python

    def detect_delimiter(file_path, num_lines=10):
        """
        Detects the delimiter used in a given CSV file by reading the first 'num_lines' lines.
        The function checks for common delimiters such as commas (,), semicolons (;), tabs (\t), 
        and spaces (' ').
        
        :param file_path: Path to the CSV file.
        :param num_lines: Number of lines to read to detect the delimiter.
        :return: The detected delimiter.
        """
        delimiters = [',', ';', '\t', '|', ' ']
        with open(file_path, 'r', newline='') as file:
            # Read the first 'num_lines' lines
            lines = [file.readline() for _ in range(num_lines)]
            
            # Count occurrences of potential delimiters
            delimiter_counts = {delim: 0 for delim in delimiters}
            for line in lines:
                counter = Counter(line)
                for delim in delimiters:
                    delimiter_counts[delim] += counter[delim]
            
            # Detect the delimiter with the highest number of occurrences
            detected_delimiter = max(delimiter_counts, key=delimiter_counts.get)
            
        return detected_delimiter

This function detects the delimiter used in a CSV file by examining the first few lines. It supports common delimiters like commas, semicolons, tabs, pipes, and spaces.

assign_window_colors
--------------------

.. autofunction:: module_name.assign_window_colors

.. code-block:: python

    def assign_window_colors(data, genotype_colors):
        """
        Assigns colors to windows based on genotype, using a dictionary that maps genotypes to colors.
        Handles NaN values by converting them to 'NA'.

        :param data: DataFrame containing the data, including a 'genotype' column for which colors will be assigned.
        :param genotype_colors: Dictionary where keys are genotypes and values are the corresponding colors.
        """
        data = data.copy()
        # Convert all 'genotype' values to strings, explicitly handling NaN values
        data['genotype'] = data['genotype'].apply(lambda x: 'NA' if pd.isnull(x) else str(x).strip())
        # Assign colors based on the 'genotype_colors' dictionary, with 'gray' as the default color
        data['color'] = data['genotype'].apply(lambda x: genotype_colors.get(x, 'gray'))
        return data

Assigns colors to data windows based on their genotype values. It ensures that NaN values are handled appropriately by assigning them a default color.

find_files_in_directories
-------------------------

.. autofunction:: module_name.find_files_in_directories

.. code-block:: python

    def find_files_in_directories(directory, filename, max_dirs=None):
        """
        Searches through subdirectories of 'directory' for 'filename'.
        Returns a dictionary with the names of directories where the file was found and the full path of the file.
        The 'max_dirs' parameter allows to limit the number of directories to be considered.
        
        :param directory: Parent directory to search for subdirectories.
        :param filename: Name of the file to search for in the subdirectories.
        :param max_dirs: Maximum number of directories to consider (None for unlimited).
        :return: Dictionary of file paths with directory names as keys.
        """
        file_paths = {}
        dir_count = 0
        # Traverse directories and subdirectories
        for subdir, dirs, files in os.walk(directory):
            # Skip the root directory
            if subdir == directory:
                continue
            # Check if the desired file is in the current directory
            if filename in files:
                # Add the full path of the file to the list
                file_paths[os.path.basename(subdir)] = os.path.join(subdir, filename)
                dir_count += 1
                # Stop searching if the maximum number of directories is reached
                if max_dirs is not None and dir_count >= max_dirs:
                    break
        return file_paths
Searches through specified directories for a particular file name, returning a dictionary of the paths where the file is found. It allows for a limit on the number of directories to search.

visualsize_genotope
-------------------

.. autofunction:: module_name.visualsize_genotope

.. code-block:: python

    def visualsize_genotope(file_paths, outdir, genotype_colors):
        """
        Generates genomic size visualizations for each chromosome 
        from specified files and saves them in an output directory.
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
            plt.xlabel('Window Index', fontsize=fontsize)
            n = len(chromosome_data_with_colors)
            step = max(1, n // 10)
            plt.xticks(range(0, n, step), chromosome_data_with_colors['chr_window'].iloc[::step], rotation=90, fontsize=fontsize)
            plt.ylim(-(len(data_list) / 2.0), (len(data_list) / 2.0))
            plt.yticks(y_positions, [name for name, _ in data_list], fontsize=fontsize)
            plt.tight_layout()
            plt.subplots_adjust(right=0.75)
            
            # Dynamically create legend elements based on 'genotype_colors'
            legend_elements = [plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=color, markersize=10, label=genotype) for genotype, color in genotype_colors.items()]
            plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=20)
            

            plt.title(f'Chromosome {chromosome}', fontsize=fontsize)
            plt.savefig(f'{outdir}/{chromosome}.png', dpi=100, bbox_inches='tight')
            #plt.show()

Generates visual representations of genomic data organized by chromosome and stores these visualizations in an output directory. It uses matplotlib for plotting and saves the plots as images.





