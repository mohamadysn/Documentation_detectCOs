Launcher_Rules_RIL.py
=======================

This Python script is a launcher for the detectCOs tool. It's designed to be run from the command line with a YAML configuration file as input.

It imports necessary modules and functions from the detectCOs package, including functions for reading files, performing sliding window analysis, and identifying crossover events.

.. code-block:: python

    import pandas as pd
    import os, sys
    import pprint
    import yaml
    from detectCOs_read_files import *
    from collections import OrderedDict
    from VisualiseGenotypesChromosome import *
    from detectCOs_identifyCOs_RIL import *

..

It provides usage instructions for running the script from the command line, specifying the path to the configuration file as an argument.

1. **Loading Configuration File and Creating Output Directory :**  
    - The script begins by loading a configuration file in YAML format specified as a command-line argument. If no argument is provided, it raises an `ImportError`.
    - The current working directory is set to the root of the `polyrec` project based on the configuration.
    - A variable `outdir` is set to save all results, and an output directory is created if it doesn't exist.


.. code-block:: python

    print("# Load config file and create output directory")
    print("-------------------------------------")

    if len(sys.argv) == 1:
        raise ImportError("Please give the config file (.yaml) for detectCOs")

    config_file = sys.argv[1]
    with open(config_file, "r") as cf:
        config = yaml.load(cf, Loader=yaml.FullLoader)


    # Set variable to save all results
    outdir = config['output'] 
    # And create output directory if not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print("Create output directory")

    # Load the Excel file
    data = pd.read_excel(config['Path_info_RIL'] + "7RV498_window_call.xlsx")

    # Set the column names to the values of the first row
    data.columns = data.iloc[0]

    # Remove the first row from the dataframe
    df = data.iloc[1:].reset_index(drop=True)

    print(df.head())

    df_indexed = df.set_index('chr_window')

    # Now create the dictionary
    RILs = df_indexed.apply(lambda row: row.tolist(), axis=1).to_dict()


    print("\n# STEP 8. Smooth offspring sliding window")
    print("-------------------------------------")

    if os.path.exists(config['Path_offspring_EMS'] + "offspring_genotype_window_normalized_smoothed_6_snps_ratio_0.2.txt") and \
        os.path.exists(config['Path_offspring_basic'] + "offspring_genotype_window_normalized_smoothed.txt"):
        print("Load data...")
        offspring_genotype_window_smoothed = load_file_in_dict(config['Path_offspring_basic'] + "offspring_genotype_window_normalized_smoothed.txt")
        offspring_genotype_window_smoothed_EMS = load_file_in_dict( config['Path_offspring_EMS'] + "offspring_genotype_window_normalized_smoothed_6_snps_ratio_0.2.txt")


    test, nb_window_chr = Rules_RIls(info_RIL = RILs,offspring_EMS = offspring_genotype_window_smoothed_EMS,offspring_basic = offspring_genotype_window_smoothed)


    export_dict_in_file(my_dict = nb_window_chr, 
                    output_file=outdir + "nb_window_chr.txt",
                    header="chr\tnb_window",
                    overwrite=False)


    export_dict_in_file(my_dict = test, 
                    output_file= outdir + "test_rules.txt",
                    header="chr_window\tstart\tstop\tADref/DP\tADalt/DP\tgenotype_base\tgenotype_EMS\tgenotype_RIL\tResponse\tWARNING\tGT_Final",
                    overwrite=False)


    file_paths = {
        # 'directory_name': 'offspring_genotype_window_normalized_smoothed.txt'
    }

    file_paths['detectCOs_scripts_modify'] = outdir + 'test_rules.txt'


    genotype_colors = {
        'A': 'blue',
        'RIL': 'red',
        'NA': 'gray'  
    }

    visualsize_genotope(file_paths, outdir , genotype_colors, column = 'GT_Final')



    print("\n# STEP 9. Identify COs")
    print("-------------------------------------")


    # Configurez logging pour écrire dans un fichier et afficher sur le terminal
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logger = logging.getLogger()
    logger.addHandler(logging.FileHandler(outdir + 'IdentifyCOs_output_modif.txt', 'a'))
    console = logging.StreamHandler()
    logger.addHandler(console)


    candidates_co, db_co = IdentifyCOs_RILs(test, nb_window_chr)


    updated_candidates_co = OrderedDict()  # modify MY
    warnings = []
    last_co = None
    last_chr_window = None
    warnings_co = OrderedDict()

    for chr_window, details in candidates_co.items(): # modify MY
        start_win_stop_win, co_start, co_stop, pre_geno, cur_geno = details # modify MY
        diff = co_stop - co_start  # Calcul de la différence # modify MY
        updated_candidates_co[chr_window] =  [start_win_stop_win, co_start, co_stop,diff, pre_geno, cur_geno] # modify MY
        
        # Extract the window number of the previous crossover's end
        if last_co is not None:
            # Extract the current chromosome number and compare it to the last processed CO's chromosome
            current_chr = chr_window.split('_')[0]
            last_chr = last_chr_window.split('_')[0]
            co_start_win = float(chr_window.split('_')[1])
            
            if current_chr == last_chr:  # Check if the CO is on the same chromosome
                window_distance = co_start_win - last_co  # Calculate the distance between the COs
                if window_distance < 40:
                    warnings.append(f"Warning: CO at {chr_window} is only {int(window_distance)} windows away from CO at {last_chr_window}. Possible double CO event.")
                    warnings_co[last_chr_window] = updated_candidates_co[last_chr_window]
                    warnings_co[chr_window] = updated_candidates_co[chr_window]
        
        last_co = float(chr_window.split('_')[1])
        last_chr_window = chr_window

    print("\nCandidateCOs: ")
    print("[chr_window] = [start_win:stop_win, co_start, co_stop, pre_geno, cur_geno]")
    pprint.pprint(list(updated_candidates_co.items()))  # modify MY


    export_dict_in_file(my_dict = updated_candidates_co,   # modify MY
                    output_file= outdir + "candidateCO.txt",
                    header="chr_mean_win\tstart_win:stop_win\tstart_co\tstop_co\t(stop_co-start_co)\tpre_geno\tcur_geno",
                    overwrite=True)



    if len(warnings_co) != 0:
        print("\nWARNING dbCOs")
        pprint.pprint(list(warnings_co.items()))
        print("\n")
        for warning in warnings:
            print(warning)

        export_dict_in_file(my_dict = warnings_co,   # modify MY
                    output_file= outdir + "Warnings_dbCO.txt",
                    header="chr_mean_win\tstart_win:stop_win\tstart_co\tstop_co\t(stop_co-start_co)\tpre_geno\tcur_geno",
                    overwrite=True)

    print("\n")
    print("-------------------------------------")
    print("\n")