Test_laucher_With_diff_SNPs_EMS.py
===================================

This Python script is a launcher for the detectCOs tool. It's designed to be run from the command line with a YAML configuration file as input.

It imports necessary modules and functions from the detectCOs package, including functions for reading files, performing sliding window analysis, and identifying crossover events.

.. code-block:: python

    import os, sys
    import pprint
    import yaml

    from detectCOs_read_files_EMS import *
    from detectCOs_sliding_window_EMS import *
    from detectCOs_identifyCOs_EMS import *
    from VisualiseGenotypesChromosomeEMS import *
    from collections import OrderedDict
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

    # Set current working directory at the root of polyrec project
    os.chdir(config['path_to_polyrec_project'])

    # Set variable to save all results
    outdir = config['path_to_polyrec_project'] + config['output_dir_in_polyrec'] + \
        config['analyze_id'] + "/"
    # And create output directory if not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print("Create output directory")
..

2. **Loading Chromosome Length**:
    - Chromosome lengths are loaded from a file specified in the configuration.
    - The loaded chromosome lengths are printed.

.. code-block:: python

    print("\n# STEP 1. Load chromosome length")
    print("-------------------------------------")

    chr_length = ReadChrLen(input_chr_len=config['chr_len'], 
                            prefix_chr=config['prefix_chr'])

    print("\nchr_length[chr] = lenght of chromosome")
    pprint.pprint(list(chr_length.items()))
..

3. **Loading Positions of the Centromeric Region**:
    - Positions of the centromeric region for each chromosome are loaded from a file specified in the configuration.
    - The loaded centromeric regions are printed.

.. code-block:: python

    print("\n# STEP 2. Load position of the centromeric region")
    print("-------------------------------------")

    centromere = ReadCentroReg(input_centro_reg=config['centromere_reg'],
                            prefix_chr=config['prefix_chr'])

    print("\ncentromere[chr] = [left_border, right_border]")
    pprint.pprint(list(centromere.items()))

..

4. **Loading Parental SNPs from VCF**:
    - Parental SNP data is loaded from a VCF (Variant Call Format) file specified in the configuration.
    - If the data has been previously processed and saved, it is loaded from files. Otherwise, the `ReadParentalVCF` function is called to read the VCF file.
    - The loaded parental SNPs and last SNP per chromosome are printed.
    - Outputs of `ReadParentalVCF` are saved to files.

.. code-block:: python

    print("\n# STEP 3. Load parental snps from vcf")
    print("-------------------------------------")

    if os.path.exists(outdir + "parental_snps.txt") and \
        os.path.exists(outdir + "parental_snps_last_pos_chr.txt"):
        print("Load data...")
        parental_snps = load_file_in_dict(outdir + "parental_snps.txt")
        parental_last_snp_chr = load_file_in_dict(outdir + "parental_snps_last_pos_chr.txt")

    else:
        print("Run ReadParentalVCF()...")
        parental_snps, parental_last_snp_chr = ReadParentalVCF(
            input_vcf=config['parental_vcf'], 
            prefix_chr=config['prefix_chr']
            )

    print("\nparental_snps[chr_pos] = [chr, pos, GT]")
    pprint.pprint(list(parental_snps.items())[:5])
    print("\nparental_last_snp_chr:")
    pprint.pprint(list(parental_last_snp_chr.items()))

    #------------------------------------------------------------------------------
    print ("\nSave outputs of ReadParentalVCF()")
    export_dict_in_file(my_dict = parental_snps, 
                    output_file=outdir + "parental_snps.txt",
                    header="chr_pos\tchr\tpos\tGT",
                    overwrite=False)

..

5. **Loading Offspring SNPs from VCF**:
    - If the data has been previously processed and saved, it is loaded from files. Otherwise, the `ReadOffspringVCF` function is called to read the VCF file and process the data.
    - The loaded offspring SNPs and last SNP per chromosome are printed.
    - Outputs of `ReadOffspringVCF` are saved to files.

**Saving Last SNPs per Chromosome**:
    - The last SNPs per chromosome, obtained from the parental data, are saved to a file.

.. code-block:: python

    print("\n# STEP 4. Load offspring snps from vcf")
    print("-------------------------------------")

    if os.path.exists(outdir + "offspring_snps.txt") and \
        os.path.exists(outdir + "offspring_snps_last_pos_chr.txt"):
        print("Load data...")
        offspring_snps = load_file_in_dict(outdir + "offspring_snps.txt")
        offspring_last_snp_chr = load_file_in_dict(outdir + "offspring_snps_last_pos_chr.txt")

    else:
        print("Run ReadOffspringVCF()...")
        offspring_snps, offspring_last_snp_chr = ReadOffspringVCF(
            input_vcf=config['sample_vcf'], 
            parental_snps=parental_snps,
            prefix_chr=config['prefix_chr'],
            geno_ref=config['genotype_ref'], 
            geno_alt=config['genotype_alt'],
            analyze_id=config['analyze_id']
            )

    print("\noffspring_snps[chr_pos] = [chr, pos,GT, ADref,ADalt, genotype]")
    pprint.pprint(list(offspring_snps.items())[:5])


    #------------------------------------------------------------------------------
    print ("\nSave outputs of ReadOffspringVCF()")
    export_dict_in_file(my_dict = offspring_snps, 
                    output_file=outdir + "offspring_snps.txt",
                    header="chr\tpos\tGT\tADref\tADalt\tgenotype",
                    overwrite=False)

    last_snps_chr = parental_last_snp_chr


    export_dict_in_file(my_dict = last_snps_chr, 
                    output_file=outdir + "last_snps_chr.txt",
                    header="chr\tpos_last_SNP",
                    overwrite=False)

..

6. **Summary Info About Parental SNPs by Sliding Window**:
    - If the data has been previously processed and saved, it is loaded from a file. Otherwise, the `ParentalSlidingWindow` function is called to summarize parental SNPs by sliding window.
    - The summarized parental SNPs by sliding window are printed.

**Saving Output of Parental Sliding Window**:
    - The output of the parental sliding window analysis is saved to a file.
.. code-block:: python

    print("\n# STEP 5. Summary info about parental SNPs by sliding window")
    print("-------------------------------------")

    print("Run ParentalSlidingWindow()...")
    parental_snps_windowsnps = ParentalSlidingWindowSNPs(
        parental_snps=parental_snps,
        snps_per_window=config['snps_per_window']
        )


    print("\nparental_snps_windowsnps[chr_window] = [start, stop, nb_snps]")
    pprint.pprint(list(parental_snps_windowsnps.items())[5:9])

    #------------------------------------------------------------------------------
    print ("\nSave output of ParentalSlidingWindow()")
    export_dict_in_file(my_dict = parental_snps_windowsnps, 
                        output_file=outdir + "parental_snps_windowsnps_"+ str(config['snps_per_window']) +"_snps.txt",
                        header="chr_window\tstart\tstop\tnb_snps",
                        overwrite=False)

..

7. **Summary Info About Offspring SNPs by Sliding Window**:
    - If the data has been previously processed and saved, it is loaded from a file. Otherwise, the `OffspringSlidingWindow` function is called to summarize offspring SNPs by sliding window.
    - The summarized offspring SNPs by sliding window are printed.

**Saving Output of Offspring Sliding Window**:
    - The output of the offspring sliding window analysis is saved to a file.

.. code-block:: python

    print("\n# STEP 6. Summary info about offspring SNPs by sliding window")
    print("-------------------------------------")

    print("Run ParentalSlidingWindow()...")
    offspring_snps_windowsnps = OffspringSlidingWindowBySNPs(
        offspring_snps = offspring_snps , 
        snps_per_window = config['snps_per_window']) 


    print("\offspring_snps_windowsnps[chr_window] = [start, stop, nb_snps]")
    pprint.pprint(list(offspring_snps_windowsnps.items())[5:9])

    #------------------------------------------------------------------------------
    print ("\nSave output of ParentalSlidingWindow()")
    export_dict_in_file(my_dict = offspring_snps_windowsnps, 
                        output_file=outdir + "offspring_snps_windowsnps_"+ str(config['snps_per_window'])+"_snps.txt",
                        header="chr_window\tstart\tstop\tADref\tADalt\tDP\tnbSNP_Aref\tnbSNP_Aalt\ttot_snps", #\tADref/DP\tADalt/DP",
                        overwrite=False)

..

8. **Identify Offspring Genotype in Normalized Sliding Window**:
    - If the data has been previously processed and saved, it is loaded from files. Otherwise, the `NormalizeOffspringSlidingWindow` function is called to identify the genotype of the offspring in the normalized sliding window.
    - The identified offspring genotype and related information are printed.

**Saving Outputs of Normalize Offspring Sliding Window**:
    - The outputs of the normalization process are saved to files.

.. code-block:: python

    print("\n# STEP 7. Identify the genotype of the offspring normalized sliding window")
    print("-------------------------------------")

    offspring_snps_window_normalized, offspring_genotype_window_normalized, \
        nb_windows_chr = NormalizeOffspringSlidingWindowsnps(
            parental_snps_window=parental_snps_windowsnps,
            offspring_snps_window=offspring_snps_windowsnps,
            geno_ref=config['EMS_ref'],
            geno_alt=config['EMS_alt'],
            min_snp_num=config['min_snp_num'],
            min_reads_num=config['min_reads_num'],
            ratio_min_homo=config['ratio_min_EMS'],
            depth_division_th = config['depth_division_th']  # modify MY
            )

    print("\noffspring_snps_window_normalized: ")
    print("[chr_window] = [start, stop, ADref, ADalt, DP]")
    pprint.pprint(list(offspring_snps_window_normalized.items())[5:9])

    print("\noffspring_genotype_window_normalized: ")
    print("[chr_window] = [start, stop, ADref/DP, ADalt/DP, genotype]")
    pprint.pprint(list(offspring_genotype_window_normalized.items())[5:9], width = 120, compact = True)

    print("\nnb_windows_chr: ")
    pprint.pprint(list(nb_windows_chr.items()))

    #------------------------------------------------------------------------------    
        
    print ("\nSave outputs of NormalizeOffspringSlidingWindow()")
    export_dict_in_file(my_dict = offspring_snps_window_normalized, 
                    output_file=outdir + "offspring_snps_window_normalized_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
                    header="chr_window\tstart\tstop\tADref\tADalt\tDP",
                    overwrite=False)

    export_dict_in_file(my_dict = offspring_genotype_window_normalized, 
                    output_file=outdir + "offspring_genotype_window_normalized_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
                    header="chr_window\tstart\tstop\tADref/DP\tADalt/DP\tgenotype",
                    overwrite=False)

    export_dict_in_file(my_dict = nb_windows_chr, 
                    output_file=outdir + "nb_window_chr_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
                    header="chr\tnb_window",
                    overwrite=False)


..



9. **Smooth Offspring Sliding Window**:
    - If the smoothed data has been previously processed and saved, it is loaded from files. Otherwise, the `SmoothNormalizedOffspringSlidingWindow` function is called to smooth the offspring sliding window data.
    - The smoothed offspring sliding window data, offspring genotype in the smoothed window, and information about the smoothing process are printed.

**Saving Outputs of Smoothed Offspring Sliding Window**:
    - The smoothed offspring sliding window data, offspring genotype in the smoothed window, and information about the smoothing process are saved to separate files.

.. code-block:: python

    print("\n# STEP 8. Smooth offspring sliding window")
    print("-------------------------------------")

    offspring_snps_window_smoothed, offspring_genotype_window_smoothed, \
        check_smoothing = SmoothNormalizedOsffspringSlidingWindowsnps(\
            offspring_snps_window=offspring_snps_window_normalized,
            nb_windows_chr=nb_windows_chr,
            geno_ref=config['EMS_ref'],
            geno_alt=config['EMS_alt'],
            ratio_min_homo=config['ratio_min_EMS'],
            depth_division_th = config['depth_division_th']  # modify MY
            )


    print("\noffspring_snps_window_smoothed: ")
    print("[chr_window] = [start, stop, ADref, ADalt, DP]")
    pprint.pprint(list(offspring_snps_window_smoothed.items())[5:9])

    print("\noffspring_genotype_window_normalized_smoothed")
    print("[chr_window] = [start, stop, ADref/DP, ADalt/DP, probHomoRef, probHetero, probHomoAlt, genotype]")
    pprint.pprint(list(offspring_genotype_window_smoothed.items())[5:9], width = 120, compact = True)

    print("\ncheck_smoothing[chr_window] = [start_window,stop_window,start_smooth,stop_smooth,window_smoothed")
    pprint.pprint(list(check_smoothing.items())[5:9])

    #------------------------------------------------------------------------------

    print("\nSave outputs of SmoothNormalizedOsffspringSlidingWindow()")

    export_dict_in_file(my_dict = offspring_snps_window_smoothed, 
                    output_file=outdir + "offspring_snps_window_normalized_smoothed_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
                    header="chr_window\tstart\tstop\tADref\tADalt\tDP",
                    overwrite=False)

    export_dict_in_file(my_dict = offspring_genotype_window_smoothed, 
                    output_file=outdir + "offspring_genotype_window_normalized_smoothed_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
                    header="chr_window\tstart\tstop\tADref/DP\tADalt/DP\tgenotype",
                    overwrite=False)

    export_dict_in_file(my_dict = check_smoothing, 
                    output_file=outdir + "check_smoothing_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
                    header="chr_window\tstart_window\tstop_window\tstart_smooth\tstop_smooth\twindow_smoothed",
                    overwrite=False)


..

9.bis. **Visualize Genotype Example**:

This section of the code is used to visualize genotype information based on the configuration settings. It maps file paths to genotype data and assigns colors to different genotypes for visual representation.

.. code-block:: python

    print("\n# STEP 8.bis. visualise Genotype ")
    print("-------------------------------------")

    # 'file_paths' est maintenant un dictionnaire avec les noms des répertoires comme clés et les chemins complets comme valeurs.


    # Exemple d'utilisation
    file_paths = {
        # 'directory_name': 'offspring_genotype_window_normalized_smoothed.txt'
    }

    file_paths[config['analyze_id']] = outdir + "offspring_genotype_window_normalized_smoothed_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt"
    #file_paths[config['analyze_id'] + "_" + str(config['new_window_size']) + "_kb"] = outdir + "offspring_genotype_window_normalized_smoothed_" + str(config['new_window_size']) + "_kb" + ".txt"

    genotype_colors = {
        config['EMS_ref']: 'blue',
        config['EMS_alt']: 'red',
        'NA': 'gray'  
    }


    visualsize_genotope(file_paths, outdir, genotype_colors, snps_per_window = config['snps_per_window'], ratio_min_EMS = config['ratio_min_EMS'])



10. **Identify COs**:
    - The `IdentifyCOs` function is called to identify COs based on the smoothed offspring genotype window data.
    - Candidate COs are printed and exported to a file named "candidateCO.txt". 

**Identify COs (Qichao version)**:
    - Another approach for identifying COs, referred to as the "Qichao version," is implemented.
    - This method uses additional information such as smoothed probabilities and sliding genotype ratios.
    - Candidate COs according to this method are printed and exported to a file named "candidates_co_qichao.txt".

.. code-block:: python

    print("\n# STEP 9. Identify COs")
    print("-------------------------------------")

    candidates_co, db_co = IdentifyCOs_snps(offspring_genotype_window_smoothed, nb_windows_chr)


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
                if window_distance < 10:
                    warnings.append(f"Warning: CO at {chr_window} is only {int(window_distance)} windows away from CO at {last_chr_window}. Possible double CO event.")
                    warnings_co[last_chr_window] = updated_candidates_co[last_chr_window]
                    warnings_co[chr_window] = updated_candidates_co[chr_window]
        
        last_co = float(chr_window.split('_')[1])
        last_chr_window = chr_window



    print("\nCandidateCOs: ")
    print("[chr_window] = [start_win:stop_win, co_start, co_stop, pre_geno, cur_geno]")
    pprint.pprint(list(updated_candidates_co.items()))  # modify MY

    export_dict_in_file(my_dict = updated_candidates_co,   # modify MY
                    output_file=outdir + "candidateCO_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
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

..

11. Sliding with new number of snps in window:

.. code-block:: python

    print("\n# STEP 5.bis. Summary info about parental SNPs by sliding window" + str(config['new_snp_per_window']) )
    print("-------------------------------------")


    if os.path.exists(outdir + "parental_snps_window_" + str(config['new_snp_per_window']) + ".txt"):
        print("Load data...")
        parental_snps_window_new = load_file_in_dict(outdir + "parental_snps_window_" + str(config['new_snp_per_window']) + ".txt")

    else:
        print("Run ParentalSlidingWindow()...")
        parental_snps_window_new = ParentalSlidingWindowSNPs(
            parental_snps=parental_snps,
            snps_per_window=config['new_snp_per_window']
        )

    print("\nparental_snps_window_" + str(config['new_snp_per_window']) + "[chr_window] = [start, stop, nb_snps]")
    pprint.pprint(list(parental_snps_window_new.items())[5:10])

    #------------------------------------------------------------------------------
    print ("\nSave output of ParentalSlidingWindow()")
    export_dict_in_file(my_dict = parental_snps_window_new, 
                    output_file=outdir + "parental_snps_window_" + str(config['new_snp_per_window']) + ".txt",
                    header="chr_window\tstart\tstop\tnb_snps",
                    overwrite=False)

    print("\n# STEP 6.bis. Summary info about offspring SNPs by sliding window" + str(config['new_snp_per_window']) )
    print("-------------------------------------")

    if os.path.exists(outdir + "offspring_snps_window_" + str(config['new_snp_per_window']) +  ".txt"):
        print("Load data...")
        offspring_snps_window_new = load_file_in_dict(outdir + "offspring_snps_window_" + str(config['new_snp_per_window']) + ".txt")

    else:
        print("Run OffspringSlidingWindow()...")
        offspring_snps_window_new = OffspringSlidingWindowBySNPs(
            offspring_snps=offspring_snps,
            snps_per_window=config['new_snp_per_window']
        )


    print("\noffspring_snps_window_" + str(config['new_snp_per_window']) +  "[chr_window] = [start,stop,ADref,ADalt,DP,nbSNP_Aref,nbSNP_Aalt,TOTsnps-window]")
    pprint.pprint(list(offspring_snps_window_new.items())[5:10])

    #------------------------------------------------------------------------------
    print ("\nSave output of OffspringSlidingWindow()")
    export_dict_in_file(my_dict = offspring_snps_window_new, 
                    output_file=outdir + "offspring_snps_window_" + str(config['new_snp_per_window']) + ".txt",
                    header="chr_window\tstart\tstop\tADref\tADalt\tDP\tnbSNP_Aref\tnbSNP_Aalt\ttot_snps",
                    overwrite=False)

    print("\n# STEP 7.bis. Identify the genotype of the offspring normalized sliding window" + str(config['new_snp_per_window']) )
    print("-------------------------------------")

    if os.path.exists(outdir + "offspring_snps_window_normalized_" + str(config['new_snp_per_window']) +  ".txt") and \
        os.path.exists(outdir + "offspring_genotype_window_normalized_" + str(config['new_snp_per_window']) + ".txt") and \
            os.path.exists(outdir + "nb_window_chr_" + str(config['new_snp_per_window']) +  ".txt"):
        print("Load data...")
        offspring_snps_window_normalized_new = load_file_in_dict(outdir + "offspring_snps_window_normalized_" + str(config['new_snp_per_window']) +  ".txt")
        offspring_genotype_window_normalized_new = load_file_in_dict(outdir + "offspring_genotype_window_normalized_" + str(config['new_snp_per_window']) + ".txt")
        nb_windows_chr_new = load_file_in_dict(outdir + "nb_window_chr_" + str(config['new_snp_per_window']) +  ".txt")

    else:
        offspring_snps_window_normalized_new, offspring_genotype_window_normalized_new, \
            nb_windows_chr_new = NormalizeOffspringSlidingWindowsnps(
                parental_snps_window=parental_snps_window_new,
                offspring_snps_window=offspring_snps_window_new,
                geno_ref=config['EMS_ref'],
                geno_alt=config['EMS_alt'],
                min_snp_num=config['min_snp_num'],
                min_reads_num=config['min_reads_num'],
                ratio_min_homo=config['ratio_min_EMS'],
                depth_division_th = config['depth_division_th']   # modify MY
                )




    print("\noffspring_snps_window_normalized_" + str(config['new_snp_per_window']))
    print("[chr_window] = [start, stop, ADref, ADalt, DP]")
    pprint.pprint(list(offspring_snps_window_normalized_new.items())[5:10])

    print("\noffspring_genotype_window_normalized_" + str(config['new_snp_per_window']) )
    print("[chr_window] = [start, stop, ADref/DP, ADalt/DP, probHomoRef, probHetero, probHomoAlt, genotype]")
    pprint.pprint(list(offspring_genotype_window_normalized_new.items())[5:10], width = 120, compact = True)

    print("\nnb_windows_chr_" + str(config['new_snp_per_window']))
    pprint.pprint(list(nb_windows_chr_new.items()))

    #------------------------------------------------------------------------------    
        
    print ("\nSave outputs of NormalizeOffspringSlidingWindow() " + str(config['new_snp_per_window']) )
    export_dict_in_file(my_dict = offspring_snps_window_normalized_new, 
                    output_file=outdir + "offspring_snps_window_normalized_" + str(config['new_snp_per_window']) + ".txt",
                    header="chr_window\tstart\tstop\tADref\tADalt\tDP",
                    overwrite=False)

    export_dict_in_file(my_dict = offspring_genotype_window_normalized_new, 
                    output_file=outdir + "offspring_genotype_window_normalized_" + str(config['new_snp_per_window'])  + ".txt",
                    header="chr_window\tstart\tstop\tADref/DP\tADalt/DP\tprobHomoRef\tprobHetero\tprobHomoAlt\tgenotype",
                    overwrite=False)

    export_dict_in_file(my_dict = nb_windows_chr_new, 
                    output_file=outdir + "nb_window_chr_" + str(config['new_snp_per_window'])+ ".txt",
                    header="chr\tnb_window",
                    overwrite=False)

    print("\n# STEP 8. Smooth offspring sliding window " + str(config['new_snp_per_window']) )
    print("-------------------------------------")

    if os.path.exists(outdir + "offspring_snps_window_normalized_smoothed_" + str(config['new_snp_per_window'])  + ".txt") and \
        os.path.exists(outdir + "offspring_genotype_window_normalized_smoothed_" + str(config['new_snp_per_window'])  + ".txt") and \
            os.path.exists(outdir + "offspring_nb_window_normalized_" + str(config['new_snp_per_window']) + ".txt"):
        print("Load data...")
        offspring_snps_window_smoothed_new = load_file_in_dict(outdir + "offspring_snps_window_normalized_smoothed_" + str(config['new_snp_per_window'])  + ".txt")
        offspring_genotype_window_smoothed_new = load_file_in_dict(outdir + "offspring_genotype_window_normalized_smoothed_" + str(config['new_snp_per_window']) + ".txt")
        check_smoothing_new = load_file_in_dict(outdir + "check_smoothing_" + str(config['new_snp_per_window'])  + ".txt")

    else:
        offspring_snps_window_smoothed_new, offspring_genotype_window_smoothed_new, \
            check_smoothing_new = SmoothNormalizedOsffspringSlidingWindowsnps(\
                offspring_snps_window=offspring_snps_window_normalized_new,
                nb_windows_chr=nb_windows_chr_new,
                geno_ref=config['EMS_ref'],
                geno_alt=config['EMS_alt'],
                ratio_min_homo=config['ratio_min_EMS'],
                depth_division_th = config['depth_division_th']  # modify MY
                )



    print("\noffspring_snps_window_smoothed_" + str(config['new_snp_per_window'])  + " : ")
    print("[chr_window] = [start, stop, ADref, ADalt, DP]")
    pprint.pprint(list(offspring_snps_window_smoothed_new.items())[5:10])

    print("\noffspring_genotype_window_normalized_smoothed_" + str(config['new_snp_per_window']) )
    print("[chr_window] = [start, stop, ADref/DP, ADalt/DP, probHomoRef, probHetero, probHomoAlt, genotype]")
    pprint.pprint(list(offspring_genotype_window_smoothed_new.items())[5:10], width = 120, compact = True)

    print("\ncheck_smoothing_" + str(config['new_snp_per_window'])  + "[chr_window] = [start_window,stop_window,start_smooth,stop_smooth,window_smoothed")
    pprint.pprint(list(check_smoothing_new.items())[5:10])

    #------------------------------------------------------------------------------

    print("\nSave outputs of SmoothNormalizedOsffspringSlidingWindow()")

    export_dict_in_file(my_dict = offspring_snps_window_smoothed_new, 
                    output_file=outdir + "offspring_snps_window_normalized_smoothed_" + str(config['new_snp_per_window'])  + ".txt",
                    header="chr_window\tstart\tstop\tADref\tADalt\tDP",
                    overwrite=False)

    export_dict_in_file(my_dict = offspring_genotype_window_smoothed_new, 
                    output_file=outdir + "offspring_genotype_window_normalized_smoothed_" + str(config['new_snp_per_window'])  + ".txt",
                    header="chr_window\tstart\tstop\tADref/DP\tADalt/DP\tprobHomoRef\tprobHetero\tprobHomoAlt\tgenotype",
                    overwrite=False)

    export_dict_in_file(my_dict = check_smoothing_new, 
                    output_file=outdir + "check_smoothing_" + str(config['new_snp_per_window']) + ".txt",
                    header="chr_window\tstart_window\tstop_window\tstart_smooth\tstop_smooth\twindow_smoothed",
                    overwrite=False)


11.bis **Pricise COs**:


    - The `IdentifyCOs` function is called to identify crossover events (COs) based on the new smoothed offspring genotype window data, which is adjusted according to the newly specified window size.
    - Precise COs are printed and exported to a file named "preciseCOs.txt" for further analysis and record-keeping.

.. code-block:: python

    print("\n# STEP 9.bis. Precise COs")
    print("-------------------------------------")


    preciseCOs = OrderedDict()
    db_co_new = OrderedDict()

    for i in range (len(list(updated_candidates_co.keys()))):
        start_co=updated_candidates_co[list(updated_candidates_co.keys())[i]][1] 
        end_co=updated_candidates_co[list(updated_candidates_co.keys())[i]][2] 
        chr = list(updated_candidates_co.keys())[i].split("_")[0]
        before =OrderedDict()
        after = OrderedDict()
        
        for key, snp_info in offspring_genotype_window_smoothed.items():
            candidate = key
            chromosome = candidate.split("_")[0]
            index = candidate.split("_")[1]
            start,	stop, ratio_ADref, ratio_ADalt, genotype = snp_info

            
            if start < start_co and chromosome == chr and genotype != 'NA' :
                before[candidate] = snp_info
            if stop > end_co and chromosome == chr and genotype != 'NA' :
                after[candidate] = snp_info
        before = dict(list(before.items())[-3:])
        after = dict(list(after.items())[:3])


        start_co_2=before[list(before.keys())[0]][0]
        if len(list(after.keys()))>2:
            end_co_2= after[list(after.keys())[2]][1]
        elif len(list(after.keys())) == 0:
            end_co_2 = end_co
        elif len(list(after.keys())) == 2:
            end_co_2 = after[list(after.keys())[1]][1]
        else:
            end_co_2= after[list(after.keys())[0]][1]

        new_offspring_genotype_window_smoothed = OrderedDict()
        for cle, valeur in offspring_genotype_window_smoothed_new.items():
            if valeur[0] >= start_co_2 + 1 and valeur[1] <= end_co_2 and cle.split("_")[0] == chr:
                new_offspring_genotype_window_smoothed[cle] = valeur
                print(cle, valeur)
        
        start_window = int(list(new_offspring_genotype_window_smoothed.keys())[0].split("_")[1])
        end_window = int(list(new_offspring_genotype_window_smoothed.keys())[-1].split("_")[1])
        print(start_window, end_window)


        candidates_co_2, db_co_2 = PreciseCOs_snps(new_offspring_genotype_window_smoothed, start_window=start_window, end_window=end_window)
        
        if candidates_co_2 :
            candidates_co_2 = candidates_co_2
        else:
            start_co=updated_candidates_co[list(updated_candidates_co.keys())[i]][1] 
            end_co=updated_candidates_co[list(updated_candidates_co.keys())[i]][2]
            pre_geno=updated_candidates_co[list(updated_candidates_co.keys())[i]][4]
            cur_geno=updated_candidates_co[list(updated_candidates_co.keys())[i]][5]
            cur_chr = list(updated_candidates_co.keys())[i].split("_")[0]
            key = cur_chr + "_" + str(round((start_window + end_window) / 2, 1))
            candidates_co_2[key] = [str(start_window) + ":" + str(end_window), start_co, end_co, pre_geno, cur_geno]
        preciseCOs.update(candidates_co_2)
        db_co_new.update(db_co_2)


    updated_preciseCOs = OrderedDict() # modify MY
    for chr_window, details in preciseCOs.items(): # modify MY
        start_win_stop_win, co_start, co_stop, pre_geno, cur_geno = details # modify MY
        diff = co_stop - co_start  # Calcul de la différence # modify MY
        updated_preciseCOs[chr_window] =  [start_win_stop_win, co_start, co_stop,diff, pre_geno, cur_geno] # modify MY


    print("\nCandidateCOs new: ")
    print("[chr_window] = [start_win:stop_win, co_start, co_stop, pre_geno, cur_geno]")
    pprint.pprint(list(updated_preciseCOs.items())) # modify MY


    export_dict_in_file(my_dict = updated_preciseCOs,  # modify MY
                    output_file=outdir + "preciseCOs_" + str(config['new_snp_per_window'])+"_snps_ratio_" + str(config['ratio_min_EMS']) + ".txt",
                    header="chr_mean_win\tstart_win:stop_win\tstart_co\tstop_co\t(stop_co-start_co)\tpre_geno\tcur_geno",
                    overwrite=True)

..

