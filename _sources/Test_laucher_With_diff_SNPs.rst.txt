Test_laucher_With_diff_SNPs.py
===============================

This Python script is a launcher for the detectCOs tool. It's designed to be run from the command line with a YAML configuration file as input.

It imports necessary modules and functions from the detectCOs package, including functions for reading files, performing sliding window analysis, and identifying crossover events.

.. code-block:: python

    import os, sys
    import pprint
    import yaml

    from detectCOs_read_files import *
    from detectCOs_sliding_window import *
    from detectCOs_identifyCOs import *
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

    if os.path.exists(outdir + "parental_snps_window.txt"):
        print("Load data...")
        parental_snps_window = load_file_in_dict(outdir + "parental_snps_window.txt")

    else:
        print("Run ParentalSlidingWindow()...")
        parental_snps_window = ParentalSlidingWindow(
            chr_len=chr_length,
            parental_snps=parental_snps,
            last_snps_chr=last_snps_chr,
            window_size=config['window_size']
        )

    print("\nparental_snps_window[chr_window] = [start, stop, nb_snps]")
    pprint.pprint(list(parental_snps_window.items())[2012:2017])

    #------------------------------------------------------------------------------
    print ("\nSave output of ParentalSlidingWindow()")
    export_dict_in_file(my_dict = parental_snps_window, 
                    output_file=outdir + "parental_snps_window.txt",
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

    if os.path.exists(outdir + "offspring_snps_window.txt"):
        print("Load data...")
        offspring_snps_window = load_file_in_dict(outdir + "offspring_snps_window.txt")

    else:
        print("Run OffspringSlidingWindow()...")
        offspring_snps_window = OffspringSlidingWindow(
            chr_len=chr_length,
            offspring_snps=offspring_snps,
            last_snps_chr=last_snps_chr,
            window_size=config['window_size']
        )

    print("\noffspring_snps_window[chr_window] = [start,stop,ADref,ADalt,DP,nbSNP_Aref,nbSNP_Aalt,TOTsnps-window]")
    pprint.pprint(list(offspring_snps_window.items())[2012:2017])

    #------------------------------------------------------------------------------
    print ("\nSave output of OffspringSlidingWindow()")
    export_dict_in_file(my_dict = offspring_snps_window, 
                    output_file=outdir + "offspring_snps_window.txt",
                    header="chr_window\tstart\tstop\tADref\tADalt\tDP\tnbSNP_Aref\tnbSNP_Aalt\ttot_snps",
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

    if os.path.exists(outdir + "offspring_snps_window_normalized.txt") and \
        os.path.exists(outdir + "offspring_genotype_window_normalized.txt") and \
            os.path.exists(outdir + "offspring_nb_window_normalized.txt"):
        print("Load data...")
        offspring_snps_window_normalized = load_file_in_dict(outdir + "offspring_snps_window_normalized.txt")
        offspring_genotype_window_normalized = load_file_in_dict(outdir + "offspring_genotype_window_normalized.txt")
        nb_windows_chr = load_file_in_dict(outdir + "nb_window_chr.txt")

    else:
        offspring_snps_window_normalized, offspring_genotype_window_normalized, \
            nb_windows_chr = NormalizeOffspringSlidingWindow(
                parental_snps_window=parental_snps_window,
                offspring_snps_window=offspring_snps_window,
                geno_ref=config['genotype_ref'],
                geno_alt=config['genotype_alt'],
                min_snp_num=config['min_snp_num'],
                min_reads_num=config['min_reads_num'],
                ratio_min_homo=config['ratio_min_homo']
                )

    print("\noffspring_snps_window_normalized: ")
    print("[chr_window] = [start, stop, ADref, ADalt, DP]")
    pprint.pprint(list(offspring_snps_window_normalized.items())[2012:2017])

    print("\noffspring_genotype_window_normalized: ")
    print("[chr_window] = [start, stop, ADref/DP, ADalt/DP, probHomoRef, probHetero, probHomoAlt, genotype]")
    pprint.pprint(list(offspring_genotype_window_normalized.items())[2012:2017], width = 120, compact = True)

    print("\nnb_windows_chr: ")
    pprint.pprint(list(nb_windows_chr.items()))

    #------------------------------------------------------------------------------    
        
    print ("\nSave outputs of NormalizeOffspringSlidingWindow()")
    export_dict_in_file(my_dict = offspring_snps_window_normalized, 
                    output_file=outdir + "offspring_snps_window_normalized.txt",
                    header="chr_window\tstart\tstop\tADref\tADalt\tDP",
                    overwrite=False)

    export_dict_in_file(my_dict = offspring_genotype_window_normalized, 
                    output_file=outdir + "offspring_genotype_window_normalized.txt",
                    header="chr_window\tstart\tstop\tADref/DP\tADalt/DP\tprobHomoRef\tprobHetero\tprobHomoAlt\tgenotype",
                    overwrite=False)

    export_dict_in_file(my_dict = nb_windows_chr, 
                    output_file=outdir + "nb_window_chr.txt",
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

    if os.path.exists(outdir + "offspring_snps_window_normalized_smoothed.txt") and \
        os.path.exists(outdir + "offspring_genotype_window_normalized_smoothed.txt") and \
            os.path.exists(outdir + "offspring_nb_window_normalized.txt"):
        print("Load data...")
        offspring_snps_window_smoothed = load_file_in_dict(outdir + "offspring_snps_window_normalized_smoothed.txt")
        offspring_genotype_window_smoothed = load_file_in_dict(outdir + "offspring_genotype_window_normalized_smoothed.txt")
        check_smoothing = load_file_in_dict(outdir + "check_smoothing.txt")

    else:
        offspring_snps_window_smoothed, offspring_genotype_window_smoothed, \
            check_smoothing = SmoothNormalizedOsffspringSlidingWindow(\
                offspring_snps_window=offspring_snps_window_normalized,
                nb_windows_chr=nb_windows_chr,
                centromere=centromere,
                geno_ref=config['genotype_ref'],
                geno_alt=config['genotype_alt'],
                ratio_min_homo=config['ratio_min_homo']
                )


    print("\noffspring_snps_window_smoothed: ")
    print("[chr_window] = [start, stop, ADref, ADalt, DP]")
    pprint.pprint(list(offspring_snps_window_smoothed.items())[2012:2017])

    print("\noffspring_genotype_window_normalized_smoothed")
    print("[chr_window] = [start, stop, ADref/DP, ADalt/DP, probHomoRef, probHetero, probHomoAlt, genotype]")
    pprint.pprint(list(offspring_genotype_window_smoothed.items())[2012:2017], width = 120, compact = True)

    print("\ncheck_smoothing[chr_window] = [start_window,stop_window,start_smooth,stop_smooth,window_smoothed")
    pprint.pprint(list(check_smoothing.items())[2012:2017])

    #------------------------------------------------------------------------------

    print("\nSave outputs of SmoothNormalizedOsffspringSlidingWindow()")

    export_dict_in_file(my_dict = offspring_snps_window_smoothed, 
                    output_file=outdir + "offspring_snps_window_normalized_smoothed.txt",
                    header="chr_window\tstart\tstop\tADref\tADalt\tDP",
                    overwrite=False)

    export_dict_in_file(my_dict = offspring_genotype_window_smoothed, 
                    output_file=outdir + "offspring_genotype_window_normalized_smoothed.txt",
                    header="chr_window\tstart\tstop\tADref/DP\tADalt/DP\tprobHomoRef\tprobHetero\tprobHomoAlt\tgenotype",
                    overwrite=False)

    export_dict_in_file(my_dict = check_smoothing, 
                    output_file=outdir + "check_smoothing.txt",
                    header="chr_window\tstart_window\tstop_window\tstart_smooth\tstop_smooth\twindow_smoothed",
                    overwrite=False)

..

9.bis. **Visualize Genotype Example**:

This section of the code is used to visualize genotype information based on the configuration settings. It maps file paths to genotype data and assigns colors to different genotypes for visual representation.

.. code-block:: python

    # STEP 8.bis. Visualize Genotype
    print("-------------------------------------")

    file_paths = {
        # 'directory_name': 'offspring_genotype_window_normalized_smoothed.txt'
    }

    file_paths[config['analyze_id']] = outdir + 'offspring_genotype_window_normalized_smoothed.txt'

    genotype_colors = {
        config['genotype_ref']: 'blue',
        config['genotype_alt']: 'red',
        'Col/Ct': 'green',
        'NA': 'gray'  
    }

    visualsize_genotope(file_paths, outdir, genotype_colors, column = 'genotype')




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

    candidates_co, db_co = IdentifyCOs(offspring_genotype_window_smoothed, nb_windows_chr)

    print("\nCandidateCOs: ")
    print("[chr_window] = [start_win:stop_win, co_start, co_stop, pre_geno, cur_geno]")
    pprint.pprint(list(candidates_co.items())[:5])

    export_dict_in_file(my_dict = candidates_co, 
                    output_file=outdir + "candidateCO.txt",
                    header="chr_mean_win\tstart_win:stop_win\tstart_co\tstop_co\tpre_geno\tcur_geno",
                    overwrite=True)


    ### IdentifyCOs:  Qichao version
    print("\nCandidateCOs Qichao: ")

    offspring_smooth_probs = OrderedDict()
    for key, value in  offspring_genotype_window_smoothed.items():
        offspring_smooth_probs[key] = value[0:2] + value[4:8]

    offspring_sliding_geno_ratios = OrderedDict()
    for key, value in  offspring_genotype_window_normalized.items():
        offspring_sliding_geno_ratios[key] = value[0:4]

    candidates_co_qichao = IdentifyCOsQichao(
        Offspring_smoothProbs= offspring_smooth_probs, 
        Offspring_smoothWinNums= nb_windows_chr, 
        Offspring_slidingGenoNums= offspring_snps_window_normalized,
        Offspring_slidingGenoRatios=offspring_sliding_geno_ratios, 
        Centromere=centromere,
        genoRef=config['genotype_ref'],
        genoAlt=config['genotype_alt'])

    pprint.pprint(list(candidates_co_qichao.items()))

    export_dict_in_file(my_dict = candidates_co_qichao, 
                    output_file=outdir + "candidates_co_qichao.txt",
                    header="chr_mean_win\tstart_win:stop_win\tstart_co\tstop_co\tpre_geno\tcur_geno",
                    overwrite=True)

..

10.bis **Pricise COs**:


    - The `IdentifyCOs` function is called to identify crossover events (COs) based on the new smoothed offspring genotype window data, which is adjusted according to the newly specified window size.
    - Precise COs are printed and exported to a file named "preciseCOs.txt" for further analysis and record-keeping.

.. code-block:: python

    print("\n# STEP 9.bis. Precise COs")
    print("-------------------------------------")


    preciseCOs = OrderedDict()
    db_co_new = OrderedDict()

    for i in range (len(list(candidates_co.keys()))):
        start_co=candidates_co[list(candidates_co.keys())[i]][1] - 100000
        end_co=candidates_co[list(candidates_co.keys())[i]][2] + 100000

        
        candidate = list(candidates_co.keys())[i].split("_")[0]
        new_offspring_genotype_window_smoothed = OrderedDict()
        for cle, valeur in offspring_genotype_window_smoothed_new.items():
            if valeur[0] >= start_co + 1 and valeur[1] <= end_co and cle.split("_")[0] == candidate:
                new_offspring_genotype_window_smoothed[cle] = valeur
                print(cle, valeur)
        
        start_window = int(list(new_offspring_genotype_window_smoothed.keys())[0].split("_")[1])
        end_window = int(list(new_offspring_genotype_window_smoothed.keys())[-1].split("_")[1])
        print(start_window, end_window)

        #print(new_offspring_genotype_window_smoothed.items())    

        candidates_co_2, db_co_2 = PreciseCOs(new_offspring_genotype_window_smoothed, start_window=start_window, end_window=end_window)
        
        if candidates_co_2 :
            candidates_co_2 = candidates_co_2
        else:
            start_co=candidates_co[list(candidates_co.keys())[i]][1] 
            end_co=candidates_co[list(candidates_co.keys())[i]][2]
            pre_geno=candidates_co[list(candidates_co.keys())[i]][3]
            cur_geno=candidates_co[list(candidates_co.keys())[i]][4]
            cur_chr = list(candidates_co.keys())[i].split("_")[0]
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
                    output_file=outdir + "preciseCOs_" + str(config['new_window_size']) + "_kb" + ".txt",
                    header="chr_mean_win\tstart_win:stop_win\tstart_co\tstop_co\t(stop_co-start_co)\tpre_geno\tcur_geno",
                    overwrite=True)

..


11. **Refine COs Border (Qichao version)**:
    - The COs identified using the Qichao version are also refined using the `RefineCOBordersQichao` function.
    - Refined COs according to this method are printed and exported to a file named "refinedCOs.txt".

.. code-block:: python

    print("\n# STEP 10. Refines COs border")
    print("-------------------------------------")

    print("\nRefineCOBorders")
    offspring_refinedCOs = RefineCOBorders(candidates_co, offspring_snps, window_size=config['window_size'])

    print("[chr_window] = [co_start, co_stop, pre_geno, cur_geno]")
    pprint.pprint(list(offspring_refinedCOs.items())[:5])

    export_dict_in_file(my_dict = offspring_refinedCOs, 
                    output_file=outdir + "refinedCOs.txt",
                    header="chr_mean_win\tstart_co\tstop_co\tpre_geno\tcur_geno",
                    overwrite=True)



    print("\nRefineCOBorders Qichao: ")
    offspring_refinedCOs_qichao = RefineCOBordersQichao(candidates_co_qichao, offspring_snps, window_size=config['window_size'])

    print("[chr_window] = [co_start, co_stop, pre_geno, cur_geno]")
    pprint.pprint(list(offspring_refinedCOs_qichao.items())[:5])

    export_dict_in_file(my_dict = offspring_refinedCOs_qichao, 
                    output_file=outdir + "refinedCOs.txt",
                    header="chr_mean_win\tstart_co\tstop_co\tpre_geno\tcur_geno",
                    overwrite=True)

..

Each step is accompanied by explanatory comments and print statements for tracking the progress of the script.
