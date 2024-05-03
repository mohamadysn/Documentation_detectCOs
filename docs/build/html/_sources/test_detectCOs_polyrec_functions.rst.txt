test_detectCOs_polyrec_functions.py
====================================

This Python script is a launcher for the detectCOs tool. It's designed to be run from the command line with a YAML configuration file as input.

It imports necessary modules and functions from the detectCOs package, including functions for reading files, performing sliding window analysis, and identifying crossover events.

.. code-block:: python

    import os, sys
    import pprint
    import yaml

    from detectCOs_read_files import *
    from detectCOs_sliding_window import *
    from detectCOs_identifyCOs import *
..


It provides usage instructions for running the script from the command line, specifying the path to the configuration file as an argument.

1. **Loading Configuration :** It loads the YAML configuration file specified as a command line argument. The configuration file contains various parameters required for the analysis.

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
..

It sets the current working directory to the root of the Polyrec project as specified in the configuration.

2. **Reading EMS Lines :** It reads EMS lines from the configuration file using the `ReadEMSline` function. EMS lines are specific sets of genetic markers used in the analysis.

.. code-block:: python

    snps_list_line1, pos_snps_A = ReadEMSline(config['EMS_line1'], config['prefix_chr'])
    snps_list_line2, pos_snps_D = ReadEMSline(config['EMS_line2'], config['prefix_chr'])

    print("\nEMS_line1:")
    print(snps_list_line1[:5])
    print("\nEMS_line2:")
    print(snps_list_line2[:5])

    snps_list = ["./.", "0/0", "0/1", "0/2", "0/3", "1/1", "1/2", "1/3", "2/2", "2/3", "3/3"]
..

3. **Loading Parental SNPs :** It loads parental SNPs from a VCF file specified in the configuration using the `ReadParentalVCF` function. This function extracts SNP information including chromosome, position, and genotype.



.. code-block:: python

    print("\n# STEP 3. Load parental snps from vcf")
    print("-------------------------------------")

    parental_snps, parental_last_snp_chr = ReadParentalVCF(
        input_vcf=config['parental_vcf'], 
        prefix_chr=config['prefix_chr']
        )

    print("\nparental_snps[chr_pos] = [chr, pos, GT]")
    pprint.pprint(list(parental_snps.items())[:5])
    print("\nparental_last_snp_chr:")
    pprint.pprint(list(parental_last_snp_chr.items()))

    count_A = 0
    for A in pos_snps_A:
        if A in parental_snps.keys():
            count_A +=1 

    count_D = 0
    for D in pos_snps_D:
        if D in parental_snps.keys():
            count_D +=1 

    print("nb snps A in golden snps:\t", count_A)
    print("nb snps D in golden snps:\t", count_D)

..

4. **Printing Information :** It prints information about the loaded parental SNPs, such as the first few entries and the count of SNPs from each EMS line that match the parental SNPs.


Overall, this script serves as a central component for initiating the detectCOs tool's analysis pipeline, loading necessary data, and providing feedback on the loaded data.