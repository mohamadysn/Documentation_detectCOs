Test_yaml_detectCOs_BC_EMS.yaml
================================

CONFIG FILE FOR detectCOs TOOLS
---------------------------------

This configuration file specifies the paths to various files required for the analysis as well as the tool's parameters.

.. note::

    ALL PARAMETERS ARE MANDATORY !!! 

**1. General path informations**

.. code-block:: yaml

    path_to_polyrec_project: path/to/folder/polyrec/
    output_dir_in_polyrec: path/to/folder/
..


**2. Sample information**

Name of your sample, a directory `analyze_id` will be created in `output_dir_in_polyrec`

.. code-block:: yaml

    analyze_id: parental_Col-Ler_offspring_Plate1_A11
..

Prefix used or to be added in CHROM column in parental_vcf and offspring_vcf (they have to be the same between the two files).

.. code-block:: yaml

    prefix_chr: "Chr"
..


Path to chromosome length file

.. code-block:: yaml

    chr_len: path/to/folder/chr_len_athaliana_qichao.txt

..

Path to centromeric region file

.. code-block:: yaml

    centromere_reg: path/to/folder/centromere_athaliana_qichao.txt

..

Path to parental vcf

.. code-block:: yaml

    parental_vcf: path/to/folder/parental_golden_snps_139_genotyped_Indiv_28_merged_FORMATTED_IDpoint_noChr_dict_FORMAT.vcf

..

Path to offspring vcf

.. code-block:: yaml

    sample_vcf: path/to/folder/111355_AA_gt_by_golden_AWKED.vcf

..

**3. Tool parameters**

Reference genotype 
.. code-block:: yaml

    genotype_ref: Col

..

Alternative genotype 

.. code-block:: yaml

    genotype_alt: Ct
..

Size of the window (kb). The step size is automatically calcultaed as half the size window. (default: 100 kb)

.. code-block:: yaml

    window_size: 100
    
..

minimum number of coverage in the window to calculate the ratio and determine the genotype (default: 10)

.. code-block:: yaml

    min_reads_num: 10

..

minimum number of SNPs in the window to calculate the ratio and determine the genotype (default: 16)

.. code-block:: yaml

    min_snp_num: 16

..

minimum frequency of AD/DP to be homozygous (float between 0 and 1 excluded, default: 0.75)

.. code-block:: yaml

    ratio_min_homo: 0.75 
..


depth division based on window size, default = 1 when window_size = 100 kb                   # modify MY

.. code-block:: yaml

    depth_division_th: 0.5 

.. 


Size of the window (kb). To refined the position of CO                     # modify MY

.. code-block:: yaml

    new_window_size: 4

..

new depth division based on window size                                  # modify MY

.. code-block:: yaml

    new_depth_division_th: 0.5

..

Number of snps needed to define a window
.. code-block:: yaml

    snps_per_window: 6 
..

Reference genotype EMS
.. code-block:: yaml

    EMS_ref: A 

..

Alternative genotype EMS
.. code-block:: yaml

    EMS_alt: "-" 

..


minimum frequency of AD/DP to be EMS (float between 0 and 1 excluded, default: 0.1) # modify MY
.. code-block:: yaml

    ratio_min_EMS: 0.2 
..

Number of snps needed to define a window
.. code-block:: yaml

    new_snp_per_window: 4
..



This configuration file allows customization of the detectCOs tool's parameters based on the specific data to be analyzed.




