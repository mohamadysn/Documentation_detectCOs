Usage
=====

.. _installation:

Installation
------------


To use detectCOs, first create a virtual environment using mamba:

.. code-block:: console

   $ mamba env create -f config_detectCOs_CG.yaml -p ~/miniconda3/envs/detectCOs_env
..

.. code-block:: console

   $ conda activate ~/miniconda3/envs/detectCOs_env
..

.. _launch:

launch
----------

.. code-block:: console

   $ python3 path/to/folder/Test_laucher_With_diff_SNPs.py config_detectCOs_CG.yaml
..


Preprocessing Step: Variant Calling with GATK HaplotypeCaller
--------------------------------------------------------------


Before moving on to the Crossing Over (COs) detection step, it's crucial to generate a Variant Call Format (VCF) file containing genetic variants for each individual in the population. This step is performed using the GATK (Genome Analysis Toolkit) tool with the HaplotypeCaller command.

Here's why this step is important:

1. **Genetic Variant Identification**: The generated VCF file contains information about genetic variants, such as Single Nucleotide Polymorphisms (SNPs) and indels (insertions and deletions of DNA), present in the analyzed DNA sample. These variants are identified by comparing the DNA data to a genomic reference.

2. **Reference SNP Detection**: The `--alleles` option allows specifying a VCF file containing reference SNPs. This helps improve the sensitivity and specificity of variant calling by informing the model about known and confirmed variants in the parental population.

3. **Region of Interest Selection**: The `-L` option specifies regions of interest for analysis. In this case, it appears to be used to limit the analysis to regions where reference SNPs are present, reducing computational burden by considering only relevant regions for CO analysis.

4. **Output Mode**: The `--output-mode EMIT_ALL_ACTIVE_SITES` option ensures that all active sites are emitted, including those that do not vary from the reference. This enables comprehensive analysis of sites, which is necessary for accurate CO detection.

In summary, this preprocessing step with GATK enables generating a comprehensive and accurate VCF file containing genetic variants, which will serve as the basis for subsequent CO detection.



.. code-block:: bash

   java -jar path/to/folder/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar HaplotypeCaller \
      --native-pair-hmm-threads 8 \
      -R path/to/folder/TAIR10_chr_all.fasta \
      -I path/to/folder/7RV498_recalibrated_sorted.bam \
      -O path/to/folder/7RV498_golden_SNPs.vcf \
      --alleles path/to/folder/parental_goldensnps_Indiv_28.vcf \
      -L path/to/folder/parental_goldensnps_Indiv_28.vcf \
      --output-mode EMIT_ALL_ACTIVE_SITES
..






