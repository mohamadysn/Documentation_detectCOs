detectCOs_required_functions_EMS.py
=====================================


.. code-block:: python

    import math
    import numpy as np
    from decimal import *
    import scipy
..


Function : GetGenoWindowsnps
-----------------------------

**Code :**

.. code-block:: python

    def GetGenoWindowsnps(cur_window: list, genoRef: str, genoAlt: str, No_EMS_freq: float = 0.1, depth_division_th: float = 1.0):
        """
        Determine the genotype based on the allele frequency ratios (ADref/DP and ADalt/DP).

        Parameters:
        - cur_window: List containing the start, stop positions, ADref, ADalt, and DP.
        - genoRef: Reference genotype.
        - genoAlt: Alternative genotype.
        - No_EMS_freq: Frequency threshold below which a genotype is considered not EMS (mutated).
        - depth_division_th: Depth threshold to consider for genotyping.

        Returns:
        - cur_geno: List with start, stop positions, ADref/DP, ADalt/DP, and inferred genotype.
        """
        # Unpack the current window details
        start, stop, ADref, ADalt, depth = cur_window

        # Calculate the allele frequency ratios
        if depth == 0:
            ratio_ADref_DP = ratio_ADalt_DP = 0.0
        else:
            ratio_ADref_DP = round(ADref / depth, 3)
            ratio_ADalt_DP = round(ADalt / depth, 3)

        # Initialize the genotype determination as 'NA'
        genotype = 'NA'

        # Determine the genotype based on the ratios and depth threshold
        if depth >= 10 * depth_division_th:  # Adjust the multiplier as needed
            if  ratio_ADalt_DP >= No_EMS_freq:
                genotype = genoRef
            else:
                genotype = genoAlt

        # Compile the current genotype information
        cur_geno = [start, stop, ratio_ADref_DP, ratio_ADalt_DP, genotype]

        return cur_geno


**Explanation:**

`GetGenoWindow`: This function determines the genotype of a genetic window based on the ratios of reads aligned to the reference and alternative alleles, as well as sequencing depth. It returns a list containing information about the genetic window, including the predicted genotype, etc.
- **Input:**
  - `cur_window`: List containing elements of the window [start, stop, ADref, ADalt, DP].
  - `genoRef`: Reference genotype.
  - `genoAlt`: Alternative genotype.
  - `No_EMS_freq`: Minimum frequency (defaulted to 0.1).
- **Functionality:**
  - The function first calculates the ADref/DP and ADalt/DP ratios.
  - Finally, it determines the genotype of the window based on the probabilities and returns it in `cur_geno`.
- **Output:**
  - `cur_geno`: List containing the elements of the window along with the determined probabilities and genotype.
