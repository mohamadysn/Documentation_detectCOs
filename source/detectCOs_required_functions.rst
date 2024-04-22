detectCOs_required_functions.py
================================


.. code-block:: python

    import math
    import numpy as np
    from decimal import *
    import scipy
..

Function : CalculateHomoProbability
------------------------------------

**Code :**

.. code-block:: python

    def CalculateHomoProbability(a:float, b:float):
        """
        Description:	calculate the probability to be homozygous a
        Input:			a = ratio ADa/DP ; b = ratio ADb/DP
        Output: 		probability to be homozygous a 
        """
        prob = scipy.special.comb(a+b, a) * pow(0.99, a) * pow(0.01, b)
        if math.isnan(prob) or math.isinf(prob):
            prob = np.nan_to_num(prob)
        return prob

..


**Explanation:**

`CalculateHomoProbability`: This function calculates the probability of being homozygous for a specific allele using the binomial law formula, where `a` is the ADa/DP ratio (number of reads supporting allele a divided by sequencing depth) and `b` is the ADb/DP ratio (number of reads supporting allele b divided by sequencing depth).
- **Input:**
  - `a`: ADa/DP ratio.
  - `b`: ADb/DP ratio.
- **Functionality:**
  - The function uses `scipy.special.comb()` to calculate the binomial coefficient (\(a + b\) choose \(a\)).
  - It multiplies the binomial coefficient by the probability of `a` successes (0.99) and `b` failures (0.01).
  - If the result is a non-finite floating number (`nan` or `inf`), it converts it to a finite number using `np.nan_to_num()`.
- **Output:**
  - `prob`: Probability of being homozygous for allele `a`.

This function is useful for estimating the probability of homozygosity of an allele based on the AD/DP ratios of both alleles in a given SNP.


Function : CalculateHeteroProbability
--------------------------------------

**Code :**

.. code-block:: python

    def CalculateHeteroProbability(a:float, b:float):
        """
        Description:	calculate the probability to be heterozygous
        Input:			a = ratio ADa/DP ; b = ratio ADb/DP
        Output:			probability to be heterozygous
        """
        prob = scipy.special.comb(a+b, a) * pow(0.5, a) * pow(0.5, b)
        if math.isnan(prob) or math.isinf(prob):
            prob = np.nan_to_num(prob)
        return prob
..

**Explanation:**

`CalculateHeteroProbability`: This function calculates the probability of being heterozygous for two different alleles. The probability is calculated using the binomial law formula, where `a` is the ADa/DP ratio (number of reads supporting allele a divided by sequencing depth) and `b` is the ADb/DP ratio (number of reads supporting allele b divided by sequencing depth).
- **Input:**
  - `a`: ADa/DP ratio.
  - `b`: ADb/DP ratio.
- **Functionality:**
  - The function uses `scipy.special.comb()` to calculate the binomial coefficient (\(a + b\) choose \(a\)).
  - It multiplies the binomial coefficient by the probability of `a` successes (0.5) and `b` failures (0.5).
  - If the result is a non-finite floating number (`nan` or `inf`), it converts it to a finite number using `np.nan_to_num()`.
- **Output:**
  - `prob`: Probability of being heterozygous.

This function is useful for estimating the probability of heterozygosity of a SNP based on the AD/DP ratios of both alleles.

Function : GetGenoWindow
-------------------------

**Code :**

.. code-block:: python

    def GetGenoWindow(cur_window:list, genoRef:str, genoAlt:str, min_homo_freq:float=0.9):
        """
        Description:	determine the genotype according to the ratio of ADref/DP 
                        and ADalt/DP
        Input:	cur_window = [start, stop, ADref, ADalt, DP]
        Output:	cur_geno = [start, stop, ADref/DP, ADalt/DP, probHomoRef, 
                            probHetero, probHomoAlt, genotype]
        """
            
        min_hetero_freq = 1 - min_homo_freq

        # Save each element of the list cur_window into specific variable
        start, stop, ADref, ADalt, depth = cur_window 

        # Calculate ratio ADref/DP and ADalt/DP
        if depth == 0 :
            ratio_ADref_DP = float(0.0)
            ratio_ADalt_DP = float(0.0)
        else : 
            ratio_ADref_DP = ADref/depth
            ratio_ADalt_DP = ADalt/depth

        # Define the probability to be homozygous Ref/Alt and heterozygous	
        prob_homoRef = prob_homoAlt = prob_hetero = 0.0

        if ratio_ADref_DP > ratio_ADalt_DP:
            if ratio_ADref_DP >= min_homo_freq or ratio_ADalt_DP <= min_hetero_freq:
                prob_homoRef = 1.0
            elif ratio_ADref_DP >= 0.3 and ratio_ADref_DP <= 0.7:
                prob_hetero = 1.0
            else:
                prob_homoRef = CalculateHomoProbability(ADref, ADalt)
                prob_hetero = CalculateHeteroProbability(ADref, ADalt)
                prob_homoAlt = CalculateHomoProbability(ADalt, ADref)
                if prob_homoRef == prob_homoAlt == prob_hetero == 0:
                    prob_hetero = 1.0

        elif ratio_ADref_DP == ratio_ADalt_DP and ratio_ADref_DP != 0:
            prob_hetero = 1.0

        elif ratio_ADref_DP < ratio_ADalt_DP:
            if ratio_ADalt_DP >= min_homo_freq or ratio_ADref_DP <= min_hetero_freq:
                prob_homoAlt = 1.0
            elif ratio_ADalt_DP >= 0.3 and ratio_ADalt_DP <= 0.7:
                prob_hetero = 1.0
            else:
                prob_homoRef = CalculateHomoProbability(ADref, ADalt)
                prob_hetero = CalculateHeteroProbability(ADref, ADalt)
                prob_homoAlt = CalculateHomoProbability(ADalt, ADref)
                if prob_homoRef == prob_homoAlt == prob_hetero == 0:
                    prob_hetero = 1.0

        else: # ratio ADref/DP == ADalt/DP == 0
            # do nothing: prob_homoRef = prob_homoAlt = prob_hetero = 0.0
            pass
        
        # Edit prob_geno
        prob_geno = [prob_homoRef, prob_hetero, prob_homoAlt, 'NA']
        if prob_homoRef > prob_hetero and prob_homoRef > prob_homoAlt:
            prob_geno[3] = genoRef
        elif prob_homoAlt > prob_hetero and prob_homoAlt > prob_homoRef:
            prob_geno[3] = genoAlt
        elif prob_hetero > prob_homoRef and prob_hetero > prob_homoAlt:
            prob_geno[3] = genoRef + "/" + genoAlt
        else: # prob_homoRef == prob_homoAlt == heteroAB == 0.0
            pass # do nothing: prob_geno = [0.0, 0.0, 0.0, "NA"]
        
        cur_geno = [start, stop, ratio_ADref_DP, ratio_ADalt_DP] + prob_geno

        return cur_geno

..



**Explanation:**

`GetGenoWindow`: This function determines the genotype of a genetic window based on the ratios of reads aligned to the reference and alternative alleles, as well as sequencing depth. It returns a list containing information about the genetic window, including probabilities of being homozygous for the reference and alternative alleles, the probability of being heterozygous, the predicted genotype, etc.
- **Input:**
  - `cur_window`: List containing elements of the window [start, stop, ADref, ADalt, DP].
  - `genoRef`: Reference genotype.
  - `genoAlt`: Alternative genotype.
  - `min_homo_freq`: Minimum homozygosity frequency (defaulted to 0.9).
- **Functionality:**
  - The function first calculates the ADref/DP and ADalt/DP ratios.
  - Then, it determines the probabilities of being homozygous for the reference genotype (`prob_homoRef`), homozygous for the alternative allele (`prob_homoAlt`), and heterozygous (`prob_hetero`) based on these ratios.
  - If the probabilities are not clear from the ratios, they are calculated using the `CalculateHomoProbability` and `CalculateHeteroProbability` functions.
  - Finally, it determines the genotype of the window based on the probabilities and returns it in `cur_geno`.
- **Output:**
  - `cur_geno`: List containing the elements of the window along with the determined probabilities and genotype.



Function : Updated GetGenoWindow
--------------------------------

**Code :**

.. code-block:: python

    def GetGenoWindow(cur_window:list, genoRef:str, genoAlt:str, min_homo_freq:float=0.75, depth_division_th:float=1.0):  # modif MY
        """
        Description:    determine the genotype according to the ratio of ADref/DP and ADalt/DP.
        Input:  cur_window = [start, stop, ADref, ADalt, DP, depth_division]
        Output: cur_geno = [start, stop, ADref/DP, ADalt/DP, genotype]
        """
        min_hetero_freq = 1 - min_homo_freq
        # Save each element of the list cur_window into specific variable
        start, stop, ADref, ADalt, depth = cur_window  

        # Calculate ratio ADref/DP and ADalt/DP
        if depth == 0:
            ratio_ADref_DP = ratio_ADalt_DP = 0.0
        else: 
            ratio_ADref_DP = round(ADref/depth, 3)  # modif MY
            ratio_ADalt_DP = round(ADalt/depth, 3) # modif MY
            
        
        # Define the probability to be homozygous Ref/Alt and heterozygous	
        prob_homoRef = prob_homoAlt = prob_hetero = 0.0

        if ratio_ADref_DP > ratio_ADalt_DP:
            if ratio_ADref_DP >= min_homo_freq or ratio_ADalt_DP <= min_hetero_freq:
                prob_homoRef = 1.0
            else:
                prob_hetero = 1.0

        elif ratio_ADref_DP == ratio_ADalt_DP and ratio_ADref_DP != 0:
            prob_hetero = 1.0

        elif ratio_ADref_DP < ratio_ADalt_DP:
            if ratio_ADalt_DP >= min_homo_freq or ratio_ADref_DP <= min_hetero_freq:
                prob_homoAlt = 1.0
            else:
                prob_hetero = 1.0

        else: # ratio ADref/DP == ADalt/DP == 0
            # do nothing: prob_homoRef = prob_homoAlt = prob_hetero = 0.0
            pass
        
        prob_geno = [prob_homoRef, prob_hetero, prob_homoAlt, 'NA']
        if depth >= 1200*depth_division_th :  # modif MY
            # Edit prob_geno
            if prob_homoRef > prob_hetero and prob_homoRef > prob_homoAlt:
                prob_geno[3] = genoRef
            elif prob_homoAlt > prob_hetero and prob_homoAlt > prob_homoRef:
                prob_geno[3] = genoAlt
            elif prob_hetero > prob_homoRef and prob_hetero > prob_homoAlt:
                prob_geno[3] = genoRef + "/" + genoAlt
            else: # prob_homoRef == prob_homoAlt == heteroAB == 0.0
                pass # do nothing: prob_geno = [0.0, 0.0, 0.0, "NA"]
        else:
            pass
        cur_geno = [start, stop, ratio_ADref_DP, ratio_ADalt_DP] + prob_geno
        return cur_geno


**Explanation:**

`GetGenoWindow`: This function determines the genotype of a genetic window based on the ratios of reads aligned to the reference and alternative alleles, as well as sequencing depth. It returns a list containing information about the genetic window, including probabilities of being homozygous for the reference and alternative alleles, the probability of being heterozygous, the predicted genotype, etc.
- **Input:**
  - `cur_window`: List containing elements of the window [start, stop, ADref, ADalt, DP].
  - `genoRef`: Reference genotype.
  - `genoAlt`: Alternative genotype.
  - `min_homo_freq`: Minimum homozygosity frequency (defaulted to 0.75).
- **Functionality:**
  - The function first calculates the ADref/DP and ADalt/DP ratios.
  - Then, it determines the probabilities of being homozygous for the reference genotype (`prob_homoRef`), homozygous for the alternative allele (`prob_homoAlt`), and heterozygous (`prob_hetero`) based on these ratios.
  - Finally, it determines the genotype of the window based on the probabilities and returns it in `cur_geno`.
- **Output:**
  - `cur_geno`: List containing the elements of the window along with the determined probabilities and genotype.
