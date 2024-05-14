#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#==============================================================================
# Program		:	additional functions for detectCOs
# Author		:	Qichao Lian [qlian@mpipz.mpg.de]
# Corrector 	:   Maëla Sémery [maela.semery@inrae.fr]
# Corrector 2	:   Mohamad Yassine [mohamad.yassine@inrae.fr]
# Date			:	06.07.2022
# Last update	:	18.03.2024
# Version		:	2.0
#==============================================================================

import math
import numpy as np
from decimal import *
import scipy
import math

#==============================================================================

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


