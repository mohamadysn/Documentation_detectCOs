#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#==============================================================================
# Program		:	additional functions for detectCOs
# Author		:	Qichao Lian [qlian@mpipz.mpg.de]
# Corrector	:   Mohamad Yassine [mohamad.yassine@inrae.fr]
# Date			:	06.07.2022
# Last update	:	26.04.2024
# Version		:	2.0
#==============================================================================

import hues
from collections import OrderedDict
from detectCOs_required_functions_EMS import *

#==============================================================================

def ParentalSlidingWindowSNPs(parental_snps: OrderedDict, snps_per_window: int = 20):
    '''
    Description: determine the number of SNPs per window (key= chr_window). 
    Input: 
        - parental_snps:	first element of the output of ReadParentalVCF()
        - snps_per_window: fixed number of SNPs

    Output: 
        - snps_window[chr_num-window] = [start, stop, nb-snps] 
    
    Determine the SNP distribution over sliding windows, each containing a fixed number of SNPs,
    with each new window starting at the midpoint of the previous window in terms of SNP count.
    '''
    snps_window = OrderedDict()
    current_snps = []  # To keep track of SNPs in the current window
    cur_chr = ""  # Current chromosome
    window_index = 1  # Window index for naming
    half_window_size = snps_per_window // 2
    initial_window = True  # To control the first two special windows
    special_window_count = 0  # To keep track of the first two windows

    for key, snp_info in parental_snps.items():
        chr, pos, _ = snp_info

        if chr != cur_chr:
            if current_snps:
                try:
                    if len(current_snps) > half_window_size:
                        snps_window[f"{cur_chr}_{window_index}"] = [current_snps[0], current_snps[half_window_size-1], half_window_size]
                        window_index += 1
                        snps_window[f"{cur_chr}_{window_index}"] = [current_snps[half_window_size], current_snps[-1], len(current_snps) - half_window_size]
                    else:
                        snps_window[f"{cur_chr}_{window_index}"] = [current_snps[0], current_snps[-1], len(current_snps)]
                except IndexError as e:
                    print(f"Error handling remaining SNPs on chromosome change: {str(e)}")
            cur_chr = chr
            current_snps = []
            window_index = 1
            initial_window = True
            special_window_count = 0

        current_snps.append(pos)

        current_window_size = half_window_size if initial_window or special_window_count < 2 else snps_per_window
        if len(current_snps) == current_window_size:
            start_pos = current_snps[0]
            stop_pos = current_snps[-1]
            snps_window[f"{chr}_{window_index}"] = [start_pos, stop_pos, len(current_snps)]
            window_index += 1
            if initial_window:
                if special_window_count < 1:
                    current_snps = current_snps[half_window_size:]
                else:
                    initial_window = False
                    current_snps = []
                special_window_count += 1
            else:
                current_snps = current_snps[half_window_size:]

    if current_snps:
        try:
            if len(current_snps) > half_window_size:
                snps_window[f"{cur_chr}_{window_index}"] = [current_snps[0], current_snps[half_window_size-1], half_window_size]
                window_index += 1
                snps_window[f"{cur_chr}_{window_index}"] = [current_snps[half_window_size], current_snps[-1], len(current_snps) - half_window_size]
            else:
                snps_window[f"{cur_chr}_{window_index}"] = [current_snps[0], current_snps[-1], len(current_snps)]
        except IndexError as e:
            print(f"Error handling final SNPs: {str(e)}")

    return snps_window

#==============================================================================

def OffspringSlidingWindowBySNPs(offspring_snps: OrderedDict, snps_per_window: int = 20):
    '''
    Description: determine the number of SNPs per window (key= chr_window). 
    Input: 
        - offspring_snps:	first element of the output of ReadOffspringVCF()
        - snps_per_window: fixed number of SNPs

    Output: 
        - snps_window[chr_num-window] = [start, stop, nb-snps] 
    
    Determine the SNP distribution over sliding windows, each containing a fixed number of SNPs,
    with each new window starting at the midpoint of the previous window in terms of SNP count.
    '''    
    snps_window = OrderedDict()
    current_snps = []
    cur_chr = ""
    window_index = 1
    initial_window = True
    special_window_count = 0

    for key, snp_info in offspring_snps.items():
        chr, pos, GT, ADref, ADalt, genotype = snp_info
        ADref, ADalt = int(ADref), int(ADalt)

        if chr != cur_chr:
            if current_snps:  # Handle the last two windows of the previous chromosome
                process_final_windows(current_snps, snps_window, cur_chr, window_index, snps_per_window)
                current_snps = []  # Ensure no carry-over of SNPs
            # Reset for new chromosome
            cur_chr = chr
            window_index = 1
            initial_window = True
            special_window_count = 0

        current_snps.append((pos, ADref, ADalt, GT))

        # Determine if current window should be processed
        if initial_window:
            if len(current_snps) == snps_per_window // 2:
                snps_window[f"{chr}_{window_index}"] = create_window_output(current_snps, chr, window_index)
                window_index += 1
                current_snps = []
                special_window_count += 1
                if special_window_count == 2:
                    initial_window = False
        elif len(current_snps) >= snps_per_window:
            snps_window[f"{chr}_{window_index}"] = create_window_output(current_snps, chr, window_index)
            window_index += 1
            current_snps = current_snps[len(current_snps) // 2:]  # Start next window at midpoint of current

    # Handle any remaining SNPs for the last window of the last chromosome
    if current_snps:
        process_final_windows(current_snps, snps_window, cur_chr, window_index, snps_per_window)
        current_snps = []  # Ensure no carry-over of SNPs

    return snps_window

def create_window_output(snps, chr, window_index):
    ADref_sum = sum(snp[1] for snp in snps)
    ADalt_sum = sum(snp[2] for snp in snps)
    DP_sum = ADref_sum + ADalt_sum
    return [
        snps[0][0],  # start_pos
        snps[-1][0],  # stop_pos
        ADref_sum, 
        ADalt_sum, 
        DP_sum,
        sum("0" in snp[3] for snp in snps),  # nbSNP_Aref
        sum("1" in snp[3] for snp in snps),  # nbSNP_Aalt
        len(snps)  # count of SNPs
    ]

def process_final_windows(current_snps, snps_window, chr, window_index, snps_per_window):
    half_window_size = snps_per_window // 2
    try:
        if len(current_snps) >= 2 * half_window_size:
            print("Creating two half-sized windows")
            snps_window[f"{chr}_{window_index}"] = create_window_output(current_snps[:half_window_size], chr, window_index)
            window_index += 1
            snps_window[f"{chr}_{window_index}"] = create_window_output(current_snps[half_window_size:2*half_window_size], chr, window_index)
        elif len(current_snps) >= half_window_size:
            print("Creating one half-sized window")
            snps_window[f"{chr}_{window_index}"] = create_window_output(current_snps[:half_window_size], chr, window_index)
            window_index += 1
            if len(current_snps) > half_window_size:
                snps_window[f"{chr}_{window_index}"] = create_window_output(current_snps[half_window_size:], chr, window_index)
                window_index += 1
        else:
            print("Creating last window with remaining SNPs")
            snps_window[f"{chr}_{window_index}"] = create_window_output(current_snps, chr, window_index)
    except IndexError as e:
        print(f"Error processing final windows: {str(e)}")

#==============================================================================

def NormalizeOffspringSlidingWindowsnps(parental_snps_window:OrderedDict, \
		offspring_snps_window:OrderedDict, geno_ref:str, geno_alt:str, \
		min_snp_num:int=16, min_reads_num:int=10, ratio_min_homo:float=0.1, depth_division_th:float=1.0):   # modify MY
	"""
	Description:	dertermine the main genotype of each window according to the 
					ratio of ADref/DP and ADalt/DP and by calculating the 
					probability to be homozygous and heterozygous.
	Input:
	- parental_snps_window:		output of ParentalSlidingWindow()
	- offspring_snps_window:	output of OffspringSlidingWindow()
	- geno_ref: 			genotype of reference parent 
	- geno_alt: 			genotype of alternative parent
	- min_snp_num: 		minimum number of SNPs in the window to calculate the 
						ratio and determine the genotype (default: 16)
	- min_reads_num:	minimum number of coverage in the window to calculate 
						the ratio and determine the genotype (default: 10)
	- ratio_min_homo:	minimum frequency of AD/DP to be homozygous (float 
						between 0 and 1	excluded, default = 0.9)
	- depth_division_th:   depth division based on window size, default = 1 when window_size = 100 kb    # modify MY
	Output:
	- snps_window[chr_num-window] = [start, stop, ADref, ADalt, DP]
	- geno_window[chr_num-window] = [start, stop, ratio_ADref/DP, 
									ratio_ADalt_DP, prob_homo_ref, prob_hetero,
									prob_homo_alt, genotype]
	- nb_window_chr[chr_num-window] = number of window by chromosome
	"""

	print()
	hues.info("Normalize offspring sliding window")

	snps_window = OrderedDict() 
	# snps_window[chr_window] = [start, stop, ADref, ADalt, DP]
	geno_window = OrderedDict() 
	# geno_window[chr_window] = [start, stop, ADref/DP, ADalt/DP, probHomoA,probHeteroAB, probHomoB, genotype]
	nb_window_chr = OrderedDict()
	# nb_window_chr[chr] = nb_window

	if ratio_min_homo <= 0 or ratio_min_homo >= 1:
		raise ValueError("Invalid homozygous frequency. This value must be between 0 and 1 excluded !")
    
	for key_window, value in offspring_snps_window.items():
		## offspring_snps_window[chr_num-window] = [start,stop,ADref,ADalt,DP,nbSNPref,nbSNPalt,TOTsnps-window] 
		start, stop, ad_ref, ad_alt, dp = value[:5]
		
		if not key_window in parental_snps_window.keys():
			hues.warn("NOT include in parental sliding window!")
			raise ValueError("All Chr_window of the offspring must be in parental")
		
		else:
			## parental_snps_window[chr_num-window] = [start,stop,TOTsnps-window] 
			nb_snps_window = parental_snps_window[key_window][2] # = nb_snps_window
			cur_window = [start, stop, 0.0, 0.0, 0.0] # = [start, stop, ADref, ADalt, DP]

			if dp >= min_reads_num and nb_snps_window >= min_snp_num :
				cur_window = [start, stop, ad_ref, ad_alt, dp]
									
			snps_window[key_window] = cur_window
			geno_window[key_window] = GetGenoWindowsnps(cur_window, geno_ref, geno_alt, ratio_min_homo, depth_division_th)

			# save the last window of each chromosome in a dictionnary nb_window_chr[chr]=num_last_window
			chr = key_window.split("_")[0]
			win_id = int(key_window.split("_")[1])
			if chr in nb_window_chr:
				if win_id > nb_window_chr[chr]:
					nb_window_chr[chr] = win_id
			else:
				nb_window_chr[chr] = win_id

	return snps_window, geno_window, nb_window_chr

#==============================================================================

def SmoothNormalizedOsffspringSlidingWindowsnps(offspring_snps_window:OrderedDict,\
		nb_windows_chr:OrderedDict, geno_ref:str,\
		geno_alt:str, ratio_min_homo:float=0.1, depth_division_th:float=1.0):           # modify MY
	"""
	Description:	determine the sum of each parameters of all SNPs in the
					window. Paramaters studied: ADref, ADalt, DP, SNP with at 
					least one reference allele, SNP with at least one reference
					allele, number of SNPs to the window size and the sliding 
					size.
	Input: 
	- offspring_snps_window:	first output of NormalizeOffspringSlidingWindow()
	- nb_windows_chr:	number of window per chromosome, determine with function
						NormalizeOffspringSlidingWindow()
	- geno_ref:	genotype of reference
	- geno_alt:	genotype of alternative
	- ratio_min_homo:	minimum ratio of AD/DP to consider SNPs as homozygous
	- depth_division_th:   depth division based on window size, default = 1 when window_size = 100 kb              # modify MY
	Output:
	- snps_window[chr_num-window] = [start, stop, ADref, ADalt, DP]
	- geno_window[chr_num-window]  = [start, stop, ratio_ADref/DP, ratio_ADalt_DP,
	prob_homo_ref, prob_hetero, prob_homo_alt, genotype]
	- nb_window_chr[chr_num-window] = number of window by chromosome
	- check_smoothed : allow to check whether smoothing takes into account the right windows
	check_smoothed[chr_window] = [start_window, stop_window, start_smooth, stop_smooth] 
	Notes: 
	- All windows overlapping or inside centromeric region are associated 
	with NA genotype. But to smooth window next to centromeric region, we use 
	window overlapping or inside centromeric region.
	- ratio_min_hetero is automatically consider as 1-ratio_min_homo
	"""

	print()
	hues.info("Smooth normalized offspring sliding window")

	snps_window = OrderedDict() 
	geno_window = OrderedDict()
	nb_window_chr = OrderedDict()
	check_smoothed = OrderedDict()

	for cur_chr, nb_win_chr in nb_windows_chr.items():
		nb_window_chr[cur_chr] = nb_win_chr # nb window for each chromosome

		for num_win in range(1, nb_win_chr + 1): 
			key_window = cur_chr + "_" + str(num_win)
			start_win_smooth = num_win - 1
			stop_win_smooth = num_win + 1 

			if start_win_smooth < 1:
				start_win_smooth = 1
			if stop_win_smooth > nb_win_chr:
				stop_win_smooth = nb_win_chr

			pos_start_cur_win = int(offspring_snps_window[key_window][0])
			pos_stop_cur_win = int(offspring_snps_window[key_window][1])
			# Note: Offspring_slidingGenoWindow[chr_num-window] = [start,stop,ADref,ADalt,DP]
						
			cur_window = [pos_start_cur_win, pos_stop_cur_win, 0.0, 0.0, 0.0]
			# cur_window = [pos_start_cur_win, pos_stop_cur_wind, ADref, ADalt, depth]]

			win = ""
			for window in range(start_win_smooth, stop_win_smooth + 1):
				key = cur_chr + "_" + str(window) 
				cur_window = [pos_start_cur_win, pos_stop_cur_win,\
				int(cur_window[2] + offspring_snps_window[key][2]),\
				int(cur_window[3] + offspring_snps_window[key][3]),\
				int(cur_window[4] + offspring_snps_window[key][4])]
				if win == "" :
					win = str(window)
				else:
					win = win + ":" + str(window)

			
			snps_window[key_window] = cur_window
			
			geno_window[key_window] = GetGenoWindowsnps(cur_window, geno_ref, geno_alt, ratio_min_homo, depth_division_th)
			
			check_smoothed[key_window] = [pos_start_cur_win, pos_stop_cur_win, \
				start_win_smooth, stop_win_smooth, win]

	return snps_window, geno_window, check_smoothed

