#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#==============================================================================
# Program		:	additional functions for detectCOs
# Author		:	Qichao Lian [qlian@mpipz.mpg.de]
# Corrector 	:   Maëla Sémery [maela.semery@inrae.fr]
# Corrector 2	:   Mohamad Yassine [mohamad.yassine@inrae.fr]
# Date			:	06.07.2022
# Last update	:	24.08.2023
# Version		:	2.0
#==============================================================================

import hues
from collections import OrderedDict
from detectCOs_required_functions import *

#==============================================================================

def ParentalSlidingWindow(chr_len:OrderedDict, parental_snps:OrderedDict, \
	last_snps_chr:OrderedDict, window_size:int=100):
	"""
	Description: determine the number of SNPs per window (key= chr_window). 
	Input: 
	- chr_len:			length of each chromosome
	- parental_snps:	first element of the output of ReadParentalVCF()
	- last_snps_chr:	position of the last SNPs per chromosome
	- window_size: 		size of the window in kb (default = 100kb, max=1Mb/1000kb),
						the sliding size is half of window size
	Output: 
	- snps_window[chr_num-window] = [start, stop, nb-snps] 
	"""
	hues.info("Parental SNPs sliding window")
	
    # Check the window size 
	if window_size > 1000 :
		raise ValueError("The maximum of window size is 1Mb (1000 kb) !")
	elif window_size % 2 != 0:
		raise ValueError("Please enter an even number to avoid half step when calculting sliding size !")
	
	# convert window and sliding in pb
	window = window_size * 1000
	sliding = int(window / 2)
	hues.log("Parental SNPs:\twindow size = " + str(window_size) + \
		" kb / sliding size = " + str(window_size/2) + " kb")

	snps_window = OrderedDict()

	for cur_chr , last_snp_pos in last_snps_chr.items():
		chr_length = int(chr_len[cur_chr])
		last_snp_pos = int(last_snp_pos)
		nb_win = int(chr_length / sliding) # Note: round to inf integer

		hues.log(
			cur_chr + "\n- length:\t\t\t" + str(chr_length) + \
			"\n- position of the last SNPs:\t" + str(last_snp_pos) + \
			"\n- number of window:\t\t" + str(nb_win) 
		)

		for num_win in range(1, nb_win + 1):
			start = 1 + sliding * (num_win - 1)
			stop = window + sliding * (num_win - 1)
			if stop > chr_length:
				stop = chr_length

			cur_window = [start,stop,0] ## [start, stop, total_snps]

			for pos in range(start, stop + 1):
				key = cur_chr  + "_" + str(pos)
				if key in parental_snps.keys():
					cur_window[2] += 1

			snps_window[cur_chr  + "_" + str(num_win)] = cur_window

	return snps_window

#==============================================================================

def OffspringSlidingWindow(chr_len:OrderedDict, offspring_snps:OrderedDict,\
	last_snps_chr:OrderedDict, window_size:int=100):
	"""
	Description:	determine the sum of each parameters of all SNPs in the 
					window. Paramaters studied: ADref, ADalt, DP, SNP with at 
					least one reference allele,	SNP with at least one reference 
					allele, number of SNPs to the window size and the sliding 
					size.
	Input: 
	- chr_len:			length of each chromosome
	- offspring_snps:	first element of the output of ReadOffspringVCF()
	- last_snps_chr:	position of the last SNPs per chromosome
	- window_size:		size of the window in kb (default = 100kb),
						the sliding size is half of window size
	Output:
	- snps_window[chr_num-window] = [start, stop, ADref, ADalt, DP, 
									nbSNP_Aref, nbSNP_Aalt, TOTsnps-window] 
	Note: 
	nbSNP_Aref/alt => nb SNPs with at least one reference or alternative allele.
	"""
	hues.info("Offspring SNPs sliding window")

	# Check the window size 
	if window_size > 1000 :
		raise ValueError("The maximum of window size is 1Mb (1000 kb) !")
	elif window_size % 2 != 0:
		raise ValueError("Please enter a even number to avoid half step when calculting sliding size !")
	
	# Convert window size in pb
	window = window_size * 1000
	sliding = int(window / 2)
	hues.log("Offspring SNPs:\twindow size = " + str(window_size) + \
		" kb / sliding size = " + str(window_size/2) + " kb")
	
	snps_window = OrderedDict()
	
	for cur_chr, last_snp_pos in last_snps_chr.items():
		chr_length = int(chr_len[cur_chr])
		last_snp_pos = int(last_snp_pos)
		nb_win = int(chr_length / sliding) # Note: round to inf integer
		
		hues.log(cur_chr + "\n- length:\t\t\t" + str(chr_length) + \
	   			"\n- position of the last SNPs:\t" + str(last_snp_pos) + \
				"\n- number of window:\t\t" + str(nb_win))
		
		for num_win in range(1, nb_win + 1):
			start = 1 + sliding * (num_win - 1)
			stop = window + sliding * (num_win - 1)
			
			if stop > chr_length:
				stop = chr_length

			cur_window = [start,stop,0,0,0,0,0,0] 
			## cur_window = [start,stop,ADref,ADalt,DP,nbSNPref,nbSNPalt,TOTsnps-window] 

			for pos in range(start, stop + 1):
				key = cur_chr  + "_" + str(pos)

				if key in offspring_snps.keys():
					snp = offspring_snps[key] ## snp = [chr,pos,GT,ADref,ADalt,genotype]

					# HomoA
					if snp[2] == "0/0":
						cur_window[2] += int(snp[3]) # ADref += snp[ADref]
						cur_window[4] += int(snp[3]) # DP += snp[ADref]

						cur_window[5] += 1 # nbSNPref += 1
						cur_window[7] += 1 # TOTsnps-window += 1

					# HeteAB
					if snp[2] == "0/1":
						cur_window[2] += int(snp[3]) # ADref += snp[ADref]
						cur_window[4] += int(snp[3]) # DP += snp[ADref]
						cur_window[5] += 1 # nbSNPref += 1

						cur_window[3] += int(snp[4]) # ADalt += snp[ADalt]
						cur_window[4] += int(snp[4]) # DP += snp[ADalt]
						cur_window[6] += 1 # nbSNPalt += 1

						cur_window[7] += 1 # TOTsnps-window += 1

					# HomoB
					if snp[2] == "1/1":
						cur_window[3] += int(snp[4]) # ADalt += snp[ADalt]
						cur_window[4] += int(snp[4]) # DP += snp[ADalt]

						cur_window[6] += 1 # nbSNPalt += 1
						cur_window[7] += 1 # TOTsnps-window += 1
			
			snps_window[cur_chr  + "_" + str(num_win)] = cur_window

	return snps_window

#==============================================================================

def NormalizeOffspringSlidingWindow(parental_snps_window:OrderedDict, \
		offspring_snps_window:OrderedDict, geno_ref:str, geno_alt:str, \
		min_snp_num:int=16, min_reads_num:int=10, ratio_min_homo:float=0.9, depth_division_th:float=1.0):   # modify MY
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
			geno_window[key_window] = GetGenoWindow(cur_window, geno_ref, geno_alt, ratio_min_homo, depth_division_th)

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

def SmoothNormalizedOsffspringSlidingWindow(offspring_snps_window:OrderedDict,\
		nb_windows_chr:OrderedDict, centromere:OrderedDict, geno_ref:str,\
		geno_alt:str, ratio_min_homo:float=0.9, depth_division_th:float=1.0):           # modify MY
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
	- centromere:	ordered dictionnay with border of centromeric region for 
					each chromosome
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
		cen_left, cen_right = centromere[cur_chr]
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

				if pos_stop_cur_win <= cen_left or pos_start_cur_win >= cen_right:
					cur_window = [pos_start_cur_win, pos_stop_cur_win,\
		       		int(cur_window[2] + offspring_snps_window[key][2]),\
					int(cur_window[3] + offspring_snps_window[key][3]),\
					int(cur_window[4] + offspring_snps_window[key][4])]
					if win == "" :
						win = str(window)
					else:
						win = win + ":" + str(window)
				else:
					cur_window = [pos_start_cur_win, pos_stop_cur_win,\
		       		int(cur_window[2] + 0), int(cur_window[3] + 0), int(cur_window[4] + 0)]
			
			snps_window[key_window] = cur_window
			
			geno_window[key_window] = GetGenoWindow(cur_window, geno_ref, geno_alt, ratio_min_homo, depth_division_th)
			
			check_smoothed[key_window] = [pos_start_cur_win, pos_stop_cur_win, \
				start_win_smooth, stop_win_smooth, win]

	return snps_window, geno_window, check_smoothed 
