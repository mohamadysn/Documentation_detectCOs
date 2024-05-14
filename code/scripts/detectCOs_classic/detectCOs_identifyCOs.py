#!/usr/bin/env python3
# -*- coding: utf-8 -*-



#==============================================================================
# Program		:	additional functions for detectCOs
# Author		:	Qichao Lian [qlian@mpipz.mpg.de]
# Corrector 	:   Maëla Sémery [maela.semery@inrae.fr]
# Corrector 2	:   Mohamad Yassine [mohamad.yassine@inrae.fr]
# Date			:	06.07.2022
# Last update	:	26.04.2024
# Version		:	2.0
#==============================================================================



import hues
from collections import OrderedDict
from detectCOs_required_functions import *

#==============================================================================


def IdentifyCOs(offspring_genotype_window_smoothed, nb_windows_chr): 
	"""
	Description:	Find location of hypothetical cross-over according to the change of genotype
	
	Input:			- offspring_genotype_window_smoothed: Ordered dictionnary 
					from SmoothNormalizedOsffspringSlidingWindow containing :
					[chr_window] = [start, stop, ADref/DP, ADalt/DP, probHomoRef, probHetero, probHomoAlt, genotype]
					
					- nb_windows_chr:	number of window per chromosome, determine with function
					NormalizeOffspringSlidingWindow()
	
	Output: 		Cross-over candidate from smoothed windows
					candidates_co[chr_mean_window] = [left_window_border + ":" + right_window_border, 
													co pos start, co pos stop, previous genotype, new genotype]

					All double cross-over found during the analysis are saved, even if they are under the threshold 
					to support the genotype (thr_support) : 
					db_co[chr_mean_window] = [cur_chr, start_co_win:stop_co_win, previous genotype, new genotype, status]

	Upgrade? : 	add input parameter : length required to support genotype (kb, default 300kb) and window size (kb)
				modify :
					if stop_win <= 10 or start_win > (nb_win - 10):
						thr_support = int((thr_support_kb/window_size) - 1) # (thr_support_kb/2)/(window_size/2) - 1
					else: 
						thr_support = int((thr_support_kb/(window_size/2)) - 1)	
	"""
	
	print()
	hues.info("Identify candidate crossovers")

	candidates_co = OrderedDict()
	db_co = OrderedDict() # save all homo to another homo genotype even if they are not enough supported
	thr_support = int() # number of windows to support genotype = threshold support

	for cur_chr, nb_win in nb_windows_chr.items(): ## iteration on each chromosome
		nb_win = int(nb_win)

		# information about the window num_win
		cur_support = 1 # number of window that support current genotype
		pre_num = 1  # number of window that support current genotype after a NA transition # modif MY
		cur_geno = str() # genotype of current window(=geno_stop)
		cur_win = int() # current window = num_win + 1 

		# information about the previous supported geonotype
		pre_geno = "" # genotype of the window before the hypothetical COs
		pre_win = int() # last window of the previous genotype (before hypothetical COs) 

		# count between previous genotype and current genotype
		count_NA = 0 # number of consecutive NA window
		support_preNA = 0 # support of genotype before starting NA region
		geno_preNA = "" # initialize geno_preNA at the start of processing each chromosome # modify MY


		for num_win in range(1, nb_win): # iteration on the each window of each chromosome cur_chr

			# Define genotype for window start and stop
			start_win = num_win # from 1 to the penultimate num_win (due to range function)
			key_start = cur_chr + "_" + str(start_win)
			geno_start = offspring_genotype_window_smoothed[key_start][7]

			stop_win = num_win + 1 # from 2 to the last num_win
			key_stop = cur_chr + "_" + str(stop_win)
			geno_stop = offspring_genotype_window_smoothed[key_stop][7]

			## Edit the thr_support = min number of window to validate a COs   # modify MY
			if start_win <= 10: 
				# for 9 first  slides   # modify MY
				if stop_win == 11 and geno_stop == geno_start != "NA":   # modify MY
					thr_support = 5 
				else:
					thr_support = 2  
			elif start_win >= (nb_win - 10):
				# for 9 last slides
				if start_win >= (nb_win - 10) and start_win <= (nb_win - 6) and geno_stop == geno_start == offspring_genotype_window_smoothed[cur_chr + "_" + str(nb_win - 10)][7] != "NA":  # modify MY
					thr_support = 5
				elif start_win == (nb_win - 10) and geno_stop != geno_start != "NA":
					thr_support = 5
				else:
					thr_support = 2
			else: # for all other windows 
				thr_support = 5

			## Initialization of the start and stop posistion of the CO
			co_start = int() # mean position of the start window of the crossover
			co_stop = int() # mean position of the stop window of the crossover
			stop_co_win = int() # number of window of the end of the crossover region


			if geno_start != geno_stop :

				if geno_start == "NA" or geno_stop == "NA":
					if geno_start != "NA" and geno_stop == "NA" : 
						# when enter in NA region:
						geno_preNA = geno_start
						win_preNA = start_win
						count_NA = 1
						support_preNA = cur_support
						continue
					
					elif geno_start == "NA" and geno_stop == geno_preNA :
						cur_support +=1
						support_preNA = cur_support  # modify MY


						if support_preNA >= thr_support and pre_geno != "" and pre_geno !=  geno_stop:	## voir cas 13   # modif MY
							cur_win = stop_win
							cur_geno = geno_stop
							stop_co_win = cur_win - cur_support - count_NA + 1 
							
							co_start = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][0] + \
							offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][1]) / 2)

							co_stop = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][0] + \
		    				offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][1]) / 2)

							key = cur_chr + "_" + str(round((pre_win + stop_co_win) / 2, 1))
							candidates_co[key] = [str(pre_win) + ":" + str(stop_co_win), co_start, co_stop, pre_geno, cur_geno]

							if ("/" not in pre_geno) and ("/" not in cur_geno) :
								db_co[key] = [cur_chr, str(pre_win) + ":" + str(stop_co_win), pre_geno, cur_geno, "supported"]
							
							pre_geno = ""                             # modify MY
							pre_win = int()                           # modify MY
							count_NA = 0                              # modify MY
							support_preNA = 0                         # modify MY
							geno_preNA = ""  # Also reset geno_preNA  # modify MY
							
						continue

					elif geno_start == "NA" and geno_stop != geno_preNA:
						if cur_support >= thr_support:
							# if geno_stop != geno_preNA, edit pre_geno and pre_win with the information before NA region
							pre_geno = geno_preNA # genotype before NA region
							pre_win = win_preNA # last window before NA region
							geno_preNA = ""
							win_preNA = ""
						else: # Added by MY
							pre_geno = pre_geno  # test MY
							pre_win = pre_win # test MY
							cur_support = pre_num # test MY
							cur_geno = geno_stop # test MY
							cur_win = stop_win # test MY
							support_preNA = 0 # test MY
							geno_preNA = "" # test MY
							continue # test MY

				elif cur_support >= thr_support:
					# if genotype change between start and stop num_win and
					# genotype supported by at least threshold support (2 or 5)
					pre_geno = geno_start # genotype of window before CO
					pre_win = start_win # window before CO
					support_preNA = 0                                        # modify MY
					pre_num = cur_support # useless # modify MY
				
				cur_support = 1 # reset cur_support

				if ("/" not in geno_start) and ("/" not in geno_stop) and geno_start != "NA" and geno_stop != "NA":
					key_db_co = cur_chr + "_" + str((start_win + stop_win) / 2)
					db_co[key_db_co] = [cur_chr, str(start_win) + ":" + str(stop_win), geno_start, geno_stop, "not_supported"]
		
				if pre_geno != "" and cur_support >= thr_support:  # Added by MY
					cur_geno = geno_stop
					cur_win = stop_win
					#log("Crossover candidate detected due to genotype change and validate thr_support in last window.")
					if support_preNA < thr_support and support_preNA != 0 : 
						# if there is some NA window during the count of cur_support
						stop_co_win = cur_win - cur_support - count_NA + 1
					else : 
						stop_co_win = cur_win - cur_support + 1

					co_start = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][0] + \
						offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][1]) / 2)

					co_stop = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][0] + \
		    			offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][1]) / 2)

					## Edit candidate
					key = cur_chr + "_" + str(round((pre_win + stop_co_win) / 2, 1))
					candidates_co[key] = [str(pre_win) + ":" + str(stop_co_win), co_start, co_stop, pre_geno, cur_geno]
					#log(f"key {key} supported by {cur_support} in candidateCO at last window.")

					if ("/" not in pre_geno) and ("/" not in cur_geno) :
						db_co[key] = [cur_chr, str(pre_win) + ":" + str(stop_co_win), pre_geno, cur_geno, "supported"]

					pre_geno = ""
					count_NA = 0
					support_preNA = 0
					pre_win = int()                           # test MY
			
			else: # geno_start == geno_stop 
				cur_geno = geno_stop
				cur_win = stop_win

				if cur_geno == "NA" :
					count_NA += 1 
					continue
				
				cur_support += 1

				if pre_geno == cur_geno:
					continue

				elif pre_geno != "" and cur_support >= thr_support: # modify MY
					if support_preNA < thr_support and support_preNA != 0 : 
						# if there is some NA window during the count of cur_support
						stop_co_win = cur_win - cur_support - count_NA + 1
					else : 
						stop_co_win = cur_win - cur_support + 1

					co_start = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][0] + \
						offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][1]) / 2)

					co_stop = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][0] + \
		    			offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][1]) / 2)

					## Edit candidate
					key = cur_chr + "_" + str(round((pre_win + stop_co_win) / 2, 1))
					candidates_co[key] = [str(pre_win) + ":" + str(stop_co_win), co_start, co_stop, pre_geno, cur_geno]

					if ("/" not in pre_geno) and ("/" not in cur_geno) :
						db_co[key] = [cur_chr, str(pre_win) + ":" + str(stop_co_win), pre_geno, cur_geno, "supported"]

					pre_geno = ""
					count_NA = 0
					support_preNA = 0

	### Raise warning for all double COs found in the data
	if len(db_co) != 0:
		hues.warn("Be careful, some double COs has been found during the analysis :")
		for key, value in db_co.items():
			# value = [cur_chr, start_co_win:stop_co_win, previous genotype, new genotype, status]
			print("- ", key, " (window ", value[1],"): ", value[2], ">", value[3], \
		 		" (", value[4], ")", sep="")

	return candidates_co, db_co

###############################################################################


def PreciseCOs(offspring_genotype_window_smoothed, start_window=None, end_window=None):   # Created by MY
	"""
	Find location of hypothetical cross-over according to the change of genotype
	
	:param offspring_genotype_window_smoothed: Ordered dictionary from SmoothNormalizedOsffspringSlidingWindow containing :
						[chr_window] = [start, stop, ADref/DP, ADalt/DP, probHomoRef, probHetero, probHomoAlt, genotype]
	:param nb_windows_chr: Number of window per chromosome, determined with function NormalizeOffspringSlidingWindow()
	:param start_window: Start window to begin analysis (inclusive)
	:param end_window: End window to stop analysis (inclusive)
	:return: Cross-over candidate from smoothed windows and all double cross-over found during the analysis
	"""
	
	# Initial setup
	candidates_co = OrderedDict()
	db_co = OrderedDict()
	thr_support = 5  
	# If start_window or end_window are not specified, use defaults
	if start_window is None:
		start_win = 1
	else:
		start_win = start_window
	if end_window is None:
		end_window = int(list(offspring_genotype_window_smoothed.keys())[-1].split("_")[1])
	
	cur_chr = list(offspring_genotype_window_smoothed.keys())[-1].split("_")[0]
	

	# Initial variables for tracking changes

	cur_support, cur_geno, cur_win = 1, '', 0
	pre_num = 1  # number of window that support current genotype after a NA transition # modif MY
	pre_geno, pre_win = '', 0
	count_NA, support_preNA, geno_preNA = 0, 0, ''

	for num_win in range(start_win, end_window):  # Adjusted to work within specified window range
		
		# Define genotype for window start and stop
		start_win = num_win # from 1 to the penultimate num_win (due to range function)
		key_start = cur_chr + "_" + str(start_win)
		geno_start = offspring_genotype_window_smoothed[key_start][7]

		stop_win = num_win + 1 # from 2 to the last num_win
		key_stop = cur_chr + "_" + str(stop_win)
		geno_stop = offspring_genotype_window_smoothed[key_stop][7]
		#(f"Window start genotype: {geno_start}, Window stop genotype: {geno_stop}")

		relative_position = num_win - start_window + 1
		relative_end_position = end_window - start_window + 1
		
		## Initialization of the start and stop posistion of the CO
		co_start = int() # mean position of the start window of the crossover
		co_stop = int() # mean position of the stop window of the crossover
		stop_co_win = int() # number of window of the end of the crossover region


		if geno_start != geno_stop :
			#print("Genotype change detected between consecutive windows.")
			
			if geno_start == "NA" or geno_stop == "NA":
				#print("Handling NA genotype scenario...")
				if geno_start != "NA" and geno_stop == "NA" :
					#print("Handling NA genotype scenario...")

					# when enter in NA region:
					geno_preNA = geno_start
					win_preNA = start_win
					count_NA = 1
					support_preNA = cur_support

					continue
				
				elif geno_start == "NA" and geno_stop == geno_preNA :
					cur_support +=1
					support_preNA = cur_support  # modify MY

					if support_preNA >= thr_support and pre_geno != "" and pre_geno !=  geno_stop:	## voir cas 13  # modif MY
						cur_win = stop_win
						cur_geno = geno_stop
						stop_co_win = cur_win - cur_support - count_NA + 1 
						
						co_start = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][0] + \
						offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][1]) / 2)

						co_stop = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][0] + \
						offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][1]) / 2)

						key = cur_chr + "_" + str(round((pre_win + stop_co_win) / 2, 1))
						candidates_co[key] = [str(pre_win) + ":" + str(stop_co_win), co_start, co_stop, pre_geno, cur_geno]

						if ("/" not in pre_geno) and ("/" not in cur_geno) :
							db_co[key] = [cur_chr, str(pre_win) + ":" + str(stop_co_win), pre_geno, cur_geno, "supported"]
						# Reset variables after detecting a CO
						pre_geno = "" # modify MY
						count_NA = 0  # modify MY
						support_preNA = 0 # modify MY
						geno_preNA = ""  # Also reset geno_preNA  # modify MY #test

					continue

				elif geno_start == "NA" and geno_stop != geno_preNA:
					if cur_support >= thr_support:
						# Logic to record crossover event
						# if geno_stop != geno_preNA, edit pre_geno and pre_win with the information before NA region
						pre_geno = geno_preNA # genotype before NA region
						pre_win = win_preNA # last window before NA region
						geno_preNA = ""
						win_preNA = ""
					else: # Added by MY
						pre_geno = pre_geno  # test MY
						pre_win = pre_win # test MY
						cur_support = pre_num # test MY
						cur_geno = geno_stop # test MY
						cur_win = stop_win # test MY
						support_preNA = 0 # test MY
						geno_preNA = "" # test MY
						continue # test MY
					

			elif cur_support >= thr_support:
				# if genotype change between start and stop num_win and
				# genotype supported by at least threshold support (2 or 5)
				pre_geno = geno_start # genotype of window before CO
				pre_win = start_win # window before CO
				support_preNA = 0 # modify MY
				pre_num = cur_support # useless
			
			cur_support = 1 # reset cur_support

			if ("/" not in geno_start) and ("/" not in geno_stop) and geno_start != "NA" and geno_stop != "NA":
				key_db_co = cur_chr + "_" + str((start_win + stop_win) / 2)
				db_co[key_db_co] = [cur_chr, str(start_win) + ":" + str(stop_win), geno_start, geno_stop, "not_supported"]

			if pre_geno != "" and cur_support >= thr_support:  # Added by MY
				cur_geno = geno_stop
				cur_win = stop_win
				#log("Crossover candidate detected due to genotype change and validate thr_support in last window.")
				if support_preNA < thr_support and support_preNA != 0 : 
					# if there is some NA window during the count of cur_support
					stop_co_win = cur_win - cur_support - count_NA + 1
				else : 
					stop_co_win = cur_win - cur_support + 1

				co_start = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][0] + \
					offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][1]) / 2)

				co_stop = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][0] + \
					offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][1]) / 2)

				## Edit candidate
				key = cur_chr + "_" + str(round((pre_win + stop_co_win) / 2, 1))
				candidates_co[key] = [str(pre_win) + ":" + str(stop_co_win), co_start, co_stop, pre_geno, cur_geno]
				#log(f"key {key} supported by {cur_support} in candidateCO at last window.")

				if ("/" not in pre_geno) and ("/" not in cur_geno) :
					db_co[key] = [cur_chr, str(pre_win) + ":" + str(stop_co_win), pre_geno, cur_geno, "supported"]

				pre_geno = ""
				count_NA = 0
				support_preNA = 0	
		
		else: # geno_start == geno_stop 
			cur_geno = geno_stop
			cur_win = stop_win

			if cur_geno == "NA" :
				count_NA += 1 
				continue
			
			cur_support += 1
			#print(f"Genotype {geno_stop} supported by {cur_support} consecutive windows.")
			
			if pre_geno == cur_geno:
				continue
			
			elif pre_geno != "" and cur_support >= thr_support:
				#print("Crossover candidate detected due to genotype change.")
				if support_preNA < thr_support and support_preNA != 0 : 
					# if there is some NA window during the count of cur_support
					stop_co_win = cur_win - cur_support - count_NA + 1
				else : 
					stop_co_win = cur_win - cur_support + 1

				co_start = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][0] + \
					offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][1]) / 2)

				co_stop = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][0] + \
					offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][1]) / 2)

				## Edit candidate
				key = cur_chr + "_" + str(round((pre_win + stop_co_win) / 2, 1))
				candidates_co[key] = [str(pre_win) + ":" + str(stop_co_win), co_start, co_stop, pre_geno, cur_geno]

				if ("/" not in pre_geno) and ("/" not in cur_geno) :
					db_co[key] = [cur_chr, str(pre_win) + ":" + str(stop_co_win), pre_geno, cur_geno, "supported"]

				pre_geno = ""
				count_NA = 0
				support_preNA = 0

	### Raise warning for all double COs found in the data
	if len(db_co) != 0:
		hues.warn("Be careful, some double COs has been found during the analysis :")
		for key, value in db_co.items():
			# value = [cur_chr, start_co_win:stop_co_win, previous genotype, new genotype, status]
			print("- ", key, " (window ", value[1],"): ", value[2], ">", value[3], " (", value[4], ")", sep="")

	
	return candidates_co, db_co


###############################################################################

def RefineCOBorders(CandidateCOs, Offspring_infoSNPs, window_size):
	print()
	hues.info("Refine crossover borders")

	RefinedCOs = OrderedDict() 

	for co_key, co_value in CandidateCOs.items():

		cur_chr, cur_pos = co_key.split("_")
		win, cur_start, cur_stop, pre_geno, cur_geno = co_value

		cur_start = int(cur_start - window_size / 2 * 1000)
		cur_stop = int(cur_stop + window_size / 2 * 1000)

		if ("/" not in pre_geno) and ("/" in cur_geno):

			cur_geno1, cur_geno2 = cur_geno.split("/")
			pre_geno_3 = pre_geno
			cur_geno_3 = ""
			comm_geno = ""

			if pre_geno == cur_geno1:
				cur_geno_3 = cur_geno2
			else:
				cur_geno_3 = cur_geno1
			comm_geno = pre_geno

			refined_infosnp = OrderedDict()
			refined_infosnp_cnt = 1
			refined_infosnp_pos = OrderedDict()
			for i in range(cur_start, cur_stop):
				key = cur_chr + "_" + str(i)

				if key in Offspring_infoSNPs.keys():
					snp_geno = Offspring_infoSNPs[key][5]

					if pre_geno_3 == snp_geno:
						refined_infosnp[refined_infosnp_cnt] = [i, 1]
						refined_infosnp_cnt += 1
						refined_infosnp_pos[i] = 1
					if cur_geno_3 == snp_geno:
						refined_infosnp[refined_infosnp_cnt] = [i, 2]
						refined_infosnp_cnt += 1
						refined_infosnp_pos[i] = 2

			refined_info_start = cur_start
			refined_info_stop = cur_stop
			for i in range(1, refined_infosnp_cnt):
				
				snp_pos, snp_type = refined_infosnp[i]
				if snp_type == 1 and snp_pos > refined_info_start:
					refined_info_start = snp_pos
				if snp_type == 2 and snp_pos < refined_info_stop:
					refined_info_stop = snp_pos
					if refined_info_start != cur_start:
						break

			if refined_info_stop < refined_info_start:
				refined_info_start = cur_start - 2000
				# hues.warn("ERROR-info-2, stop < start!")
				
			cur_start = int(cur_start + window_size / 2 * 1000)
			cur_stop = int(cur_stop - window_size / 2 * 1000)
			# hues.info("#CandidateCOs-2: " + str(cur_start) + ", " + str(cur_stop) + ", " + pre_geno + ", " + cur_geno)
			# hues.info("#RefinedInfoCOs-2: " + str(refined_info_start) + ", " + str(refined_info_stop) + ", " + pre_geno + ", " + cur_geno)
			RefinedCOs[co_key] = [refined_info_start, refined_info_stop, pre_geno, cur_geno]			

		elif ("/" in pre_geno) and ("/" not in cur_geno):

			pre_geno1, pre_geno2 = pre_geno.split("/")
			pre_geno_3 = ""
			cur_geno_3 = cur_geno
			comm_geno = ""

			if pre_geno1 == cur_geno:
				pre_geno_3 = pre_geno2
			else:
				pre_geno_3 = pre_geno1
			comm_geno = pre_geno1

			refined_infosnp = OrderedDict()
			refined_infosnp_cnt = 1
			refined_infosnp_pos = OrderedDict()
			for i in range(cur_start, cur_stop):
				key = cur_chr + "_" + str(i)

				if key in Offspring_infoSNPs.keys():
					snp_geno = Offspring_infoSNPs[key][5]

					if pre_geno_3 == snp_geno:
						refined_infosnp[refined_infosnp_cnt] = [i, 1]
						refined_infosnp_cnt += 1
						refined_infosnp_pos[i] = 1
					if cur_geno_3 == snp_geno:
						refined_infosnp[refined_infosnp_cnt] = [i, 2]
						refined_infosnp_cnt += 1
						refined_infosnp_pos[i] = 2

			refined_info_start = cur_start
			refined_info_stop = cur_stop
			for i in range(refined_infosnp_cnt - 1, 0, -1):
				
				snp_pos, snp_type = refined_infosnp[i]
				if snp_type == 1 and snp_pos > refined_info_start:
					refined_info_start = snp_pos

					if refined_info_stop != cur_stop:
						break

				if snp_type == 2 and snp_pos < refined_info_stop:
					refined_info_stop = snp_pos

			if refined_info_stop < refined_info_start:
				refined_info_stop = cur_stop + 2000
				hues.warn("ERROR-info-3, stop < start!")
				
			cur_start = int(cur_start + window_size / 2 * 1000)
			cur_stop = int(cur_stop - window_size / 2 * 1000)
			# hues.info("#CandidateCOs-3: " + str(cur_start) + ", " + str(cur_stop) + ", " + pre_geno + ", " + cur_geno)
			# hues.info("#RefinedInfoCOs-3: " + str(refined_info_start) + ", " + str(refined_info_stop) + ", " + pre_geno + ", " + cur_geno)
			RefinedCOs[co_key] = [refined_info_start, refined_info_stop, pre_geno, cur_geno]			

		else:

			pre_geno_3 = pre_geno
			cur_geno_3 = cur_geno
			comm_geno = ""

			refined_infosnp = OrderedDict()
			refined_infosnp_cnt = 1
			refined_infosnp_pos = OrderedDict()
			refined_totalsnp = OrderedDict()
			refined_totalsnp_cnt = 1
			for i in range(cur_start, cur_stop):
				key = cur_chr + "_" + str(i)

				if key in Offspring_infoSNPs.keys():
					snp_geno = Offspring_infoSNPs[key][5]

					if pre_geno_3 == snp_geno:
						refined_infosnp[refined_infosnp_cnt] = [i, 1]
						refined_infosnp_cnt += 1
						refined_infosnp_pos[i] = 1
					if cur_geno_3 == snp_geno:
						refined_infosnp[refined_infosnp_cnt] = [i, 2]
						refined_infosnp_cnt += 1
						refined_infosnp_pos[i] = 2

			refined_info_start = cur_start
			refined_info_stop = cur_stop
			for i in range(1, refined_infosnp_cnt):
				
				snp_pos, snp_type = refined_infosnp[i]
				if snp_type == 1 and snp_pos > refined_info_start:
					refined_info_start = snp_pos

				if snp_type == 2 and snp_pos < refined_info_stop:
					refined_info_stop = snp_pos

					if refined_info_start != cur_start:
						break

			if refined_info_stop < refined_info_start:
				refined_info_start = int(cur_start + window_size / 2 * 1000)
				refined_info_stop = int(cur_stop - window_size / 2 * 1000)
				hues.warn("ERROR-info-4, stop < start!")

			cur_start = int(cur_start + window_size / 2 * 1000)
			cur_stop = int(cur_stop - window_size / 2 * 1000)
			hues.info("#CandidateCOs-4: " + str(cur_start) + ", " + str(cur_stop) + ", " + pre_geno + ", " + cur_geno)
			hues.info("#RefinedInfoCOs-4: " + str(refined_info_start) + ", " + str(refined_info_stop) + ", " + pre_geno + ", " + cur_geno)
			RefinedCOs[co_key] = [refined_info_start, refined_info_stop, pre_geno, cur_geno]			

	hues.info("Re-refined crossovers")
	Re_RefinedCOs = OrderedDict()
	pre_co_chr = ""
	pre_co_pos = ""
	pre_co_start = ""
	pre_co_stop = ""
	pre_co_geno1 = ""
	pre_co_geno2 = ""
	cur_item = 0
	for co_key, co_value in RefinedCOs.items():

		co_chr, co_pos = co_key.split("_")
		co_start, co_stop, co_geno1, co_geno2 = co_value

		if co_chr != pre_co_chr:
			if pre_co_chr == "":
				pass
			else:
				Re_RefinedCOs[pre_co_chr + "_" + str(pre_co_pos)] = [pre_co_start, pre_co_stop, pre_co_geno1, pre_co_geno2]

			pre_co_chr = co_chr
			pre_co_pos = float(co_pos)
			pre_co_start = int(co_start)
			pre_co_stop = int(co_stop)
			pre_co_geno1 = co_geno1
			pre_co_geno2 = co_geno2
		else:
			if co_start <= pre_co_stop:
				if (pre_co_geno2 == co_geno1):
					if pre_co_geno1 == co_geno2:
						pre_co_chr = ""
						pre_co_pos = ""
						pre_co_start = ""
						pre_co_stop = ""
						pre_co_geno1 = ""
						pre_co_geno2 = ""
					else:
						pre_co_chr = co_chr
						pre_co_pos = round((pre_co_pos + float(co_pos)) / 2, 1)
						pre_co_start = int(pre_co_start)
						pre_co_stop = int(co_stop)
						pre_co_geno1 = pre_co_geno1
						pre_co_geno2 = co_geno2
						Re_RefinedCOs[pre_co_chr + "_" + str(pre_co_pos)] = [pre_co_start, pre_co_stop, pre_co_geno1, pre_co_geno2]
				else:
					hues.error("#3#" + "\t" + "Diff GENO!")
			else:
				if (co_start % 1000 == 0) and (pre_co_geno1 in co_geno1) and (pre_co_geno2 in co_geno1) and (pre_co_geno2 == co_geno2):
					pre_co_chr = co_chr
					pre_co_pos = round((pre_co_pos + float(co_pos)) / 2, 1)
					pre_co_start = int(pre_co_start)
					pre_co_stop = int(pre_co_stop)
					pre_co_geno1 = pre_co_geno1
					pre_co_geno2 = co_geno2
					Re_RefinedCOs[pre_co_chr + "_" + str(pre_co_pos)] = [pre_co_start, pre_co_stop, pre_co_geno1, pre_co_geno2]
				else:
					Re_RefinedCOs[pre_co_chr + "_" + str(pre_co_pos)] = [pre_co_start, pre_co_stop, pre_co_geno1, pre_co_geno2]
					pre_co_chr = co_chr
					pre_co_pos = round(float(co_pos), 1)
					pre_co_start = int(co_start)
					pre_co_stop = int(co_stop)
					pre_co_geno1 = co_geno1
					pre_co_geno2 = co_geno2

		cur_item += 1
		if cur_item == len(RefinedCOs.keys()):
			Re_RefinedCOs[pre_co_chr + "_" + str(pre_co_pos)] = [pre_co_start, pre_co_stop, pre_co_geno1, pre_co_geno2]
		
	return Re_RefinedCOs


###############################################################################
###############################################################################
def IdentifyCOsQichao(Offspring_smoothProbs, Offspring_smoothWinNums, \
					Offspring_slidingGenoNums, Offspring_slidingGenoRatios, \
					Centromere, genoRef, genoAlt):
	print()
	hues.info("Identify candidate crossovers")

	CandidateCOs = OrderedDict()

	window = 2

	for cur_chr, win_num in Offspring_smoothWinNums.items():
		min_win_num = 5
		max_win_num = int(win_num)

		cur_geno = ""
		cur_pos = 0
		cur_num = 1
		pre_geno = ""
		pre_pos = 0
		pre_num = 0

		cur1_geno = ""
		cur1_pos1 = 0
		cur1_pos2 = 0
		cur1_num = 1

		for x in range(1, max_win_num):

			start = x
			stop = start + window - 1 # start+1!

			if stop > win_num:
				stop = win_num

			if stop < 10 or start > (win_num - 10):
				min_win_num = 2
			else:
				min_win_num = 5

			key1 = cur_chr + "_" + str(start)
			geno1 = Offspring_smoothProbs[key1][5]
			key2 = cur_chr + "_" + str(stop)
			geno2 = Offspring_smoothProbs[key2][5]

			if geno1 == geno2: #geno
				if geno1 == "NA" or geno2 == "NA": # should not happen ?
					continue

				cur_geno = geno2
				cur_pos = stop
				cur_num += 1

				if pre_geno == "": # pre_geno only defined when geno1!=geno2 and cur_num > (min_win_num - 1)
					continue
				else:
					if pre_geno == "NA" or cur_geno == "NA": # pre_geno only defined when geno1!=geno2 and cur_num > (min_win_num - 1)
						continue
					else: # when a change in genotype is finally detected between geno1
						if cur_num < min_win_num or pre_geno == cur_geno:
							if cur_num > 1 and ("/" not in pre_geno) and ("/" in cur_geno):
								cur1_geno = cur_geno
								cur1_pos1 = cur_pos - cur_num + 1
								cur1_pos2 = cur_pos
								cur1_num = cur_num
							if cur1_geno != "" and pre_geno == cur_geno and cur_num >= min_win_num:
								hues.warn("Close double COs! -- same homo")
							else:
								continue

						if ("/" not in pre_geno) and ("/" not in cur_geno) and cur1_geno != "":
							# hues.warn("Close double COs!")

							if cur1_pos1 < pre_pos or cur1_pos2 > cur_pos:
								co_start = int( (Offspring_smoothProbs[cur_chr + "_" + str(pre_pos)][0] + Offspring_smoothProbs[cur_chr + "_" + str(pre_pos)][1]) / 2 )
								co_stop = int( (Offspring_smoothProbs[cur_chr + "_" + str(cur_pos - min_win_num + 1)][0] + Offspring_smoothProbs[cur_chr + "_" + str(cur_pos - min_win_num + 1)][1]) / 2 )
							else:
								if cur1_num < 4:
									if pre_geno != cur_geno:
										co_start = int( (Offspring_smoothProbs[cur_chr + "_" + str(pre_pos)][0] + Offspring_smoothProbs[cur_chr + "_" + str(pre_pos)][1]) / 2 )
										co_stop = int( (Offspring_smoothProbs[cur_chr + "_" + str(cur1_pos1)][0] + Offspring_smoothProbs[cur_chr + "_" + str(cur1_pos1)][1]) / 2 )
										key = cur_chr + "_" + str(round((pre_pos + cur1_pos1) / 2, 1))
										if key not in CandidateCOs:
											CandidateCOs[key] = [co_start, co_stop, pre_geno, cur1_geno]

										co_start = int( (Offspring_smoothProbs[cur_chr + "_" + str(cur1_pos2)][0] + Offspring_smoothProbs[cur_chr + "_" + str(cur1_pos2)][1]) / 2 )
										co_stop = int( (Offspring_smoothProbs[cur_chr + "_" + str(cur_pos - min_win_num + 1)][0] + Offspring_smoothProbs[cur_chr + "_" + str(cur_pos - min_win_num + 1)][1]) / 2 )
										key = cur_chr + "_" + str(round((cur1_pos2 + cur_pos - min_win_num + 1) / 2, 1))
										if key not in CandidateCOs:
											CandidateCOs[key] = [co_start, co_stop, cur1_geno, cur_geno]
									else:
										co_start = int( (Offspring_smoothProbs[cur_chr + "_" + str(pre_pos)][0] + Offspring_smoothProbs[cur_chr + "_" + str(pre_pos)][1]) / 2 )
										co_stop = int( (Offspring_smoothProbs[cur_chr + "_" + str(cur_pos - min_win_num + 1)][0] + Offspring_smoothProbs[cur_chr + "_" + str(cur_pos - min_win_num + 1)][1]) / 2 )
								else: # cur1_num = 4
									co_start = int( (Offspring_smoothProbs[cur_chr + "_" + str(pre_pos)][0] + Offspring_smoothProbs[cur_chr + "_" + str(pre_pos)][1]) / 2 )
									co_stop = int( (Offspring_smoothProbs[cur_chr + "_" + str(cur1_pos1)][0] + Offspring_smoothProbs[cur_chr + "_" + str(cur1_pos1)][1]) / 2 )
									key = cur_chr + "_" + str(round((pre_pos + cur1_pos1) / 2, 1))
									if key not in CandidateCOs:
										CandidateCOs[key] = [co_start, co_stop, pre_geno, cur1_geno]

									co_start = int( (Offspring_smoothProbs[cur_chr + "_" + str(cur1_pos2)][0] + Offspring_smoothProbs[cur_chr + "_" + str(cur1_pos2)][1]) / 2 )
									co_stop = int( (Offspring_smoothProbs[cur_chr + "_" + str(cur_pos - min_win_num + 1)][0] + Offspring_smoothProbs[cur_chr + "_" + str(cur_pos - min_win_num + 1)][1]) / 2 )
									key = cur_chr + "_" + str(round((cur1_pos2 + cur_pos - min_win_num + 1) / 2, 1))
									if key not in CandidateCOs:
										CandidateCOs[key] = [co_start, co_stop, cur1_geno, cur_geno]

							pre_geno = ""
							pre_pos = 0
							pre_num = 0

							continue
						else:						
							co_start = int( (Offspring_smoothProbs[cur_chr + "_" + str(pre_pos)][0] + Offspring_smoothProbs[cur_chr + "_" + str(pre_pos)][1]) / 2 )
							co_stop = int( (Offspring_smoothProbs[cur_chr + "_" + str(cur_pos - min_win_num + 1)][0] + Offspring_smoothProbs[cur_chr + "_" + str(cur_pos - min_win_num + 1)][1]) / 2 )

						cur1_geno = ""


						if co_stop - co_start > 300000:
							hues.warn("large breakpoint interval, correct it!")
							
							cor_pre_pos = pre_pos - 10
							if cor_pre_pos < 1:
								cor_pre_pos = 1
							
							cor_cur_pos = cur_pos + 8
							if cor_cur_pos > win_num - 1:
								cor_cur_pos = win_num - 1
							
							rawCO = [cur_chr, cor_pre_pos, cor_cur_pos, pre_geno, cur_geno]
							correctCO = CorrectLargeCOsQichao(Offspring_slidingGenoNums,\
										Offspring_slidingGenoRatios, rawCO, Centromere, genoRef, genoAlt)

							if len(correctCO.keys()) != 0:
								for key, value in correctCO.items():
									CandidateCOs[key] = value
							else:
								key = cur_chr + "_" + str(round((pre_pos + cur_pos - 2) / 2, 1))
								CandidateCOs[key] = [co_start, co_stop, pre_geno, cur_geno]
						else:
							key = cur_chr + "_" + str(round((pre_pos + cur_pos - 2) / 2, 1))
							CandidateCOs[key] = [co_start, co_stop, pre_geno, cur_geno]

						pre_geno = ""
						pre_pos = 0
						pre_num = 0
			else:
				if cur_num > (min_win_num - 1):
					pre_geno = geno1
					pre_pos = start
					pre_num = cur_num
				cur_num = 1

	return CandidateCOs


###############################################################################

def CorrectLargeCOsQichao(Offspring_slidingGenoNums, Offspring_slidingGenoRatios, rawCO, Centromere, genoRef,genoAlt, min_homo_freq=0.9):
	hues.info("Correct large crossovers (centromere regions)")

	CandidateCOs = OrderedDict()
	CorrectedCO = OrderedDict()
	cur_chr, rawCO_start, rawCO_stop, rawCO_geno1, rawCO_geno2 = rawCO

	cur_chr_centro_left = int(Centromere[cur_chr][0]) - 150 * 1000
	cur_chr_centro_right = int(Centromere[cur_chr][1]) + 150 * 1000

	window = 2
	sliding = 1

	cur_geno = ""
	cur_pos = 0
	cur_num = 1
	pre_geno = ""
	pre_pos = 0
	pre_num = 0

	pre_na_geno = ""
	pre_na_pos = 0
	pre_na_num = 0

	for x in range(rawCO_start, rawCO_stop + 1):

		start = x
		stop = start + window # x+2

		key1 = cur_chr + "_" + str(start) # window = x
		geno1 =	GetGenoWindow(Offspring_slidingGenoNums[key1], genoRef=genoRef, genoAlt=genoAlt)[7]

		key2 = cur_chr + "_" + str(stop - 1) # window = x +1
		geno2 = GetGenoWindow(Offspring_slidingGenoNums[key2], genoRef=genoRef, genoAlt=genoAlt)[7]
		
		if geno1 == geno2:

			if geno1 == "NA": # or geno2 == "NA":
				continue

			cur_geno = geno2
			cur_pos = stop - 1
			cur_num += 1

			if pre_geno == "":
				continue
			else:
				if pre_geno == "NA" or cur_geno == "NA":
					continue
				else:
					if cur_num < 3 or pre_geno == cur_geno:
						continue

					if pre_na_geno != "": # elif
						co_start = int( (Offspring_slidingGenoRatios[cur_chr + "_" + str(pre_pos)][0] + Offspring_slidingGenoRatios[cur_chr + "_" + str(pre_pos)][1]) / 2 )
						co_stop = int( (Offspring_slidingGenoRatios[cur_chr + "_" + str(pre_na_pos)][0] + Offspring_slidingGenoRatios[cur_chr + "_" + str(pre_na_pos)][1]) / 2 )
						key = cur_chr + "_" + str(round((pre_pos + pre_na_pos) / 2, 1))
						CandidateCOs[key] = [co_start, co_stop, pre_geno, pre_na_geno]

						pre_na_geno = ""
						pre_na_pos = 0
						pre_na_num = 0
					else:						
						co_start = int( (Offspring_slidingGenoRatios[cur_chr + "_" + str(pre_pos)][0] + Offspring_slidingGenoRatios[cur_chr + "_" + str(pre_pos)][1]) / 2 )
						co_stop = int( (Offspring_slidingGenoRatios[cur_chr + "_" + str(cur_pos - 2)][0] + Offspring_slidingGenoRatios[cur_chr + "_" + str(cur_pos - 2)][1]) / 2 )
						key = cur_chr + "_" + str(round((pre_pos + cur_pos - 2) / 2, 1))
						CandidateCOs[key] = [co_start, co_stop, pre_geno, cur_geno]

					pre_geno = ""
					pre_pos = 0
					pre_num = 0
		else:
			if cur_num > 2 and geno1 != "NA":
				pre_geno = geno1
				pre_pos = start
				pre_num = cur_num

			if geno2 == "NA" and geno1 != pre_geno and pre_geno != "":
				pre_na_geno = geno1
				pre_na_pos = start - cur_num + 1
				pre_na_num = cur_num

			if pre_pos > pre_na_pos:
				pre_na_geno = ""
				pre_na_pos = 0
				pre_na_num = 0

			if geno2 == "NA":
				continue
			elif geno1 == "NA" and geno2 == pre_geno:
				cur_num += 1
			else:
				cur_num = 1

	keepRight = False
	if "/" in rawCO_geno1:
		keepRight = True
	elif "/" in rawCO_geno2:
		keepRight = False
	else:
		pass
	needOneMore = False

	for co_key, co_value in CandidateCOs.items(): # not created yet !!!
		co_start = int(co_value[0])
		co_stop = int(co_value[1])
		pre_geno = co_value[2]
		cur_geno = co_value[3]

		if pre_geno == rawCO_geno1 and cur_geno == rawCO_geno2:
			if co_stop <= cur_chr_centro_left or co_start >= cur_chr_centro_right:
				if CorrectedCO and needOneMore == False:
					break
				else:
					CorrectedCO[co_key] = co_value
			else:
				if needOneMore:
					CorrectedCO[co_key] = co_value
					needOneMore = False
					break
				elif len(CandidateCOs.keys()) > 2 and keepRight == True and needOneMore == False:
					CorrectedCO = OrderedDict()
					CorrectedCO[co_key] = co_value
				elif len(CandidateCOs.keys()) > 2 and keepRight == False and needOneMore == False:
					CorrectedCO[co_key] = co_value
					break
				else:
					CorrectedCO[co_key] = co_value
		elif pre_geno == rawCO_geno1 and cur_geno != rawCO_geno2:
			CorrectedCO[co_key] = co_value
		elif pre_geno != rawCO_geno1 and cur_geno == rawCO_geno2:
			CorrectedCO[co_key] = co_value
		else:
			if len(CandidateCOs.keys()) > 2 and pre_geno == rawCO_geno2 and cur_geno == rawCO_geno1:
				if co_stop <= cur_chr_centro_left or co_start >= cur_chr_centro_right:
					CorrectedCO[co_key] = co_value
					needOneMore = True
				else:
					continue
			else:
				continue

	return CorrectedCO


###############################################################################

def RefineCOBordersQichao(CandidateCOs, Offspring_infoSNPs, window_size):
	print()
	hues.info("Refine crossover borders")

	RefinedCOs = OrderedDict() 

	for co_key, co_value in CandidateCOs.items():

		cur_chr, cur_pos = co_key.split("_")
		cur_start, cur_stop, pre_geno, cur_geno = co_value

		cur_start = int(cur_start - window_size / 2 * 1000)
		cur_stop = int(cur_stop + window_size / 2 * 1000)

		if ("/" not in pre_geno) and ("/" in cur_geno):

			cur_geno1, cur_geno2 = cur_geno.split("/")
			pre_geno_3 = pre_geno
			cur_geno_3 = ""
			comm_geno = ""

			if pre_geno == cur_geno1:
				cur_geno_3 = cur_geno2
			else:
				cur_geno_3 = cur_geno1
			comm_geno = pre_geno

			refined_infosnp = OrderedDict()
			refined_infosnp_cnt = 1
			refined_infosnp_pos = OrderedDict()
			for i in range(cur_start, cur_stop):
				key = cur_chr + "_" + str(i)

				if key in Offspring_infoSNPs.keys():
					snp_geno = Offspring_infoSNPs[key][5]

					if pre_geno_3 == snp_geno:
						refined_infosnp[refined_infosnp_cnt] = [i, 1]
						refined_infosnp_cnt += 1
						refined_infosnp_pos[i] = 1
					if cur_geno_3 == snp_geno:
						refined_infosnp[refined_infosnp_cnt] = [i, 2]
						refined_infosnp_cnt += 1
						refined_infosnp_pos[i] = 2

			refined_info_start = cur_start
			refined_info_stop = cur_stop
			for i in range(1, refined_infosnp_cnt):
				
				snp_pos, snp_type = refined_infosnp[i]
				if snp_type == 1 and snp_pos > refined_info_start:
					refined_info_start = snp_pos
				if snp_type == 2 and snp_pos < refined_info_stop:
					refined_info_stop = snp_pos
					if refined_info_start != cur_start:
						break

			if refined_info_stop < refined_info_start:
				refined_info_start = cur_start - 2000
				# hues.warn("ERROR-info-2, stop < start!")
				
			cur_start = int(cur_start + window_size / 2 * 1000)
			cur_stop = int(cur_stop - window_size / 2 * 1000)
			# hues.info("#CandidateCOs-2: " + str(cur_start) + ", " + str(cur_stop) + ", " + pre_geno + ", " + cur_geno)
			# hues.info("#RefinedInfoCOs-2: " + str(refined_info_start) + ", " + str(refined_info_stop) + ", " + pre_geno + ", " + cur_geno)
			RefinedCOs[co_key] = [refined_info_start, refined_info_stop, pre_geno, cur_geno]			

		elif ("/" in pre_geno) and ("/" not in cur_geno):

			pre_geno1, pre_geno2 = pre_geno.split("/")
			pre_geno_3 = ""
			cur_geno_3 = cur_geno
			comm_geno = ""

			if pre_geno1 == cur_geno:
				pre_geno_3 = pre_geno2
			else:
				pre_geno_3 = pre_geno1
			comm_geno = pre_geno1

			refined_infosnp = OrderedDict()
			refined_infosnp_cnt = 1
			refined_infosnp_pos = OrderedDict()
			for i in range(cur_start, cur_stop):
				key = cur_chr + "_" + str(i)

				if key in Offspring_infoSNPs.keys():
					snp_geno = Offspring_infoSNPs[key][5]

					if pre_geno_3 == snp_geno:
						refined_infosnp[refined_infosnp_cnt] = [i, 1]
						refined_infosnp_cnt += 1
						refined_infosnp_pos[i] = 1
					if cur_geno_3 == snp_geno:
						refined_infosnp[refined_infosnp_cnt] = [i, 2]
						refined_infosnp_cnt += 1
						refined_infosnp_pos[i] = 2

			refined_info_start = cur_start
			refined_info_stop = cur_stop
			for i in range(refined_infosnp_cnt - 1, 0, -1):
				
				snp_pos, snp_type = refined_infosnp[i]
				if snp_type == 1 and snp_pos > refined_info_start:
					refined_info_start = snp_pos

					if refined_info_stop != cur_stop:
						break

				if snp_type == 2 and snp_pos < refined_info_stop:
					refined_info_stop = snp_pos

			if refined_info_stop < refined_info_start:
				refined_info_stop = cur_stop + 2000
				hues.warn("ERROR-info-3, stop < start!")
				
			cur_start = int(cur_start + window_size / 2 * 1000)
			cur_stop = int(cur_stop - window_size / 2 * 1000)
			# hues.info("#CandidateCOs-3: " + str(cur_start) + ", " + str(cur_stop) + ", " + pre_geno + ", " + cur_geno)
			# hues.info("#RefinedInfoCOs-3: " + str(refined_info_start) + ", " + str(refined_info_stop) + ", " + pre_geno + ", " + cur_geno)
			RefinedCOs[co_key] = [refined_info_start, refined_info_stop, pre_geno, cur_geno]			

		else:

			pre_geno_3 = pre_geno
			cur_geno_3 = cur_geno
			comm_geno = ""

			refined_infosnp = OrderedDict()
			refined_infosnp_cnt = 1
			refined_infosnp_pos = OrderedDict()
			refined_totalsnp = OrderedDict()
			refined_totalsnp_cnt = 1
			for i in range(cur_start, cur_stop):
				key = cur_chr + "_" + str(i)

				if key in Offspring_infoSNPs.keys():
					snp_geno = Offspring_infoSNPs[key][5]

					if pre_geno_3 == snp_geno:
						refined_infosnp[refined_infosnp_cnt] = [i, 1]
						refined_infosnp_cnt += 1
						refined_infosnp_pos[i] = 1
					if cur_geno_3 == snp_geno:
						refined_infosnp[refined_infosnp_cnt] = [i, 2]
						refined_infosnp_cnt += 1
						refined_infosnp_pos[i] = 2

			refined_info_start = cur_start
			refined_info_stop = cur_stop
			for i in range(1, refined_infosnp_cnt):
				
				snp_pos, snp_type = refined_infosnp[i]
				if snp_type == 1 and snp_pos > refined_info_start:
					refined_info_start = snp_pos

				if snp_type == 2 and snp_pos < refined_info_stop:
					refined_info_stop = snp_pos

					if refined_info_start != cur_start:
						break

			if refined_info_stop < refined_info_start:
				refined_info_start = int(cur_start + window_size / 2 * 1000)
				refined_info_stop = int(cur_stop - window_size / 2 * 1000)
				hues.warn("ERROR-info-4, stop < start!")

			cur_start = int(cur_start + window_size / 2 * 1000)
			cur_stop = int(cur_stop - window_size / 2 * 1000)
			hues.info("#CandidateCOs-4: " + str(cur_start) + ", " + str(cur_stop) + ", " + pre_geno + ", " + cur_geno)
			hues.info("#RefinedInfoCOs-4: " + str(refined_info_start) + ", " + str(refined_info_stop) + ", " + pre_geno + ", " + cur_geno)
			RefinedCOs[co_key] = [refined_info_start, refined_info_stop, pre_geno, cur_geno]			

	hues.info("Re-refined crossovers")
	Re_RefinedCOs = OrderedDict()
	pre_co_chr = ""
	pre_co_pos = ""
	pre_co_start = ""
	pre_co_stop = ""
	pre_co_geno1 = ""
	pre_co_geno2 = ""
	cur_item = 0
	for co_key, co_value in RefinedCOs.items():

		co_chr, co_pos = co_key.split("_")
		co_start, co_stop, co_geno1, co_geno2 = co_value

		if co_chr != pre_co_chr:
			if pre_co_chr == "":
				pass
			else:
				Re_RefinedCOs[pre_co_chr + "_" + str(pre_co_pos)] = [pre_co_start, pre_co_stop, pre_co_geno1, pre_co_geno2]

			pre_co_chr = co_chr
			pre_co_pos = float(co_pos)
			pre_co_start = int(co_start)
			pre_co_stop = int(co_stop)
			pre_co_geno1 = co_geno1
			pre_co_geno2 = co_geno2
		else:
			if co_start <= pre_co_stop:
				if (pre_co_geno2 == co_geno1):
					if pre_co_geno1 == co_geno2:
						pre_co_chr = ""
						pre_co_pos = ""
						pre_co_start = ""
						pre_co_stop = ""
						pre_co_geno1 = ""
						pre_co_geno2 = ""
					else:
						pre_co_chr = co_chr
						pre_co_pos = round((pre_co_pos + float(co_pos)) / 2, 1)
						pre_co_start = int(pre_co_start)
						pre_co_stop = int(co_stop)
						pre_co_geno1 = pre_co_geno1
						pre_co_geno2 = co_geno2
						Re_RefinedCOs[pre_co_chr + "_" + str(pre_co_pos)] = [pre_co_start, pre_co_stop, pre_co_geno1, pre_co_geno2]
				else:
					hues.error("#3#" + "\t" + "Diff GENO!")
			else:
				if (co_start % 1000 == 0) and (pre_co_geno1 in co_geno1) and (pre_co_geno2 in co_geno1) and (pre_co_geno2 == co_geno2):
					pre_co_chr = co_chr
					pre_co_pos = round((pre_co_pos + float(co_pos)) / 2, 1)
					pre_co_start = int(pre_co_start)
					pre_co_stop = int(pre_co_stop)
					pre_co_geno1 = pre_co_geno1
					pre_co_geno2 = co_geno2
					Re_RefinedCOs[pre_co_chr + "_" + str(pre_co_pos)] = [pre_co_start, pre_co_stop, pre_co_geno1, pre_co_geno2]
				else:
					Re_RefinedCOs[pre_co_chr + "_" + str(pre_co_pos)] = [pre_co_start, pre_co_stop, pre_co_geno1, pre_co_geno2]
					pre_co_chr = co_chr
					pre_co_pos = round(float(co_pos), 1)
					pre_co_start = int(co_start)
					pre_co_stop = int(co_stop)
					pre_co_geno1 = co_geno1
					pre_co_geno2 = co_geno2

		cur_item += 1
		if cur_item == len(RefinedCOs.keys()):
			Re_RefinedCOs[pre_co_chr + "_" + str(pre_co_pos)] = [pre_co_start, pre_co_stop, pre_co_geno1, pre_co_geno2]
		
	return Re_RefinedCOs
 
 
