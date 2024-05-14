#!/usr/bin/env python3
# -*- coding: utf-8 -*-



#==============================================================================
# Program		:	additional functions for detectCOs for RIL
# Author		:   Mohamad Yassine [mohamad.yassine@inrae.fr]
# Date			:	30.04.2024
# Last update	:	30.04.2024
# Version		:	2.0
#==============================================================================



import hues
from collections import OrderedDict
from detectCOs_required_functions import *

#==============================================================================

import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()

def log(message):
    '''
    Écrit un message à la fois dans le fichier de log et sur la console.
    '''
    logger.info(message)

def IdentifyCOs_RILs(offspring_genotype_window_smoothed, nb_windows_chr): 
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
	log("Identify candidate crossovers")

	candidates_co = OrderedDict()
	db_co = OrderedDict() # save all homo to another homo genotype even if they are not enough supported
	thr_support = int() # number of windows to support genotype = threshold support

	for cur_chr, nb_win in nb_windows_chr.items(): ## iteration on each chromosome
		log(f"\nProcessing chromosome: {cur_chr} with {nb_win} windows.")

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
		win_preNA = int() # test MY

		
		chrom = []
		for key, value in offspring_genotype_window_smoothed.items():
			if cur_chr == key.split('_')[0]:
				chrom.append(key.split('_')[1])
		
		#nb_win = int(len(chrom))
		nb_win = int(chrom[-1])
		#print(nb_win)
		

		for num_win in range(len(chrom) -1): # iteration on the each window of each chromosome cur_chr
			
			#print(num_win)
			

			# Define genotype for window start and stop
			start_win = int(chrom[num_win]) # from 1 to the penultimate num_win (due to range function)
			log(f"\nAnalyzing window {start_win}/{nb_win} on chromosome {cur_chr}.")
			key_start = cur_chr + "_" + str(start_win)
			geno_start = offspring_genotype_window_smoothed[key_start][9]

			stop_win = int(chrom[num_win + 1])# from 2 to the last num_win
			key_stop = cur_chr + "_" + str(stop_win)
			geno_stop = offspring_genotype_window_smoothed[key_stop][9]

			## Edit the thr_support = min number of window to validate a COs   # modify MY
			if start_win <= 10: 
				# for 9 first  slides   # modify MY
				if stop_win == 11 and geno_stop == geno_start != "NA":   # modify MY
					thr_support = 5 
				else:
					thr_support = 2  
			elif start_win >= (nb_win - 10):
				# for 9 last slides
				if start_win >= (nb_win - 10) and start_win <= (nb_win - 6) and geno_stop == geno_start == offspring_genotype_window_smoothed[cur_chr + "_" + str(nb_win - 10)][9] != "NA":  # modify MY
					thr_support = 5
				elif start_win == (nb_win - 10) and geno_stop != geno_start != "NA":
					thr_support = 5
				else:
					thr_support = 2
			else: # for all other windows 
				thr_support = 5

			log(f"Window start genotype: {geno_start}, Window stop genotype: {geno_stop}")

			log(f"Current threshold support: {thr_support}.")
			## Initialization of the start and stop posistion of the CO
			co_start = int() # mean position of the start window of the crossover
			co_stop = int() # mean position of the stop window of the crossover
			stop_co_win = int() # number of window of the end of the crossover region


			if geno_start != geno_stop :
				log("Genotype change detected between consecutive windows.")

				if geno_start == "NA" or geno_stop == "NA":
					log("Handling NA genotype scenario...")
					if geno_start != "NA" and geno_stop == "NA" : 
						log("Handling geno_stop == 'NA' and  geno_start != 'NA' scenario...")
						log(f'geno_start : {geno_start}  , geno_stop :  {geno_stop} ')
						# when enter in NA region:
						geno_preNA = geno_start
						win_preNA = start_win
						count_NA = 1
						support_preNA = cur_support
						log(f'geno_preNA = {geno_preNA}, win_preNA = {win_preNA}, count_NA = {count_NA}, support_preNA = {support_preNA}')
						continue
					
					elif geno_start == "NA" and geno_stop == geno_preNA :
						log(f"Handling geno_start == 'NA' and  geno_stop ({geno_stop}) == geno_preNA ({geno_preNA}) scenario...")
						log(f'pre_win = {pre_win}, pre_geno = {pre_geno}')
						cur_support +=1
						log(f'cur_support = {cur_support}')
						support_preNA = cur_support  # modify MY
						log(f'support_preNA = {support_preNA}')


						if support_preNA >= thr_support and pre_geno != "" and pre_geno !=  geno_stop:	## voir cas 13   # modif MY
							cur_win = stop_win
							cur_geno = geno_stop
							stop_co_win = cur_win - cur_support - count_NA + 1 
							
							co_start = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][0] + \
							offspring_genotype_window_smoothed[cur_chr + "_" + str(pre_win)][1]) / 2)

							co_stop = int((offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][0] + \
		    				offspring_genotype_window_smoothed[cur_chr + "_" + str(stop_co_win)][1]) / 2)

							key = cur_chr + "_" + str(round((pre_win + stop_co_win) / 2, 1))
							log(f"key {key} supported by {cur_support} in candidateCO with NA.")
							candidates_co[key] = [str(pre_win) + ":" + str(stop_co_win), co_start, co_stop, pre_geno, cur_geno]

							if ("/" not in pre_geno) and ("/" not in cur_geno) :
								db_co[key] = [cur_chr, str(pre_win) + ":" + str(stop_co_win), pre_geno, cur_geno, "supported"]
							
							pre_geno = ""                             # modify MY
							pre_win = int()                           # test MY
							count_NA = 0                              # modify MY
							support_preNA = 0                         # modify MY
							geno_preNA = ""  # Also reset geno_preNA  # modify MY
							
						continue

					elif geno_start == "NA" and geno_stop != geno_preNA:
						log(f"Handling geno_start == 'NA' and  geno_stop ({geno_stop}) != geno_preNA ({geno_preNA}) scenario...")
						log(f'pre_win = {pre_win}, pre_geno = {pre_geno}')
						log(f'cur_support = {cur_support}, thr_support = {thr_support}')
						log(f'win_preNA = {win_preNA}, count_NA = {count_NA}, pre_num = {pre_num}')
						if cur_support >= thr_support:
							# if geno_stop != geno_preNA, edit pre_geno and pre_win with the information before NA region
							pre_geno = geno_preNA # genotype before NA region
							pre_win = win_preNA # last window before NA region
							geno_preNA = ""
							win_preNA = int()
						else:
							log(f' test test') # test MY
							log(f'cur_geno : {cur_geno}, cur_win : {cur_win}') # test MY
							pre_geno = pre_geno  # test MY
							pre_win = pre_win # test MY
							cur_support = pre_num # test MY
							cur_geno = geno_stop # test MY
							cur_win = stop_win # test MY
							support_preNA = 0 # test MY
							geno_preNA = "" # test MY
							continue # test MY



								


				elif cur_support >= thr_support:
					log(f"cur_support = {cur_support}, in the first side.")
					# if genotype change between start and stop num_win and
					# genotype supported by at least threshold support (2 or 5)
					pre_geno = geno_start # genotype of window before CO
					pre_win = start_win # window before CO
					support_preNA = 0                                        # modify MY
					pre_num = cur_support # useless  # modify MY
				
				cur_support = 1 # reset cur_support

				if ("/" not in geno_start) and ("/" not in geno_stop) and geno_start != "NA" and geno_stop != "NA":
					key_db_co = cur_chr + "_" + str((start_win + stop_win) / 2)
					db_co[key_db_co] = [cur_chr, str(start_win) + ":" + str(stop_win), geno_start, geno_stop, "not_supported"]
					log(f"We are here {key_db_co}.")
					log(f'pre_win = {pre_win}, pre_geno = {pre_geno}')
					log(f'geno_start = {geno_start}= pre_geno, geno_stop = {geno_stop}')
		
				if pre_geno != "" and cur_support >= thr_support:  # Added by MY
					cur_geno = geno_stop
					cur_win = stop_win
					log("Crossover candidate detected due to genotype change and validate thr_support in last window.")
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
					log(f"key {key} supported by {cur_support} in candidateCO at last window.")
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
				log(f"Genotype {geno_stop} supported by {cur_support} consecutive windows.")
				log(f"pre_geno {pre_geno}.")
				log(f"cur_geno {cur_geno}.")
				log(f"thr_support {thr_support}.")

				if pre_geno == cur_geno:
					continue

				elif pre_geno != "" and cur_support >= thr_support: # modify MY
					log("Crossover candidate detected due to genotype change and validate thr_support.")
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
					log(f"key {key} supported by {cur_support} in candidateCO without NA.")
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
			log("- " + key + " (window " + value[1] + "): " + value[2] + ">" + value[3] + " (" + value[4] + ")")

	return candidates_co, db_co


#==============================================================================



def Rules_RIls(info_RIL,offspring_EMS,offspring_basic ):
	"""
    Analyzes genotype data across multiple datasets to determine genetic rules and discrepancies.
    
    This function cross-references genotype data from basic and EMS offspring datasets with RIL information,
    identifies overlapping windows, and evaluates genotype discrepancies to establish a final genotype
    decision with associated warnings based on predefined rules.
    
    Parameters:
        info_RIL (dict): Dictionary containing RIL genotype information with window identifiers as keys.
        offspring_EMS (dict): Dictionary with EMS window data, where keys are window identifiers and
                              values include start and stop positions and genotypes.
        offspring_basic (dict): Dictionary similar to offspring_EMS but for basic genotype data.
    
    Returns:
        tuple: A tuple containing two dictionaries:
            1. RULEs (OrderedDict): Keys are window identifiers from offspring_basic and values are
                lists containing window start, stop, additional SNP info, and final genotype decisions.
            2. nb_window_chr (OrderedDict): Records the highest window number for each chromosome.
    """ 	
	RULEs = OrderedDict()   
	nb_window_chr = OrderedDict()
	for EMS_key, EMS_value in offspring_EMS.items():
		start_EMS, stop_EMS, *_ = EMS_value
		genotype_EMS = EMS_value[4] 
		
		for key, snp_info in offspring_basic.items():
			candidate = key
			chromosome= candidate.split("_")[0]
			start, stop, *_, genotype = snp_info
			if chromosome == EMS_key.split("_")[0]:
				if start_EMS <= stop and start <= stop_EMS :
					
					if candidate in info_RIL and info_RIL[candidate][:2] == [start, stop]:
						genotype_RIL = info_RIL[candidate][18]
					
					genotype_EMS = EMS_value[4]
					
					response = "Yes" if genotype_EMS == "A" else "No"
					
					if genotype_RIL == "Col":
						if genotype == "Col":
							if response == "Yes":
								WARNING = "A"
								GT_Final = "A"
							else:
								WARNING = "RIL"
								GT_Final = "RIL"						
						elif genotype == "Col/Ct":
							if response == "Yes":
								WARNING = "IMPO"
								GT_Final = "NA"
							else:
								WARNING = "IMPO"
								GT_Final = "NA"
						elif genotype == "Ct":
							if response == "Yes":
								WARNING = "IMPO"
								GT_Final = "NA"
							else:
								WARNING = "IMPO"
								GT_Final = "NA"
						else:
							WARNING = "GT_NA"
							GT_Final = "NA"
					elif genotype_RIL == "Ct":
						if genotype == "Col":
							if response == "Yes":
								WARNING = "A"
								GT_Final = "A"
							else:
								WARNING = "NEG_A"
								GT_Final = "A"			# modif			
						elif genotype == "Col/Ct":
							if response == "Yes":
								WARNING = "FALSE_A"
								GT_Final = "RIL"   # modif
							else:
								WARNING = "RIL"
								GT_Final = "RIL"
						elif genotype == "Ct":
							if response == "Yes":
								WARNING = "IMPO"
								GT_Final = "NA"
							else:
								WARNING = "IMPO"
								GT_Final = "NA"
						else:
							WARNING = "GT_NA"
							GT_Final = "NA"
					else:
						WARNING = "RIL_NA"
						GT_Final = "NA"

					RULEs[candidate] = [start,stop,snp_info[2],snp_info[3],genotype, genotype_EMS, genotype_RIL,response,WARNING,GT_Final]
				
				elif (start <= offspring_EMS[chromosome + "_" + str(1)][0] and stop <= offspring_EMS[chromosome + "_" + str(1)][0]):
					if candidate in info_RIL and info_RIL[candidate][:2] == [start, stop]:
							genotype_RIL = info_RIL[candidate][18]
					genotype_EMS = offspring_EMS[chromosome + "_" + str(1)][4]
					response = "Yes" if genotype_EMS == "A" else "No"
					
					if genotype_RIL == "Col":
						if genotype == "Col":
							if response == "Yes":
								WARNING = "A"
								GT_Final = "A"
							else:
								WARNING = "RIL"
								GT_Final = "RIL"						
						elif genotype == "Col/Ct":
							if response == "Yes":
								WARNING = "IMPO"
								GT_Final = "NA"
							else:
								WARNING = "IMPO"
								GT_Final = "NA"
						elif genotype == "Ct":
							if response == "Yes":
								WARNING = "IMPO"
								GT_Final = "NA"
							else:
								WARNING = "IMPO"
								GT_Final = "NA"
						else:
							WARNING = "GT_NA"
							GT_Final = "NA"
					elif genotype_RIL == "Ct":
						if genotype == "Col":
							if response == "Yes":
								WARNING = "A"
								GT_Final = "A"
							else:
								WARNING = "NEG_A"
								GT_Final = "A"			# modif			
						elif genotype == "Col/Ct":
							if response == "Yes":
								WARNING = "FALSE_A"
								GT_Final = "RIL"   # modif
							else:
								WARNING = "RIL"
								GT_Final = "RIL"
						elif genotype == "Ct":
							if response == "Yes":
								WARNING = "IMPO"
								GT_Final = "NA"
							else:
								WARNING = "IMPO"
								GT_Final = "NA"
						else:
							WARNING = "GT_NA"
							GT_Final = "NA"
					else:
						WARNING = "RIL_NA"
						GT_Final = "NA"

					RULEs[candidate] = [start,stop,snp_info[2],snp_info[3],genotype, genotype_EMS, genotype_RIL,response,WARNING,GT_Final]

				elif (start > offspring_EMS[EMS_key][1] and stop > offspring_EMS[EMS_key][1]):
					if candidate in info_RIL and info_RIL[candidate][:2] == [start, stop]:
						genotype_RIL = info_RIL[candidate][18]
					genotype_EMS = offspring_EMS[EMS_key][4]
					response = "Yes" if genotype_EMS == "A" else "No"
					
					if genotype_RIL == "Col":
						if genotype == "Col":
							if response == "Yes":
								WARNING = "A"
								GT_Final = "A"
							else:
								WARNING = "RIL"
								GT_Final = "RIL"						
						elif genotype == "Col/Ct":
							if response == "Yes":
								WARNING = "IMPO"
								GT_Final = "NA"
							else:
								WARNING = "IMPO"
								GT_Final = "NA"
						elif genotype == "Ct":
							if response == "Yes":
								WARNING = "IMPO"
								GT_Final = "NA"
							else:
								WARNING = "IMPO"
								GT_Final = "NA"
						else:
							WARNING = "GT_NA"
							GT_Final = "NA"
					elif genotype_RIL == "Ct":
						if genotype == "Col":
							if response == "Yes":
								WARNING = "A"
								GT_Final = "A"
							else:
								WARNING = "NEG_A"
								GT_Final = "A"			# modif			
						elif genotype == "Col/Ct":
							if response == "Yes":
								WARNING = "FALSE_A"
								GT_Final = "RIL"   # modif
							else:
								WARNING = "RIL"
								GT_Final = "RIL"
						elif genotype == "Ct":
							if response == "Yes":
								WARNING = "IMPO"
								GT_Final = "NA"
							else:
								WARNING = "IMPO"
								GT_Final = "NA"
						else:
							WARNING = "GT_NA"
							GT_Final = "NA"
					else:
						WARNING = "RIL_NA"
						GT_Final = "NA"

					RULEs[candidate] = [start,stop,snp_info[2],snp_info[3],genotype, genotype_EMS, genotype_RIL,response,WARNING,GT_Final]

				chr = candidate.split("_")[0]
				win_id = int(candidate.split("_")[1])
				if chr in nb_window_chr:
					if win_id > nb_window_chr[chr]:
						nb_window_chr[chr] = win_id
				else:
					nb_window_chr[chr] = win_id
	return RULEs, nb_window_chr






