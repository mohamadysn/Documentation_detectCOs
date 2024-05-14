#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#==============================================================================
# Program		:	additional functions for detectCOs
# Author		:	Qichao Lian [qlian@mpipz.mpg.de]
# Corrector 	:   Maëla Sémery [maela.semery@inrae.fr]
# Corrector 2	:   Mohamad Yassine [mohamad.yassine@inrae.fr]
# Date			:	06.07.2022
# Last update	:	25.08.2023
# Version		:	2.0
#==============================================================================

import hues,os,sys,re
from collections import OrderedDict
from detectCOs_required_functions_EMS import *

#==============================================================================

def CheckInput(input_file:str):
	"""
	Description:	Check if the file exist and if it's not empty
	Input: 			path-to-file/input_file
	"""
	if os.path.exists(input_file):
		hues.log("Found input file:\t" + input_file)

		if os.path.getsize(input_file):
			hues.log("Input file:\t" + input_file + " is not empty")
		else:
			hues.error("Input file:\t" + input_file + " is empty")
			sys.exit()
	else:
		hues.error("Not found input file:\t" + input_file)
		sys.exit()

#==============================================================================

def RemoveFile(output_file:str): #old name CheckOuput
	"""
	Description:	Remove output file if it exist
	Input:			path-to-file/output_file
	"""
	if os.path.exists(output_file):
		hues.warn("Found output file:\t" + output_file)
		hues.warn("Remove\t" + output_file + "\t...")
		os.remove(output_file)

#==============================================================================

def load_file_in_dict(input_file:str, prefix_chr:str="Chr"):
	"""
	Description:	convert file in dictionnary
	Input:			Tab-separated file
	Output: 		my_dict[value_col0] = [value_col1, ..., value_coln]
	"""
	my_dict = OrderedDict()
	values_list = []

	with open(input_file) as input:
		for line in input:
			value = line.strip().split("\t")
			# .strip() remove space at the beginning and at the end of the string
			if value[0].startswith(prefix_chr):
				values_list = value[1:]
				for i in range(len(values_list)):
					# convert all element in value_list if it match with
					# regex_integer and regex_float
					regex_int = "\d+" # equal to "[0-9]+"
					regex_float = "\d+\.\d+" # equal to "[0-9]+\.[0-9]+"

					if re.fullmatch(regex_int, values_list[i]):
						values_list[i] = int(values_list[i])
					elif re.fullmatch(regex_float, values_list[i]):
						values_list[i] = float(values_list[i])

				if len(values_list) == 1 :
					# if the list has only 1 element, return the element and not the list
					my_dict[value[0]] = values_list[0]
				else:	
					my_dict[value[0]] = values_list

	return my_dict

#==============================================================================

def export_dict_in_file(my_dict:OrderedDict, output_file:str, header:str, overwrite:bool=False):
	"""
	Description: 	save dictionnary into file
	Input: 
	- my_dict: 		the dictionnary to convert
	- output_file: 	path to output file
	- header: 		Each column must be separated by tabulation "\t".
					If header = "", no header will be added.
	- overwrite:	Overwrite file if it already exist
	"""

	if os.path.exists(output_file) and overwrite is False:
		print(output_file, " already exists.")
		exit
	
	else:
		if os.path.exists(output_file) and overwrite is True:
			RemoveFile(output_file)

		with open(output_file, "w") as output:
			# write header
			if header != "" :
				output.write(header + "\n")

			# write one item key, value per line
			for key, value in my_dict.items():
				line = str(key)

				if isinstance(value, list):
					for i in range(len(value)):
						line = line + "\t" + str(value[i])
				else:
					line = line + "\t" + str(value)
				output.write(line + "\n")

			print(output_file, "created.")

#==============================================================================

def ReadChrLen(input_chr_len:str, prefix_chr:str="Chr"):
	"""
	Description:	create a dictionnary containing the length of each 
					chromosome
	Input:	
	- input_chr_len:	Tab-separated file containing the name of each 
						chromosome and its length.
	- prefix_str:	Prefix of the chromosome number (e.g. "Chr", "chr", etc.)
	Output: 
	- chr_len[chr_x] = length-of-chromosome-x
	"""
	CheckInput(input_chr_len)
	chr_len = OrderedDict()

	if prefix_chr == "":
		raise ValueError("Invalid prefix chromosome")

	with open(input_chr_len, 'r') as input_file:
		for line in input_file:
			lines = line.strip().split("\t")
			if lines[0].startswith(prefix_chr):
				chr_len[lines[0]] = lines[1]
			else:
				chr_len[prefix_chr + lines[0]] = lines[1]


	hues.log(str(len(chr_len)) + " Chromosome length loaded!")
	return chr_len

#==============================================================================

def ReadCentroReg(input_centro_reg:str, prefix_chr="Chr"):
	"""
	Description:	create a dictionnary containing the coordonate of the 
					pericentromeric region of each chromosome
	Input:	
	- input_centro_reg: Tab-separated file containing the name of each chromosome 
						(line) and the position of the beginning (left border) 
						and the end (right border) of the pericentromeric region.
	- prefix_str:	Prefix of the chromosome number (e.g. "Chr", "chr", etc.)
	Output:	
	- ordered dictionnary: centromere[chr_x] = [left-border, right-border]
	"""
	CheckInput(input_centro_reg)
	centromere = OrderedDict()

	if prefix_chr == "":
		raise ValueError("Invalid prefix chromosome")

	with open(input_centro_reg, 'r') as input_file:
		for line in input_file:
			lines = line.strip().split("\t")

			if len(lines) == 3:
				if lines[0] == prefix_chr: # ignore header if it exists
					continue
				region = [int(lines[1]), int(lines[2])]
				if lines[0].startswith(prefix_chr):
					chr = lines[0]
				else:
					chr = prefix_chr + lines[0]
			else:
				raise ValueError ("Invalid input file !")

			centromere[chr] = region

	hues.log(str(len(centromere)) + " Chromosome centromere regions loaded!")
	return centromere


#==============================================================================

def ReadParentalVCF(input_vcf:str, prefix_chr:str="Chr"):
	"""
	Description:	create 2 dictionnaries containing all SNPs per chromosome 
					and the last SNPs per chromosome
	Input: 
	- input_vcf:	SNP markers between parental lines (vcf format).
	- prefix_str:	Prefix of the chromosome number (e.g. "Chr", "chr", etc.)
	Output:	
	- snps[chr_pos] = [chr,pos,GT]
	- last_snps[chr]  = position-last-snp-per-chr 
	"""
	hues.info("Reading Parental SNPs file")
	
	CheckInput(input_vcf)
	snps = OrderedDict()
	last_snps = OrderedDict()

	if prefix_chr == "":
		raise ValueError("Invalid prefix chromosome")

	with open(input_vcf, 'r') as input_file:
		for line in input_file:
			if line.startswith("#"):
				continue
			lines = line.strip("\n").split("\t")
			
			chr = lines[0]
			pos = int(lines[1])

			# correct chr if needed
			if not chr.startswith(prefix_chr):
				chr = prefix_chr + chr
			key = chr + "_" + str(pos)

			# snps contains all SNPs
			snps[key] = [chr, pos]

			# last_snps contains the last SNPs per chromosome
			if chr in last_snps.keys():
				if pos > last_snps[chr]:
					last_snps[chr] = pos
			else:
				last_snps[chr] = pos

			info = lines[9].split(":")
			geno = info[0]
			snps[key].append(geno)

	hues.log(str(len(snps)) + " Parental SNP markers loaded!")
	return snps, last_snps


#==============================================================================

def ReadOffspringVCF(input_vcf:str, parental_snps:OrderedDict, \
	geno_ref:str, geno_alt:str, analyze_id:str, prefix_chr:str="Chr"):
	"""
	Description: create 2 dictionnaries containing all SNPs per chromosome and 
	the last SNPs per chromosome
	Input: 
	- input_vcf:		SNP markers between offspring lines (vcf format).
	- parental_snps:	first element in output of ReadParentalVCF()
	- geno_ref: 			reference genotype 
	- geno_alt: 			alternative genotype
	- prefix_chr: 		Prefix of the chromosome number (e.g. "Chr", "chr", etc.)
	Output:
	- info_snps[chr_pos] = [chr,pos,GT,ADref,ADalt,genotype]
	- info_last_snps[chr] = position-last-snp-per-chr 
	"""
	hues.info("Reading Offspring SNPs file")

	if prefix_chr == "":
		raise ValueError("Invalid prefix chromosome")
	
	# Prepare log directory and remove files if already exist
	log_path = os.path.dirname(__file__) + "/log/" + analyze_id + "/"

	RemoveFile(log_path + "ReadOffsrpingVCF_new_snps.log")
	RemoveFile(log_path + "ReadOffsrpingVCF_weird_snps.log")
	CheckInput(input_vcf)

	info_snps = OrderedDict()
	info_last_snps = OrderedDict()
	total_snp_num = 0
	weird_snp = 0
	new_snp = 0
	count_dp0 = 0

	with open(input_vcf, 'r') as input_file:
		for line in input_file:
			if line.startswith("#"):
				# ignore all comments and header from vcf file
				continue

			if "DP=0" in line:
				count_dp0 +=1
				# ignore all snps with depth coverage = 0

			total_snp_num += 1
			lines = line.strip("\n").split("\t")
			
			chr = lines[0]
			pos = int(lines[1])

			# correct chr if needed
			if not chr.startswith(prefix_chr):
				chr = prefix_chr + chr
			
			#set key
			key = chr + "_" + str(pos)
			
			#get GT value from the info field
			geno = lines[9].split(":")[0]
			
			# get AD ref and AD alt from info field 
			ref_supp = lines[9].split(":")[1].split(",")[0]
			alt_supp = lines[9].split(":")[1].split(",")[1]

			if key in parental_snps.keys():
				if geno == "0/0":
					info_snps[key] = [chr, pos, geno, ref_supp, alt_supp, geno_ref]
				elif geno == "0/1":
					info_snps[key] = [chr, pos, geno, ref_supp, alt_supp, str(geno_ref + "/" + geno_alt)]
				elif geno == "1/1":
					info_snps[key] = [chr, pos, geno, ref_supp, alt_supp, geno_alt]
				else:
					weird_snp += 1
					try:
						os.makedirs(log_path)
					except FileExistsError: ## pass if directory already exists
						pass
					with open(log_path + "ReadOffsrpingVCF_weird_snps.log","a+") as logfile:
						logfile.write(line) 		
			else:
				new_snp +=1
				try:
					os.makedirs(log_path)
				except FileExistsError: ## pass if directory already exists
					pass
				with open(log_path + "ReadOffsrpingVCF_new_snps.log","a+") as logfile:
					logfile.write(line) 

			if chr in info_last_snps.keys():
				if pos > info_last_snps[chr]:
					info_last_snps[chr] = pos
			else:
				info_last_snps[chr] = pos

	hues.log(str(total_snp_num) + " Offspring genotyped SNP markers loaded!")
	hues.log(str(len(info_snps)) + " Offspring genotyped informative SNP markers kept!")
	hues.warn(str(new_snp) + " Offspring specific SNPs (not found in parent)")
	hues.warn(str(weird_snp) + " Weird genotype (1/2): heterozygous genotype composed of two different ALT alleles")
	hues.warn(str(count_dp0) + " snps with a depth coverage = 0")
	
	return info_snps, info_last_snps

#==============================================================================

def ReadEMSline(input:str, prefix_chr:str):
	"""
	Description: 
	Input:
	Output:
	"""
	hues.info("Reading Offspring SNPs file")
	CheckInput(input)

	snps_list = list()
	snps_list2 = list()
	if prefix_chr == "":
		raise ValueError("Invalid prefix chromosome")
	
	with open(input, 'r') as input_file:
		for lines in input_file:
			if lines.startswith("#"):
				# ignore all comments and header from vcf file
				continue

			line = lines.strip("\n").split("\t")		
			chr = line[0]
			pos = int(line[1])
			alt = line[3]

			# correct chr if needed
			if not chr.startswith(prefix_chr):
				chr = prefix_chr + chr

			snp_key = chr + "_" + str(pos) + "_" + alt
			snp_key2 = chr + "_" + str(pos)
			snps_list.append(snp_key)
			snps_list2.append(snp_key2)
	

	return snps_list, snps_list2

###############################################################################

def IdentifyGeno(chr, pos, alt_list, geno, snps_list_line1, snps_list_line2, 
				 geno_ref, geno_alt, geno_line1, geno_line2):
	
	if geno == "./.":
		final_geno = "NA"
	
	else :
		genotype=list()
		gt_list = geno.split("/")

		for gt in gt_list:
			if gt == "0":
				genotype.append(geno_ref)
			else :
				# define snp_key
				if gt =="1":
					snp_key = chr + "_" + str(pos) + "_" + alt_list[0]
				elif gt =="2":
					snp_key = chr + "_" + str(pos) + "_" + alt_list[1]
				elif gt =="3":
					snp_key = chr + "_" + str(pos) + "_" + alt_list[2]
				
				# define associated genotype
				if snp_key in snps_list_line1:
						genotype.append(geno_line1)
				elif snp_key in snps_list_line2:
					genotype.append(geno_line2)
				else:
					genotype.append(geno_alt)
		
		# if Homo, give genotype once only
		if genotype[0] == genotype[1]:
			final_geno = genotype[0]
		else:
			final_geno = "/".join(genotype)
	
	return final_geno

#==============================================================================

def ReadRecombinedOffspringVCF(input_vcf:str, input_ems_line1: str, input_ems_line2: str, 
	parental_snps:OrderedDict, geno_ref:str, geno_alt:str, geno_line1:str, geno_line2:str, 
	analyze_id:str, prefix_chr:str="Chr"):
	"""
	Description: create 2 dictionnaries containing all SNPs per chromosome and 
	the last SNPs per chromosome
	Input: 
	- input_vcf:		SNP markers between offspring lines (vcf format).
	- parental_snps:	first element in output of ReadParentalVCF()
	- geno_ref: 		reference genotype 
	- geno_alt: 		alternative genotype
	- geno_line1:		genotype of EMS line 1
	- geno_line2:		genotype of EMS line 2
	- prefix_chr: 		Prefix of the chromosome number (e.g. "Chr", "chr", etc.)
	Output:
	- info_snps[chr_pos] = [chr,pos,GT,ADref,ADalt,genotype]
	- info_last_snps[chr] = position-last-snp-per-chr 
	"""
	hues.info("Reading Offspring SNPs file")

	if prefix_chr == "":
		raise ValueError("Invalid prefix chromosome")
	
	# Prepare log directory and remove files if already exist
	log_path = os.path.dirname(__file__) + "/log/" + analyze_id + "/"

	RemoveFile(log_path + "ReadOffsrpingVCF_new_snps.log")
	RemoveFile(log_path + "ReadOffsrpingVCF_weird_snps.log")
	CheckInput(input_vcf)

	info_snps = OrderedDict()
	info_last_snps = OrderedDict()
	total_snp_num = 0
	new_snp = 0
	unknown_snps = 0

	snps_list_line1 = ReadEMSline(input_ems_line1, prefix_chr)
	snps_list_line2 = ReadEMSline(input_ems_line2, prefix_chr)

	with open(input_vcf, 'r') as input_file:
		for line in input_file:
			if line.startswith("#"):
				# ignore all comments and header from vcf file
				continue

			if "./." in line:
				unknown_snps += 1

			total_snp_num += 1
			lines = line.strip("\n").split("\t")
			
			chr = lines[0]
			pos = int(lines[1])

			# correct chr if needed
			if not chr.startswith(prefix_chr):
				chr = prefix_chr + chr

			key = chr + "_" + str(pos)
			geno = lines[9].split(":")[0]

			ref_supp = lines[9].split(":")[1].split(",")[0]
			alt_list = lines[9].split(":")[1].split(",")[1:]

			if key in parental_snps.keys():
				info_snps[key] = [chr, pos, geno, ref_supp, ",".join(alt_list), 
						IdentifyGeno(chr, pos, alt_list, geno, snps_list_line1, snps_list_line2, 
				   		geno_ref, geno_alt, geno_line1, geno_line2)]
			else:
				new_snp +=1
				try:
					os.makedirs(log_path)
				except FileExistsError: ## pass if directory already exists
					pass
				with open(log_path + "ReadOffsrpingVCF_new_snps.log","a+") as logfile:
					logfile.write(line) 

			if chr in info_last_snps.keys():
				if pos > info_last_snps[chr]:
					info_last_snps[chr] = pos
			else:
				info_last_snps[chr] = pos

	hues.log(str(total_snp_num) + " Offspring genotyped SNP markers loaded!")
	hues.log(str(len(info_snps)) + " Offspring genotyped informative SNP markers kept!")
	hues.warn(str(new_snp) + " Offspring specific SNPs (not found in parent)")
	hues.warn(str(unknown_snps) + " not genotyped snps")
	
	return info_snps, info_last_snps
