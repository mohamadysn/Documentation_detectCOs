import os, sys
import pprint
import yaml


#==============================================================================
# Program		:	Launcher of for detectCOs EMS
# Author		:	Mohamad Yassine [mohamad.yassine@inrae.fr]
# Date			:	25.02.2024
# Last update	:	26.04.2024
# Version		:	2.0
#==============================================================================


from detectCOs_read_files_EMS import *
from detectCOs_sliding_window_EMS import *
from detectCOs_identifyCOs_EMS import *
from VisualiseGenotypesChromosomeEMS import *
from collections import OrderedDict

# Usage: python3 path/to/folder/launcher_detectCOs_v1_backcross.py config_detectCOs.yaml

###############################################################################

print("# Load config file and create output directory")
print("-------------------------------------")

if len(sys.argv) == 1:
    raise ImportError("Please give the config file (.yaml) for detectCOs")

config_file = sys.argv[1]
with open(config_file, "r") as cf:
    config = yaml.load(cf, Loader=yaml.FullLoader)

# Set current working directory at the root of polyrec project
os.chdir(config['path_to_polyrec_project'])

# Set variable to save all results
outdir = config['path_to_polyrec_project'] + config['output_dir_in_polyrec'] + \
    config['analyze_id'] + "/"
# And create output directory if not exist
if not os.path.exists(outdir):
    os.makedirs(outdir)
    print("Create output directory")


###############################################################################

print("\n# STEP 1. Load chromosome length")
print("-------------------------------------")

chr_length = ReadChrLen(input_chr_len=config['chr_len'], 
                        prefix_chr=config['prefix_chr'])

print("\nchr_length[chr] = lenght of chromosome")
pprint.pprint(list(chr_length.items()))


###############################################################################

print("\n# STEP 2. Load position of the centromeric region")
print("-------------------------------------")

centromere = ReadCentroReg(input_centro_reg=config['centromere_reg'],
                           prefix_chr=config['prefix_chr'])

print("\ncentromere[chr] = [left_border, right_border]")
pprint.pprint(list(centromere.items()))


###############################################################################

print("\n# STEP 3. Load parental snps from vcf")
print("-------------------------------------")

if os.path.exists(outdir + "parental_snps.txt") and \
    os.path.exists(outdir + "parental_snps_last_pos_chr.txt"):
    print("Load data...")
    parental_snps = load_file_in_dict(outdir + "parental_snps.txt")
    parental_last_snp_chr = load_file_in_dict(outdir + "parental_snps_last_pos_chr.txt")

else:
    print("Run ReadParentalVCF()...")
    parental_snps, parental_last_snp_chr = ReadParentalVCF(
        input_vcf=config['parental_vcf'], 
        prefix_chr=config['prefix_chr']
        )

print("\nparental_snps[chr_pos] = [chr, pos, GT]")
pprint.pprint(list(parental_snps.items())[:5])
print("\nparental_last_snp_chr:")
pprint.pprint(list(parental_last_snp_chr.items()))

#------------------------------------------------------------------------------
print ("\nSave outputs of ReadParentalVCF()")
export_dict_in_file(my_dict = parental_snps, 
                   output_file=outdir + "parental_snps.txt",
                   header="chr_pos\tchr\tpos\tGT",
                   overwrite=False)

        
###############################################################################

print("\n# STEP 4. Load offspring snps from vcf")
print("-------------------------------------")

if os.path.exists(outdir + "offspring_snps.txt") and \
    os.path.exists(outdir + "offspring_snps_last_pos_chr.txt"):
    print("Load data...")
    offspring_snps = load_file_in_dict(outdir + "offspring_snps.txt")
    offspring_last_snp_chr = load_file_in_dict(outdir + "offspring_snps_last_pos_chr.txt")

else:
    print("Run ReadOffspringVCF()...")
    offspring_snps, offspring_last_snp_chr = ReadOffspringVCF(
        input_vcf=config['sample_vcf'], 
        parental_snps=parental_snps,
        prefix_chr=config['prefix_chr'],
        geno_ref=config['genotype_ref'], 
        geno_alt=config['genotype_alt'],
        analyze_id=config['analyze_id']
        )

print("\noffspring_snps[chr_pos] = [chr, pos,GT, ADref,ADalt, genotype]")
pprint.pprint(list(offspring_snps.items())[:5])


#------------------------------------------------------------------------------
print ("\nSave outputs of ReadOffspringVCF()")
export_dict_in_file(my_dict = offspring_snps, 
                   output_file=outdir + "offspring_snps.txt",
                   header="chr\tpos\tGT\tADref\tADalt\tgenotype",
                   overwrite=False)

last_snps_chr = parental_last_snp_chr


export_dict_in_file(my_dict = last_snps_chr, 
                   output_file=outdir + "last_snps_chr.txt",
                   header="chr\tpos_last_SNP",
                   overwrite=False)

###############################################################################

print("\n# STEP 5. Summary info about parental SNPs by sliding window")
print("-------------------------------------")

print("Run ParentalSlidingWindow()...")
parental_snps_windowsnps = ParentalSlidingWindowSNPs(
    parental_snps=parental_snps,
    snps_per_window=config['snps_per_window']
    )


print("\nparental_snps_windowsnps[chr_window] = [start, stop, nb_snps]")
pprint.pprint(list(parental_snps_windowsnps.items())[5:9])

#------------------------------------------------------------------------------
print ("\nSave output of ParentalSlidingWindow()")
export_dict_in_file(my_dict = parental_snps_windowsnps, 
                    output_file=outdir + "parental_snps_windowsnps_"+ str(config['snps_per_window']) +"_snps.txt",
                    header="chr_window\tstart\tstop\tnb_snps",
                    overwrite=False)



###############################################################################

print("\n# STEP 6. Summary info about offspring SNPs by sliding window")
print("-------------------------------------")




print("Run ParentalSlidingWindow()...")
offspring_snps_windowsnps = OffspringSlidingWindowBySNPs(
    offspring_snps = offspring_snps , 
    snps_per_window = config['snps_per_window']) 


print("\offspring_snps_windowsnps[chr_window] = [start, stop, nb_snps]")
pprint.pprint(list(offspring_snps_windowsnps.items())[5:9])

#------------------------------------------------------------------------------
print ("\nSave output of ParentalSlidingWindow()")
export_dict_in_file(my_dict = offspring_snps_windowsnps, 
                    output_file=outdir + "offspring_snps_windowsnps_"+ str(config['snps_per_window'])+"_snps.txt",
                    header="chr_window\tstart\tstop\tADref\tADalt\tDP\tnbSNP_Aref\tnbSNP_Aalt\ttot_snps", #\tADref/DP\tADalt/DP",
                    overwrite=False)

#######################################################################################################



offspring_snps_window_normalized, offspring_genotype_window_normalized, \
    nb_windows_chr = NormalizeOffspringSlidingWindowsnps(
        parental_snps_window=parental_snps_windowsnps,
        offspring_snps_window=offspring_snps_windowsnps,
        geno_ref=config['EMS_ref'],
        geno_alt=config['EMS_alt'],
        min_snp_num=config['min_snp_num'],
        min_reads_num=config['min_reads_num'],
        ratio_min_homo=config['ratio_min_EMS'],
        depth_division_th = config['depth_division_th']  # modify MY
        )

print("\noffspring_snps_window_normalized: ")
print("[chr_window] = [start, stop, ADref, ADalt, DP]")
pprint.pprint(list(offspring_snps_window_normalized.items())[5:9])

print("\noffspring_genotype_window_normalized: ")
print("[chr_window] = [start, stop, ADref/DP, ADalt/DP, genotype]")
pprint.pprint(list(offspring_genotype_window_normalized.items())[5:9], width = 120, compact = True)

print("\nnb_windows_chr: ")
pprint.pprint(list(nb_windows_chr.items()))

#------------------------------------------------------------------------------    
    
print ("\nSave outputs of NormalizeOffspringSlidingWindow()")
export_dict_in_file(my_dict = offspring_snps_window_normalized, 
                   output_file=outdir + "offspring_snps_window_normalized_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
                   header="chr_window\tstart\tstop\tADref\tADalt\tDP",
                   overwrite=False)

export_dict_in_file(my_dict = offspring_genotype_window_normalized, 
                   output_file=outdir + "offspring_genotype_window_normalized_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
                   header="chr_window\tstart\tstop\tADref/DP\tADalt/DP\tgenotype",
                   overwrite=False)

export_dict_in_file(my_dict = nb_windows_chr, 
                   output_file=outdir + "nb_window_chr_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
                   header="chr\tnb_window",
                   overwrite=False)






###############################################################################
"""
print("\n# STEP 8. Smooth offspring sliding window")
print("-------------------------------------")


offspring_snps_window_smoothed, offspring_genotype_window_smoothed, \
    check_smoothing = SmoothNormalizedOsffspringSlidingWindowsnps(\
        offspring_snps_window=offspring_snps_window_normalized,
        nb_windows_chr=nb_windows_chr,
        geno_ref=config['EMS_ref'],
        geno_alt=config['EMS_alt'],
        ratio_min_homo=config['ratio_min_EMS'],
        depth_division_th = config['depth_division_th']  # modify MY
        )


print("\noffspring_snps_window_smoothed: ")
print("[chr_window] = [start, stop, ADref, ADalt, DP]")
pprint.pprint(list(offspring_snps_window_smoothed.items())[5:9])

print("\noffspring_genotype_window_normalized_smoothed")
print("[chr_window] = [start, stop, ADref/DP, ADalt/DP, probHomoRef, probHetero, probHomoAlt, genotype]")
pprint.pprint(list(offspring_genotype_window_smoothed.items())[5:9], width = 120, compact = True)

print("\ncheck_smoothing[chr_window] = [start_window,stop_window,start_smooth,stop_smooth,window_smoothed")
pprint.pprint(list(check_smoothing.items())[5:9])

#------------------------------------------------------------------------------

print("\nSave outputs of SmoothNormalizedOsffspringSlidingWindow()")

export_dict_in_file(my_dict = offspring_snps_window_smoothed, 
                   output_file=outdir + "offspring_snps_window_normalized_smoothed_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
                   header="chr_window\tstart\tstop\tADref\tADalt\tDP",
                   overwrite=False)

export_dict_in_file(my_dict = offspring_genotype_window_smoothed, 
                   output_file=outdir + "offspring_genotype_window_normalized_smoothed_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
                   header="chr_window\tstart\tstop\tADref/DP\tADalt/DP\tgenotype",
                   overwrite=False)

export_dict_in_file(my_dict = check_smoothing, 
                   output_file=outdir + "check_smoothing_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
                   header="chr_window\tstart_window\tstop_window\tstart_smooth\tstop_smooth\twindow_smoothed",
                   overwrite=False)

"""

##############################################   # modify MY

print("\n# STEP 8.bis. visualise Genotype ")
print("-------------------------------------")

# 'file_paths' est maintenant un dictionnaire avec les noms des répertoires comme clés et les chemins complets comme valeurs.


# Exemple d'utilisation
file_paths = {
    # 'directory_name': 'offspring_genotype_window_normalized_smoothed.txt'
}

file_paths[config['analyze_id']] = outdir + "offspring_genotype_window_normalized_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt"
#file_paths[config['analyze_id'] + "_" + str(config['new_window_size']) + "_kb"] = outdir + "offspring_genotype_window_normalized_smoothed_" + str(config['new_window_size']) + "_kb" + ".txt"

genotype_colors = {
    config['EMS_ref']: 'blue',
    config['EMS_alt']: 'red',
    'NA': 'gray'  
}


visualsize_genotope(file_paths, outdir, genotype_colors, snps_per_window = config['snps_per_window'], ratio_min_EMS = config['ratio_min_EMS'])


###############################################################################


print("\n# STEP 9. Identify COs")
print("-------------------------------------")


candidates_co, db_co = IdentifyCOs_snps(offspring_genotype_window_normalized, nb_windows_chr)


updated_candidates_co = OrderedDict()  # modify MY


for chr_window, details in candidates_co.items(): # modify MY
    start_win_stop_win, co_start, co_stop, pre_geno, cur_geno = details # modify MY
    diff = co_stop - co_start  # Calcul de la différence # modify MY
    updated_candidates_co[chr_window] =  [start_win_stop_win, co_start, co_stop,diff, pre_geno, cur_geno] # modify MY

print("\nCandidateCOs: ")
print("[chr_window] = [start_win:stop_win, co_start, co_stop, pre_geno, cur_geno]")
pprint.pprint(list(updated_candidates_co.items()))  # modify MY

export_dict_in_file(my_dict = updated_candidates_co,   # modify MY
                   output_file=outdir + "candidateCO_"+ str(config['snps_per_window'])+"_snps_ratio_"+str(config['ratio_min_EMS']) +".txt",
                   header="chr_mean_win\tstart_win:stop_win\tstart_co\tstop_co\t(stop_co-start_co)\tpre_geno\tcur_geno",
                   overwrite=True)


############################################################################### # modify MY

print("\n# STEP 5.bis. Summary info about parental SNPs by sliding window" + str(config['new_snp_per_window']) )
print("-------------------------------------")


if os.path.exists(outdir + "parental_snps_window_" + str(config['new_snp_per_window']) + ".txt"):
    print("Load data...")
    parental_snps_window_new = load_file_in_dict(outdir + "parental_snps_window_" + str(config['new_snp_per_window']) + ".txt")

else:
    print("Run ParentalSlidingWindow()...")
    parental_snps_window_new = ParentalSlidingWindowSNPs(
        parental_snps=parental_snps,
        snps_per_window=config['new_snp_per_window']
    )

print("\nparental_snps_window_" + str(config['new_snp_per_window']) + "[chr_window] = [start, stop, nb_snps]")
pprint.pprint(list(parental_snps_window_new.items())[5:10])

#------------------------------------------------------------------------------
print ("\nSave output of ParentalSlidingWindow()")
export_dict_in_file(my_dict = parental_snps_window_new, 
                   output_file=outdir + "parental_snps_window_" + str(config['new_snp_per_window']) + ".txt",
                   header="chr_window\tstart\tstop\tnb_snps",
                   overwrite=False)


###############################################################################  # modify MY

print("\n# STEP 6.bis. Summary info about offspring SNPs by sliding window" + str(config['new_snp_per_window']) )
print("-------------------------------------")

if os.path.exists(outdir + "offspring_snps_window_" + str(config['new_snp_per_window']) +  ".txt"):
    print("Load data...")
    offspring_snps_window_new = load_file_in_dict(outdir + "offspring_snps_window_" + str(config['new_snp_per_window']) + ".txt")

else:
    print("Run OffspringSlidingWindow()...")
    offspring_snps_window_new = OffspringSlidingWindowBySNPs(
        offspring_snps=offspring_snps,
        snps_per_window=config['new_snp_per_window']
    )


print("\noffspring_snps_window_" + str(config['new_snp_per_window']) +  "[chr_window] = [start,stop,ADref,ADalt,DP,nbSNP_Aref,nbSNP_Aalt,TOTsnps-window]")
pprint.pprint(list(offspring_snps_window_new.items())[5:10])

#------------------------------------------------------------------------------
print ("\nSave output of OffspringSlidingWindow()")
export_dict_in_file(my_dict = offspring_snps_window_new, 
                   output_file=outdir + "offspring_snps_window_" + str(config['new_snp_per_window']) + ".txt",
                   header="chr_window\tstart\tstop\tADref\tADalt\tDP\tnbSNP_Aref\tnbSNP_Aalt\ttot_snps",
                   overwrite=False)



###############################################################################   # modify MY

print("\n# STEP 7.bis. Identify the genotype of the offspring normalized sliding window" + str(config['new_snp_per_window']) )
print("-------------------------------------")

if os.path.exists(outdir + "offspring_snps_window_normalized_" + str(config['new_snp_per_window']) +  ".txt") and \
    os.path.exists(outdir + "offspring_genotype_window_normalized_" + str(config['new_snp_per_window']) + ".txt") and \
        os.path.exists(outdir + "nb_window_chr_" + str(config['new_snp_per_window']) +  ".txt"):
    print("Load data...")
    offspring_snps_window_normalized_new = load_file_in_dict(outdir + "offspring_snps_window_normalized_" + str(config['new_snp_per_window']) +  ".txt")
    offspring_genotype_window_normalized_new = load_file_in_dict(outdir + "offspring_genotype_window_normalized_" + str(config['new_snp_per_window']) + ".txt")
    nb_windows_chr_new = load_file_in_dict(outdir + "nb_window_chr_" + str(config['new_snp_per_window']) +  ".txt")

else:
    offspring_snps_window_normalized_new, offspring_genotype_window_normalized_new, \
        nb_windows_chr_new = NormalizeOffspringSlidingWindowsnps(
            parental_snps_window=parental_snps_window_new,
            offspring_snps_window=offspring_snps_window_new,
            geno_ref=config['EMS_ref'],
            geno_alt=config['EMS_alt'],
            min_snp_num=config['min_snp_num'],
            min_reads_num=config['min_reads_num'],
            ratio_min_homo=config['ratio_min_EMS'],
            depth_division_th = config['depth_division_th']   # modify MY
            )




print("\noffspring_snps_window_normalized_" + str(config['new_snp_per_window']))
print("[chr_window] = [start, stop, ADref, ADalt, DP]")
pprint.pprint(list(offspring_snps_window_normalized_new.items())[5:10])

print("\noffspring_genotype_window_normalized_" + str(config['new_snp_per_window']) )
print("[chr_window] = [start, stop, ADref/DP, ADalt/DP, probHomoRef, probHetero, probHomoAlt, genotype]")
pprint.pprint(list(offspring_genotype_window_normalized_new.items())[5:10], width = 120, compact = True)

print("\nnb_windows_chr_" + str(config['new_snp_per_window']))
pprint.pprint(list(nb_windows_chr_new.items()))

#------------------------------------------------------------------------------    
    
print ("\nSave outputs of NormalizeOffspringSlidingWindow() " + str(config['new_snp_per_window']) )
export_dict_in_file(my_dict = offspring_snps_window_normalized_new, 
                   output_file=outdir + "offspring_snps_window_normalized_" + str(config['new_snp_per_window']) + ".txt",
                   header="chr_window\tstart\tstop\tADref\tADalt\tDP",
                   overwrite=False)

export_dict_in_file(my_dict = offspring_genotype_window_normalized_new, 
                   output_file=outdir + "offspring_genotype_window_normalized_" + str(config['new_snp_per_window'])  + ".txt",
                   header="chr_window\tstart\tstop\tADref/DP\tADalt/DP\tprobHomoRef\tprobHetero\tprobHomoAlt\tgenotype",
                   overwrite=False)

export_dict_in_file(my_dict = nb_windows_chr_new, 
                   output_file=outdir + "nb_window_chr_" + str(config['new_snp_per_window'])+ ".txt",
                   header="chr\tnb_window",
                   overwrite=False)



###############################################################################  # modify MY
"""
print("\n# STEP 8. Smooth offspring sliding window " + str(config['new_snp_per_window']) )
print("-------------------------------------")

if os.path.exists(outdir + "offspring_snps_window_normalized_smoothed_" + str(config['new_snp_per_window'])  + ".txt") and \
    os.path.exists(outdir + "offspring_genotype_window_normalized_smoothed_" + str(config['new_snp_per_window'])  + ".txt") and \
        os.path.exists(outdir + "offspring_nb_window_normalized_" + str(config['new_snp_per_window']) + ".txt"):
    print("Load data...")
    offspring_snps_window_smoothed_new = load_file_in_dict(outdir + "offspring_snps_window_normalized_smoothed_" + str(config['new_snp_per_window'])  + ".txt")
    offspring_genotype_window_smoothed_new = load_file_in_dict(outdir + "offspring_genotype_window_normalized_smoothed_" + str(config['new_snp_per_window']) + ".txt")
    check_smoothing_new = load_file_in_dict(outdir + "check_smoothing_" + str(config['new_snp_per_window'])  + ".txt")

else:
    offspring_snps_window_smoothed_new, offspring_genotype_window_smoothed_new, \
        check_smoothing_new = SmoothNormalizedOsffspringSlidingWindowsnps(\
            offspring_snps_window=offspring_snps_window_normalized_new,
            nb_windows_chr=nb_windows_chr_new,
            geno_ref=config['EMS_ref'],
            geno_alt=config['EMS_alt'],
            ratio_min_homo=config['ratio_min_EMS'],
            depth_division_th = config['depth_division_th']  # modify MY
            )



print("\noffspring_snps_window_smoothed_" + str(config['new_snp_per_window'])  + " : ")
print("[chr_window] = [start, stop, ADref, ADalt, DP]")
pprint.pprint(list(offspring_snps_window_smoothed_new.items())[5:10])

print("\noffspring_genotype_window_normalized_smoothed_" + str(config['new_snp_per_window']) )
print("[chr_window] = [start, stop, ADref/DP, ADalt/DP, probHomoRef, probHetero, probHomoAlt, genotype]")
pprint.pprint(list(offspring_genotype_window_smoothed_new.items())[5:10], width = 120, compact = True)

print("\ncheck_smoothing_" + str(config['new_snp_per_window'])  + "[chr_window] = [start_window,stop_window,start_smooth,stop_smooth,window_smoothed")
pprint.pprint(list(check_smoothing_new.items())[5:10])

#------------------------------------------------------------------------------

print("\nSave outputs of SmoothNormalizedOsffspringSlidingWindow()")

export_dict_in_file(my_dict = offspring_snps_window_smoothed_new, 
                   output_file=outdir + "offspring_snps_window_normalized_smoothed_" + str(config['new_snp_per_window'])  + ".txt",
                   header="chr_window\tstart\tstop\tADref\tADalt\tDP",
                   overwrite=False)

export_dict_in_file(my_dict = offspring_genotype_window_smoothed_new, 
                   output_file=outdir + "offspring_genotype_window_normalized_smoothed_" + str(config['new_snp_per_window'])  + ".txt",
                   header="chr_window\tstart\tstop\tADref/DP\tADalt/DP\tprobHomoRef\tprobHetero\tprobHomoAlt\tgenotype",
                   overwrite=False)

export_dict_in_file(my_dict = check_smoothing_new, 
                   output_file=outdir + "check_smoothing_" + str(config['new_snp_per_window']) + ".txt",
                   header="chr_window\tstart_window\tstop_window\tstart_smooth\tstop_smooth\twindow_smoothed",
                   overwrite=False)

"""


print("\n# STEP 9.bis. Precise COs")
print("-------------------------------------")


preciseCOs = OrderedDict()
db_co_new = OrderedDict()

for i in range (len(list(updated_candidates_co.keys()))):
    start_co=updated_candidates_co[list(updated_candidates_co.keys())[i]][1] 
    end_co=updated_candidates_co[list(updated_candidates_co.keys())[i]][2] 
    chr = list(updated_candidates_co.keys())[i].split("_")[0]
    before =OrderedDict()
    after = OrderedDict()
    
    for key, snp_info in offspring_genotype_window_normalized.items():
        candidate = key
        chromosome = candidate.split("_")[0]
        index = candidate.split("_")[1]
        start,	stop, ratio_ADref, ratio_ADalt, genotype = snp_info

        
        if start < start_co and chromosome == chr and genotype != 'NA' :
            before[candidate] = snp_info
        if stop > end_co and chromosome == chr and genotype != 'NA' :
            after[candidate] = snp_info
    before = dict(list(before.items())[-3:])
    after = dict(list(after.items())[:3])


    start_co_2=before[list(before.keys())[0]][0]
    if len(list(after.keys()))>2:
        end_co_2= after[list(after.keys())[2]][1]
    elif len(list(after.keys())) == 0:
        end_co_2 = end_co
    elif len(list(after.keys())) == 2:
        end_co_2 = after[list(after.keys())[1]][1]
    else:
        end_co_2= after[list(after.keys())[0]][1]

    new_offspring_genotype_window = OrderedDict()
    for cle, valeur in offspring_genotype_window_normalized_new.items():
        if valeur[0] >= start_co_2 + 1 and valeur[1] <= end_co_2 and cle.split("_")[0] == chr:
            new_offspring_genotype_window[cle] = valeur
            print(cle, valeur)
    
    start_window = int(list(new_offspring_genotype_window.keys())[0].split("_")[1])
    end_window = int(list(new_offspring_genotype_window.keys())[-1].split("_")[1])
    print(start_window, end_window)


    candidates_co_2, db_co_2 = PreciseCOs_snps(new_offspring_genotype_window, start_window=start_window, end_window=end_window)
    
    if candidates_co_2 :
        candidates_co_2 = candidates_co_2
    else:
        start_co=updated_candidates_co[list(updated_candidates_co.keys())[i]][1] 
        end_co=updated_candidates_co[list(updated_candidates_co.keys())[i]][2]
        pre_geno=updated_candidates_co[list(updated_candidates_co.keys())[i]][4]
        cur_geno=updated_candidates_co[list(updated_candidates_co.keys())[i]][5]
        cur_chr = list(updated_candidates_co.keys())[i].split("_")[0]
        key = cur_chr + "_" + str(round((start_window + end_window) / 2, 1))
        candidates_co_2[key] = [str(start_window) + ":" + str(end_window), start_co, end_co, pre_geno, cur_geno]
    preciseCOs.update(candidates_co_2)
    db_co_new.update(db_co_2)


updated_preciseCOs = OrderedDict() # modify MY
for chr_window, details in preciseCOs.items(): # modify MY
    start_win_stop_win, co_start, co_stop, pre_geno, cur_geno = details # modify MY
    diff = co_stop - co_start  # Calcul de la différence # modify MY
    updated_preciseCOs[chr_window] =  [start_win_stop_win, co_start, co_stop,diff, pre_geno, cur_geno] # modify MY


print("\nCandidateCOs new: ")
print("[chr_window] = [start_win:stop_win, co_start, co_stop, pre_geno, cur_geno]")
pprint.pprint(list(updated_preciseCOs.items())) # modify MY


export_dict_in_file(my_dict = updated_preciseCOs,  # modify MY
                   output_file=outdir + "preciseCOs_" + str(config['new_snp_per_window'])+"_snps_ratio_" + str(config['ratio_min_EMS']) + ".txt",
                   header="chr_mean_win\tstart_win:stop_win\tstart_co\tstop_co\t(stop_co-start_co)\tpre_geno\tcur_geno",
                   overwrite=True)