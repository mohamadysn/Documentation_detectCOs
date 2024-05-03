detectCOs_sliding_window.py
=============================



.. code-block:: python

    import hues
    from collections import OrderedDict
    from detectCOs_required_functions import *
..


Function : ParentalSlidingWindow
---------------------------------

**Code :**

.. code-block:: python

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

    
..

**Explanation:**

`ParentalSlidingWindow`: This function determines the number of SNPs per window along parental chromosomes. It considers the length of each chromosome, parental SNPs, and the position of the last SNPs per chromosome. The window size is specified in kilobases (kb), with a step size equal to half the window size.
- **Input:**
  - `chr_len`: A dictionary containing the length of each chromosome.
  - `parental_snps`: A dictionary containing the SNPs of the parent.
  - `last_snps_chr`: A dictionary containing the position of the last SNP of each chromosome.
  - `window_size`: The size of the window in kilobases (default: 100 kb, maximum: 1 Mb / 1000 kb). The step size is half of the window size.
- **Functionality:**
  - The function starts by verifying if the window size is valid. It must be less than or equal to 1000 kb (1 Mb) and must be an even number to avoid having a half-step when calculating the step size.
  - It then converts the window size and step size into base pairs (bp) for ease of calculations.
  - For each chromosome, it determines the number of windows by dividing the chromosome length by the step size. It rounds this number up to the nearest integer.
  - For each window on each chromosome, it determines the start and end positions of the window. If the end position exceeds the length of the chromosome, it adjusts the end position to the length of the chromosome.
  - It then iterates over each position within the window and checks if it corresponds to a SNP in the `parental_snps` dictionary. If so, it increments the SNP counter for that window.
  - Finally, it records the total number of SNPs in each window in a dictionary `snps_window` with a unique key for each window.
- **Output:**
  - `snps_window`: A dictionary containing the number of SNPs per sliding window for each chromosome. Each entry in the dictionary has a key consisting of the chromosome number and window number, and the value is a list containing the start position, end position, and total number of SNPs in that window.

Function : OffspringSlidingWindow
----------------------------------

**Code :**

.. code-block:: python

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
    
..


**Explanation:**

`OffspringSlidingWindow`: This function is similar to the previous one but processes SNPs from offspring instead of parents. It also calculates various statistics for each window, such as the number of SNPs with at least one reference or alternative allele.
- **Input:**
  - `chr_len`: A dictionary containing the length of each chromosome.
  - `offspring_snps`: A dictionary containing the SNPs of the offspring.
  - `last_snps_chr`: A dictionary containing the position of the last SNP of each chromosome.
  - `window_size`: The size of the window in kilobases (default: 100 kb). The step size is half the window size.
- **Functionality:**
  - The function starts by verifying if the window size is valid. It must be less than or equal to 1000 kb (1 Mb) and must be an even number to avoid having a half-step when calculating the step size.
  - It then converts the window size and step size into base pairs (bp) for ease of calculations.
  - For each chromosome, it determines the number of windows by dividing the chromosome length by the step size. It rounds this number up to the nearest integer.
  - For each window on each chromosome, it determines the start and end positions of the window. If the end position exceeds the length of the chromosome, it adjusts the end position to the length of the chromosome.
  - It then iterates through each position within the window and checks if it corresponds to a SNP in the `offspring_snps` dictionary. If so, it updates the window parameters based on the SNP's genotype.
  - The parameters analyzed for each window include:
     - `ADref`: The sum of the reference allele (HomoA).
     - `ADalt`: The sum of the alternative allele (HomoB).
     - `DP`: The total depth.
     - `nbSNP_Aref`: The number of SNPs with at least one reference allele.
     - `nbSNP_Aalt`: The number of SNPs with at least one alternative allele.
     - `TOTsnps-window`: The total number of SNPs in the window.
  - Finally, it records these parameters for each window in a dictionary `snps_window` with a unique key for each window.
- **Output:**
  - `snps_window`: A dictionary containing the sum of each parameter for all SNPs in each sliding window for each chromosome of the offspring. Each entry in the dictionary has a key consisting of the chromosome number and window number, and the value is a list containing the start position, end position, `ADref`, `ADalt`, `DP`, `nbSNP_Aref`, `nbSNP_Aalt`, and `TOTsnps-window`.


Function : NormalizeOffspringSlidingWindow
--------------------------------------------

**Code :**

.. code-block:: python

    def NormalizeOffspringSlidingWindow(parental_snps_window:OrderedDict, \
            offspring_snps_window:OrderedDict, geno_ref:str, geno_alt:str, \
            min_snp_num:int=16, min_reads_num:int=10, ratio_min_homo:float=0.75):  # modify MY
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
                geno_window[key_window] = GetGenoWindow(cur_window, geno_ref, geno_alt, ratio_min_homo)

                # save the last window of each chromosome in a dictionnary nb_window_chr[chr]=num_last_window
                chr = key_window.split("_")[0]
                win_id = int(key_window.split("_")[1])
                if chr in nb_window_chr:
                    if win_id > nb_window_chr[chr]:
                        nb_window_chr[chr] = win_id
                else:
                    nb_window_chr[chr] = win_id

        return snps_window, geno_window, nb_window_chr
    
..


**Explanation:**

`NormalizeOffspringSlidingWindow`: This function normalizes the SNP data of the offspring based on the parent's data. It also determines the predominant genotype of each window based on the `ADref/DP` and `ADalt/DP` ratios, using a minimum frequency to consider a SNP as homozygous.
- **Input:**
  - `parental_snps_window`: The output from the `ParentalSlidingWindow()` function.
  - `offspring_snps_window`: The output from the `OffspringSlidingWindow()` function.
  - `geno_ref`: The reference genotype of the parent.
  - `geno_alt`: The alternative genotype of the parent.
  - `min_snp_num`: The minimum number of SNPs in the window required to calculate the ratio and determine the genotype (default: 16).
  - `min_reads_num`: The minimum number of reads in the window required to calculate the ratio and determine the genotype (default: 10).
  - `ratio_min_homo`: The minimum frequency of `AD/DP` to be considered homozygous (floating between 0 and 1, default: 0.75).
- **Functionality:**
  - The function starts by initializing three ordered dictionaries to store the data:
    - `snps_window`: Contains information on each SNP window of the offspring.
    - `geno_window`: Contains information on the genotype of each SNP window of the offspring.
    - `nb_window_chr`: Contains the number of windows for each chromosome.
  - It then verifies if `ratio_min_homo` is valid. It must be between 0 and 1, exclusive.
  - For each SNP window of the offspring, it extracts information such as start, end, `ADref`, `ADalt`, and `DP`.
  - If the offspring SNP window is not present in the parental windows, it raises an error because all offspring SNP windows must be present in the parental windows.
  - Otherwise, it retrieves the total number of SNPs in the corresponding parental window.
  - If the number of reads (`DP`) is greater than or equal to `min_reads_num` and the total number of SNPs in the parental window is greater than or equal to `min_snp_num`, it updates the information for the offspring SNP window in `snps_window`.
  - It then determines the genotype of each SNP window of the offspring by calling the `GetGenoWindow` function.
  - Finally, it records the number of the last window of each chromosome in `nb_window_chr`.
- **Output:**
  Three elements are returned:
  - `snps_window`: A dictionary containing information on each SNP window of the offspring.
  - `geno_window`: A dictionary containing information on the genotype of each SNP window of the offspring.
  - `nb_window_chr`: A dictionary containing the number of windows for each chromosome.



Function : SmoothNormalizedOsffspringSlidingWindow
---------------------------------------------------

**Code :**

.. code-block:: python

    def SmoothNormalizedOsffspringSlidingWindow(offspring_snps_window:OrderedDict,\
            nb_windows_chr:OrderedDict, centromere:OrderedDict, geno_ref:str,\
            geno_alt:str, ratio_min_homo:float=0.75):   # Modify MY
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
                
                geno_window[key_window] = GetGenoWindow(cur_window, geno_ref, geno_alt, ratio_min_homo)
                
                check_smoothed[key_window] = [pos_start_cur_win, pos_stop_cur_win, \
                    start_win_smooth, stop_win_smooth, win]

        return snps_window, geno_window, check_smoothed 
    
..


**Explanation:**

`SmoothNormalizedOffspringSlidingWindow`: This function smooths the normalized SNP data of the offspring by considering neighboring windows. It accounts for centromeric regions and assigns a "NA" genotype to windows that overlap or are within these regions.
- **Input:**
  - `offspring_snps_window`: The first element from the output of the `NormalizeOffspringSlidingWindow()` function.
  - `nb_windows_chr`: The number of windows per chromosome, determined with the `NormalizeOffspringSlidingWindow()` function.
  - `centromere`: An ordered dictionary with the boundaries of the centromeric region for each chromosome.
  - `geno_ref`: The genotype of the reference parent.
  - `geno_alt`: The genotype of the alternative parent.
  - `ratio_min_homo`: The minimum `AD/DP` ratio to consider SNPs as homozygous.

- **Functionality:**
  - The function begins by initializing three ordered dictionaries to store the data:
    - `snps_window`: Contains information on each SNP window of the offspring.
    - `geno_window`: Contains information on the genotype of each SNP window of the offspring.
    - `check_smoothed`: Used to verify if the smoothing takes into account the correct windows.
  - For each chromosome and each SNP window of the offspring:
    - It retrieves the boundaries of the centromeric region for the current chromosome.
    - It defines the smoothing windows as the current window minus one and the current window plus one, ensuring they do not exceed the bounds.
    - It calculates the sum of `ADref`, `ADalt`, and `DP` for all smoothing windows that do not intersect or are inside the centromeric region.
    - Updates the information in the corresponding dictionaries.
    - Also stores information about the smoothing windows for verification.

- **Output:**
  Three elements are returned:
  - `snps_window`: A dictionary containing information on each SNP window of the offspring.
  - `geno_window`: A dictionary containing information on the genotype of each SNP window of the offspring.
  - `check_smoothed`: A dictionary to verify if the smoothing considers the correct windows.