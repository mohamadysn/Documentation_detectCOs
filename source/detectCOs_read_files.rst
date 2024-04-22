detectCOs_read_files.py
========================



.. code-block:: python

    import hues,os,sys,re
    from collections import OrderedDict
    from detectCOs_required_functions import *
..

Function : CheckInput
-----------------------

**Code :**

.. code-block:: python
    
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
..

**Explanation:**

`CheckInput`: This function checks if an input file exists and is not empty.
- **Input:**
  - `input_file`: The path to the input file.
- **Functionality:**
  - Checks if the specified file exists using `os.path.exists()`.
  - If the file exists:
    - Checks if its size is greater than zero using `os.path.getsize()`.
    - If the file is not empty, displays a message indicating that the input file has been found and is not empty.
    - Otherwise, displays an error message indicating that the input file is empty and exits the program using `sys.exit()`.
  - If the file does not exist, displays an error message indicating that the input file has not been found and exits the program using `sys.exit()`.
- **Output:**
  - No explicit output, but the function displays messages to the console based on the result of the verification. If the input file does not exist or is empty, the program terminates with an exit code of 1.

Function : RemoveFile
-----------------------

**Code :**

.. code-block:: python

    def RemoveFile(output_file:str): #old name CheckOuput
        """
        Description:	Remove output file if it exist
        Input:			path-to-file/output_file
        """
        if os.path.exists(output_file):
            hues.warn("Found output file:\t" + output_file)
            hues.warn("Remove\t" + output_file + "\t...")
            os.remove(output_file)
..

**Explanation:**

`RemoveFile`: This function deletes an output file if it exists.
- **Input:**
  - `output_file`: The path of the output file to be deleted.
- **Functionality:**
  - Checks if the specified file exists using `os.path.exists()`.
  - If the file exists, displays a warning message indicating that the output file has been found.
  - Then deletes the file using `os.remove(output_file)`.
- **Output:**
  - No explicit output, but the function displays messages to the console based on the result of the file deletion.

Function : load_file_in_dict
-----------------------------

**Code :**

.. code-block:: python

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
    
..

**Explanation:**

`load_file_in_dict`: This function converts a tab-delimited file into a dictionary.
- **Input:**
  - `input_file`: The path of the input file to be loaded.
  - `prefix_chr`: The prefix used to filter lines from the input file. It defaults to `"Chr"`.
- **Functionality:**
  - Opens the specified file using the syntax `with open(input_file) as input`.
  - Iterates through each line in the file.
  - Splits each line into elements using `line.strip().split("\t")` to separate values by tabs.
  - If the first element of the line starts with the specified prefix (`prefix_chr`), then the following values in the line are added to a list.
  - Each element in the list is converted into an integer or a float if it matches integer or floating-point number patterns, respectively.
  - If the resulting list contains only one element, that element is added to the dictionary with the key being the first element of the line. Otherwise, the entire list is added to the dictionary.
  - Once all lines are processed, the resulting dictionary is returned.
- **Output:**
  - A dictionary where the keys are the values of the first column (after filtering with `prefix_chr`) and the values are either lists of the values from the other columns or directly the value if only one column is present after filtering.

Function : export_dict_in_file
-------------------------------

**Code :**

.. code-block:: python

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
    
..

**Explanation:**

`export_dict_in_file`: This function saves a dictionary into a file.
- **Input:**
  - `my_dict`: The dictionary to save.
  - `output_file`: The path of the output file where the dictionary will be recorded.
  - `header`: A string representing the header of the output file. Each column of the header must be separated by a tabulation (`"\t"`). If `header` is an empty string (`""`), no header will be added to the output file.
  - `overwrite`: A boolean indicating whether to overwrite the output file if it already exists. By default, it is set to `False`, which means the file will not be overwritten if it already exists.
- **Functionality:**
  - Checks if the output file already exists and if the `overwrite` option is enabled. If the file exists and `overwrite` is `False`, the function displays a message and exits.
  - If the output file exists and `overwrite` is `True`, the function calls `RemoveFile(output_file)` to delete the existing file.
  - Opens the output file in write mode (`"w"`).
  - If the header is not an empty string, writes the header into the file followed by a newline.
  - Iterates through each item in the dictionary. For each item, writes the key followed by its values in one line of the file, separated by tabs. If the value is a list, each element of the list is added to the line separated by a tab.
  - Once all items are written, the function displays a message indicating that the output file has been created.
- **Output:**
  - A text file containing the data from the dictionary, with the option of a header if specified.


Function : ReadChrLen
-----------------------

**Code :**


.. code-block:: python

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
    
..

**Explanation:**

`ReadChrLen`: This function reads a file containing the length of each chromosome and returns a dictionary.
- **Input:**
  - `input_chr_len`: The path to the tab-delimited file containing the lengths of the chromosomes.
  - `prefix_chr`: The prefix for the chromosome number, such as `"Chr"`, `"chr"`, etc. This allows concatenating the prefix with the chromosome number in the output dictionary. By default, the prefix is `"Chr"`.
- **Functionality:**
  - The function begins by checking if the input file exists and is not empty by calling `CheckInput(input_chr_len)`.
  - Then, a dictionary `chr_len` is created to store the lengths of the chromosomes.
  - The function opens the input file and iterates through each line. For each line, it separates the elements by tabs (`\t`) and checks if the chromosome name starts with the specified prefix. If so, it adds an entry to the `chr_len` dictionary where the key is the chromosome name (with or without the prefix, as applicable) and the value is its length.
  - After processing all the lines, the function displays a message indicating the number of chromosome lengths loaded.
  - Finally, it returns the `chr_len` dictionary containing the lengths of each chromosome.
- **Output:**
  - A dictionary where the keys are the names of the chromosomes and the values are their respective lengths.
        

Function : ReadCentroReg
-------------------------

**Code :**


.. code-block:: python

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
    
..

**Explanation:**

`ReadCentroReg`: This function reads a file containing the coordinates of the pericentromeric region for each chromosome and returns a dictionary.
- **Input:**
  - `input_centro_reg`: The path to the tab-delimited file containing the coordinates of the pericentromeric region.
  - `prefix_chr`: The prefix for the chromosome number, such as `"Chr"`, `"chr"`, etc. This allows concatenating the prefix with the chromosome number in the output dictionary. By default, the prefix is `"Chr"`.
- **Functionality:**
  - The function begins by checking if the input file exists and is not empty by calling `CheckInput(input_centro_reg)`.
  - Then, a dictionary `centromere` is created to store the coordinates of the pericentromeric region of each chromosome.
  - The function opens the input file and iterates through each line. For each line, it separates the elements by tabs (`\t`) and checks if there are three elements in the line (chromosome name, start, and end of the pericentromeric region). If so, it extracts the chromosome name and the coordinates of the pericentromeric region and adds them to the `centromere` dictionary.
  - After processing all the lines, the function displays a message indicating the number of pericentromeric regions of chromosomes loaded.
  - Finally, it returns the `centromere` dictionary containing the coordinates of the pericentromeric region for each chromosome.
- **Output:**
  - A dictionary where the keys are the names of the chromosomes and the values are lists containing the start and end coordinates of the pericentromeric region for each chromosome.

Function : ReadParentalVCF
---------------------------

**Code :**

.. code-block:: python

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
    
..


**Explanation:**

`ReadParentalVCF`: This function reads a VCF (Variant Call Format) file containing SNP markers between parental lines and returns two dictionaries containing all SNPs by chromosome and the last SNP by chromosome.
- **Input:**
  - `input_vcf`: The path to the VCF file containing SNP markers between the parental lines.
  - `prefix_chr`: The prefix for the chromosome number, such as `"Chr"`, `"chr"`, etc. This allows concatenating the prefix with the chromosome number in the output dictionary. By default, the prefix is `"Chr"`.
- **Functionality:**
  - The function starts by displaying a message informing that it is reading the VCF file.
  - It then checks if the input file exists and is not empty by calling `CheckInput(input_vcf)`.
  - Two empty dictionaries, `snps` and `last_snps`, are initialized to store the SNPs and the last SNP by chromosome, respectively.
  - The function opens the VCF file and iterates through each line. If the line starts with `#`, it is ignored as it represents comments in the VCF file.
  - For each line, it splits the elements by tabs (`\t`) and extracts the chromosome name (`chr`) and the position (`pos`) of the SNP.
  - If the chromosome name does not start with the specified prefix, it corrects it by adding the prefix.
  - It creates a key for the `snps` dictionary combining the chromosome name and the SNP position.
  - It updates the `last_snps` dictionary to store the position of the last SNP for each chromosome.
  - It also extracts the genotype of the SNP from the ninth column of the VCF line and adds it to the corresponding entry in the `snps` dictionary.
  - Once all lines are processed, the function displays a message indicating the number of parental SNP markers loaded.
  - Finally, it returns the dictionaries `snps` and `last_snps` containing the SNP information.
- **Output:**
  - Two dictionaries:
    - `snps`: A dictionary where keys are strings of the form `chr_pos` (chromosome name + SNP position) and values are lists containing the chromosome name, SNP position, and genotype.
    - `last_snps`: A dictionary where the keys are the names of the chromosomes and the values are the position of the last SNP for each chromosome.

Function : ReadOffspringVCF
----------------------------

**Code :**

.. code-block:: python

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
    
..


**Explanation:**

`ReadOffspringVCF`: This function reads a VCF (Variant Call Format) file containing SNP markers among offspring lines and returns two dictionaries. The first dictionary, `info_snps`, contains all SNPs by chromosome with details such as genotype, number of reference reads (ADref), and number of alternative reads (ADalt). The second dictionary, `info_last_snps`, stores the position of the last SNP for each chromosome.
- **Input:**
  - `input_vcf`: The path to the VCF file containing the SNP markers of the offspring lines.
  - `parental_snps`: The first output of the `ReadParentalVCF()` function, containing parental SNPs.
  - `geno_ref`: The reference genotype.
  - `geno_alt`: The alternative genotype.
  - `analyze_id`: The analysis ID, used to create the path for the log file.
  - `prefix_chr`: The prefix for the chromosome number, such as `"Chr"`, `"chr"`, etc. The default is `"Chr"`.
- **Functionality:**
  - The function starts by displaying a message informing that it is reading the VCF file.
  - It checks if the chromosome number prefix is valid. If the prefix is empty, it raises an exception.
  - It initializes three counters: `total_snp_num` for the total number of SNPs, `weird_snp` for SNPs with strange genotypes, and `new_snp` for new offspring-specific SNPs.
  - It iterates through each line of the VCF file. If the line starts with `#`, it is ignored as it represents comments in the VCF file.
  - It checks if the sequencing depth (`DP`) is zero. If so, it skips the SNP due to zero coverage.
  - It extracts the chromosome name (`chr`) and the position (`pos`) of the SNP from the line.
  - It corrects the chromosome name if necessary by adding the specified prefix.
  - It creates a key for the `info_snps` dictionary by combining the chromosome name and SNP position.
  - It extracts the genotype (`geno`), the number of reference reads (`ref_supp`), and the number of alternative reads (`alt_supp`) from the ninth column of the VCF line.
  - It checks if the key exists in the `parental_snps` dictionary (parental SNPs). If it does, it updates the `info_snps` dictionary with the appropriate information.
  - If the genotype does not conform (`0/0`, `0/1`, or `1/1`), it increments the `weird_snp` counter and logs the corresponding line in a file named `ReadOffspringVCF_weird_snps.log`.
  - If the key does not exist in the `parental_snps` dictionary, it increments the `new_snp` counter and logs the corresponding line in a file named `ReadOffspringVCF_new_snps.log`.
  - It updates the `info_last_snps` dictionary to store the position of the last SNP for each chromosome.
  - Once all lines are processed, the function displays information about the total number of SNPs loaded, the number of SNPs retained, the number of new offspring-specific SNPs, the number of SNPs with strange genotypes, and the number of SNPs with zero coverage.
  - Finally, it returns the `info_snps` and `info_last_snps` dictionaries containing information about the offspring SNPs.
- **Output:**
  - Two dictionaries:
    - `info_snps`: A dictionary where the keys are strings of the form `chr_pos` (chromosome name + SNP position) and the values are lists containing the chromosome name, SNP position, genotype, number of reference reads (ADref), and number of alternative reads (ADalt).
    - `info_last_snps`: A dictionary where the keys are the names of the chromosomes and the values are the position of the last SNP for each chromosome.

Function : ReadEMSline
-----------------------

**Code :**

.. code-block:: python

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
    
..


**Explanation:**

`ReadEMSline`: This function reads a file and extracts the necessary information for EMS (Ethyl Methane Sulfonate) analysis. It extracts the chromosome, position, and alternative allele from each line of the file. Then, it combines these pieces of information to form a unique key for each SNP.
- **Input:**
  - `input`: The path to the file to be read.
  - `prefix_chr`: The prefix for the chromosome number, such as `"Chr"`, `"chr"`, etc.
- **Functionality:**
  - The function starts by displaying a message indicating that it is reading the file.
  - It checks if the chromosome number prefix is valid. If the prefix is empty, it raises an exception.
  - It initializes two empty lists to store the SNP keys.
  - It iterates through each line of the file.
  - If the line starts with `#`, it is ignored as it represents comments.
  - It splits the line into elements using the tabulator as a separator.
  - It extracts the chromosome name (`chr`), position (`pos`), and alternative allele (`alt`) from the line.
  - It corrects the chromosome name if necessary by adding the specified prefix.
  - It forms a unique key for each SNP by combining the chromosome name, position, and alternative allele.
  - It forms a second key with only the chromosome name and position.
  - It adds both keys to their respective lists.
  - Once all lines are processed, it returns the two lists containing the SNP keys.
- **Output:**
  - Two lists:
    - `snps_list`: A list containing unique keys for each SNP, formed by combining the chromosome name, position, and alternative allele.
    - `snps_list2`: A list containing unique keys for each SNP, formed by combining only the chromosome name and position.
        

Function : IdentifyGeno
-------------------------

**Code :**

.. code-block:: python

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
    
..


**Explanation:**

`IdentifyGeno`: This function determines the genotype of a SNP using information from the VCF file and EMS SNP markers. It considers the raw genotype (`geno`) and information about alternative alleles, reference and alternative genotypes for each SNP in two lines (`geno_line1` and `geno_line2`), as well as SNP keys for these lines (`snps_list_line1` and `snps_list_line2`).

- **Input:**
  - `chr`: The chromosome number of the SNP.
  - `pos`: The position of the SNP on the chromosome.
  - `alt_list`: A list of alternative alleles for the SNP.
  - `geno`: The raw genotype of the SNP for the considered line.
  - `snps_list_line1`: The list of SNP keys for the first line.
  - `snps_list_line2`: The list of SNP keys for the second line.
  - `geno_ref`: The reference genotype.
  - `geno_alt`: The alternative genotype.
  - `geno_line1`: The genotype of the first line.
  - `geno_line2`: The genotype of the second line.

- **Functionality:**
  - If the raw genotype is "./.", indicating it is missing, the function returns "NA" to signify that the final genotype cannot be determined.
  - Otherwise, the function starts by initializing an empty list called `genotype` to store genotypes associated with raw alleles.
  - It splits the raw genotype into a list `gt_list` using the separator "/".
  - For each element `gt` in `gt_list`:
    - If `gt` is "0", indicating it is a reference allele, the function adds `geno_ref` to the `genotype` list.
    - Otherwise, the function creates a unique key `snp_key` by combining the chromosome number, SNP position, and corresponding alternative allele.
    - Then, it determines the genotype associated with this key in the SNP lists of both lines (`snps_list_line1` and `snps_list_line2`).
    - If the key is present in `snps_list_line1`, it adds `geno_line1` to the `genotype` list.
    - Alternatively, if the key is present in `snps_list_line2`, it adds `geno_line2` to the `genotype` list.
    - If the key is not present in any of the lists, it adds `geno_alt` to the `genotype` list.
  - Then, it checks if both genotypes in the `genotype` list are identical. If so, this indicates the line is homozygous for this SNP, and the final genotype is simply the genotype of the first position in `genotype`. Otherwise, it creates the final genotype by joining the genotypes in `genotype` with the separator "/".
  - Finally, it returns the final genotype.

- **Output:**
  - `final_geno`: The final genotype of the SNP for the considered line.
    

Function : ReadRecombinedOffspringVCF
--------------------------------------

**Code :**


.. code-block:: python

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
    
..

**Explanation:**

`ReadRecombinedOffspringVCF`: This function reads a VCF file containing SNP markers between recombined offspring lines and returns two dictionaries containing all SNPs by chromosome and the last SNPs by chromosome.
- **Input:**
  - `input_vcf`: Path to the VCF file containing SNP markers between offspring lines.
  - `input_ems_line1`: Path to the file containing SNPs from the first EMS line.
  - `input_ems_line2`: Path to the file containing SNPs from the second EMS line.
  - `parental_snps`: The first output of the `ReadParentalVCF` function.
  - `geno_ref`: Reference genotype.
  - `geno_alt`: Alternative genotype.
  - `geno_line1`: Genotype of the first EMS line.
  - `geno_line2`: Genotype of the second EMS line.
  - `analyze_id`: Analysis identifier.
  - `prefix_chr`: Prefix for the chromosome number (e.g., "Chr", "chr", etc.).
- **Functionality:**
  - The function begins by validating the chromosome prefix and checking if the VCF input file exists.
  - It initializes empty dictionaries to store SNPs and the last SNPs.
  - It reads the SNPs from the EMS lines from the specified files.
  - Then, it processes the VCF file line by line.
  - For each line, it checks if it starts with "#" (comment or header); if so, the line is skipped.
  - It also checks if the genotype is missing ("./."). If so, it increments the counter for unknown SNPs.
  - For each line containing a SNP:
    - It extracts the chromosome number and the SNP position.
    - It corrects the chromosome number if necessary.
    - It checks if the SNP is present in the parental SNPs.
    - If the SNP is found in the parental SNPs, it retrieves the raw genotype, AD support, and alternative alleles.
    - It determines the final genotype from the information about alternative alleles and genotypes of the EMS lines using the function `IdentifyGeno`.
    - If the SNP is not found in the parental SNPs, it is logged to a journal file and the counter for new SNPs is incremented.
  - Finally, it returns the dictionaries containing SNPs and the last SNPs.
- **Output:**
  - `info_snps`: Dictionary containing information about all SNPs by chromosome.
  - `info_last_snps`: Dictionary containing the positions of the last SNPs by chromosome.

This function is used to process VCF files from recombined offspring lines to identify descendant-specific SNPs and potentially impute missing genotypes.

