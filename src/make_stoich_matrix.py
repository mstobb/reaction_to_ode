#!/usr/bin/python
# Python code that generates a Stoichiometric Matrix for Law of Mass action Reactions
# 
# Use: make_stoich_matrix.py foo.txt foo.csv
#
# foo.txt must have the following structure:
# 
#	SINGLE_REACTION , RATE_VARIABLE
# ex:
#	A + 2 * B -> C , k_1
#
# In the case of reversible reactions, two rates must be specified, the forward reaction
#	is always FIRST, eg:
#
#	A + 2 * B <-> C , k_1 , k_2
#
# Output File: foo.csv
#
# General Requirements:
#
#	1) Only the operators '*', '+', '=', '<->', and '->' are allowed
#	2) Reactions should be pre-simplified, this code will NOT reduce algebra
#	3) Species names can only contain [0-9a-zA-Z_:]
#	4) Species names MUST begin with a letter - no leading numbers allowed
#	5) While not technically required, given rate parameters are perferred
#
#
# Created by Michael T. Stobb
# Created on 1/30/15
# Update for Python3 8/3/2023

try:
    from tabulate import tabulate
except ImportError:
    print(" ************************************************ \n")
    print(" * Python tabulate module not found....         * \n")
    print(" * Consider running: sudo pip install tabulate  * \n")
    print(" *            THIS IS OPTIONAL                  * \n")
    print(" ************************************************ \n")
    from pprint import pprint as pp

import re , sys, types, csv, textwrap

if len(sys.argv) == 1:
    print(textwrap.dedent("""

                          Python code that generates a Stoichiometric Matrix for Law of Mass action Reactions

                          Use: make_stoich_matrix.py foo.txt foo.csv

                          foo.txt must have the following structure:

                          SINGLE_REACTION , RATE_VARIABLE
                          ex:
                          A + 2 * B -> C , k_1

                          In the case of reversible reactions, two rates must be specified, the forward reaction
                          is always FIRST, eg:

                          A + 2 * B <-> C , k_1 , k_2

                          Output File: foo.csv

                          General Requirements:

                          1) Only the operators '*', '+', '=', '<->', and '->' are allowed
          2) Reactions should be pre-simplified, this code will NOT reduce algebra
    3) Species names can only contain [0-9a-zA-Z_:]
       4) Species names MUST begin with a letter - no leading numbers allowed
       5) While not technically required, given rate parameters are perferred

"""))
else:
    input_file = sys.argv[1]
    num_reactions = sum(1 for line in open(input_file))
    num_reactions_orig = num_reactions
    double_react = 0
    print(" Opening Input File: ", input_file)
    reactions = open(input_file)
    operators = ['+', '=', '<->', '->', '*']
    species_list = []
    species_set  = set()
    line_number = 0
    a = re.compile(r"^\s*([1-9]+\d*)?\s*\*?\s*([\w:]+)\s*$")

    # First find the species and count the number of double reactions
    for line in reactions:
        line = line.rstrip() 
        line = line.replace(" ", "")
        #x = line.split() # Breakup line by whitespace
        react_type = line.count("<->") + line.count("=")
        if react_type > 1:
            print("Too many main operators on line:", line_number+1)
            raise NameError('Input file syntax')
        elif react_type == 1:
            # Reversible Reaction
            double_react += 1
            if line.count("=") != 0:
                x = line.split("=")
            else: 
                x = line.split("<->")
            react = x[0].split("+")
            for i in react:
                terms = re.match(a,i)
                spec = terms.groups()[1]
                if spec not in species_set:
                    species_set.add(spec)
                    species_list.append(spec)
            try:
                prods_full = x[1].split(",")
            except:
                prods_full = x[1]
            prods = prods_full[0]
            prods = prods.split("+")
            for i in prods:
                terms = re.match(a,i)
                spec = terms.groups()[1]
                if spec not in species_set:
                    species_set.add(spec)
                    species_list.append(spec)

        else:
            # Single Reaction
            x = line.split("->")
            react = x[0].split("+")
            for i in react:
                terms = re.match(a,i)
                spec = terms.groups()[1]
                if spec not in species_set:
                    species_set.add(spec)
                    species_list.append(spec)
            prods_full = x[1].split(",")
            try:
                prods = prods_full[0]
            except:
                prods_full = x[1]
            prods = prods.split("+")
            for i in prods:
                terms = re.match(a,i)
                spec = terms.groups()[1]
                if spec not in species_set:
                    species_set.add(spec)
                    species_list.append(spec)

    num_reactions += double_react 
    num_species = len(species_list)
    stoich_matrix = [[0 for x in range(num_reactions)] for x in range(num_species)] 
    numrows = len(stoich_matrix)
    numcols = len(stoich_matrix[0])
    rates = [0 for x in range(num_reactions)]
    print("\n Species found: ", species_list)
    print("\n Number of species:", num_species)
    print("\n Number of reactions:", num_reactions_orig)
    print("\n Numer of total reactions (counting reversible): ",num_reactions )


    # Generate the Stoichiometry Matrix


    reactions = open(input_file)
    rc = 0 # reaction counter
    for line in reactions:
        line = line.rstrip() 
        line = line.replace(" ", "")
        #x = line.split() # Breakup line by whitespace
        react_type = line.count("<->") + line.count("=")
        if react_type > 1:
            print("Too many main operators on line:", line_number)
        elif react_type == 1:
            # Reversible Reaction
            double_react += 1
            if line.count("=") != 0:
                x = line.split("=")
            else: 
                x = line.split("<->")
            react = x[0].split("+")
            for i in react:
                terms = re.match(a,i)
                fact = terms.groups()[0]
                if fact == None:
                    fact = 1
                spec = terms.groups()[1]
                r = species_list.index(spec)
                stoich_matrix[r][rc] = -int(fact)
                stoich_matrix[r][rc+1] = int(fact)    

            prods_full = x[1].split(",")
            if len(prods_full) > 3:
                print("\n\n Error: Too many reaction variables found on line ", line_number+1)
                break
            elif len(prods_full) == 1:
                print("\n\n No reaction variable found for both reactions on line ", line_number+1)
                arb_rates = [ "arb_var_"+str(line_number+1)+"_1" , "arb_var_"+str(line_number+1)+"_2"]
                print(" Will use arbitrary variable names:", arb_rates, "\n\n")
                prods_full += arb_rates
            elif len(prods_full) == 2:
                print("\n\n No reaction variable found for 2nd reaction on line ", line_number+1)
                arb_rates = [ "arb_var_"+str(line_number+1)+"_1"]
                print(" Will use arbitary variable name:", arb_rates,"\n\n")
                prods_full += arb_rates
            prods = prods_full[0]
            rates[rc] = prods_full[1]
            rates[rc+1] = prods_full[2]
            prods = prods.split("+")
            for i in prods:
                terms = re.match(a,i)
                fact = terms.groups()[0]
                spec = terms.groups()[1]
                if fact == None:
                    fact = 1
                spec = terms.groups()[1]
                r = species_list.index(spec)
                stoich_matrix[r][rc] = int(fact)
                stoich_matrix[r][rc+1] = - int(fact)   

            rc+=2 # two reactions finished
            line_number += 1
        else:
            # Single Reaction
            x = line.split("->")
            react = x[0].split("+")
            for i in react:
                terms = re.match(a,i)
                fact = terms.groups()[0]
                spec = terms.groups()[1]
                if fact == None:
                    fact = 1
                spec = terms.groups()[1]
                r = species_list.index(spec)
                stoich_matrix[r][rc] = - int(fact)
            prods_full = x[1].split(",")
            if len(prods_full) > 2:
                print("\n\n Error: Too many reaction variables found on line ", line_number+1)
                break
            elif len(prods_full) == 1:
                print("\n\n No reaction variable found for single reaction on line ", line_number+1)
                arb_rates = [ "arb_var_"+str(line_number+1)+"_1"]
                print(" Will use arbitrary variable name:", arb_rates,"\n\n")
                prods_full += arb_rates
            prods = prods_full[0]
            rates[rc] = prods_full[1]
            prods = prods.split("+")
            for i in prods:
                terms = re.match(a,i)
                fact = terms.groups()[0]
                spec = terms.groups()[1]
                if fact == None:
                    fact = 1
                spec = terms.groups()[1]
                r = species_list.index(spec)
                stoich_matrix[r][rc] = int(fact)
            rc+=1
            line_number += 1



    # Format the output file

    stoich_matrix_output = [[0 for x in range(num_reactions+1)] for x in range(num_species+1)]

    for i in range(1,num_reactions+1):
        stoich_matrix_output[0][i] = rates[i-1]
    for i in range(1,num_species+1):
        stoich_matrix_output[i][0] = species_list[i-1]
    for i in range(1,num_species+1):
        for j in range(1,num_reactions+1):
            stoich_matrix_output[i][j] = str(stoich_matrix[i-1][j-1])

    # Print stoichiometry matrix to screen
    if len(stoich_matrix_output) < 30:
        print("\n Stoichiometry Matrix:")
        try:
            print(tabulate(stoich_matrix_output))
        except:
            pp(stoich_matrix_output)

    # Write stoichiometry matrix to csv file
    try:
        output_file = sys.argv[2]
    except:
        output_file = input_file.split('.')
        output_file = output_file[0]
        output_file+=".csv"
    print('\n Output File: \n', output_file)
    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerows(stoich_matrix_output)
