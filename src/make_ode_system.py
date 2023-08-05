#!/bin/python3
"""Script reimplementing the reaction_to_ode code.

Motivations to rewrite it instead of just refactor:
    - Better organization (which could help others extend)
    - Fix some systematic issues regarding duplicate species
    - Update the code to use modern Python3 features (like f-strings!)

Written by Michael Stobb

Originally written on 8/4/2023 (as procrastination for other projects)
"""

# Imports
import os, re, json, argparse # Yeah, that's it!


# Create class definitions for species
class Species:
    def __init__(self, name: str, multiplicity: int) -> None:
        self.name = name
        self.multiplicity = multiplicity

    # Define species equality to be their names match, not multiplicity
    def __eq__(self, other):
        if not isinstance(other, Species):
            return NotImplemented
        return self.name == other.name

    # Define the spec for hashing species (just by name)
    def __hash__(self):
        return hash(self.name)
    
    # Define way to print the species
    def __str__(self):
        return f"Name:{self.name}, Multiplicity: {self.multiplicity}"
    
    def __repr__(self):
        return f"Species:\nName:{self.name}, Multiplicity: {self.multiplicity}"

    
# Create class for reactions, each is only one direction
class Reaction:
    def __init__(self, raw_text: str, reactants: list, products: list, rate: str) -> None:
        self.raw_text = raw_text
        self.reactants = reactants
        self.products = products
        self.rate = rate
    
    def get_species(self):
        all_species = self.reactants + self.products
        all_species = list(set(all_species))
        return all_species
    
    # Define way to print the reaction
    def __str__(self):
        return f"Raw: {self.raw_text}\nReactants: {self.reactants}\nProducts: {self.products}"
    
    def __repr__(self):
        return f"\nReaction:\nRaw: {self.raw_text}\nReactants: {self.reactants}\nProducts: {self.products}\nRate: {self.rate}"
    
    
# Create class for the stoichiometry matrix 
class Stoich:
    def __init__(self, species: list, rates: list) -> None:
        self.species = dict()
        c = 0
        for s in species:
            if s not in self.species:
                self.species[s] = c
                c += 1
       
        self.rates = dict()
        c = 0
        for r in rates:
            if r not in self.rates:
                self.rates[r] = c
                c += 1
    
        self.N = len(self.species)
        self.M = len(self.rates)
        self.stoich = [[0 for i in range(self.M)] for j in range(self.N)]
    
    def update(self, reaction: Reaction):
        # Find the index of the rate
        rate_ind = self.rates[reaction.rate]
       
        # Deal with reactants first
        for spec in reaction.reactants:
            spec_ind = self.species[spec]

            # Subtract species multiplicities from stoich matrix
            self.stoich[spec_ind][rate_ind] += - spec.multiplicity

        # Deal with products next
        for spec in reaction.products:
            spec_ind = self.species[spec]

            # Add species multiplicities from stoich matrix
            self.stoich[spec_ind][rate_ind] += spec.multiplicity
            
    # Define how to print stoich objects
    def __str__(self):
        return f"Stoich Info:\nSpecies: {list(self.species.keys())}\nRates: {list(self.rates.keys())}\nMatrix: {self.stoich}\n"
    def __repr__(self):
        return f"Stoich Info:\nSpecies: {list(self.species.keys())}\nRates: {list(self.rates.keys())}\nMatrix: {self.stoich}\n"


# Define helper functions for parsing file
def species_parser(text: str):
    # Species should only be additive for mass action reactions
    split_text = re.split(r"\W*\+\W*", text)
    
    species = []

    for s in split_text:
        sm = re.findall(r"(\d+)",s)
        ss = re.findall(r"([a-zA-Z]+)\W*$",s)

        if len(sm) == 0:
            species += [Species(ss[0],1)]
        else:
            species += [Species(ss[0],int(sm[0]))]
            
    return species

def line_parser(text: str, param_start: int = 1):
    split_text = re.split("\W*,\W*",  text)
    reactants = split_text[0]
    if len(split_text) > 1:
        rates = split_text[1:]
    else:
        rates = []
    
    sr = re.split(r"(=|->|<->)", reactants)
    assert len(sr) % 2 == 1, "Error: Line has wrong number of reaction operators"
    
    reactions = []
    
    c = param_start
    
    for i in range(1,len(sr),2):
        if (sr[i] == "=") or (sr[i] == "<->"):
            if len(rates) < 2:
                rates += [f"q{c}", f"q{c+1}"]
                c += 2
            reactions += [Reaction(sr[i-1]+"->"+sr[i+1], species_parser(sr[i-1]), species_parser(sr[i+1]), rates.pop(0).strip())]
            reactions += [Reaction(sr[i-1]+"<-"+sr[i+1], species_parser(sr[i+1]), species_parser(sr[i-1]), rates.pop(0).strip())]
        else:
            if len(rates) < 1:
                rates += [f"q{c}"]
                c += 1
            reactions += [Reaction(sr[i-1]+"->"+sr[i+1], species_parser(sr[i-1]), species_parser(sr[i+1]), rates.pop(0).strip())]
            
    return reactions

def get_all_species_and_rates(reactions: list):
    spec = []
    rates = []
    for r in reactions:
        spec += r.get_species()
        rates += [r.rate]
    return spec, rates



#####################
# Output Generators
#####################


def create_matlab_ouput(input_file: str, output_file: str, s: Stoich, v: bool = False):

    Ns = len(s.species)
    Nr = len(s.rates)
    species = [i.name for i in s.species.keys()]
    rates = [i for i in s.rates.keys()]
    if v: print('\nOutput File: \n', output_file)
    if v: print("\n")
    f = open(output_file,'w')

    f.write("function [time,y] = ")
    f.write(output_file.strip(".m"))
    f.write("(t_final,t_start)\n")
    f.write("% Solves a system of ODEs from t=t_start to t=t_final \n")
    f.write("% If no start time is given, then t_start = 0 \n")
    f.write("% If no start or final time is given, then t_start = 0, t_final = 1 \n")
    f.write("%\n")
    f.write("%\n")
    f.write("% This file was created by issuing command: \n")
    f.write("% 	python make_ode_system.py ")
    f.write(input_file)
    f.write("\n")
    f.write("%\n")
    f.write("\n")
    f.write("\nif nargin == 1\n")
    f.write("     t_start = 0;  % Default start time is 0 \n")
    f.write("elseif nargin == 0\n")
    f.write("     t_start = 0; % Default start time is 0\n")
    f.write("     t_final = 1; % Default final time is 1\n")
    f.write("end\n\n\n")

    f.write("% Kinetic Parameters \n")
    for i in range(Nr):
        f.write(rates[i])
        f.write(" = 1; \n")

    f.write("\np = [ ")
    f.write(rates[0])
    for i in range(1, Nr):
        f.write(", ")
        f.write(rates[i])
    f.write("];\n\n\n")

    
    f.write("% Initial Conditions \n")
    for i in range(Ns):
        f.write(species[i])
        f.write("_IC")
        f.write(" = 0; \n")

    f.write("\ninit_cond = [ ")
    f.write(species[0])
    f.write("_IC")
    for i in range(1,Ns):
        f.write(", ")
        f.write(species[i])
        f.write("_IC")
    f.write("];\n\n\n")

    f.write("options = odeset('RelTol',1e-12,'AbsTol',1e-23);\n\n\n")

    f.write("%------------------------------ Main Solve ---------------------------%\n")
    f.write("[time,y] = ode15s(@(t,y)RHS(t,y,p), [t_start t_final], init_cond, options);\n")
    f.write("%---------------------------------------------------------------------%\n\n\n")

    f.write("% Rename solution components\n")
    for i in range(Ns):
        f.write(species[i])
        f.write(" = y(:,")
        f.write(str(i+1))
        f.write("); \n")

    f.write("\n\n\n")

    f.write("%  \n")
    f.write("% Place plots or other calculations here\n")
    f.write("%   \n")
    f.write("% Example: \n")
    f.write("% plot(time, ")
    f.write(str(species[0]))
    f.write(", 'k-o', 'LineWidth', 4, 'MarkerSize', 4); legend('")
    f.write(str(species[0]))
    f.write("');\n\n\n")

    f.write("end\n\n\n\n")

    f.write("%-----------------------------------------------------%\n")
    f.write("%-------------------- RHS Function -------------------%\n")
    f.write("%-----------------------------------------------------%\n\n")


    f.write("function dy = RHS(t,y,p)\n\n")
    f.write("dy = zeros(") 
    f.write(str(Ns))
    f.write(",1);\n")
    f.write("\n\n")

    f.write("% Rename Variables \n\n")

    for i in range(Ns):
        f.write(str(species[i]))
        f.write("   = y(")
        f.write(str(i+1))
        f.write("); \n")
    f.write("\n\n")

    f.write("% Rename Kinetic Parameters \n")
    for i in range(Nr):
        f.write(str(rates[i]))
        f.write(" = p(")
        f.write(str(i+1))
        f.write(");  \n")

    f.write("\n\n")


    f.write("\n\n")

    f.write("% ODEs from reaction equations \n\n")

    if v: print("Writing ODEs now....\n")

    for i in range(Ns):
        f.write("% ")
        f.write(str(species[i]))
        f.write("\n dy(")
        f.write(str(i+1))
        f.write(")  =")
        for j in range(Nr):
            if int(s.stoich[i][j]) < 0:
                f.write("  -  ")
                f.write(str(rates[j]))
                for k in range(Ns):
                    if int(s.stoich[k][j]) < 0:
                        f.write(" * ")
                        f.write(species[k])
                        if abs(int(s.stoich[k][j])) != 1:
                            f.write("^")
                            f.write(str(abs(int(s.stoich[k][j]))))
            elif int(s.stoich[i][j]) > 0:
                f.write("  +  ")
                f.write(str(rates[j]))
                for k in range(Ns):
                    if int(s.stoich[k][j]) < 0:
                        f.write(" * ")
                        if int(s.stoich[i][j]) > 1:
                            f.write(str(int(s.stoich[i][j])))
                            f.write(" * ")
                        f.write(species[k])
                        if abs(int(s.stoich[k][j])) != 1:
                            f.write("^")
                            f.write(str(abs(int(s.stoich[k][j]))))
        if v: print(species[i]," complete")
        f.write(";\n\n")


    f.write("\n\n\n\n")
    f.write("end")
    

    

###############################################
# Main Function
###############################################


# Main function
def main(filename: str, verbose: bool = False, log: bool = True, output: str = 'm', output_file: str = ''):
    file_prefix = filename.split(".")[0]
    if len(output_file) == 0:
        output_file = file_prefix
    reactions = []
    lines = ""

    # Read file and parse reactions
    if verbose: print("Reading file and parsing reactions\n")
    with open(filename,"r") as fid:
        for line in fid:
            lines += line
            reactions += line_parser(line)
    
    # Split out species and rates
    species, rates = get_all_species_and_rates(reactions)
    
    # Initialize the stoichiometry matrix
    if verbose: print("Creating Stoichiometry Matrix")
    s = Stoich(species, rates)
    
    # Fill in the matrix
    for r in reactions:
        s.update(r)
        
    # Create log if desired
    if log:
        log_dict = {"input":line, "reactions":[repr(i) for i in reactions], "species":[repr(i) for i in s.species], "rates":[repr(i) for i in s.rates], "stoich":repr(s)}
        json.dump(log_dict, open(file_prefix+".log","w"), indent=4)
        
    # Generate output file(s)
    if output == 'm':
        create_matlab_ouput(filename, output_file+".m", s, verbose)
    else:
        print("That output type has not yet been implemented.\n\nExiting without output.")
    
    
    if verbose: print("\nFinished successfully\n")


########################################################
# Script Code
########################################################

if __name__ == "__main__":
    desc_text = f"""

This small script converts a set of mass action reactions to
a set of numerical differential equations that are ready for
further simulation.

"""
    
    parser = argparse.ArgumentParser(description = desc_text, epilog="Written by: Michael Stobb.  Last updated: 8/5/2023")
    parser.add_argument('filename', help="Name of input file with reactions")
    parser.add_argument('-t','-type', default='m', help="Type of output file to generate. 'm' for Matlab")
    parser.add_argument('-o','-output', default="", help="Output file name.  If none, will use input filename prefix")
    parser.add_argument('-l','-log', action="store_true", help="Turn on extra logging, will produce json file of calculations")
    parser.add_argument('-v','-verbose', action="store_true", help="Print more comments to screen")
    
    args = parser.parse_args()
    
    print(args.filename)
    
    input_file = args.filename      
              
    # Check if file exists
    if not os.path.exists(input_file):
        raise Exception("File does not exist!")

    main(input_file, args.v, args.l, args.t, args.o)