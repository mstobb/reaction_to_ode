"""Script reimplementing the reaction_to_ode code.

Motivations to rewrite it instead of just refactor:
    - Better organization (which could help others extend)
    - Fix some systematic issues regarding duplicate species
    - Update the code to use modern Python3 features (like f-strings!)

Written by Michael Stobb

Originally written on 8/4/2023 (as procrastination for other projects)
"""

# Imports
import re, sys, json, numpy as np  # Yeah, that's it!


# Create class definitions for reactions and species
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

class Reaction:
    def __init__(self, raw_text: str, reactants: list, products: list, rate: str) -> None:
        self.raw_text = raw_text
        self.reactants = reactants
        self.products = products
        self.rate = rate
    
    def get_all_species(self):
        all_species = self.reactants + self.products
        return all_species

class Stoich:
    def __init__(self, species: list, rates: list) -> None:
        self.species = list(set(species)) # Remove all duplicates
        self.rates = list(set(rates))
        self.N = len(self.species)
        self.M = len(self.rates)
        self.stoich = np.zeros((self.N, self.M), dtype=int)
    
    def update(self, reaction: Reaction):
        # Find the index of the rate
        rate_ind = [i for i,r in enumerate(self.rates) if r == reaction.rate]
        assert len(rate_ind) == 1, "Code Error: Multilpe rates with same name found"
       
        # Deal with reactants first
        for spec in reaction.reactants:
            spec_ind = [i for i,s in enumerate(self.species) if s == spec]
            assert len(spec_ind) == 1 , "Code Error: Multiple species with same name found"

            # Subtract species multiplicities from stoich matrix
            self.stoich[spec_ind[0],rate_ind[0]] += - spec.multiplicity

        # Deal with products next
        for spec in reaction.products:
            spec_ind = [i for i,s in enumerate(self.species) if s == spec]
            assert len(spec_ind) == 1 , "Code Error: Multiple species with same name found"

            # Add species multiplicities from stoich matrix
            self.stoich[spec_ind[0],rate_ind[0]] += spec.multiplicity
            

