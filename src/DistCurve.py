import numpy as np
import cantera as ct
import time as chrono
import Global
import logging
import GC_Python.distCurve as DC


class DistCurve(object):

    def __init__(self,spec_list):
       self.valid_species_list = spec_list

       return

    def storeMultiplier(self, species_index_exclude, species_index_damp, multipliers):
        self.species_index_exclude = species_index_exclude
        self.species_index_damp = species_index_damp
        self.multipliers = multipliers
        return

    def putMultiplier(self, curgas):
        reactants = curgas.reactant_stoich_coeffs()
        products = curgas.product_stoich_coeffs()

        # Damping of reactions
        multicoeffs = np.ones(curgas.n_reactions)
        for k in range(curgas.n_reactions):
            for kk, index in enumerate(self.species_index_damp):
                if (not(reactants[index, k] == 0)) or (not(products[index, k] == 0)):
                    multicoeffs[k] *= self.multipliers[kk]
        for k in range(curgas.n_reactions):
            curgas.set_multiplier(multicoeffs[k], k)

        # Setting to zero for excluded species
        for k in range(curgas.n_reactions):
            for index in self.species_index_exclude:
                if (not(reactants[index, k] == 0)) or (not(products[index, k] == 0)):
                    curgas.set_multiplier(0.0, k)
        return

    def getDistCurveError(self,x):
        MAX = 0.0
        # Return 0 if invalid composition
        if (x[-1] < -Global.sum_relax):
          return MAX
        elif (x[-1] > -Global.sum_relax and x[-1] < 0.0):
          x[-1] = 0.0

        x_val = 0.0
        # Set the species vector if unset
        if (all(x == 0)):
          # Corresponds to palette setting
          x_val = self.set_using_palette()
        else:
          x_val = x

        # Calculate objective function
        obj = DC.distMoleFrac(self.valid_species_list,x_val,'posf4658_simdist','fractional')

        return obj

    # Generic functions
    def getSpeciesIndex(self,species,specList):
      for idx, sp in enumerate(specList):
        if (species == sp):
          return idx
      return -1

    def set_using_palette(self):
      x = np.zeros(len(Global.valid_species))
      for idx, species in enumerate(Global.palette):
        idx_spec = self.getSpeciesIndex(species,Global.valid_species)
        if (idx_spec != -1):
          x[idx_spec] = Global.test_comp[idx]
      return x

    def getQuantity(self,x):
        return self.getDistCurveError(x)
