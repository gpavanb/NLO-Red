import numpy as np
import cantera as ct
import time as chrono
import Global
import logging

class ReactorCstVolume(object):

    def __init__(self, gas, Tin, P, phi, fuel, n2_o2_ratio):
        self.Tin = Tin
        self.P = P
        self.phi = phi
        self.fuel = fuel
        self.n2_o2_ratio = n2_o2_ratio

        ifuel = gas.species_index(self.fuel)
        io2 = gas.species_index('O2')
        in2 = gas.species_index('N2')

        x = np.zeros(gas.n_species)
        x[ifuel] = self.phi
        x[io2] = gas.n_atoms(self.fuel, 'C') + 0.25 * \
            gas.n_atoms(self.fuel, 'H')
        x[in2] = x[io2] * self.n2_o2_ratio

        gas.TPX = self.Tin, self.P, x

        gas.equilibrate('UV')
        # Ignition time = 90% of the final temperature
        self.Ttarget = 0.9 * gas.T
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

    # Make private
    def create_comp_string(self,specList,x):
        gas = Global.gas
        final_string = ''

        n_c_atoms = 0.0
        n_h_atoms = 0.0
        for idx, compound in enumerate(specList):
            final_string += compound + ":" + str(x[idx])
            final_string += ","
            n_c_atoms += gas.n_atoms(compound,'C') * x[idx]
            n_h_atoms += gas.n_atoms(compound,'H') * x[idx]

        # Add oxidizer and nitrogen
        xo2 = (n_c_atoms + 0.25 * n_h_atoms)/self.phi 
        xn2 = xo2 * self.n2_o2_ratio
        final_string += "O2:" + str(xo2) + ",N2:" + str(xn2)
  
        #print final_string
        return final_string

    def set_gas_using_palette(self):
        gas = Global.gas
        comp_string = self.create_comp_string(Global.palette,Global.test_comp)
        gas.TPX = self.Tin, self.P, comp_string
        print comp_string
        return comp_string

    # Make private
    def set_gas_using_valid_species(self,x):
        gas = Global.gas
        comp_string = self.create_comp_string(Global.valid_species,x)
        gas.TPX = self.Tin, self.P, comp_string
        return comp_string

    def getAI(self,x):
        # Return 0 if invalid composition
        if (x[-1] < 0.0):
          return 0.0
         
        gas = Global.gas
        ifuel = gas.species_index(self.fuel)
        io2 = gas.species_index('O2')
        in2 = gas.species_index('N2')

        #self.putMultiplier(gas)

        # Set the species vector if unset
        if (all(x == 0)):
          self.set_gas_using_palette()
        else:
          self.set_gas_using_valid_species(x) 

        # Create reactor network
        r = ct.Reactor(gas)
        sim = ct.ReactorNet([r])

        # Advance the reactor
        time = 0.0
        time0 = 0.0
        time1 = 0.0
        curT0 = r.T
        curT1 = r.T

        #sim.atol = 1e-7
        #sim.rtol = 1e-4

        while (curT1 < self.Ttarget):
            time = sim.step()
            curT0 = curT1
            curT1 = r.T
            time0 = time1
            time1 = time
 
            if (time > Global.tMax):
              return time

        weight = (self.Ttarget - curT0)/(curT1-curT0)

        AItime = (1.-weight) * time0 + weight * time1

        logging.debug('Done Ignition time computation')

        return AItime

    def getQuantity(self,x):
        return self.getAI(x)
