import numpy as np

def create_comp_string(gas, options, x):
    final_string = ''

    for idx, compound in enumerate(options.palette):
        final_string += compound + ":" + str(x[idx])
        if idx != len(x)-1:
            final_string += ","

    return final_string

def set_gas_using_palette(gas, t, p, x):
    comp_string = create_comp_string(gas, options, x)
    gas.TPX = t, p, comp_string
    return comp_string

def vec_from_palette(valid_species,palette,test_comp):
  num_valid = len(valid_species)
  ret = np.zeros(num_valid)

  for idx_palette,comp in enumerate(palette):
    idx_valid = valid_species.index(comp)
    ret[idx_valid] = test_comp[idx_palette]

  return ret
