import numpy as np
from pemfc_model import pemfc_model
from pemfc_init import ptr

# Here we define a function 'polarization' which runs the pemfc model for a 
# range of current densities to make a polarization curve.
#
# The funcion accepts one input - the temperature at which it is to run. It 
# passes this temperature, plus the current density, to the fuel cell model. 
# The `T=None` specification provides the default.  If no value is provided, it # sets T=None and sends this to the pemfc model. The model will handle that as 
# it sees fit - either throw an error or assume some default value.
def polarization(T=None):
    print("Temperature = ", T, " K")
    i_array = np.linspace(0.01,20000,25)
    V_cell = np.zeros_like(i_array)

    for j, current in enumerate(i_array):
        print('    i_ext = ',current)
        solution = pemfc_model(current,T)
        V_cell[j] = solution.y[ptr.phi_dl_ca,-1] - solution.y[ptr.phi_dl_an[0],-1]
        print('        V_cell =', V_cell[j])

    return i_array, V_cell

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    
    i, V = polarization()

    plt.plot(i/10000, V)
    plt.show()
