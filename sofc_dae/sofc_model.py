"""To run this model you will need to have scikits.odes installed, which has a solver capable of integrating differential-algebraic equation systems (DAEs).

If your conda environment for this class is called 'echem', then activate it and then install scikits.odes, using the conda-forge chanel.  In the command line this looks like:

    conda activate echem
    conda install -c conda-forge scikits.odes
    type 'y' and then hit return, when prompted.
    
This model solves two types of algebraic variables:
    1. The cathode electric potential is set as an input (potentiostatic boundary condition). This is in contrast to the galvanostatic boudnary condition (i_ext is an input) of the other simulations.
    2. We discretize the electrolyte, and therefore the electric potential profile as a function of electrolyte depth must be such that i_io = i_ext at each depth.
    
The code structure at this level is meant to be very simple, and delegates most functions to lower-level modules, which can be updated with new capabilties, over time.  The code:

    1 - Initializes the model by calling sofc_init.py
        a - sofc_init.py reads in the user inputs from sofc_inputs.py, 
        b - sofc_init.py then creates an initial solution vector 'SV_0'
        c - sofc_init.py returns SV_0 and a class 'pars' holding all simulation parameters.
    2 - Calls the function sofc_function.py and integrate over the 
        user-defined time span from sofc_intputs.py.  For this simulation, we want steady-state behvaior, so we simply integrate over a sufficiently long time span (100 s), with a fixed boundary conition.
    3 - The simulation then returns the solution vector at the end of the 
        integration, which can be processed as needed to plot or analyze any quantities of interest.
"""

def sofc_model(V_cell=None, Temp=None):
    # Import necessary modules:
    from scikits.odes.dae import dae #integration function for DAE system.
    from sofc_function import residual # point to the residual function

    # Import the initial solution, parameters, and pointer structures:
    from sofc_init import SV_0, SVdot_0, ptr, pars

    # Parse and overwrite any variables passed to the function call:
    if V_cell:
        pars.V_cell = V_cell
    if Temp:
        # Adjust ideal gas concentrations in SV_0:
        # SL stands for "support layer"
        # AL stands for "active layer"
        SV_0[ptr.C_k_an_SL] *= pars.T/Temp
        SV_0[ptr.C_k_an_AL] *= pars.T/Temp
        # Overwrite temperature:
        pars.T = Temp

     # Set up the differential algebraic equation (dae) solver:
    options =  {'user_data':(pars, ptr), 'rtol':1e-4, 'atol':1e-6,  
        'algebraic_vars_idx':pars.algvars, 'compute_initcond':'yp0'}
    
    solver = dae('ida', residual, **options)

    # solver.init_step(0.0, SV_0, SVdot_0)
    solution = solver.solve(pars.time_span, SV_0, SVdot_0)

    # Return the solution results to whatever routine called the function:
    return solution


# If you want to run this as a script, doing one-off simulations, this will 
#    call the pemfc_model function, above
if __name__ == '__main__':
    sofc_model()