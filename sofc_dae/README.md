# SOFC Model - WORK IN PROGRESS

This SOFC Model illustrates use of a Differential-Algegraic Equation (DAE) solver.  
At present, the system of equations has one algebraic variable - the cathode electric potential, which is set by the user (i.e. this model has a potentiostatic boundary condition, rather than the galvanostatic i_ext condition).

To run this model you need to add `scikits.odes` to your software environment.  This software is installed from the popular `conda-forge` conda channel.  Via the command line, assuming you have named your conda environment `echem`, you can install via:

``` 
conda activate echem
conda install -c conda-forge scikits.odes
````
There are a few key differences from our use of `solve_ivp`:
- The setup and call to the solver (in sofc_model.py) looks a little different.  You send the time, the solution, the time derivative, and the residual to the solver. Any extra parameters you want to use can be sent via `options`, using the `user-data` field.
- Your residual function specifies a system of equations that _all_ must equal zero, always. More on this below.
- Because the residual is sent as an input to the residual function, the residual function does not actually need to `return` anything. It simply modifies the residual in place.
- You also need to tell the solver the indices of the algebraic variables.  These are set in the `options` dictionary, using the 'algebraic_vars_idx` key.

### Residual function:
The residual specifies a system of equations that equals zero.
-  For a differential solution variable, assuming the variable `SVdot` is the array of time-derivatives that you sent to the funciton:  
$$ {\tt resid}[{\tt variable\,index}]  = {\tt SVdot}[{\tt variable\,index}] - [{\tt your\, equation\, for\, dSV/dt}]$$
- This is equivalent to saying:
$$  {\tt SVdot}[{\tt variable\,index}] = [{\tt your\, equation\, for\, dSV/dt}]$$
- For an algebraic variable, you just write the equation that equals zero. In this code, for example, the cathode electric potential is set by the user, which is stored as `pars.V_cell`.  So our algebraic equation (it is algebraic because it does not involve `SVdot`) is:
$$ {\tt resid}[{\tt index\,of\,\phi_{ca}}] = SV[{\tt index\,of\,\phi_{ca}}] - {\tt pars.V\_cell}$$


## SOFC Model Details:
NOTE: this is adopted and reconfigured from the PEMFC model.  Many of the variable names have not yet been converted, which may cause some confusion.

This model simulates constant-voltage operation of a solid oxide fuel cell.  At present, the domain includes:
- Anode support layer (SL/GDL), with constant gas-phase composition (`C_k`) and electric potential(`phi`) assumed equal to zero.
- Anode active layer (AL/CL), which consists of a Ni-YSZ cermet and gas phase.
-Cathode active layer, which is assumed to have constant gas phase, for now.

The gas-phase composition and double-layer electric potential (`dPhi_dl`) in the active layers are simulated via physically-based partial differential equations.  Modeled phenomena include:
- Global half-cell electrochemical reactions, with associated Faradaic current (`i_Far`).
- Double layer charging current (`i_dl`) to maintain charge neutrality.
- Ionic (`i_io`) and electronic (`i_el`) current into and out of the AL, calculated via the electrolyte ionic current $i_{io} = \sigma_{\rm elyte}\frac{\phi_{\rm elyte,\,an} - \phi_{\rm elyte,\,ca}}{\Delta y_{\rm elyte}}$
- Gas-phase transport (diffusion + convection) between the SL and AL.

Setting parameters and conditions:
The full range of parmeters and input values can be set in `sofc_inputs.py`. 
