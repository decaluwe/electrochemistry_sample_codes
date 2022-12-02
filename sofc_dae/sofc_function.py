""" Residual for 1D sofc model.
        Inputs:
            t: current simulation time (s).  Required by the integrator, but not used.
            SV: current state of the SOFC simulation domain.
            SVdot: time derivative of solution vector variables SV.
        Returns:   
            residual: array of function values that should at all times equal zero. 
"""
import numpy as np
from math import exp

# Constants:
F = 96485    # Faraday's constant, C/mol of equivalent
R = 8.3145   # Universal gas constant, J/mol-K

def residual(t, SV, SVdot, resid, inputs):
    # Read out parameters and ptr class:
    pars, ptr = inputs

    # Calculate the ionic current in the electrolyte:
    phi_an = 0 # anode potential is the reference potential.
    phi_elyte_an = phi_an - SV[ptr.phi_dl_an]
    phi_ca = pars.V_cell
    phi_elyte_ca = phi_ca - SV[ptr.phi_dl_ca]
    
    i_io = pars.sigma_dy_elyte*(phi_elyte_an - phi_elyte_ca)

    # GDL gas phase:
    j = 0
    # Support layer, save as state s1:
    s1, _  = read_state(SV, ptr, pars, 'sl_an', j)
    
    # Active layer, save as state s2:
    """Note: 'read_state' is a funciton I made up, that reads the gas phase state variables from the solution vector for node j."""
    s2, eta_2 = read_state(SV, ptr, pars, 'al_an', j)
    
    gas_props = {'T':pars.T, 'D_k':pars.D_k_g_an, 'mu':pars.mu_g_an}

    # Calculate the flux between states 1 (GDL) and 2 (first CL node)
    #   The subscript 'o' refers to the flux 'out' of the s2 node.
    N_k_o = sofc_gas_flux(s1, s2, gas_props)

    "========ANODE ACTIVE LAYER==========="
    # Re-assign state 2 's2' as 's1', previously 'out' flux as 'in' flux.
    s1 = s2
    N_k_i = N_k_o
    eta = eta_2

    # Calculate currents and chemical production rates:
    """chem_calcs is another function I made up. It calculates i_Far from the BV equation, and then calculates the double layer current."""
    i_dl, sdot_k = chem_calcs(eta, pars, i_io, 'an')
    
    # Change in double layer potential per time:
    resid[ptr.phi_dl_an[j]] = SVdot[ptr.phi_dl_an[j]] - i_dl/pars.C_dl_an

    # Gas properties used in flux calculations:
    gas_props = {'T':pars.T, 'D_k':pars.D_k_g_an, 'mu':pars.mu_g_an}\
    
    # Zero gas flux at electrolyte boundary:
    N_k_o = np.zeros_like(N_k_i)
    
    # Change in catalyst layer gas phase mole fractions:
    dCk_dt = (N_k_i - N_k_o + sdot_k*pars.A_fac_Ni)*pars.eps_g_dy_Inv_AL
    resid[ptr.C_k_an_CL[j,:]] =  SVdot[ptr.C_k_an_CL[j,:]] - dCk_dt
    resid[ptr.C_k_an_GDL] = SVdot[ptr.C_k_an_GDL]
   
    
    "========CATHODE==========="
    " THIS SHOULD EVENTUALLY DEPEND ON SPECIES ACTIVITIES. "
    eta = SV[ptr.phi_dl_ca] - pars.delta_Phi_eq_ca
    i_dl, _ = chem_calcs(eta, pars, i_io, 'ca')
    
    resid[ptr.phi_dl_ca] = SVdot[ptr.phi_dl_ca] - i_dl/pars.C_dl_ca
    
    resid[ptr.phi_ca] = SV[ptr.phi_ca] - pars.V_cell

def read_state(SV, ptr, pars, domain, j):
    if domain=='sl_an':
        # No ionomer, so no eta:
        eta = 0

        # Read out gas phase molar concentrations:
        C_k = SV[ptr.C_k_an_GDL]

        # Load the GDL state as s1
        state = {'C_k': C_k, 'dy':pars.dy_GDL, 'eps_g':pars.eps_g_GDL, 
            'n_Brugg':pars.n_Brugg_GDL, 'd_solid':pars.d_solid_GDL}
        
    elif domain=='al_an':
        # Calculate overpotential.
        # SHOULD BE A FUNCTION OF SPECIES ACTIVITIES:
        eta = SV[ptr.phi_dl_an[j]] - pars.delta_Phi_eq_an 

        # Read out gas phase molar concentrations:
        C_k = SV[ptr.C_k_an_CL[j][:]]
        # Load the GDL state as s1
        state = {'C_k': C_k, 'dy':pars.dy_CL, 'eps_g':pars.eps_g_CL, 
            'n_Brugg':pars.n_Brugg_CL, 'd_solid':pars.d_solid_CL}

    return state, eta

def chem_calcs(eta, pars, i_io, domain):
    
    if domain == 'an':
        # Butler-Vollmer equation:
        # i_o SHOULD VARY WITH SPECIES ACTIVITIES:
        i_Far = pars.i_o_an*(exp(-pars.n_an*F*pars.beta_an*eta/R/pars.T)
                        - exp(pars.n_an*F*(1-pars.beta_an)*eta/R/pars.T))

        # Molar production rates due to Faradaic current:
        sdot_k = i_Far*pars.nu_k_an/pars.n_an/F

        # Double layer current density (per unit area surface)
        i_dl = -i_io*pars.A_fac_dl - i_Far*pars.f_Pt


    elif domain=='ca':
        i_Far = -pars.i_o_ca*(exp(-pars.n_ca*F*pars.beta_ca*eta/R/pars.T)
                          - exp(pars.n_ca*F*(1-pars.beta_ca)*eta/R/pars.T))
        # Molar production rates due to Faradaic current:
        sdot_k = i_Far*pars.nu_k_ca/pars.n_ca/F

        # Double layer current density (per unit area surface)
        i_dl = i_io*pars.A_fac_dl - i_Far*pars.f_Pt
        
    return i_dl, sdot_k

def sofc_gas_flux(node1, node2, gas_props):
    # Initialize molar fluxes:
    N_k  = np.zeros_like(node1['C_k'])

    # Weighting fractions between current (GDL) and next (CL) node:
    f1 = node1['dy']/(node1['dy'] + node2['dy'])
    f2 = 1-f1

    C_int = f1*node1['C_k'] + f2*node2['C_k']

    X_k_1 = node1['C_k']/np.sum(node1['C_k'])
    X_k_2 = node2['C_k']/np.sum(node2['C_k'])
    X_k_int = f1*X_k_1 + f2*X_k_2

    P_1 = np.sum(node1['C_k'])*R*gas_props['T']
    P_2 = np.sum(node2['C_k'])*R*gas_props['T']
    # print(P_1, P_2)

    eps_g = f1*node1['eps_g'] + f2*node2['eps_g']
    tau_fac = (f1*node1['eps_g']**node1['n_Brugg'] 
        + f2*node2['eps_g']**node2['n_Brugg'])
    D_k_eff = eps_g*gas_props['D_k']/tau_fac
    
    # Distance between node centers:
    dY =0.5*(node1['dy'] + node2['dy'])

    # Note: this is not correct.  It is a temporary placeholder until after HW 
    # 6 is due :)
    N_k = D_k_eff*(node1['C_k'] - node2['C_k'])/dY

    return N_k