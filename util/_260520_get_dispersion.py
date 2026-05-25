import os
import numpy as np
from sympy import symbols, sin, asin, diff

def _get_dispersion(setting):

    c_value = 2.9979E8
    wl0_value = setting['center_wavelength (nm)'] * 1E-9
    d_value = 1/setting['grating_line_density (mm-1)'] * 1E-3
    gamma_value = setting['incident_angle (deg)']/180*np.pi
    Lg_value = setting['grating_separation (mm)'] * 1E-3
    w0_value = 2*np.pi*c_value/wl0_value

    littrow_angle = asin(wl0_value/2/d_value)/np.pi*180
    diffration_angle = asin(wl0_value/d_value - sin(gamma_value))/np.pi*180

    w, c, wl, d, gamma, Lg= symbols('\\omega, c, \\lambda, d, \\gamma, L_g')
    phi_expr = 2*w*Lg/c*(1-(2*np.pi*c/w/d - sin(gamma))**2)**0.5
    phi_expr2 = phi_expr.subs([
        (d, d_value),
        (Lg, Lg_value),
        (c, c_value),
        (gamma, gamma_value),
        ])

    gvd_value = (diff(phi_expr2, w, 2)).subs(w, w0_value)
    tod_value = (diff(phi_expr2, w, 3)).subs(w, w0_value)
    fod_value = (diff(phi_expr2, w, 4)).subs(w, w0_value)

    return w0_value, littrow_angle, diffration_angle, gvd_value, tod_value, fod_value