import numpy as np

def _get_full_duration(
    center_wavelength,
    wavelength_fwhm,
    GVD=0,
    TOD=0,
    FOD=0,
    ):

    c_value = 2.9979E8

    w0 = 2*np.pi*c_value/center_wavelength
    w_min = 2*np.pi*c_value/(center_wavelength + wavelength_fwhm)
    w_max = 2*np.pi*c_value/(center_wavelength - wavelength_fwhm)

    t_min = GVD*(w_min - w0) + 1/2.0*TOD*(w_min - w0)**2 + 1/6.0*FOD*(w_min - w0)**3
    t_max = GVD*(w_max - w0) + 1/2.0*TOD*(w_max - w0)**2 + 1/6.0*FOD*(w_max - w0)**3

    return np.abs(t_max - t_min)

def _print_full_duration(full_duration):
    print("")
    print(f"Full_pulse_duration_with_10nm_FWHM: {full_duration:.2E} s")
    print("")

"""
Reference:
Kapteyn, H. (1998). High power ultrafast lasers. Review of Scientific Instruments. 
https://doi.org/10.1063/1.1148795
"""