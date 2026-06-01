"""
Environment: numpy, sympy，LaTeX compiler (e.g., pdflatex)
"""

"""Part I: Parameter Setting"""

setting = {

'center_wavelength (nm)':2060,
'grating_line_density (mm-1)': 890,
'incident_angle (deg)': 62.77,
'grating_separation (mm)': 1000,
'author': 'Hua Lin',

}

"""Part II: Run"""
from util._260520_get_dispersion import _get_dispersion, _print_dispersion
dispersion = _get_dispersion(
    center_wavelength=setting['center_wavelength (nm)'] * 1E-9,
    grating_line_density=setting['grating_line_density (mm-1)'] * 1E3,
    incident_angle=setting['incident_angle (deg)'],
    grating_separation=setting['grating_separation (mm)'] * 1E-3
    )
_print_dispersion(dispersion)

from util._260521_get_full_duration import _get_full_duration, _print_full_duration
full_duration = _get_full_duration(
    center_wavelength=setting['center_wavelength (nm)'] * 1E-9,
    wavelength_fwhm=10 * 1E-9,
    GVD=dispersion[3], TOD=dispersion[4], FOD=dispersion[5]
    )
_print_full_duration(full_duration)

from util._260522_latex_comp_results import _latex_comp_results
_latex_comp_results(setting, dispersion)