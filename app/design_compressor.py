"""
Environment: numpy, sympy
"""

"""Part I: Parameter Setting"""

setting = {

'center_wavelength (nm)':2060,
'grating_line_density (mm-1)': 890,
'incident_angle (deg)': 62.77,
'grating_separation (mm)': 1000,
'author': 'Hua Lin',
'date': '2026-5-20',

}

"""Part II: Run"""
from util._260520_get_dispersion import _get_dispersion
dispersion = _get_dispersion(setting)
