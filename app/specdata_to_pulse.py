"""
Environment: numpy, scipy, pandas, matplotlib
"""

"""Part I: Parameter Setting"""

setting = {

'fname': 'data//specdata.xlsx',
'GVD (s^2)': 134E-25,
'TOD (s^3)': -3.87E-37,
'FOD (s^4)': 0,

}

"""Part II: Run"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

specdata = pd.read_excel(setting['fname'])
wavelength = specdata['wavelength (nm)'].values
wavelength_intensity = specdata['intensity (a.u.)'].values
center_wavelength = np.sum(wavelength*wavelength_intensity)/np.sum(wavelength_intensity)
wavelength_fwhm = wavelength[np.where(wavelength_intensity >= np.max(wavelength_intensity)/2)[0][-1]] \
    - wavelength[np.where(wavelength_intensity >= np.max(wavelength_intensity)/2)[0][0]]

print(f"Center wavelength: {center_wavelength:.2f} nm")
print(f"Wavelength FWHM: {wavelength_fwhm:.2f} nm")


# 1. 根据色散参数计算展宽后的脉宽全宽
from util._260521_get_full_duration import _get_full_duration, _print_full_duration
full_duration = _get_full_duration(
    center_wavelength=center_wavelength*1E-9,
    wavelength_fwhm=wavelength_fwhm*1E-9,
    GVD=setting['GVD (s^2)'], TOD=setting['TOD (s^3)'], FOD=setting['FOD (s^4)']
    )

_print_full_duration(full_duration, wavelength_fwhm)

# 2. 计算频谱中心角频率、带宽、以及最优采样点数
from util._260523_ifft import _get_angular_frequency, _get_bandwidth, \
    _1d_gaus_fwhm, _1d_sech2_fwhm, _get_optimal_N, _ifft
w0, w_fwhm = _get_angular_frequency(
    center_wavelength=center_wavelength*1E-9,
    wavelength_fwhm=wavelength_fwhm*1E-9
    )
bandwidth = _get_bandwidth(
    center_wavelength=center_wavelength*1E-9,
    wavelength_fwhm=wavelength_fwhm*1E-9
    )
time, w = _get_optimal_N(full_duration, bandwidth)[4:6]

# 3. 构造频谱强度和相位，计算时域电场和强度
from scipy.interpolate import CubicSpline, UnivariateSpline
from scipy.interpolate import interp1d

w = w + w0

wl = wavelength
I_wl = wavelength_intensity

pad_total = len(w) - len(wl)
pad_left = pad_total // 2
pad_right = pad_total - pad_left


f = interp1d(np.arange(len(wl)), wl, kind='linear', fill_value='extrapolate')
x_new = np.arange(-pad_left, len(wl) + pad_right)
wl_padded = f(x_new)

I_wl_padded = np.pad(I_wl, (pad_left, pad_right), 'constant', constant_values=0)

spline = CubicSpline(np.flip(2*np.pi*3E8/(wl_padded*1E-9)), np.flip(I_wl_padded), bc_type='natural')
spectrum_intensity = spline(w)

wavelength = np.flip(2*np.pi*3E8/w)

phi_omega = setting['GVD (s^2)']*(w-w0)**2/2 + setting['TOD (s^3)']*(w-w0)**3/6 + \
    setting['FOD (s^4)']*(w-w0)**4/24
E_time, I_time, I_time_normalized = _ifft(spectrum_intensity, phi_omega)

# 4. 绘制光谱数据和对应的时域脉冲
fig, axes = plt.subplots(2, 1, figsize=(10, 6), tight_layout=True)

ax1 = axes[0]
ax1.plot(wavelength*1E9, np.flip(spectrum_intensity), 'b', label='Spectrum Intensity')
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Intensity (a.u.)')
ax1.set_title('Input Spectrum')
ax1.grid()

ax1_phase = ax1.twinx()
ax1_phase.plot(wavelength*1E9, np.flip(phi_omega), 'r--', label='Spectral Phase')
ax1_phase.set_ylabel('Phase (rad)', color='r')
ax1_phase.tick_params(axis='y', color='r')

# 计算峰值位置并将时间轴调整为以峰值为中心
peak_I_time_index = np.argmax(I_time_normalized)
peak_I_time = time[peak_I_time_index]
offset_time = time - peak_I_time

ax2 = axes[1]
ax2.plot(offset_time, I_time_normalized, 'k')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Intensity (a.u.)')
ax2.set_title('Output Pulse')
ax2.grid()

plt.suptitle('Spectrum to Pulse Conversion', fontsize=16)
plt.show()