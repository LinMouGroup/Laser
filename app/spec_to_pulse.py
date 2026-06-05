"""
Environment: numpy, scipy, matplotlib
"""

"""Part I: Parameter Setting"""

setting = {

'spectrum_shape': 'gaus', #gaus, sech2
'center_wavelength (nm)': 2000,
'wavelength_fwhm (nm)': 20,
'GVD (s^2)': 134E-25,
'TOD (s^3)': -3.87E-37,
'FOD (s^4)': 0,

}

"""Part II: Run"""
import numpy as np
import matplotlib.pyplot as plt

# 1. 根据色散参数计算展宽后的脉宽全宽
from util._260521_get_full_duration import _get_full_duration, _print_full_duration
full_duration = _get_full_duration(
    center_wavelength=setting['center_wavelength (nm)']*1E-9,
    wavelength_fwhm=setting['wavelength_fwhm (nm)']*1E-9,
    GVD=setting['GVD (s^2)'], TOD=setting['TOD (s^3)'], FOD=setting['FOD (s^4)']
    )

# 2. 计算频谱中心角频率、带宽、以及最优采样点数
from util._260523_ifft import _get_angular_frequency, _get_bandwidth, \
    _1d_gaus_fwhm, _1d_sech2_fwhm, _get_optimal_N, _ifft
w0, w_fwhm = _get_angular_frequency(
    center_wavelength=setting['center_wavelength (nm)']*1E-9,
    wavelength_fwhm=setting['wavelength_fwhm (nm)']*1E-9
    )
bandwidth = _get_bandwidth(
    center_wavelength=setting['center_wavelength (nm)']*1E-9,
    wavelength_fwhm=setting['wavelength_fwhm (nm)']*1E-9
    )
time, w = _get_optimal_N(full_duration, bandwidth)[4:6]

# 3. 构造频谱强度和相位，计算时域电场和强度
w = w + w0

if setting['spectrum_shape'] == 'sech2':
    spectrum_intensity = _1d_sech2_fwhm(w, 1, w0, w_fwhm)
else: # 默认使用高斯光谱
    spectrum_intensity = _1d_gaus_fwhm(w, 1, w0, w_fwhm)

wavelength = np.flip(2*np.pi*3E8/w)

phi_omega = setting['GVD (s^2)']*(w-w0)**2/2 + setting['TOD (s^3)']*(w-w0)**3/6 + \
    setting['FOD (s^4)']*(w-w0)**4/24
E_time, I_time, I_time_normalized = _ifft(spectrum_intensity, phi_omega)

# 4. 绘制光谱数据和对应的时域脉冲
from util._260524_mpl import mark_pulse_fwhm, mark_spectrum_fwhm

fig, axes = plt.subplots(2, 1, figsize=(10, 6), tight_layout=True)

ax1 = axes[0]
ax1.plot(wavelength*1E9, np.flip(spectrum_intensity), 'b', label='Spectrum Intensity')
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Intensity (a.u.)')
ax1.set_title('Input Spectrum')
ax1.grid()
mark_spectrum_fwhm(ax1, wavelength*1E9, np.flip(spectrum_intensity))

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
mark_pulse_fwhm(ax2, offset_time, I_time_normalized)

plt.suptitle('Spectrum to Pulse Conversion', fontsize=16)
print("")
print("Results as shown in matplotlib figure.")
plt.show()