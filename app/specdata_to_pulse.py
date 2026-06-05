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

# 1. 读取光谱数据，并计算中心波长和光谱半高宽度，打印结果

# 1.1. 读取光谱数据，并把波长单位从nm转换为m，把强度归一化到0-1范围
specdata = pd.read_excel(setting['fname'])
wavelength = specdata['wavelength (nm)'].values * 1E-9
wavelength_intensity = specdata['intensity (a.u.)'].values
wavelength_intensity = (wavelength_intensity - wavelength_intensity.min()) \
    / (wavelength_intensity.max() - wavelength_intensity.min())

# 1.2. 计算中心波长和光谱半高宽度
center_wavelength = np.sum(wavelength*wavelength_intensity)/np.sum(wavelength_intensity)
wavelength_fwhm = wavelength[np.where(wavelength_intensity >= np.max(wavelength_intensity)/2)[0][-1]] \
    - wavelength[np.where(wavelength_intensity >= np.max(wavelength_intensity)/2)[0][0]]
print(f"Center wavelength: {center_wavelength*1E9:.2f} nm")
print(f"Wavelength FWHM: {wavelength_fwhm*1E9:.2f} nm")

# 2. 根据色散参数计算展宽后的脉宽全宽，打印结果
from util._260521_get_full_duration import _get_full_duration, _print_full_duration
full_duration = _get_full_duration(
    center_wavelength=center_wavelength,
    wavelength_fwhm=wavelength_fwhm,
    GVD=setting['GVD (s^2)'], TOD=setting['TOD (s^3)'], FOD=setting['FOD (s^4)']
    )

_print_full_duration(full_duration, wavelength_fwhm)

# 3. 计算频谱中心角频率、带宽、以及最优采样点数
from util._260523_ifft import _get_angular_frequency, _get_bandwidth, \
    _1d_gaus_fwhm, _1d_sech2_fwhm, _get_optimal_N, _ifft
w0, w_fwhm = _get_angular_frequency(
    center_wavelength=center_wavelength,
    wavelength_fwhm=wavelength_fwhm
    )
bandwidth = _get_bandwidth(
    center_wavelength=center_wavelength,
    wavelength_fwhm=wavelength_fwhm
    )
time, w = _get_optimal_N(full_duration, bandwidth)[4:6]
w = w + w0

# 4. 构造频谱强度和相位，计算时域电场和强度
from scipy.interpolate import CubicSpline, interp1d

wl = wavelength
I_wl = wavelength_intensity

# 4.1. 如果w的长度大于原有数据的长度，则需要对原数据进行插值延长
if len(w) >= len(wl):
    pad_total = len(w) - len(wl)
    pad_left = pad_total // 2
    pad_right = pad_total - pad_left

    # 4.1.1. 把原有wavelength数据两边线性插值延长
    f = interp1d(np.arange(len(wl)), wl, kind='linear', fill_value='extrapolate')
    wl_new = np.arange(-pad_left, len(wl) + pad_right)
    wl_padded = f(wl_new)

    # 4.1.2. 把原有wavelength_intensity数据两边加零延长
    I_wl_padded = np.pad(I_wl, (pad_left, pad_right), 'constant', constant_values=0)

    # 4.1.3. 按照w的长度生成新的光谱数据
    spline = CubicSpline(np.flip(2*np.pi*3E8/(wl_padded)), np.flip(I_wl_padded), bc_type='natural')
    spectrum_intensity = spline(w)

# 4.2. 如果w的长度小于原有数据的长度，则直接按照w的长度从原有数据中选取对应的强度值    
else:
    indices = np.linspace(0, len(wavelength_intensity) - 1, len(w)).astype(int)
    spectrum_intensity = wavelength_intensity[indices]

# 5. 进行傅里叶逆变换
phi_omega = setting['GVD (s^2)']*(w-w0)**2/2 + setting['TOD (s^3)']*(w-w0)**3/6 + \
    setting['FOD (s^4)']*(w-w0)**4/24
E_time, I_time, I_time_normalized = _ifft(spectrum_intensity, phi_omega)

# 6. 绘制光谱数据和对应的时域脉冲
from util._260524_mpl import mark_pulse_fwhm, mark_spectrum_fwhm
fig, axes = plt.subplots(2, 1, figsize=(10, 6), tight_layout=True)

ax1 = axes[0]
ax1.plot(np.flip(2*np.pi*3E8/w)*1E9, np.flip(spectrum_intensity), 'b', label='Spectrum Intensity')
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Intensity (a.u.)')
ax1.set_title('Input Spectrum')
ax1.grid()
mark_spectrum_fwhm(ax1, np.flip(2*np.pi*3E8/w)*1E9, np.flip(spectrum_intensity))

ax1_phase = ax1.twinx()
ax1_phase.plot(np.flip(2*np.pi*3E8/w)*1E9, np.flip(phi_omega), 'r--', label='Spectral Phase')
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