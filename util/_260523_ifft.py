import numpy as np
from scipy.fft import ifft, fftshift

def _1d_gaus(x, a, b, c):
    return a*np.exp(-(x-b)**2/(2*c**2))

def _1d_gaus_fwhm(x, a, b, fwhm):
    return a*np.exp(-(x-b)**2/(2*(fwhm/2.35482)**2))

# 计算中心角频率
def _get_angular_frequency(center_wavelength, wavelength_fwhm):
    c_value = 2.9979E8
    freq_center = c_value / center_wavelength
    omega_center = 2 * np.pi * freq_center
    omega_fwhm = 2 * np.pi * (c_value / (center_wavelength - wavelength_fwhm) \
        - c_value / (center_wavelength + wavelength_fwhm))
    return omega_center, omega_fwhm

# 计算带宽
def _get_bandwidth(center_wavelength, wavelength_fwhm):
    c_value = 2.9979E8
    freq_center = c_value / center_wavelength
    freq_min = c_value / (center_wavelength + wavelength_fwhm)
    freq_max = c_value / (center_wavelength - wavelength_fwhm)
    bandwidth = freq_max - freq_min
    return bandwidth

# 选择最优的采样点数N以满足IFFT的性能要求
def _get_optimal_N(pulse_duration, bandwidth, safety_factor=5):

    # 时间窗口应该足够宽以包含整个脉冲
    time_window = safety_factor * pulse_duration
    
    # 根据带宽确定所需时间分辨率（奈奎斯特采样）
    dt_needed = 1 / (2 * bandwidth)
    
    # 计算最小采样点数
    N_min = int(time_window / dt_needed)
    
    # 取2的幂次以优化FFT
    N_optimal = 2 ** int(np.ceil(np.log2(N_min)))
    
    # 输出其他参数
    dt_optimal = time_window / N_optimal
    df_optimal = 1 / time_window
    dw_optimal = 2 * np.pi * df_optimal
    time = np.linspace(-N_optimal//2, N_optimal//2, N_optimal)*dt_optimal
    w = np.linspace(-N_optimal//2, N_optimal//2, N_optimal)*dw_optimal
    
    return N_optimal, dt_optimal, df_optimal, dw_optimal, time, w

# 根据频谱强度和相位计算时域电场和强度
def _ifft(spectrum_intensity, phi_omega):

    # 执行IFFT以获取时域信号
    E_time = fftshift(ifft(fftshift(spectrum_intensity*np.exp(1j*phi_omega))))
    I_time = np.abs(E_time)**2
    I_time_normalized = I_time / np.max(I_time)

    return E_time, I_time, I_time_normalized