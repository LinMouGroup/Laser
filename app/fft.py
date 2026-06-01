import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import ifft, fftshift, ifftshift

# ========== 参数设置 ==========
# 频率网格参数
N = 2048                 # 采样点数（频域和时域相同）
lambda_center = 800e-9   # 中心波长 (m)
c = 3e8                  # 光速 (m/s)
freq_center = c / lambda_center  # 中心频率 (Hz)
bandwidth = 20e12        # 带宽 (Hz) 约 20 THz

# 构建频率数组
df = bandwidth / N       # 频率分辨率
freq = np.linspace(-bandwidth/2, bandwidth/2, N) + freq_center

# ========== 1. 定义频谱强度 I(ω) = |E(ω)|^2 ==========
# 例如高斯型光谱
spectrum_intensity = np.exp(-((freq - freq_center)**2 / (2*(bandwidth/4)**2)))

# ========== 2. 定义频谱相位 φ(ω) ==========
# 二阶色散（啁啾）项： φ(ω) = a * (ω - ω0)^2
chirp_coeff = 2e-28      # 啁啾系数 (s^2/rad^2)
phi_omega = chirp_coeff * (freq - freq_center)**2

# 可选：添加线性啁啾（频域线性相位，等效时域平移）
# phi_omega += 0.0 * (freq - freq_center)   # 线性项影响时间轴位置

# ========== 3. 构建复频谱 E(ω) = sqrt(I(ω)) * exp(iφ(ω)) ==========
E_omega = np.sqrt(spectrum_intensity) * np.exp(1j * phi_omega)

# ========== 4. 傅里叶逆变换得到时域电场 E(t) ==========
# 注意：频域需进行 ifftshift 使零频位于两端，以适应 ifft 的输入要求
E_time = ifft(ifftshift(E_omega))

# 时域坐标
dt = 1 / (N * df)        # 时间分辨率
time = np.linspace(-N//2, N//2, N) * dt

# 电场包络（振幅）
envelope = np.abs(E_time)
# 相位（注意unwrap避免跳变）
instant_phase = np.unwrap(np.angle(E_time))

# 时域强度（正比于电场平方）
intensity_time = envelope**2

# ========== 5. 绘图 ==========
fig, axes = plt.subplots(3, 1, figsize=(10, 8), tight_layout=True)

# 图1: 频谱强度 I(ω) 和相位 φ(ω)
ax1 = axes[0]
ax1.plot(freq*1e-12, spectrum_intensity, 'b', label='Spectrum Intensity')
ax1.set_xlabel('Frequency (THz)')
ax1.set_ylabel('Intensity (a.u.)')
ax1.set_title('Spectral Intensity')
ax1.grid(True)

ax1_phase = ax1.twinx()
ax1_phase.plot(freq*1e-12, phi_omega, 'r--', label='Spectral Phase')
ax1_phase.set_ylabel('Phase (rad)', color='r')
ax1_phase.tick_params(axis='y', labelcolor='r')
ax1_phase.legend(loc='upper right')
ax1.legend(loc='upper left')

# 图2: 时域脉冲包络 |E(t)|
ax2 = axes[1]
ax2.plot(time*1e15, envelope, 'g')
ax2.set_xlabel('Time (fs)')
ax2.set_ylabel('Amplitude (a.u.)')
ax2.set_title('Pulse Envelope |E(t)|')
ax2.grid(True)

# 图3: 时域强度（可选）及瞬时频率（相位时间导数）
ax3 = axes[2]
ax3.plot(time*1e15, intensity_time, 'm', label='Intensity I(t)')

# 瞬时频率 ω(t) = dφ(t)/dt （可显示啁啾特性）
inst_freq = np.gradient(instant_phase, dt) / (2*np.pi)  # 转换为 Hz
ax3_f = ax3.twinx()
ax3_f.plot(time*1e15, inst_freq*1e-12, 'k--', label='Instantaneous frequency (THz)')
ax3_f.set_ylabel('Freq. (THz)', color='k')
ax3_f.tick_params(axis='y', labelcolor='k')
ax3_f.legend(loc='upper right')

ax3.set_xlabel('Time (fs)')
ax3.set_ylabel('Intensity (a.u.)')
ax3.set_title('Temporal Intensity and Instantaneous Frequency')
ax3.grid(True)

plt.suptitle('Chirped Laser Pulse Shape from Spectral Amplitude & Phase')
plt.show()