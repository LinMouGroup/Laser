import numpy as np

def mark_pulse_fwhm(ax, x, y):
    # 找到最大值和半高值
    max_val = np.max(y)
    half_max = max_val / 2
    
    # 找到半高处的左右索引
    above_half = y >= half_max
    left_idx = np.argmax(above_half)
    right_idx = len(y) - 1 - np.argmax(above_half[::-1])
    
    # 计算FWHM
    fwhm = x[right_idx] - x[left_idx]
    
    # 左箭头：从左侧外部指向半高位置
    ax.annotate('', xy=(x[left_idx], half_max), xytext=(x[left_idx] - fwhm*0.3, half_max),
                arrowprops=dict(arrowstyle='->', lw=2, alpha=0.8))
    
    # 右箭头：从右侧外部指向半高位置
    ax.annotate('', xy=(x[right_idx], half_max), xytext=(x[right_idx] + fwhm*0.3, half_max),
                arrowprops=dict(arrowstyle='->', lw=2, alpha=0.8))
    
    # 添加文字标签（放在中间）
    ax.text(x[right_idx], half_max + max_val*0.05, 
            f'{fwhm:.2E} s', ha='left', fontsize=11)
    
    return fwhm


def mark_spectrum_fwhm(ax, x, y):
    # 找到最大值和半高值
    max_val = np.max(y)
    half_max = max_val / 2
    
    # 找到半高处的左右索引
    above_half = y >= half_max
    left_idx = np.argmax(above_half)
    right_idx = len(y) - 1 - np.argmax(above_half[::-1])
    
    # 计算FWHM
    fwhm = x[right_idx] - x[left_idx]
    
    # 左箭头：从左侧外部指向半高位置
    ax.annotate('', xy=(x[left_idx], half_max), xytext=(x[left_idx] - fwhm*0.3, half_max),
                arrowprops=dict(arrowstyle='->', color='b', lw=2, alpha=0.8))
    
    # 右箭头：从右侧外部指向半高位置
    ax.annotate('', xy=(x[right_idx], half_max), xytext=(x[right_idx] + fwhm*0.3, half_max),
                arrowprops=dict(arrowstyle='->', color='b', lw=2, alpha=0.8))
    
    # 添加文字标签（放在中间）
    ax.text(x[right_idx], half_max + max_val*0.05,
            f'{fwhm:.2f} nm', ha='left', fontsize=11, color='b')
    
    return fwhm