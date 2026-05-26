import os
import numpy as np
from datetime import date, datetime
from sympy import symbols, sin, latex, Eq, pi
from util._260521_get_full_duration import _get_full_duration

# 利用光谱宽度计算频谱宽度
def _dwl_to_dv(center_wavelength, wavelength_fwhm):
    return 2.998E8/(center_wavelength-0.5*wavelength_fwhm)-\
            2.998E8/(center_wavelength+0.5*wavelength_fwhm)

# 将计算结果写入LaTeX文件，并编译成PDF
def _latex_comp_results(setting, dispersion):

    # 代入初始参数
    c_value = 2.998E8
    wl0_value = setting['center_wavelength (nm)'] * 1E-9
    d_value = 1/setting['grating_line_density (mm-1)'] * 1E-3
    gamma_value = setting['incident_angle (deg)']/180*np.pi
    Lg_value = setting['grating_separation (mm)'] * 1E-3

    # 代入色散的计算结果
    w0_value = dispersion[0]
    littrow_angle = dispersion[1]
    diffration_angle = dispersion[2]
    gvd_value = dispersion[3]
    tod_value = dispersion[4]
    fod_value = dispersion[5]

    # 计算不同光谱宽度下的频谱宽度、变换极限脉宽和拉伸后的脉宽
    dwl_val = np.array([1E-9, 2E-9, 5E-9, 10E-9, 20E-9, 50E-9, 100E-9])
    dv_val = _dwl_to_dv(wl0_value, dwl_val)
    dt_val = 0.441 / dv_val
    T_val = [_get_full_duration(wl0_value, dwl, gvd_value, tod_value,
        fod_value) for dwl in dwl_val]

    # 生成LaTeX文件并编译成PDF
    fname = str(datetime.now().strftime("%Y%m%d%H%M%S")) + '_comp_results'
    output_root=os.path.join('data', fname)
    tex_file = open(output_root + '.tex', 'w')

    document_beginings = r"""
    \documentclass{article}
    \usepackage{amsmath}
    \usepackage{amssymb}
    \title{Compressor Design}
    \author{%s}
    \date{%s}
    \begin{document}
    \maketitle
    """ % (setting['author'], str(date.today().strftime("%Y-%m-%d")))

    section1 = '\\section{Grating pair compressor equations}'

    PHI, wl, gamma, Lg, d, c, w = \
        symbols('PHI, \\lambda, \\gamma, L_g, d, c, \\omega')
    phi_expr = 2*w*Lg/c*(1-(2*pi*c/w/d - sin(gamma))**2)**0.5
    PHI_eq = Eq(PHI, phi_expr)

    section2 = '\section{Grating parameters}'

    text2_1 = '\[ \lambda_0 = %d \: nm \]' % \
        setting['center_wavelength (nm)']

    text2_2 = r'\[ d = \frac{1}{%d} \: mm \]' % \
        setting['grating_line_density (mm-1)']

    text2_3 = '\[ L_g = %d \: mm\]' % \
        setting['grating_separation (mm)']

    text2_4 = '\[ Incident \: angle \: \gamma = %.2f \: ^{\circ} \]' % \
        setting['incident_angle (deg)']

    text2_5 = '\[ c = %.0f \cdot 10^{11} \: mm/s \]' % \
        (c_value / 1E8)

    section3 = '\section{Results}'

    table3_1 = '\\begin{center}\\begin{tabular}[t]{|l|r|} \
    \hline \
    $\omega_0$ (Hz) & $%.3E$ \\\ \
    \hline \
    Littrow angle $(^{\circ})$ & $%.2f$ \\\ \
    \hline \
    Diffration angle $(^{\circ})$ & $%.2f$ \\\ \
    \hline \
    GVD $(s^2)$ & $%.3E$ \\\ \
    \hline \
    TOD $(s^3)$ & $%.3E$ \\\ \
    \hline \
    FOD $(s^4)$ & $%.3E$ \\\ \
    \hline \
    \end{tabular}\end{center} \
    ' % (w0_value,
        littrow_angle,
        diffration_angle,
        gvd_value, tod_value, fod_value)

    table3_2 = '\\begin{center}\\begin{tabular}[t]{|r|r|r|} \
    \hline \
    Bandwidth (nm) & Transform-lim-Gau (s) & Stretched full width (s) \\\ \
    \hline \
    %d & $%.2E$ & $%.2E$ \\\ \
    \hline \
    %d & $%.2E$ & $%.2E$ \\\ \
    \hline \
    %d & $%.2E$ & $%.2E$ \\\ \
    \hline \
    %d & $%.2E$ & $%.2E$ \\\ \
    \hline \
    %d & $%.2E$ & $%.2E$ \\\ \
    \hline \
    %d & $%.2E$ & $%.2E$ \\\ \
    \hline \
    %d & $%.2E$ & $%.2E$ \\\ \
    \hline \
    \end{tabular}\end{center} \
    ' % (
        dwl_val[0]*1E9, dt_val[0], T_val[0],
        dwl_val[1]*1E9, dt_val[1], T_val[1],
        dwl_val[2]*1E9, dt_val[2], T_val[2],
        dwl_val[3]*1E9, dt_val[3], T_val[3],
        dwl_val[4]*1E9, dt_val[4], T_val[4],
        dwl_val[5]*1E9, dt_val[5], T_val[5],
        dwl_val[6]*1E9, dt_val[6], T_val[6],
        )

    strings = [
    document_beginings,

    section1,
    latex(PHI_eq, mode='equation', fold_func_brackets=True),

    section2,
    text2_1,
    text2_2,
    text2_3,
    text2_4,
    text2_5,

    section3,
    table3_1,
    table3_2,

    '\end{document}',
    ]

    latex_str = ' '.join(strings)

    tex_file.write(latex_str)
    tex_file.close()

    # 编译LaTeX文件成PDF，并删除中间文件
    # 注意：>nul 2>&1 是为了在Windows系统上隐藏编译过程中的输出，不适配Linux系统和MacOS系统
    os.system('pdflatex -interaction=batchmode -output-directory=' + \
        'data\\' + ' ' + os.path.split(output_root)[1] + '.tex > nul 2>&1')
    os.remove(output_root + '.log')
    os.remove(output_root + '.aux')
    os.remove(output_root + '.tex')

    print("")
    print("Results saved as: %s.pdf" % fname)

