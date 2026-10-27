import os
import math
import numpy as np
import matplotlib.pyplot as plt

input_path = 'data/lab2/input_2_2.txt'
eps = 1e-4
x1_0, x2_0 = 1.2, 1.5

if os.path.exists(input_path):
    with open(input_path, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]
        if len(lines) >= 2:
            eps = float(lines[0])
            parts = lines[1].split()
            x1_0, x2_0 = float(parts[0]), float(parts[1])
        elif len(lines) == 1:
            parts = lines[0].split()
            eps = float(parts[0])
            x1_0, x2_0 = float(parts[1]), float(parts[2])

def sys_f1(x1, x2):
    return x1**2 - 2.0 * np.log10(x2) - 1.0

def sys_f2(x1, x2):
    return x1**2 - 2.0 * x1 * x2 + 2.0

def phi1(x2):
    return np.sqrt(2.0 * np.log10(x2) + 1.0)

def phi2(x1):
    return (x1**2 + 2.0) / (2.0 * x1)

mpi_traj = [(x1_0, x2_0)]
mpi_diffs = []
x1_c, x2_c = x1_0, x2_0

for k in range(1, 1000):
    x1_p, x2_p = x1_c, x2_c
    x1_c = float(phi1(x2_p))
    x2_c = float(phi2(x1_p))
    diff = max(abs(x1_c - x1_p), abs(x2_c - x2_p))
    mpi_traj.append((x1_c, x2_c))
    mpi_diffs.append(diff)
    if diff <= eps:
        break

newton_traj = [(x1_0, x2_0)]
newton_diffs = []
x1_c, x2_c = x1_0, x2_0

for k in range(1, 1000):
    x1_p, x2_p = x1_c, x2_c

    j00 = 2.0 * x1_p
    j01 = -2.0 / (x2_p * np.log(10.0))
    j10 = 2.0 * x1_p - 2.0 * x2_p
    j11 = -2.0 * x1_p
    
    f1 = sys_f1(x1_p, x2_p)
    f2 = sys_f2(x1_p, x2_p)
    
    det = j00 * j11 - j01 * j10
    dx1 = (-f1 * j11 - (-f2) * j01) / det
    dx2 = (j00 * (-f2) - j10 * (-f1)) / det
    
    x1_c = x1_p + dx1
    x2_c = x2_p + dx2
    diff = max(abs(dx1), abs(dx2))
    
    newton_traj.append((x1_c, x2_c))
    newton_diffs.append(diff)
    if diff <= eps:
        break

root_x1, root_x2 = newton_traj[-1]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

x1_grid = np.linspace(0.8, 1.6, 400)
curve1_x2 = 10.0 ** ((x1_grid**2 - 1.0) / 2.0)
curve2_x2 = (x1_grid**2 + 2.0) / (2.0 * x1_grid)

ax1.plot(x1_grid, curve1_x2, label=r'$x_1^2 - 2\lg x_2 - 1 = 0$', color='#1f77b4', linewidth=2.2)
ax1.plot(x1_grid, curve2_x2, label=r'$x_1^2 - 2 x_1 x_2 + 2 = 0$', color='#ff7f0e', linewidth=2.2)

nt_x1 = [p[0] for p in newton_traj]
nt_x2 = [p[1] for p in newton_traj]
ax1.plot(nt_x1, nt_x2, 'o--', color='#2ca02c', label='Шаги метода Ньютона', linewidth=1.5, markersize=5)

ax1.scatter([x1_0], [x2_0], color='purple', s=90, zorder=5, label=f'Старт ({x1_0}, {x2_0})')
ax1.scatter([root_x1], [root_x2], color='red', marker='*', s=220, zorder=6,
            label=f'Решение: ({root_x1:.4f}, {root_x2:.4f})')

ax1.set_title('Графическое отделение корней системы (Задача 2.2)', fontsize=13, fontweight='bold')
ax1.set_xlabel('$x_1$', fontsize=12)
ax1.set_ylabel('$x_2$', fontsize=12)
ax1.set_xlim(0.8, 1.6)
ax1.set_ylim(0.8, 2.5)
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(fontsize=10, loc='upper left')

k_mpi = range(1, len(mpi_diffs) + 1)
k_newton = range(1, len(newton_diffs) + 1)

ax2.plot(k_mpi, mpi_diffs, marker='o', markersize=6, label=f'МПИ (итераций: {len(mpi_diffs)})', color='#d62728', linewidth=1.8)
ax2.plot(k_newton, newton_diffs, marker='s', markersize=7, label=f'Метод Ньютона (итераций: {len(newton_diffs)})', color='#2ca02c', linewidth=2.0)
ax2.axhline(eps, color='gray', linestyle='--', linewidth=1.2, label=f'Порог $\\varepsilon = {eps}$')

ax2.set_yscale('log')
ax2.set_title('Зависимость погрешности от номера итерации', fontsize=13, fontweight='bold')
ax2.set_xlabel('Номер итерации (k)', fontsize=11)
ax2.set_ylabel(r'Погрешность $\max |x_i^{(k)} - x_i^{(k-1)}|$ (log-scale)', fontsize=11)
ax2.grid(True, which='both', linestyle=':', alpha=0.6)
ax2.legend(fontsize=10)

plt.tight_layout()

os.makedirs('data/lab2', exist_ok=True)
out_png = 'data/lab2/plot_2_2.png'
plt.savefig(out_png, dpi=300)
print(f"График сохранён в: {out_png}")

plt.show()
