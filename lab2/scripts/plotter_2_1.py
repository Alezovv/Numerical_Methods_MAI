import os
import math
import numpy as np
import matplotlib.pyplot as plt

input_path = 'data/lab2/input_2_1.txt'
eps = 1e-4
x0 = 0.5

if os.path.exists(input_path):
    with open(input_path, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]
        if len(lines) >= 2:
            eps = float(lines[0])
            x0 = float(lines[1])
        elif len(lines) == 1:
            parts = lines[0].split()
            eps = float(parts[0])
            x0 = float(parts[1])

def f(x):
    return 2.0**x + x**2 - 2.0

def df(x):
    return 2.0**x * np.log(2.0) + 2.0 * x

def phi(x):
    return np.sqrt(np.maximum(0.0, 2.0 - 2.0**x))

mpi_x = [x0]
mpi_diffs = []
x_curr = x0
for k in range(1, 1000):
    x_prev = x_curr
    x_curr = float(phi(x_prev))
    diff = abs(x_curr - x_prev)
    mpi_x.append(x_curr)
    mpi_diffs.append(diff)
    if diff <= eps:
        break

newton_x = [x0]
newton_diffs = []
x_curr = x0
for k in range(1, 1000):
    x_prev = x_curr
    x_curr = x_prev - f(x_prev) / df(x_prev)
    diff = abs(x_curr - x_prev)
    newton_x.append(x_curr)
    newton_diffs.append(diff)
    if diff <= eps:
        break

root = newton_x[-1]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

x_vals = np.linspace(-0.2, 1.2, 500)
ax1.plot(x_vals, f(x_vals), label=r'$f(x) = 2^x + x^2 - 2$', color='#1f77b4', linewidth=2.5)
ax1.plot(x_vals, x_vals**2, label=r'$y_1 = x^2$', color='#2ca02c', linestyle='--', linewidth=1.5)
ax1.plot(x_vals, 2.0 - 2.0**x_vals, label=r'$y_2 = 2 - 2^x$', color='#ff7f0e', linestyle='--', linewidth=1.5)

ax1.axhline(0, color='black', linewidth=0.8, linestyle='-')
ax1.axvline(0, color='black', linewidth=0.8, linestyle='-')

ax1.scatter([x0], [f(x0)], color='purple', s=80, zorder=5, label=f'Начальное $x_0 = {x0}$')
ax1.scatter([root], [0], color='red', marker='*', s=200, zorder=6, label=f'Корень $x^* \\approx {root:.4f}$')

ax1.set_title('Графическое отделение корня (Задача 2.1)', fontsize=13, fontweight='bold')
ax1.set_xlabel('x', fontsize=11)
ax1.set_ylabel('y', fontsize=11)
ax1.set_xlim(-0.2, 1.2)
ax1.set_ylim(-1.5, 2.5)
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(fontsize=10, loc='upper left')

k_mpi = range(1, len(mpi_diffs) + 1)
k_newton = range(1, len(newton_diffs) + 1)

ax2.plot(k_mpi, mpi_diffs, marker='o', markersize=4, label=f'МПИ (итераций: {len(mpi_diffs)})', color='#d62728', linewidth=1.8)
ax2.plot(k_newton, newton_diffs, marker='s', markersize=6, label=f'Метод Ньютона (итераций: {len(newton_diffs)})', color='#2ca02c', linewidth=2.0)
ax2.axhline(eps, color='gray', linestyle='--', linewidth=1.2, label=f'Порог $\\varepsilon = {eps}$')

ax2.set_yscale('log')
ax2.set_title('Зависимость погрешности от номера итерации', fontsize=13, fontweight='bold')
ax2.set_xlabel('Номер итерации (k)', fontsize=11)
ax2.set_ylabel(r'Погрешность $|x_k - x_{k-1}|$ (log-scale)', fontsize=11)
ax2.grid(True, which='both', linestyle=':', alpha=0.6)
ax2.legend(fontsize=10)

plt.tight_layout()

os.makedirs('data/lab2', exist_ok=True)
out_png = 'data/lab2/plot_2_1.png'
plt.savefig(out_png, dpi=300)
print(f"График сохранён в: {out_png}")

plt.show()
