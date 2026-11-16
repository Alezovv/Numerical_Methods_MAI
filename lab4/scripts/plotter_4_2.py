import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_h_path = 'data/lab4/plot_4_2.csv'
csv_h2_path = 'data/lab4/plot_4_2_h2.csv'

if not os.path.exists(csv_h_path) or not os.path.exists(csv_h2_path):
    print("Файлы данных не найдены. Сначала запустите C++ программу (lab4.exe -> пункт 2).")
    exit(1)

df_h = pd.read_csv(csv_h_path)
df_h2 = pd.read_csv(csv_h2_path)

# =============================================================================
# График 1: В точности как на образце (фото 2)
# =============================================================================
plt.figure(figsize=(9, 6.5))

# Плотная сетка для точного решения: y = x + exp(-2x)
x_dense = np.linspace(df_h['x'].min(), df_h['x'].max(), 300)
y_dense = x_dense + np.exp(-2.0 * x_dense)

# 1. Метод стрельбы при h = 0.1
plt.scatter(df_h['x'], df_h['shoot'],
            label='shooting method, step=0.1',
            color='#1f77b4', s=55, zorder=5)

# 2. Метод стрельбы при h = 0.05
plt.scatter(df_h2['x'], df_h2['shoot'],
            label='shooting method, step=0.05',
            color='#ff7f0e', s=45, zorder=6)

# 3. Конечно-разностный метод при h = 0.1
plt.scatter(df_h['x'], df_h['fd'],
            label='finite difference method, step=0.1',
            color='#2ca02c', s=55, zorder=4)

# 4. Конечно-разностный метод при h = 0.05
plt.scatter(df_h2['x'], df_h2['fd'],
            label='finite difference method, step=0.05',
            color='#d62728', s=45, zorder=5)

# 5. Точное аналитическое решение
plt.plot(x_dense, y_dense,
         label='exact solution',
         color='#1f77b4', linewidth=2.0, zorder=2)

plt.xlabel('x', fontsize=11)
plt.ylabel('y', fontsize=11)
plt.title('Boundary Value Problem 4.2: Shooting vs Finite Difference Method', fontsize=12)
plt.grid(True, linestyle='-', alpha=0.5, color='gray')
plt.legend(loc='upper left', fontsize=10, framealpha=0.95)
plt.tight_layout()

out_png = 'data/lab4/plot_4_2.png'
plt.savefig(out_png, dpi=300)
print(f"Основной график 4.2 сохранен в {out_png}")

# =============================================================================
# График 2: График погрешностей (log scale)
# =============================================================================
plt.figure(figsize=(9, 5.5))
shoot_err = np.maximum(df_h['shoot_err'], 1e-16)

plt.plot(df_h['x'], shoot_err, 's-', label='Погрешность стрельбы ($p=4$)', color='#2ca02c', linewidth=1.6)
plt.plot(df_h['x'], df_h['fd_err'], 'o--', label='Погрешность разностного ($p=2$)', color='#1f77b4', linewidth=1.6)

plt.yscale('log')
plt.xlabel('x', fontsize=11)
plt.ylabel(r'Абсолютная погрешность $|y_{num} - y_{exact}|$', fontsize=11)
plt.title('Задача 4.2: Погрешности методов (логарифмическая шкала)', fontsize=12)
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend(loc='upper left', fontsize=10, framealpha=0.95)
plt.tight_layout()

out_err_png = 'data/lab4/plot_4_2_errors.png'
plt.savefig(out_err_png, dpi=300)
print(f"График погрешностей 4.2 сохранен в {out_err_png}")

plt.show()
