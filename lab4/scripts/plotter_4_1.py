import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_h_path = 'data/lab4/plot_4_1.csv'
csv_h2_path = 'data/lab4/plot_4_1_h2.csv'

if not os.path.exists(csv_h_path) or not os.path.exists(csv_h2_path):
    print("Файлы данных не найдены. Сначала запустите C++ программу (lab4.exe -> пункт 1).")
    exit(1)

df_h = pd.read_csv(csv_h_path)
df_h2 = pd.read_csv(csv_h2_path)

# =============================================================================
# График 1: В точности как на образце (фото 1)
# =============================================================================
plt.figure(figsize=(9, 6.5))

# Плотная сетка для точного решения
x_dense = np.linspace(df_h['x'].min(), df_h['x'].max(), 300)
y_dense = (1.0 + x_dense) * np.exp(x_dense**2)

# 1. Эйлер: точки при шаге h, линия при шаге h=0.05
plt.scatter(df_h['x'], df_h['euler'], label='Метод Эйлера',
            color='#1f77b4', s=50, zorder=5)
plt.plot(df_h2['x'], df_h2['euler'], label=r'Эйлер, $h = 0.05$',
         color='#1f77b4', linewidth=1.6, zorder=3)

# 2. Рунге-Кутта 4: точки при шаге h, линия при шаге h=0.05
plt.scatter(df_h['x'], df_h['rk4'], label='Метод Рунге-Кутты 4-ого порядка',
            color='#ff7f0e', s=50, zorder=5)
plt.plot(df_h2['x'], df_h2['rk4'], label=r'Рунге-Кутта, $h = 0.05$',
         color='#ff7f0e', linewidth=1.6, zorder=3)

# 3. Адамс 4: точки при шаге h, линия при шаге h=0.05
plt.scatter(df_h['x'], df_h['adams'], label='Метод Адамса 4-ого порядка',
            color='#2ca02c', s=50, zorder=5)
plt.plot(df_h2['x'], df_h2['adams'], label=r'Адамс, $h = 0.05$',
         color='#2ca02c', linewidth=1.6, zorder=3)

# 4. Точное решение: сплошная красная линия
plt.plot(x_dense, y_dense, label='Точное решение',
         color='#d62728', linewidth=2.0, zorder=2)

plt.xlabel('x', fontsize=11)
plt.ylabel('y', fontsize=11)
plt.title('Задача Коши 4.1: Сравнение численных методов и точного решения', fontsize=12)
plt.grid(True, linestyle='-', alpha=0.5, color='gray')
plt.legend(loc='upper left', fontsize=10, framealpha=0.95)
plt.tight_layout()

out_png = 'data/lab4/plot_4_1.png'
plt.savefig(out_png, dpi=300)
print(f"Основной график 4.1 сохранен в {out_png}")

# =============================================================================
# График 2: График погрешностей (log scale)
# =============================================================================
plt.figure(figsize=(9, 5.5))
df_err = df_h.iloc[1:].copy()

plt.plot(df_err['x'], df_err['euler_err'], 'o--', label='Погрешность Эйлера ($p=1$)', color='#1f77b4', linewidth=1.5)
plt.plot(df_err['x'], df_err['rk4_err'], 's-', label='Погрешность РК4 ($p=4$)', color='#ff7f0e', linewidth=1.5)
plt.plot(df_err['x'], df_err['adams_err'], '^:', label='Погрешность Адамса ($p=4$)', color='#2ca02c', linewidth=1.5)

plt.yscale('log')
plt.xlabel('x', fontsize=11)
plt.ylabel(r'Абсолютная погрешность $|y_{num} - y_{exact}|$', fontsize=11)
plt.title('Задача 4.1: Погрешности методов (логарифмическая шкала)', fontsize=12)
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend(loc='upper left', fontsize=10, framealpha=0.95)
plt.tight_layout()

out_err_png = 'data/lab4/plot_4_1_errors.png'
plt.savefig(out_err_png, dpi=300)
print(f"График погрешностей 4.1 сохранен в {out_err_png}")

plt.show()
