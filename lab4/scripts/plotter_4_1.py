import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_path = 'data/lab4/plot_4_1.csv'
if not os.path.exists(csv_path):
    print(f"Файл {csv_path} не найден. Сначала запустите C++ программу.")
    exit(1)

df = pd.read_csv(csv_path)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

x_dense = np.linspace(df['x'].min(), df['x'].max(), 300)
y_dense = (1.0 + x_dense) * np.exp(x_dense**2)

ax1.plot(x_dense, y_dense, label=r'Точное: $y = (1+x)e^{x^2}$', color='black', linewidth=2)
ax1.plot(df['x'], df['euler'], 'o--', label='Эйлер', color='red', markersize=5)
ax1.plot(df['x'], df['rk4'], 's-', label='Рунге-Кутта 4', color='green', markersize=5)
ax1.plot(df['x'], df['adams'], '^:', label='Адамс 4', color='blue', markersize=5)

ax1.set_title('Задача 4.1: Сравнение методов')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.legend()

df_err = df.iloc[1:].copy()

ax2.plot(df_err['x'], df_err['euler_err'], 'o--', label='Погрешность Эйлера', color='red')
ax2.plot(df_err['x'], df_err['rk4_err'], 's-', label='Погрешность РК4', color='green')
ax2.plot(df_err['x'], df_err['adams_err'], '^:', label='Погрешность Адамса', color='blue')

ax2.set_yscale('log')
ax2.set_title('Погрешности методов (log scale)')
ax2.set_xlabel('x')
ax2.set_ylabel(r'$|y_{num} - y_{exact}|$')
ax2.grid(True, which='both', linestyle='--', alpha=0.6)
ax2.legend()

plt.tight_layout()

out_png = 'data/lab4/plot_4_1.png'
plt.savefig(out_png, dpi=300)
print(f"График сохранен в {out_png}")
plt.show()
