import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_path = 'data/lab4/plot_4_2.csv'
if not os.path.exists(csv_path):
    print(f"Файл {csv_path} не найден. Сначала запустите C++ программу.")
    exit(1)

df = pd.read_csv(csv_path)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

x_dense = np.linspace(df['x'].min(), df['x'].max(), 300)
y_dense = x_dense + np.exp(-2.0 * x_dense)

ax1.plot(x_dense, y_dense, label=r'Точное: $y = x + e^{-2x}$', color='black', linewidth=2)
ax1.plot(df['x'], df['shoot'], 's-', label='Метод стрельбы', color='green', markersize=5)
ax1.plot(df['x'], df['fd'], 'o--', label='Конечно-разностный', color='blue', markersize=5)

ax1.set_title('Задача 4.2: Сравнение методов')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.legend()

shoot_err = np.maximum(df['shoot_err'], 1e-16)

ax2.plot(df['x'], shoot_err, 's-', label='Погрешность стрельбы', color='green')
ax2.plot(df['x'], df['fd_err'], 'o--', label='Погрешность разностного', color='blue')

ax2.set_yscale('log')
ax2.set_title('Погрешности методов (log scale)')
ax2.set_xlabel('x')
ax2.set_ylabel(r'$|y_{num} - y_{exact}|$')
ax2.grid(True, which='both', linestyle='--', alpha=0.6)
ax2.legend()

plt.tight_layout()

out_png = 'data/lab4/plot_4_2.png'
plt.savefig(out_png, dpi=300)
print(f"График сохранен в {out_png}")
plt.show()
