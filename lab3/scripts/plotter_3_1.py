import matplotlib.pyplot as plt
import pandas as pd
import os

# Пути к файлам (предполагается запуск из корня проекта)
csv_path = 'data/lab3/plot_3_1.csv'
meta_path = 'data/lab3/meta_3_1.txt'

if not os.path.exists(csv_path) or not os.path.exists(meta_path):
    print("Файлы данных не найдены. Сначала запустите C++ программу.")
    exit()

# Чтение плотных точек для плавных кривых
df = pd.read_csv(csv_path)

# Чтение узлов и целевой точки x*
nodes_x, nodes_y = [], []
star_x, star_y = 0.0, 0.0

with open(meta_path, 'r') as f:
    lines = f.readlines()
    is_star = False
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line == 'STAR':
            is_star = True
            continue

        parts = line.split()
        if not is_star:
            nodes_x.append(float(parts[0]))
            nodes_y.append(float(parts[1]))
        else:
            star_x = float(parts[0])
            star_y = float(parts[1])

# Отрисовка
plt.figure(figsize=(10, 6))

# Ограничиваем отрисовку точной функции только для x >= 0 (т.к. y = sqrt(x))
df_exact = df[df['x'] >= 0]
plt.plot(df_exact['x'], df_exact['f(x)'],
         label='f(x) = √x', color='black', linewidth=2)

# Полиномы
plt.plot(df['x'], df['lagrange'], label='Многочлен Лагранжа',
         linestyle='--', color='blue', alpha=0.7)
plt.plot(df['x'], df['newton'], label='Многочлен Ньютона',
         linestyle=':', color='red', alpha=0.7)

# Точки
plt.scatter(nodes_x, nodes_y, color='green',
            zorder=5, s=60, label='Узлы интерполяции')
plt.scatter([star_x], [star_y], color='purple', marker='*',
            zorder=6, s=200, label=f'Точка x* = {star_x}')

# Оформление
plt.title('Интерполяция (Задача 3.1)')
plt.xlabel('x')
plt.ylabel('y')
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()
