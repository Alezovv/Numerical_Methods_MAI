import matplotlib.pyplot as plt
import pandas as pd
import os

csv_path = 'data/lab3/plot_3_2.csv'
meta_path = 'data/lab3/meta_3_2.txt'

if not os.path.exists(csv_path) or not os.path.exists(meta_path):
    print("Файлы данных не найдены. Сначала запустите C++ программу.")
    exit()

df = pd.read_csv(csv_path)

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
            star_x, star_y = float(parts[0]), float(parts[1])

plt.figure(figsize=(10, 6))
plt.plot(df['x'], df['S'], label='Кубический сплайн',
         color='blue', linewidth=2)
plt.scatter(nodes_x, nodes_y, color='red',
            zorder=5, s=60, label='Исходные узлы')
plt.scatter([star_x], [star_y], color='purple', marker='*',
            zorder=6, s=200, label=f'x* = {star_x}, S(x*) = {star_y:.4f}')

plt.title('Кубический сплайн (Задача 3.2)')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()
