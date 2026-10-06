import matplotlib.pyplot as plt
import pandas as pd
import os

csv_path = 'data/lab3/plot_3_3.csv'
meta_path = 'data/lab3/meta_3_3.txt'

if not os.path.exists(csv_path) or not os.path.exists(meta_path):
    print("Файлы данных не найдены.")
    exit()

df = pd.read_csv(csv_path)
nodes_x, nodes_y = [], []
with open(meta_path, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            nodes_x.append(float(parts[0]))
            nodes_y.append(float(parts[1]))

plt.figure(figsize=(10, 6))
plt.scatter(nodes_x, nodes_y, color='black',
            zorder=5, s=60, label='Исходные данные')
plt.plot(df['x'], df['P1'], label='P1(x) (1-я степень)', color='blue')
plt.plot(df['x'], df['P2'], label='P2(x) (2-я степень)', color='green')
plt.plot(df['x'], df['P3'], label='P3(x) (3-я степень)', color='red')

plt.title('Метод наименьших квадратов (Задача 3.3)')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()
