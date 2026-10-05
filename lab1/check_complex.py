import numpy as np

print("--- Поиск комплексно-сопряженных собственных значений (NumPy) ---")

A = np.array([
    [9.0, 0.0, 2.0],
    [-6.0, 4.0, 4.0],
    [-2.0, -7.0, 5.0]])

print("\nИсходная матрица A:")
print(A)

# Нахождение собственных значений с помощью встроенной функции
eigenvalues, eigenvectors = np.linalg.eig(A)

print("\nНайденные собственные значения:")
for i, val in enumerate(eigenvalues):
    print(f"Lambda_{i+1} = {val:.4f}")

print("\nКак мы видим, присутствуют комплексно-сопряженные корни:")
print(f"{eigenvalues[0]:.4f} и {eigenvalues[1]:.4f}")
