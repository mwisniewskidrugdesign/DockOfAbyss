import numpy as np

a = [[1, 2], [3, 4]]
macierz = np.zeros((2, 2, 3, 8))

# Konwersja listy na tablicę NumPy
a_array = np.array(a[1])
print(a_array)
# Przypisanie tablicy do odpowiedniego zakresu w macierzy
macierz[0, 1, 2, 2:4] = a_array

print("Macierz po przypisaniu:")
print(macierz)