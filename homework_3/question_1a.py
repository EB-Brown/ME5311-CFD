"""
The purpose of this script is to solve Ax = b as presented in the scanned
homework.
"""

import numpy as np

A = np.array([
    [-1, -1, -1, -1],
    [1, -1, -3, -5],
    [-1, -1, -9, -25],
    [1, -1, -27, -125],
])
b = np.array([
    [0],
    [-2],
    [0],
    [0]
])

x = np.linalg.inv(A).dot(b)
print(f"x' = [a,b,c,d]' = {[round(y[0], 3) for y in x]}'")

fourth_column = -1 * np.array([1, 1, 81, 625]) / (16 * 24)

error_term_constant = fourth_column.dot(x)[0]
print(f'The exact error term constant is {round(error_term_constant, 6)}')
print("Which is approximately 1/24")
print()
