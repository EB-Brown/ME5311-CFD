import numpy as np
from scipy import fftpack

COS_AMPLITUDE = -4 * np.pi ** 2
COS_FREQUENCY = 2 * np.pi


def d_x(x_array: np.array) -> np.array:
    return COS_AMPLITUDE * np.cos(COS_FREQUENCY * x_array)


def get_dct_solution(start, end, dx):
    x_array = np.arange(start, end, dx)
    x_len = len(x_array)
    d = d_x(x_array)
    cos_trans = fftpack.dct(d, norm='ortho')  # ortho = MATLAB func
    denominator = (2 * np.cos(np.pi * (np.arange(x_len) + 1) / x_len) - 2)
    p_k = (dx ** 2) * cos_trans / denominator
    idx = ~np.isinf(p_k)
    keep = 0.4
    answer = fftpack.idct(p_k[idx][:int(x_len * keep)], norm='ortho')
    print(max(answer))
    plt.plot(x_array[idx][:int(x_len * keep)] / keep, answer)
    plt.show()
    print()
    # return x_array, fftpack.idct(dx ** 2 * fftpack.dct(d_x(x_array)) / (2 * np.cos(np.pi * np.arange(x_len) / x_len) - 2))


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    x, p = get_dct_solution(0, 2, 0.01)
    dct_solution = {'x': x, 'p': p}

    plt.plot(x, p)
    plt.show()
    print()
