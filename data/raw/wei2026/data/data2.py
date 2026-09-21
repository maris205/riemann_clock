import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import mpmath as mp

mp.mp.dps = 20

def free_sub(d, real_part, imag_val):
    d = int(d)
    real_part = mp.mpf(real_part)
    imag_val = mp.mpf(imag_val)

    N = 2 ** d

    de = mp.mpf('0')
    for n in range(1, N + 1):
        de += 1 / (n ** real_part)

    s1 = mp.mpc(real_part, imag_val)

    sum_expr = mp.mpc(0)
    for n in range(1, N + 1):
        sign = (n % 2) * 2 - 1
        term = sign * (n ** (-s1))
        sum_expr += term

    app = -sum_expr

    abs_val = abs(app / de)
    result = -1 / real_part * mp.log(abs_val)

    return result

def compute_result_for_d(d_val):
    f = np.vectorize(lambda x: float(free_sub(d_val, real_part, x)))
    return f(data)

mat_in = loadmat('tt_type11_full.mat')

data = mat_in['tt_type11_full']
data = np.array(data, dtype=float)

real_part = mp.mpf('0.3')
results_d4_type2 = compute_result_for_d(4)
results_d8_type2 = compute_result_for_d(8)

real_part = mp.mpf('0.5')

results_d4_type1 = compute_result_for_d(4)
results_d8_type1 = compute_result_for_d(8)

np.savez(
    'data2.npz',
    results_d4_type2=results_d4_type2,
    results_d8_type2=results_d8_type2,
    results_d4_type1=results_d4_type1,
    results_d8_type1=results_d8_type1
)