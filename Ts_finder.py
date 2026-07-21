import Variables as var
import numpy as np
import matplotlib.pyplot as plt


def newton_raphson(T0, f, tol=1e-6, max_iter=100, dx = 1e-12):

    for _ in range(max_iter):

        for _ in range(max_iter):

            f_prime = (f(T0 + dx) - f(T0 - dx)) / (2 * dx)

            Ts = T0 - f(T0) / f_prime

            delta = Ts - T0

            if abs(delta) < tol:
                return Ts

            T0 = Ts
            
    raise RuntimeError(
        f"Newton-Raphson did not converge after {max_iter} iterations."
    )
