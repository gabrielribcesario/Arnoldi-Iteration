import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy

n_largest = 20
plot = True

A = pd.read_csv("../data//A.dat", header=None).values
# b = pd.read_csv("../data/b.dat", header=None).values[:,0]
# Q = pd.read_csv("../data/Q.dat", header=None).values
H = pd.read_csv("../data/H.dat", header=None).values

n_largest = min(H.shape[1], n_largest)

print(f"Comparing the top {n_largest} largest eigenvalues.")

eigenval_spy, _ = scipy.linalg.eig(A)
df_spy = pd.DataFrame({"Re": np.real(eigenval_spy), 
                       "Im": np.imag(eigenval_spy),
                       "R": np.abs(eigenval_spy)})\
           .sort_values(["R", "Im"], ascending=False, ignore_index=True)\
           .iloc[: n_largest]
df_spy.to_csv("../data/eigenvalues_scipy.dat", index=False)

eigenval_arn, _ = scipy.linalg.eig(H[:-1])
df_arn = pd.DataFrame({"Re": np.real(eigenval_arn), 
                       "Im": np.imag(eigenval_arn),
                       "R": np.abs(eigenval_arn)})\
           .sort_values(["R", "Im"], ascending=False, ignore_index=True)\
           .iloc[: n_largest]
df_arn.to_csv("../data/eigenvalues_arnoldi.dat", index=False)

residual = (df_spy.Re - df_arn.Re).values + (df_spy.Im - df_arn.Im).values * 1.0j
abs_residual = np.abs(residual)
np.savetxt("../data/residual.dat", abs_residual)

print(f"Norm of residuals: \n{np.linalg.norm(abs_residual)}")

if plot:
    plt.scatter(df_spy.Re.values, df_spy.Im.values, label="Eigenvalues of A", c="r", s=np.sqrt(n_largest) * 2.0)
    plt.scatter(df_arn.Re.values, df_arn.Im.values, label="Ritz values of H", c="k", s=np.sqrt(n_largest) * 0.8)
    plt.legend(loc="upper right")
    plt.xlabel("Re")
    plt.ylabel("Im")
    plt.tight_layout()
    plt.show()