import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
import os

os.chdir("build")

plot = True
n_largest = 500 # min(n_largest, # of evals)

# Input matrix
A = pd.read_csv("data/A.dat", header=None).values

# Arnoldi iterations + Francis Algorithm eigenvalues
df_arn = pd.read_csv("data/evals.dat", header=None, names=["Re", "Im"])
df_arn["R"] = np.sqrt(df_arn.Re**2 + df_arn.Im**2)
# Sort by magnitude and then imaginary part
df_arn = df_arn.sort_values(["R", "Im"], ascending=False, ignore_index=True) 

n_largest = min(df_arn.shape[0], n_largest)
print(f"Comparing the {n_largest} largest magnitude eigenvalues")

# Slice
df_arn = df_arn.iloc[:n_largest]

# Compute the eigenvalues using scipy (wrapper for LAPACK)
eval_spy, _ = scipy.linalg.eig(A)#scipy.sparse.linalg.eigs(A, k=n_largest)

# Scipy eigenvalues sorted by magnitude and then imaginary part
df_spy = pd.DataFrame({"Re": np.real(eval_spy),
                       "Im": np.imag(eval_spy),
                       "R": np.abs(eval_spy)}) \
           .sort_values(["R", "Im"], ascending=False, ignore_index=True)\
           .iloc[:n_largest]

# Write data out
df_spy.to_csv("data/eval_spy.dat", index=False, float_format="%+10.6g")
df_arn.to_csv("data/eval_arn.dat", index=False, float_format="%+10.6g")

residual = np.sqrt((df_spy.Re - df_arn.Re).values**2 + (df_spy.Im - df_arn.Im).values**2)
np.savetxt("data/residual.dat", residual)
print(f"RMSE: \n{np.linalg.norm(residual).mean()}")

if plot:
    plt.scatter(df_spy.Re.values, df_spy.Im.values, label="Scipy", c="r", s=np.sqrt(n_largest) * 2.0)
    plt.scatter(df_arn.Re.values, df_arn.Im.values, label="Arnoldi", c="k", s=np.sqrt(n_largest) * 0.8)
    plt.legend(loc="upper right")
    plt.xlabel("Re")
    plt.ylabel("Im")
    plt.tight_layout()
    plt.show()