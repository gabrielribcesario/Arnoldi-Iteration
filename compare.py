import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

os.chdir("build")

# A = scipy.io.FortranFile("data/A.bin", "r")
# esols = scipy.io.FortranFile("data/esols.bin", "r")
# evals = scipy.io.FortranFile("data/evals.bin", "r")

df_sol = pd.read_csv("data/esols.dat", header=None, names=["Re", "Im"])
df_sol["R"] = np.sqrt(df_sol.Re**2 + df_sol.Im**2)
df_sol = df_sol.sort_values(["R", "Im"], ascending=False, ignore_index=True) 
df_sol.to_csv("data/eval_sol.dat", index=False, float_format="%+10.6g")

df_arn = pd.read_csv("data/evals.dat", header=None, names=["Re", "Im"])
df_arn["R"] = np.sqrt(df_arn.Re**2 + df_arn.Im**2)
df_arn = df_arn.sort_values(["R", "Im"], ascending=False, ignore_index=True) 
df_arn.to_csv("data/eval_arn.dat", index=False, float_format="%+10.6g")

n_largest = min(len(df_arn), 1000)
df_sol = df_sol.iloc[:n_largest]
print(f"Comparing the {n_largest} largest magnitude eigenvalues")

residual = np.sqrt((df_sol.Re - df_arn.Re).values**2 + (df_sol.Im - df_arn.Im).values**2)
np.savetxt("data/residual.dat", residual)
print(f"RMSE: \n{np.linalg.norm(residual).mean()}")

plt.suptitle(f"Comparing the\n{n_largest} largest magnitude eigenvalues")
plt.scatter(df_sol.Re.values, df_sol.Im.values, label="Ground truth", c="r", s=np.sqrt(n_largest) * 2.0)
plt.scatter(df_arn.Re.values, df_arn.Im.values, label="Output", c="k", s=np.sqrt(n_largest) * 0.8)
plt.legend(loc="upper right")
plt.xlabel("Re")
plt.ylabel("Im")
plt.tight_layout()
plt.show()