import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

os.chdir("../build")

plt.rcParams['font.family'] = 'serif' 

df_sol = pd.read_csv("data/eval.dat")
df_sol["R"] = np.sqrt(df_sol.Re**2 + df_sol.Im**2)
df_sol = df_sol.sort_values(["R", "Im"], ascending=False, ignore_index=True) 

df_arn = pd.read_csv("data/ritz.dat")
df_arn["R"] = np.sqrt(df_arn.Re**2 + df_arn.Im**2)
df_arn = df_arn.sort_values(["R", "Im"], ascending=False, ignore_index=True) 

n_largest = min(len(df_arn), 1000)
df_sol = df_sol.iloc[:n_largest]
print(f"Comparing the {n_largest} largest magnitude eigenvalues")

residual = np.sqrt((df_sol.Re - df_arn.Re).values**2 + (df_sol.Im - df_arn.Im).values**2)
# np.savetxt("data/residual.dat", residual)
print(f"RMSE: \n{np.linalg.norm(residual).mean()}")

fig, ax = plt.subplots(figsize=(8,8))
fig.suptitle(f"{n_largest} largest magnitude eigenvalues\n|λ|~U(0,1)")
ax.grid(True)
ax.set_axisbelow(True)
# ax lim
ax.set_xlim(left=-1.1, right=1.1)
ax.set_ylim(bottom=-1.1, top=1.1)
circle = patches.Circle((0., 0.), radius=1., color="k", fill=False)
ax.add_patch(circle)
# label
ax.set_xlabel("Re(λ)")
ax.set_ylabel("Im(λ)")
# evals
ax.scatter(df_sol.Re.values, df_sol.Im.values, label="Eigenvalues", c="k", s=16.)
ax.scatter(df_arn.Re.values, df_arn.Im.values, label="Ritz values", c="r", s=2.)
fig.legend(loc="upper right")
fig.tight_layout()
plt.show()