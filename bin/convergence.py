import matplotlib.animation as animation
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

os.chdir("..")
os.system("mkdir -p figures")

plt.rcParams['font.family'] = 'serif' 

input_prefix = "build/data/convergence"

# Eigenvalues of A
df_sol = pd.read_csv(f"{input_prefix}/eval.csv")
nelem = len(df_sol)

# History of the Ritz values of A during the Arnoldi Iteration
history = np.fromfile(f"{input_prefix}/history.bin", dtype=np.complex128).reshape(nelem, nelem)

# Zero-pad string
zpad = int(np.log10(nelem)) + 1 

title_str = "Eigenvalue Approximation\nArnoldi Iteration #{:0{}}\n|λ|~U(0,1)"

# Plot approximation steps

fig, axes = plt.subplots(1, 2, figsize=(16,8))
figtitle = fig.suptitle(title_str.format(nelem, zpad))

ax1 = axes[0]
ax1.grid(True)
ax1.set_axisbelow(True)
ax1.set_aspect("equal")

# ax lim
ax1.set_xlim(left=-1.1, right=1.1)
ax1.set_ylim(bottom=-1.1, top=1.1)
circle = patches.Circle((0., 0.), radius=1., color="k", fill=False)
ax1.add_patch(circle)

# label
ax1.set_xlabel("Re(λ)")
ax1.set_ylabel("Im(λ)")

# evals
ax1.scatter(df_sol.Re, df_sol.Im, label="Eigenvalues", c="k", s=10.)
ritz_scatter = ax1.scatter(history[-1].real, history[-1].imag, label="Ritz values", c="r", s=10.)
ax1.legend(loc="upper right")

# Mean and max residual norm plots

ax2 = axes[1]
ax2.grid(True)
ax2.set_axisbelow(True)

t = np.arange(1, nelem + 1)
# Sort first on norm and then on imaginary component
idx1 = np.lexsort((df_sol.Im.values, -np.sqrt( df_sol.Re.values**2 + df_sol.Im.values**2 )))

max_residual_norm = []
mean_residual_norm = []
for i, iter in enumerate(history):
    iter = iter.copy()[:i+1]
    # Sort first on norm and then on imaginary component
    idx2 = np.lexsort((iter.imag, -np.abs(iter)))
    # Residual norm of the largest eigenvalues
    residual_norm = np.sqrt( (iter[idx2].real - df_sol.Re.values[idx1][:iter.size])**2 +
                             (iter[idx2].imag - df_sol.Im.values[idx1][:iter.size])**2 )
    max_residual_norm.append(np.max(residual_norm))
    mean_residual_norm.append(np.mean(residual_norm))

# ax label
ax2.set_xlabel("Iteration")
ax2.set_ylabel("Residual norm")

# Plots
max_norm_plot, = ax2.plot(t, max_residual_norm, label="Maximum residual norm", c="k")
mean_norm_plot, = ax2.plot(t, mean_residual_norm, label="Mean residual norm", c="r")
ax2.legend(loc="upper right")

fig.tight_layout()

def update_frame(frame):
    figtitle.set_text(title_str.format(frame+1, zpad))
    slice_i = history[frame][:frame+1]
    ritz_scatter.set_offsets(np.c_[slice_i.real, slice_i.imag])
    max_norm_plot.set_data(t[:frame+1], max_residual_norm[:frame+1])
    mean_norm_plot.set_data(t[:frame+1], mean_residual_norm[:frame+1])
    return figtitle, ritz_scatter, max_norm_plot, mean_norm_plot


fps = 30
repeat_delay = 500 # [ms]
extra_args = ["-final_delay", f"{repeat_delay/10}", "-loop", "0"]

ani = animation.FuncAnimation(fig=fig, func=update_frame, frames=nelem, interval=1000/fps, repeat_delay=repeat_delay)
ani.save("figures/convergence.gif", writer=animation.writers["ffmpeg"](fps=fps, extra_args=extra_args))
fig.savefig("figures/comparison.png")
plt.show()