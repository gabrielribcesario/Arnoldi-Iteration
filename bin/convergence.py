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

title_str = "Eigenvalues vs Ritz values\nArnoldi Iteration #{:0{}}\n|λ|~U(0,1)"

# Plot approximation steps
fig, ax = plt.subplots(figsize=(8,8))
figtitle = fig.suptitle(title_str.format(nelem, zpad))
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
ax.scatter(df_sol.Re, df_sol.Im, label="Eigenvalues", c="k", s=10.)
ritz_scatter = ax.scatter(history[-1].real, history[-1].imag, label="Ritz values", c="r", s=10.)
fig.legend(loc="upper right")
fig.tight_layout()

def update_frame(frame):
    figtitle.set_text(title_str.format(frame+1, zpad))
    slice_i = history[frame][:frame+1]
    ritz_scatter.set_offsets(np.c_[slice_i.real, slice_i.imag])
    return figtitle, ritz_scatter

fps = 30
repeat_delay = 500 # [ms]
extra_args = ["-final_delay", f"{repeat_delay/10}", "-loop", "0"]

ani = animation.FuncAnimation(fig=fig, func=update_frame, frames=nelem, interval=1000/fps, repeat_delay=repeat_delay)
ani.save("figures/convergence.gif", writer=animation.writers["ffmpeg"](fps=fps, extra_args=extra_args))
fig.savefig("figures/comparison.png")
plt.show()