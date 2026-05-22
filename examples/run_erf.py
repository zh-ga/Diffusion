"""Example: Run the original analytical erf solver."""

import matplotlib.pyplot as plt
from diffusion import ML_CVD_Model

model = ML_CVD_Model()
model("para_3layer.yaml", dcal_type=0)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(model.x_position, model.c1_res, "b-", lw=1.5)
ax1.set_xlabel("Depth (um)")
ax1.set_ylabel("Concentration (cm^-3)")
ax1.set_title("Species 1 (erf solution)")
ax1.set_yscale("log")
ax1.grid(True, alpha=0.3)

ax2.plot(model.x_position, model.c2_res, "r-", lw=1.5)
ax2.set_xlabel("Depth (um)")
ax2.set_ylabel("Concentration (cm^-3)")
ax2.set_title("Species 2 (erf solution)")
ax2.set_yscale("log")
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("result_erf.png", dpi=150)
print("Saved result_erf.png")
