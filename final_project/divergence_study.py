import os
from pathlib import Path

import numpy as np

from cfd_tools import FluidProfile, LinearMesh, plot_contour, Field


def get_velocity_field(component: np.ndarray, array: np.ndarray) -> np.ndarray:
    rand = np.random.random((len(array), len(component)))
    rand[:, :2] = rand[:, -2:] = 0
    rand[0] = rand[-1] = 0
    rand[1:-1, 2:-2] -= rand[1:-1, 2:-2].mean(axis=0)
    return rand


x_start, x_end, x_step = 0, 2, .01
x_half_s = x_step / 2
y_start, y_end, y_step = 0, 2, .01
y_half_s = y_step / 2

x_velocity_grid = LinearMesh(
    x_start=x_start - x_step, x_end=x_end + x_step, x_step=x_step,
    y_start=y_start - y_half_s, y_end=y_end + y_half_s, y_step=y_step,
)

y_velocity_grid = LinearMesh(
    x_start=x_start - x_half_s, x_end=x_end + x_half_s, x_step=x_step,
    y_start=y_start - y_step, y_end=y_end + y_step, y_step=y_step,
)

x_velocity = Field(
    values=get_velocity_field(x_velocity_grid.x, x_velocity_grid.y),
    mesh=x_velocity_grid,
)

y_velocity = Field(
    values=get_velocity_field(y_velocity_grid.y, y_velocity_grid.x).T,
    mesh=y_velocity_grid,
)

print(f"Max U average along y_axis: {x_velocity.values.mean(axis=0).max()}")
print(f"Max V average along x_axis: {y_velocity.values.mean(axis=1).max()}\n")

deliverable = FluidProfile(x_velocity, y_velocity)

prior_divergence = deliverable.velocity_divergence
print(f"Maximum divergence before subtracting pressure {prior_divergence.max()}")
deliverable.remove_divergence()
post_divergence = deliverable.velocity_divergence
print(f"Maximum divergence after subtracting pressure {post_divergence.max()}")
print()

output_dir = Path(__file__).parent / "plots/divergence_study"
output_dir.mkdir(parents=True, exist_ok=True)
for png in output_dir.rglob("*.png"):
    os.remove(png)

fig, ax = plot_contour(
    x_array=deliverable.pressure_x_domain,
    y_array=deliverable.pressure_y_domain,
    z_array=prior_divergence,
    title="Divergence Field Prior to Subtracting Pressure",
    x_label="X Domain",
    y_label="Y Domain",
)
fig.savefig(output_dir / "divergence_before_pressure.png")

fig, ax = plot_contour(
    x_array=deliverable.pressure_x_domain,
    y_array=deliverable.pressure_y_domain,
    z_array=post_divergence,
    title="Divergence Field After Subtracting Pressure",
    x_label="X Domain",
    y_label="Y Domain",
)
fig.savefig(output_dir / "divergence_after_pressure.png")

average = post_divergence.mean()
abs_average = abs(post_divergence).mean()
stndrd = np.std(post_divergence)
max_val = post_divergence.max()
min_val = post_divergence.min()
comp_precision = np.isclose(post_divergence, 0).sum() \
                 / post_divergence.shape[0] ** 2  # Square grid

print(f"""
Unfortunately the python contour plot is not able to generate a plot indicating 
regions where the divergence is less than computer precision but below is a list 
of statistics regarding the divergence after subtracting the pressure gradient:
\n
Average Value = {average}\n
Absolute Average = {abs_average}\n
Standard Deviation = {stndrd}\n
Max Value =  {max_val}\n
Min Value =  {min_val}\n
Percentage of points close to computer precision = 
{round(100 * comp_precision, 3)}%\n
""")
