import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


max_n_proc = 6
files = [f'benchmark_parallel_nproc_{i+1}.json' for i in range(max_n_proc)]

timings = []
serial_time = {}
for i, f in enumerate(files):
    with open(f, 'r') as file:
        t_per_L = json.load(file)

        for L, t in t_per_L.items():
            if i == 0:
                serial_time[L] = t
            timings.append({
                'n_proc': i+1,
                'L': int(L),
                'time': t,
                'speedup': serial_time[L] / t,
            })

timings = pd.DataFrame(timings).sort_values(['n_proc', 'L'])
timings.reset_index(drop=True, inplace=True)
print(timings)

# Usual computer screens have an aspect ratio of 16:9.
# Use this aspect ratio for presentations to fill the screen optimally!
fig, ax = plt.subplots(figsize=(8, 4.5))

# Plot the speedup curves.
sns.lineplot(
    data=timings,
    x='n_proc',
    y='speedup',
    hue='L',
    marker='o',
    ax=ax
)

# Plot the ideal scaling as a reference.
x = range(1, 7)
ax.plot(x, x, '--', label='Perfect scaling')

# Style the plot.
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend()
fig.tight_layout()

fig.savefig('benchmark_scaling.pdf')
