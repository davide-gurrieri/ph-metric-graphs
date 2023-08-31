import matplotlib.pyplot as plt
import gudhi as gd
gd.persistence_graphical_tools._gudhi_matplotlib_use_tex=False
file = "../output/alpha_persistence"
dim = input("Insert dimension (0 -> connected components, 1 -> loops, 2 -> voids):  ")
dim = int(dim)

birth_death = gd.read_persistence_intervals_in_dimension(persistence_file=file, only_this_dim=dim)

ax = gd.plot_persistence_barcode(persistence=birth_death, max_intervals=1000, alpha = 1)
ax.set_title(f"Persistence barcode for dimension {dim}")
plt.savefig(f'../output/plot/barcode_cpp.pdf')

ax = gd.plot_persistence_diagram(persistence=birth_death, max_intervals=1000, fontsize=12)
ax.set_title(f"Persistence diagram for dimension {dim}")
ax.set_aspect("equal")
plt.savefig(f'../output/plot/persistence_diagram_cpp.pdf')