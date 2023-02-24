import matplotlib.pyplot as plt
import gudhi as gd
gd.persistence_graphical_tools._gudhi_matplotlib_use_tex=False
file = "../output/alpha_persistence"
gd.plot_persistence_barcode(persistence_file=file, max_intervals=1000, alpha=1, legend=True)
plt.savefig(f'../output/plot/barcode.svg')
plt.savefig(f'../output/plot/barcode.png')

gd.plot_persistence_diagram(persistence_file=file, alpha=1,legend=True)
plt.savefig(f'../output/plot/persistence_diagram.svg')
plt.savefig(f'../output/plot/persistence_diagram.png')
