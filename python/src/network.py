import numpy as np
from pandas import read_csv

# graphics libraries
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Network:
    def __init__(self, vertices = None, edges = None):
        self.vertices = vertices
        self.edges = edges
        self.n_vertices = None
        self.n_coord = None
        self.n_edges = None
        if vertices is not None:
            self.n_vertices = vertices.shape[0]
            self.n_coord = vertices.shape[1]
        if edges is not None:
            self.n_edges = edges.shape[0]

    ############################################################################

    def read_vertices(self, name, sep, header=None):
        self.vertices = read_csv(name, delimiter=sep, header=header).to_numpy()
        self.n_vertices = self.vertices.shape[0]
        self.n_coord = self.vertices.shape[1]

    ############################################################################

    def read_edges(self, name, sep, header=None):
        edges = read_csv(name, delimiter=sep, header=header)
        if edges.isna().sum().sum() != 0:
            print("ATTENTION: missing values in the edges data")
        edges = edges.to_numpy()
        self.n_edges = edges.shape[0]
        # if there are three column the first one is useless (edges' indexes)
        if edges.shape[1] > 2:
            edges = edges[:,1:3]
        # shift self.vertices' indexes if they start from 1
        if min(edges[:,0]) > 0 and min(edges[:,1]) > 0:
            edges -= 1
        self.edges = edges

    ############################################################################

    def plot(self,
            segments_thickness = 1,
            segments_color = 'black',
            vertices_thickness = 0.5,
            vertices_color = 'red',
            show = True):

        # vertices
        vertices = self.vertices
        if(self.n_coord == 2):
            # adding a zero column for z axis
            z_col = np.zeros((self.n_vertices,1))
            vertices = np.hstack((self.vertices, z_col))

        # Create a trace with all the vertices
        f0 = go.Scatter3d(
                x=vertices[:,0],
                y=vertices[:,1],
                z=vertices[:,2],
                mode="markers",
                marker=dict(
                    size=vertices_thickness,
                    color=vertices_color
                )
            )

        data = [f0]

        # EDGES
        if self.edges is not None:
            x_lines = list()
            y_lines = list()
            z_lines = list()
            for p in self.edges:
                for i in range(2):
                    x_lines.append(vertices[:,0][p[i]])
                    y_lines.append(vertices[:,1][p[i]])
                    z_lines.append(vertices[:,2][p[i]])
                x_lines.append(None)
                y_lines.append(None)
                z_lines.append(None)
            f1 = go.Scatter3d(
                    x=x_lines,
                    y=y_lines,
                    z=z_lines,
                    mode='lines',
                    name='lines',
                    line=dict(
                        width=segments_thickness,
                        color=segments_color
                    )
                )

            data.append(f1)

        # Create the figure
        fig = go.Figure(data=data,
                        layout = go.Layout(showlegend=False))

        if show:
            fig.show()
        return([fig, data])

    ############################################################################

    # low performance in 3D with large networks
    def plot_static(self,
                    segments_thickness = 1,
                    segments_color = 'black',
                    vertices_thickness = 0.5,
                    vertices_color = 'red',
                    show = True):

        if(self.n_coord==2):
             # estraggo le coordinate x e y dei punti
            x = [p[0] for p in self.vertices]
            y = [p[1] for p in self.vertices]

            # rappresento i punti e i lati con la funzione plot
            plt.plot(x, y, '.', color = vertices_color, markersize=vertices_thickness)  # rappresenta i punti come cerchi rossi
            for i, j in self.edges:
                plt.plot([x[i], x[j]], [y[i], y[j]], color = segments_color, linewidth=segments_thickness )  # rappresenta il lato con una linea nera

        if(self.n_coord==3):
            # Creazione del grafico tridimensionale
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # Rappresentazione dei punti tridimensionali
            #for point in self.vertices:
            ax.scatter(self.vertices[:,0], self.vertices[:,1], self.vertices[:,2], c=vertices_color,s=vertices_thickness)

            # Rappresentazione dei lati del poligono
            for edge in self.edges:
                x = [self.vertices[edge[0]][0], self.vertices[edge[1]][0]]
                y = [self.vertices[edge[0]][1], self.vertices[edge[1]][1]]
                z = [self.vertices[edge[0]][2], self.vertices[edge[1]][2]]
                ax.plot(x, y, z, c=segments_color, linewidth=segments_thickness)

        if(show):
            plt.show()

    ############################################################################

    def __str__(self):
        return f"{self.n_coord}-D network:\nNumber of vertices: {self.n_vertices}\nNumber of edges: {self.n_edges}"
