"""Module providing the Network class and the NetworkHomology class."""
# Modules
import warnings
import numpy as np
import gudhi as gd
from pandas import read_csv

# graphics modules
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import matplotlib

# settings
warnings.filterwarnings("ignore")
gd.persistence_graphical_tools._gudhi_matplotlib_use_tex = False

############################################################################


class Network:
    """
    A class representing a network in a d-dimensional space.

    Attributes:
        vertices (numpy.ndarray): An array of vertex coordinates.
        edges (numpy.ndarray): An array of edge connections.
        name (str): The name or identifier of the network.
        n_vertices (int): The number of vertices in the network.
        n_edges (int): The number of edges in the network.
        n_coord (int): The space dimension of the network.

    Methods:
        read_vertices(name, sep, header=None): Read and parse vertex data from a file.
        read_edges(name, sep, header=None): Read and parse edge data from a file.
        plot(segments_thickness=1, segments_color="black", vertices_thickness=0.5, vertices_color="red", show=True):
            Plot a dynamic 3D visualization of the network.
        plot_static(segments_thickness=1, segments_color="black", vertices_thickness=0.5, vertices_color="red", show=True):
            Plot a static 2D or 3D visualization of the network.
        __str__(): Return a string representation of the network.

    Note:
        - Make sure to import the required libraries (numpy, matplotlib, plotly, pandas) before using this class.
    """

    def __init__(self, name="", vertices=None, edges=None):
        self.vertices = vertices
        self.edges = edges
        self.name = name
        if vertices is not None:
            self.n_vertices = vertices.shape[0]
            self.n_coord = vertices.shape[1]
        if edges is not None:
            self.n_edges = edges.shape[0]

    def read_vertices(self, name, sep, header=None):
        """
        Parse vertex data from a file.

        Args:
            name (str): The name of the file.
            sep (str): The delimiter used in the CSV file.
            header (int, optional): The row index to use as the header. Defaults to None.
        """
        self.vertices = read_csv(name, delimiter=sep, header=header).to_numpy()
        self.n_vertices = self.vertices.shape[0]
        self.n_coord = self.vertices.shape[1]

    def read_edges(self, name, sep, header=None):
        """
        Read and parse edge data from a file.

        Args:
            name (str): The name of the file.
            sep (str): The delimiter used in the file.
            header (int, optional): The row index to use as the header. Defaults to None.
        """
        edges = read_csv(name, delimiter=sep, header=header)
        if edges.isna().sum().sum() != 0:
            print("ATTENTION: missing values in the edges data")
        edges = edges.to_numpy()
        self.n_edges = edges.shape[0]
        # if there are three column the first one is useless (edges' indexes)
        if edges.shape[1] > 2:
            edges = edges[:, 1:3]
        # shift self.vertices' indexes if they start from 1
        if min(edges[:, 0]) > 0 and min(edges[:, 1]) > 0:
            edges -= 1
        self.edges = edges

    def plot(
        self,
        segments_thickness=1,
        segments_color="black",
        vertices_thickness=0.5,
        vertices_color="red",
        show=False,
        save=True,
    ):
        """
        Plot a dynamic 3D visualization of the network.

        Args:
            segments_thickness (int, optional). Defaults to 1.
            segments_color (str, optional). Defaults to "black".
            vertices_thickness (float, optional). Defaults to 0.5.
            vertices_color (str, optional). Defaults to "red".
            show (bool, optional). Defaults to False.
            save (bool, optional): Save the plot to an html file. Defaults to True.
        """
        # vertices
        vertices = self.vertices
        if self.n_coord == 2:
            # adding a zero column for z axis
            z_col = np.zeros((self.n_vertices, 1))
            vertices = np.hstack((self.vertices, z_col))

        # Create a trace with all the vertices
        f_0 = go.Scatter3d(
            x=vertices[:, 0],
            y=vertices[:, 1],
            z=vertices[:, 2],
            mode="markers",
            marker=dict(size=vertices_thickness, color=vertices_color),
        )
        data = [f_0]

        # EDGES
        if self.edges is not None:
            x_lines = list()
            y_lines = list()
            z_lines = list()
            for edge in self.edges:
                for i in range(2):
                    x_lines.append(vertices[:, 0][edge[i]])
                    y_lines.append(vertices[:, 1][edge[i]])
                    z_lines.append(vertices[:, 2][edge[i]])
                x_lines.append(None)
                y_lines.append(None)
                z_lines.append(None)
            f_1 = go.Scatter3d(
                x=x_lines,
                y=y_lines,
                z=z_lines,
                mode="lines",
                name="lines",
                line=dict(width=segments_thickness, color=segments_color),
            )
            data.append(f_1)

        # Create the figure
        fig = go.Figure(data=data, layout=go.Layout(showlegend=False))
        if save:
            fig.write_html(f"../output/plot_python/{self.name}_network.html")
        if show:
            fig.show()

    def plot_static(
        self,
        segments_thickness=1,
        segments_color="black",
        vertices_thickness=0.5,
        vertices_color="red",
        show=False,
        save=True,
    ):
        """
        Plot a static 2D or 3D visualization of the network.

        Args:
            segments_thickness (int, optional). Defaults to 1.
            segments_color (str, optional). Defaults to "black".
            vertices_thickness (float, optional). Defaults to 0.5.
            vertices_color (str, optional). Defaults to "red".
            show (bool, optional). Defaults to False.
            save (bool, optional): Save the plot to file. Defaults to True.
        """
        if self.n_coord == 2:
            # vertices
            x = self.vertices[:, 0]
            y = self.vertices[:, 1]
            plt.plot(x, y, ".", color=vertices_color, markersize=vertices_thickness)
            # edges
            for i, j in self.edges:
                plt.plot(
                    [x[i], x[j]],
                    [y[i], y[j]],
                    color=segments_color,
                    linewidth=segments_thickness,
                )
            plt.axis("equal")
            # Remove axis and box
            plt.axis("off")
            plt.xticks([])
            plt.yticks([])

        if self.n_coord == 3:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
            # vertices
            ax.scatter(
                self.vertices[:, 0],
                self.vertices[:, 1],
                self.vertices[:, 2],
                c=vertices_color,
                s=vertices_thickness,
            )
            # edges
            for edge in self.edges:
                x = [self.vertices[edge[0]][0], self.vertices[edge[1]][0]]
                y = [self.vertices[edge[0]][1], self.vertices[edge[1]][1]]
                z = [self.vertices[edge[0]][2], self.vertices[edge[1]][2]]
                ax.plot(x, y, z, c=segments_color, linewidth=segments_thickness)
        if save:
            plt.savefig(f"../output/plot_python/{self.name}_network.pdf")
        if show:
            plt.show()
        plt.close()

    def __str__(self):
        return (
            f"{self.n_coord}-D network:\n"
            f"Number of vertices: {self.n_vertices}\n"
            f"Number of edges: {self.n_edges}"
        )


############################################################################


class Filtration:
    """
    Represents a simplicial complex filtration and provides methods for analysis and visualization.

    Attributes:
        simplex_tree (gd.SimplexTree): The underlying simplicial complex.
        dim (int): The dimension of the simplicial complex.
        filtration_name (str): The name of the filtration.

    Methods:
        __init__(self, simplex_tree): Initialize the filtration with a given simplex tree.
        get_simplices(self, order, filtration_value=float("inf")): Get simplices of a given dimension and filtration value.
        plot_simplicial_complex(self, points, filtration_value=float("inf"), max_skeleton_order=None, points_thickness=1,
                                points_colour="red", segments_thickness=0.5, segments_colour="black", other_color="mediumaquamarine",
                                show=True): Plot the filtered simplicial complex.
        plot_persistence(self, max_intervals=100, save=True, show=False, network_name=""): Plot persistence barcode and diagram.
        __str__(self): Get a string representation of the simplicial complex information.
    """

    simplex_tree = gd.SimplexTree()
    dim = simplex_tree.dimension()
    filtration_name = ""

    def __init__(self, simplex_tree):
        self.simplex_tree = simplex_tree
        self.simplex_tree.compute_persistence()

    def get_simplices(self, order, filtration_value=float("inf")):
        """
        Get the simplices of a given dimension and with filtration value less or equal than a specified value.

        Args:
            order (int): dimension of the simplices.
            filtration_value (float, optional): threshold for the filtration value. Defaults to float("inf").
        """
        return np.array(
            [
                s[0]
                for s in self.simplex_tree.get_skeleton(order)
                if len(s[0]) == order + 1 and s[1] <= filtration_value
            ]
        )

    def plot_simplicial_complex(
        self,
        points,
        filtration_value=float("inf"),
        max_skeleton_order=None,
        points_thickness=1,
        points_colour="red",
        segments_thickness=0.5,
        segments_colour="black",
        other_color="mediumaquamarine",
        show=False,
        save=True,
        network_name="",
    ):
        """
        Plot a filtered simplicial complex.

        Args:
            points (numpy.ndarray): An array of vertex coordinates.
            filtration_value (float, optional): Plot only simplices with filtration value less or equal filtration_value. Defaults to float("inf").
            max_skeleton_order (int, optional): Plot only simplices up to dimension max_skeleton_order. Defaults to None.
            points_thickness (int, optional). Defaults to 1.
            points_colour (str, optional). Defaults to "red".
            segments_thickness (float, optional). Defaults to 0.5.
            segments_colour (str, optional). Defaults to "black".
            other_color (str, optional): Color of triangles and tetrahedra. Defaults to "mediumaquamarine".
            show (bool, optional): Show the plot. Defaults to False.
            save (bool, optional): Save the plot to an html file. Defaults to True.
            network_name (str, optional): Name of the network (if save=True). Defaults to "".
        """
        if points.shape[1] > 3 or points.shape[1] < 2:
            print("ERROR: points must be 2D or 3D")
            return

        if max_skeleton_order is None:
            max_skeleton_order = self.simplex_tree.dimension()
        if max_skeleton_order > self.simplex_tree.dimension():
            max_skeleton_order = self.simplex_tree.dimension()
        if max_skeleton_order < 0:
            print("ERROR: skeleton order must be positive")
            return

        if points.shape[1] == 2:
            # adding a zero column for z axis
            n_points = points.shape[0]
            z_col = np.zeros((n_points, 1))
            points = np.hstack((points, z_col))

        # Create a trace with all the vertices
        vertices = self.get_simplices(0, filtration_value)
        indexes = [v[0] for v in vertices]
        selected_points = points[indexes, :]
        f_0 = go.Scatter3d(
            x=selected_points[:, 0],
            y=selected_points[:, 1],
            z=selected_points[:, 2],
            mode="markers",
            marker=dict(size=points_thickness, color=points_colour),
        )
        data = [f_0]

        if max_skeleton_order >= 1:
            edges = self.get_simplices(1, filtration_value)
            if np.size(np.shape(edges)) != 1:
                x_lines = list()
                y_lines = list()
                z_lines = list()
                for edge in edges:
                    for i in range(2):
                        x_lines.append(points[:, 0][edge[i]])
                        y_lines.append(points[:, 1][edge[i]])
                        z_lines.append(points[:, 2][edge[i]])
                    x_lines.append(None)
                    y_lines.append(None)
                    z_lines.append(None)
                f_1 = go.Scatter3d(
                    x=x_lines,
                    y=y_lines,
                    z=z_lines,
                    mode="lines",
                    name="lines",
                    line=dict(width=segments_thickness, color=segments_colour),
                )
                data.append(f_1)

        if max_skeleton_order >= 2:
            triangles = self.get_simplices(2, filtration_value)
            if np.size(np.shape(triangles)) != 1:
                f_2 = go.Mesh3d(
                    x=points[:, 0],
                    y=points[:, 1],
                    z=points[:, 2],
                    i=triangles[:, 0],
                    j=triangles[:, 1],
                    k=triangles[:, 2],
                    color=other_color,
                )
                data.append(f_2)

        if max_skeleton_order >= 3:
            tetrahedra = self.get_simplices(3, filtration_value)
            if np.size(np.shape(tetrahedra)) != 1:
                f3 = go.Mesh3d(
                    x=points[:, 0],
                    y=points[:, 1],
                    z=points[:, 2],
                    i=tetrahedra[:, 0],
                    j=tetrahedra[:, 1],
                    k=tetrahedra[:, 2],
                    color=other_color,
                )
                data.append(f3)

        fig = go.Figure(data=data, layout=go.Layout(showlegend=False))
        if save:
            fig.write_html(
                f"../output/plot_python/{network_name}_{self.filtration_name}_complex_filtration{filtration_value}.html"
            )
        if show:
            fig.show()
        fig.write_html("prova.html")

    def plot_persistence(
        self, max_intervals=100, save=True, show=False, network_name=""
    ):
        """Plot the persistence barcode and diagram of the simplicial complex.

        Args:
            max_intervals (int, optional): Number of most persistent intervals to show. Defaults to 100.
            save (bool, optional): Save plots to file. Defaults to True.
            show (bool, optional): Show plots. Defaults to False.
            network_name (str, optional): Name of the network (if save=True). Defaults to "".
        """
        # persistence barcode plots
        for i in range(self.simplex_tree.dimension()):
            barcode = self.simplex_tree.persistence_intervals_in_dimension(i)
            col = matplotlib.cm.Set1.colors[i:]
            ax = gd.plot_persistence_barcode(
                barcode, max_intervals=max_intervals, colormap=col, alpha=1
            )
            ax.set_title(
                f"Persistence barcode dim {i}: {self.filtration_name} filtration"
            )
            if save:
                plt.savefig(
                    f"../output/plot_python/{network_name}_{self.filtration_name}_barcode_dim{i}.pdf"
                )
            if show:
                plt.show()
            plt.close()
        # persistence diagram plot
        ax = gd.plot_persistence_diagram(
            self.simplex_tree.persistence(), alpha=1, legend=True, max_intervals=1000
        )
        ax.set_aspect("equal")
        ax.set_title(f"Persistence diagram: {self.filtration_name} filtration")
        if save:
            plt.savefig(
                f"../output/plot_python/{network_name}_{self.filtration_name}_persistence_diagram.pdf"
            )
        if show:
            plt.show()
        plt.close()

    def __str__(self):
        complex_info_str = f"Simplicial complex of dimension {self.dim} composed by:\n"
        for i in range(self.dim + 1):
            n_i = self.get_simplices(i).shape[0]
            complex_info_str += f"- {n_i} simplices of dimension {i}\n"

        betti_numbers_str = "Betti numbers of the simplicial complex are:\n"
        for i in range(self.simplex_tree.dimension()):
            betti_numbers_str += f"- b{i} = {self.simplex_tree.betti_numbers()[i]}\n"

        return complex_info_str + betti_numbers_str


############################################################################


class RadialFiltration(Filtration):
    """
    Represents a radial filtration.

    Attributes:
        simplex_tree (gd.SimplexTree): The underlying simplicial complex.
        dim (int): The dimension of the simplicial complex.
        filtration_name (str): The name of the filtration.

    Methods:
        __init__(self, vertices, edges): Initialize the radial filtration with given vertices and edges.
        __str__(self): Get a string representation of the radial filtration information.
    """

    def __init__(self, vertices, edges):
        center = np.mean(vertices, axis=0)

        for i, vertex in enumerate(vertices):
            self.simplex_tree.insert([i], filtration=np.linalg.norm(vertex - center))

        for edge in edges:
            self.simplex_tree.insert(
                edge,
                filtration=max(
                    np.linalg.norm(vertices[edge[0]] - center),
                    np.linalg.norm(vertices[edge[1]] - center),
                ),
            )
        self.dim = self.simplex_tree.dimension()
        self.simplex_tree.set_dimension(2)
        self.filtration_name = "radial"
        self.simplex_tree.compute_persistence()

    def __str__(self):
        return "Radial Filtration\n" + super().__str__()


############################################################################


class AlphaFiltration(Filtration):
    """
    Represents an alpha filtration.

    Attributes:
        simplex_tree (gd.SimplexTree): The underlying simplicial complex.
        dim (int): The dimension of the simplicial complex.
        filtration_name (str): The name of the filtration.

    Methods:
        __init__(self, vertices): Initialize the radial filtration with given vertices.
        __str__(self): Get a string representation of the alpha filtration information.
    """

    def __init__(self, vertices):
        alpha_complex = gd.AlphaComplex(points=vertices, precision="safe")
        self.simplex_tree = alpha_complex.create_simplex_tree()
        self.dim = self.simplex_tree.dimension()
        self.filtration_name = "alpha"
        self.simplex_tree.compute_persistence()

    def __str__(self):
        return "Alpha Filtration\n" + super().__str__()


############################################################################


class NetworkAnalysis:
    """
    Represents the topological analysis of a network using radial and alpha filtrations.

    Attributes:
        network (Network): The network to be analyzed.
        radial_filtration (RadialFiltration): The radial filtration of the network.
        alpha_filtration (AlphaFiltration): The alpha filtration of the network.

    Methods:
        __init__(self, network, name): Initialize the NetworkAnalysis object with a network and a name.
        plot_persistence(self, max_intervals=100, save=True, show=False): Plot persistence diagrams for both filtrations.
        __str__(self): Get a string representation of the network analysis results.
    """

    def __init__(self, network):
        self.network = network
        self.radial_filtration = RadialFiltration(network.vertices, network.edges)
        self.alpha_filtration = AlphaFiltration(network.vertices)

    def plot_persistence(self, max_intervals=100, save=True, show=False):
        """
        Plot persistence diagrams for both radial and alpha filtrations.

        Args:
            max_intervals (int, optional): Number of most persistent intervals to show. Defaults to 100.
            save (bool, optional): Save plots to file. Defaults to True.
            show (bool, optional): Show plots. Defaults to False.
        """
        self.radial_filtration.plot_persistence(
            max_intervals, save, show, self.network.name
        )
        self.alpha_filtration.plot_persistence(
            max_intervals, save, show, self.network.name
        )

    def plot_filtration(
        self,
        filtration="radial",
        filtration_value=float("inf"),
        max_skeleton_order=None,
        points_thickness=1,
        points_colour="red",
        segments_thickness=0.5,
        segments_colour="black",
        other_color="mediumaquamarine",
        show=False,
        save=True,
    ):
        """
        Wrapper function to plot the simplicial complex of a given filtration.

        Args:
            filtration (str, optional): Filtration to plot, 'radial' or 'alpha'. Defaults to "radial".
            filtration_value (float, optional): Plot only simplices with filtration value less or equal filtration_value. Defaults to float("inf").
            max_skeleton_order (int, optional): Plot only simplices up to dimension max_skeleton_order. Defaults to None.
            points_thickness (int, optional). Defaults to 1.
            points_colour (str, optional). Defaults to "red".
            segments_thickness (float, optional). Defaults to 0.5.
            segments_colour (str, optional). Defaults to "black".
            other_color (str, optional): Color of triangles and tetrahedra. Defaults to "mediumaquamarine".
            show (bool, optional): Show the plot. Defaults to False.
            save (bool, optional): Save the plot to an html file. Defaults to True.
        """
        if filtration == "radial":
            self.radial_filtration.plot_simplicial_complex(
                self.network.vertices,
                filtration_value,
                max_skeleton_order,
                points_thickness,
                points_colour,
                segments_thickness,
                segments_colour,
                other_color,
                show,
                save,
                self.network.name,
            )
        elif filtration == "alpha":
            self.alpha_filtration.plot_simplicial_complex(
                self.network.vertices,
                filtration_value,
                max_skeleton_order,
                points_thickness,
                points_colour,
                segments_thickness,
                segments_colour,
                other_color,
                show,
                save,
                self.network.name,
            )
        else:
            print("ERROR: filtration must be 'radial' or 'alpha'")

    def __str__(self):
        return (
            f"Topological analysis of {self.network.name} network:\n\n"
            + str(self.radial_filtration)
            + "\n"
            + str(self.alpha_filtration)
        )


############################################################################
