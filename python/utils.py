from plotly.graph_objs import graph_objs as go
import numpy as np
import gudhi as gd


def get_simplices(simplex_tree, order, alpha):
    return np.array([s[0] for s in simplex_tree.get_skeleton(order) if len(s[0]) == order+1 and s[1] <= alpha])

def plot_simplicial_complex(
    points,
    simplex_tree,
    filtration_value=0.0,
    max_skeleton_order=None,
    points_thickness = 1,
    points_colour = 'red',
    segments_thickness = 0.5,
    segments_colour = 'black',
    other_color = 'mediumaquamarine',
    show = True):

    if max_skeleton_order == None:
        max_skeleton_order = simplex_tree.dimension()
    if max_skeleton_order > simplex_tree.dimension():
        max_skeleton_order = simplex_tree.dimension()
    if max_skeleton_order < 0:
        print("ERROR: skeleton order must be positive")
        return

    if(points.shape[1]==2):
        # adding a zero column for z axis
        n_points = points.shape[0]
        z_col = np.zeros((n_points,1))
        points = np.hstack((points, z_col))

    # Create a trace with all the vertices
    f0 = go.Scatter3d(
        x=points[:,0],
        y=points[:,1],
        z=points[:,2],
        mode="markers",
        marker=dict(
            size=points_thickness,
            color=points_colour
        )
    )
    data = [f0]

    if max_skeleton_order >= 1:
        edges = get_simplices(simplex_tree, 1, filtration_value)
        if np.size(np.shape(edges)) != 1:
            x_lines = list()
            y_lines = list()
            z_lines = list()
            for p in edges:
                for i in range(2):
                    x_lines.append(points[:,0][p[i]])
                    y_lines.append(points[:,1][p[i]])
                    z_lines.append(points[:,2][p[i]])
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
                    color=segments_colour
                )
            )
            data.append(f1)

    if max_skeleton_order >= 2:
        triangles = get_simplices(simplex_tree, 2, filtration_value)
        if np.size(np.shape(triangles)) != 1:
            f2 = go.Mesh3d(
                    x=points[:,0],
                    y=points[:,1],
                    z=points[:,2],
                    i = triangles[:,0],
                    j = triangles[:,1],
                    k = triangles[:,2],
                    color=other_color
            )
            data.append(f2)

    if max_skeleton_order >= 3:
        tetrahedra = get_simplices(simplex_tree, 3, filtration_value)
        if np.size(np.shape(tetrahedra)) != 1:
            f3 = go.Mesh3d(
                    x=points[:,0],
                    y=points[:,1],
                    z=points[:,2],
                    i = tetrahedra[:,0],
                    j = tetrahedra[:,1],
                    k = tetrahedra[:,2],
                    color=other_color
            )
            data.append(f3)

    fig = go.Figure(data=data,
                    layout = go.Layout(showlegend=False))
    if(show):
        fig.show()

    return([data,fig])
