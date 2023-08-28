import numpy as np
import funcs as fu
import globals as const
import plotly.graph_objs as go


test_center = 8
test_vertex = 8

# Generate cylinder
gvertices, gcenters, gaxes = fu.def_graphs()
fu.normals(gvertices)
norm_avg=0
for vertex in gvertices.nodes():
    norm_avg+=np.linalg.norm(gvertices.nodes[vertex]['normal'])
print('average of normals : ', norm_avg/len(gvertices.nodes()))


fu.simulate(gvertices, gcenters, gaxes, const.steps)
print(gcenters.nodes.data())

# test for area and perimeter of test_center
area, perim = fu.center_area_perim(test_center, gcenters, gvertices)
print('perim: ', perim ,' area: ', area)


# test for vertex normal
fu.normals(gvertices)



# Plot with Plotly
# Create Edges
edge_x, edge_y, edge_z = [], [], []
for idx0, idx1 in gvertices.edges():
    x0, y0, z0 = gvertices.nodes[idx0]["coords"]
    x1, y1, z1 = gvertices.nodes[idx1]["coords"]
    edge_x.extend([x0, x1, None])
    edge_y.extend([y0, y1, None])
    edge_z.extend([z0, z1, None])

vertx, verty, vertz = fu.coords(gvertices)
centx, centy, centz = fu.coords(gcenters)
axesx, axesy, axesz = fu.coords(gaxes)

norm_force_x, norm_force_y, norm_force_z = [], [], []
bend_force_x, bend_force_y, bend_force_z = [], [], []
for k in range(len(gvertices.nodes)):
    x0, y0, z0 = vertx[k], verty[k], vertz[k]
    x1, y1, z1 = gvertices.nodes[k]['force'][0]+x0, gvertices.nodes[k]['force'][1]+y0, gvertices.nodes[k]['force'][2]+z0
    x2, y2, z2 = -gvertices.nodes[k]['normal'][0]*const.c_alpha*const.nu+x0, -gvertices.nodes[k]['normal'][1]*const.c_alpha*const.nu+y0, -gvertices.nodes[k]['normal'][2]*const.c_alpha*const.nu+z0
    norm_force_x.extend([x0, x1, None])
    norm_force_y.extend([y0, y1, None])
    norm_force_z.extend([z0, z1, None])
    bend_force_x.extend([x0, x2, None])
    bend_force_y.extend([y0, y2, None])
    bend_force_z.extend([z0, z2, None])


edge_trace = go.Scatter3d(name="Cell edges",
    x=edge_x, y=edge_y, z=edge_z, line=dict(width=10, color='#888'), mode='lines')

force_trace = go.Scatter3d(name="Total forces",
    x=norm_force_x, y=norm_force_y, z=norm_force_z, line=dict(width=10, color='#555'), mode='lines') 

bend_trace = go.Scatter3d(name="Bending force",
    x=bend_force_x, y=bend_force_y, z=bend_force_z, line=dict(width=5, color='pink'), mode='lines') 

vert_trace=go.Scatter3d(name="Cell Vertices",
    x=vertx, 
    y=verty, 
    z=vertz, 
    mode='markers',
    marker=dict(
        size=8,
        color='blue',                # set color to an array/list of desired values
        opacity=0.8
    ),
    )



center_trace = go.Scatter3d(name="Cell centers",
    x=centx, 
    y=centy, 
    z=centz, 
    mode='markers',
    marker=dict(
        size=15,
        color='red',                # set color to an array/list of desired values
        opacity=0.8
    ))

axis_trace = go.Scatter3d(name="Axis of tube",
    x=axesx, 
    y=axesy, 
    z=axesz, 
    mode='markers',
    marker=dict(
        size=10,
        color='green',                # set color to an array/list of desired values
        opacity=0.8
    ))


# fig = go.Figure(data=[edge_trace, norm_trace, vert_trace, center_trace])
fig = go.Figure(data=[edge_trace, force_trace, bend_trace, vert_trace, center_trace, axis_trace])

# Show the figure
fig.show()

