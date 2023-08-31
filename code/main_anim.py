import numpy as np
import funcs as fu
import globals as const
import plotly.graph_objs as go


#used to plot all of the development of the tube in a plotly animation

test_center = 8
test_vertex = 8

# Generate cylinder
gvertices, gcenters, gaxes = fu.def_graphs()
fu.normals(gvertices)


all_data = fu.simulate(gvertices, gcenters, gaxes, const.steps)
steps = range(const.steps)

# Create figure

# make data

step=0

gaxes, gcenters, gvertices = all_data[step]
gaxes0, gcenters0, gvertices0=gaxes.copy(), gcenters.copy(), gvertices.copy()

edge_x, edge_y, edge_z = [], [], []
for idx0, idx1 in gvertices0.edges():
    x0, y0, z0 = gvertices0.nodes[idx0]["coords"]
    x1, y1, z1 = gvertices0.nodes[idx1]["coords"]
    edge_x.extend([x0, x1, None])
    edge_y.extend([y0, y1, None])
    edge_z.extend([z0, z1, None])

vertx, verty, vertz = fu.coords(gvertices0)
centx, centy, centz = fu.coords(gcenters0)
axesx, axesy, axesz = fu.coords(gaxes0)



edge_trace = go.Scatter3d(name="Cell edges",
    x=edge_x.copy(), y=edge_y.copy(), z=edge_z.copy(), line=dict(width=10, color='#888'), mode='lines')

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
        size=12,
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

locdata = [edge_trace, vert_trace, center_trace, axis_trace]
fig = go.Figure(data=locdata.copy())
    
# Frames

frames = []
for stepid in steps:
    gaxesx, gcentersx, gverticesx = all_data[stepid]
    gaxesloc, gcentersloc, gverticesloc = gaxesx.copy(), gcentersx.copy(), gverticesx.copy()
    edge_x, edge_y, edge_z = [], [], []
    for idx0, idx1 in gverticesloc.edges():
        x0, y0, z0 = gverticesloc.nodes[idx0]["coords"]
        x1, y1, z1 = gverticesloc.nodes[idx1]["coords"]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        edge_z.extend([z0, z1, None])

    vertx, verty, vertz = fu.coords(gverticesloc)
    centx, centy, centz = fu.coords(gcentersloc)
    axesx, axesy, axesz = fu.coords(gaxesloc)



    edge_trace = go.Scatter3d(name="Cell edges",
        x=edge_x.copy(), y=edge_y.copy(), z=edge_z.copy(), line=dict(width=10, color='#888'), mode='lines')

    vert_trace=go.Scatter3d(name="Cell Vertices",
        x=vertx.copy(), 
        y=verty.copy(), 
        z=vertz.copy(), 
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
            size=12,
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
    
    if stepid%50==0:
        print("rendering step ", stepid, " / ", const.steps)
    locdata = [edge_trace, vert_trace, center_trace, axis_trace].copy()
    frames.append(go.Frame(data=locdata.copy(), name=f'step {stepid}'))

fig.update(frames=frames)


def frame_args(duration):
    return {
            "frame": {"duration": duration, "redraw": True},
            "mode": "immediate",
            "fromcurrent": True,
            "transition": {"duration": duration, "easing": "linear"},
            }


sliders = [
    {"pad": {"b": 10, "t": 60},
     "len": 0.9,
     "x": 0.1,
     "y": 0,
     
     "steps": [
                 {"args": [[f.name], frame_args(0)],
                  "label": str(k),
                  "method": "animate",
                  } for k, f in enumerate(fig.frames)
              ]
     }
        ]

fig.update_layout(

    updatemenus = [{"buttons":[
                    {
                        "args": [None, frame_args(0)],
                        "label": "Play", 
                        "method": "animate",
                    },
                    {
                        "args": [[None], frame_args(0)],
                        "label": "Pause", 
                        "method": "animate",
                  }],
                    
                "direction": "left",
                "pad": {"r": 10, "t": 70},
                "type": "buttons",
                "x": 0.1,
                "y": 0,
            }
         ],
         sliders=sliders
    )


fig.update_layout(sliders=sliders)
fig.show()

