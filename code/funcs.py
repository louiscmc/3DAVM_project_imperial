import numpy as np
import globals as const
import networkx as nx
import pandas as pd
import plotly.express as px
from random import random, randint
from math import atan2, sqrt, cos, sin
from operator import itemgetter


np.random.seed(123) # to repeat experiences

def coords(graph):
    x=[]
    y=[]
    z=[]
    for k in range(len(graph.nodes)):
        x.append(graph.nodes[k]['coords'][0])
        y.append(graph.nodes[k]['coords'][1])
        z.append(graph.nodes[k]['coords'][2])
    return x,y,z


def hexagon(coords):
    all_points = [coords]
    all_points.append([coords[0]+const.width/(const.n_points_hor), coords[1]+const.length/const.n_points_vert])
    all_points.append([coords[0]+const.width/(const.n_points_hor), coords[1]+3*const.length/const.n_points_vert])
    all_points.append([coords[0], coords[1]+4*const.length/const.n_points_vert])
    all_points.append([coords[0]-const.width/(const.n_points_hor), coords[1]+3*const.length/const.n_points_vert])    
    all_points.append([coords[0]-const.width/(const.n_points_hor), coords[1]+const.length/const.n_points_vert])
    return all_points

def coords_center(loc_vertices, graph): # in 3D !
    l=len(loc_vertices)
    return [sum([graph.nodes[loc_vertices[i]]['coords'][0] for i in range(l)])/l, sum([graph.nodes[loc_vertices[i]]['coords'][1] for i in range(l)])/l, sum([graph.nodes[loc_vertices[i]]['coords'][2] for i in range(l)])/l]

def def_graphs():
    """
    Generate the graph of vertices, with connecting edges.
    """
    gvertices = nx.Graph()
    gcenters = nx.Graph()
    gaxes = nx.Graph()
    for axis_point in range(const.n_points_vert//3-1):
        gaxes.add_node(axis_point, **{"coords" : [], 'centers': []})
        for cell in range(const.n_points_hor//2):
            point = [const.width*(2*cell+1+(axis_point%2))/const.n_points_hor, const.length*3*axis_point/const.n_points_vert]
            loc_vertices = []
            hex_points = hexagon(point)
            center = gcenters.__len__()
            gaxes.nodes[axis_point]['centers'].append(center)
            gcenters.add_node(center, **{"axis":axis_point, "divisions":0})
            for k in hex_points:
                tmp= [x for x,y in gvertices.nodes(data=True) if np.abs(y['coords'][0]-k[0])%const.width<const.epsilon and np.abs(y['coords'][1]-k[1])<const.epsilon] # to account for float errors
                if not tmp: # if vertex not already in the graph, empty list is False
                    gvertices.add_node(gvertices.__len__(), **{"coords" : k, 'centers': [center]})
                    loc_vertices.append(gvertices.__len__()-1)
                else:
                    loc_vertices.append(tmp[0])
                    gvertices.nodes[tmp[0]]['centers'].append(center)
            for vert in range(6):
                gvertices.add_edge(loc_vertices[vert], loc_vertices[(vert+1)%6])  
            gcenters.nodes[center]['vertices']=loc_vertices
    
    cylinder_radius = const.width/(2*np.pi)  
    for vertex in gvertices.nodes:
        att=gvertices.nodes[vertex]
        theta = 2 * np.pi * att['coords'][0] / const.width
        z = att['coords'][1]
        x = cylinder_radius * np.cos(theta)
        y = cylinder_radius * np.sin(theta)
        att['coords']=[x+const.random_noise*random()+z*const.tilting/const.length,y+const.random_noise*random(),z+const.random_noise*random()]

        #boundary vertices
        if (z < const.length*2.5/const.n_points_vert) or (const.length-z < const.length*3.5/const.n_points_vert):
            att['fixed']=True
        else:
            att['fixed']=False
    for center in gcenters.nodes:
        att=gcenters.nodes[center]
        gcenters.nodes[center]['coords']=coords_center(att['vertices'], gvertices)
        if (att['coords'][0]-att['coords'][2]*const.tilting/const.length<0 and att['coords'][2]>const.length/2):
            gcenters.nodes[center]['scale']=scale_size(att['coords'][2])*const.scale_artif
        elif (att['coords'][0]-att['coords'][2]*const.tilting/const.length>0 and att['coords'][2]<const.length/2):
            gcenters.nodes[center]['scale']=scale_size(att['coords'][2])*const.scale_artif
        else:
            gcenters.nodes[center]['scale']=scale_size(att['coords'][2])
    add_area_perim(gcenters, gvertices)
    for axis_point in gaxes.nodes:
        att = gaxes.nodes[axis_point]
        gaxes.nodes[axis_point]['coords']=coords_center(att['centers'], gcenters)

    return gvertices, gcenters, gaxes


def dist(p1, p2):
    # works in 2D and 3D
    return np.sqrt(np.sum((np.array(p1) - np.array(p2))**2))

def triangle_area(p1, p2, p3):
    # heron formula
    a = dist(p1, p2)
    b = dist(p2, p3)
    c = dist(p3, p1)
    s = (a+b+c)/2
    return np.sqrt(s*(s-a)*(s-b)*(s-c))

def center_area_perim(center, gcenters, gvertices):

    # the area and perimeter of a given point / region
    vert_idx = gcenters.nodes[center]['vertices']
    center = gcenters.nodes[center]['coords']
    vertices = []
    num_vertices = len(vert_idx)
    for i in vert_idx:
        vertices.append(gvertices.nodes[i]['coords'])
    perim = 0
    total_area = 0
    for i in range(len(vertices)):
        #area 
        total_area += triangle_area(center, vertices[i], vertices[(i+1)%num_vertices])
        #perim
        perim += dist(vertices[i],vertices[(i+1)%num_vertices])
    return total_area, perim

def add_area_perim(gcenters, gvertices):
    for center in gcenters.nodes():
        c_area, c_perim = center_area_perim(center, gcenters, gvertices)
        gcenters.nodes[center]['area']=c_area
        gcenters.nodes[center]['perimeter']=c_perim

def vertex_normal(G, vertex):
    # Get the neighbors of the vertex in the graph
    neighbors = list(G.neighbors(vertex))

    # Check if vertex has exactly three neighbors
    if len(neighbors) != 3:
        return np.array([0.,0.,0.])
    # Calculate the position of points at a fixed distance along each ridge
    points_on_ridges = [np.array(G.nodes[vertex]['coords']) + (np.array(G.nodes[neighbor]['coords']) - np.array(G.nodes[vertex]['coords'])) / np.linalg.norm(np.array(G.nodes[neighbor]['coords']) - np.array(G.nodes[vertex]['coords'])) for neighbor in neighbors]

    # Compute the center of mass of the triangle
    centroid = np.mean(points_on_ridges, axis=0)
    # Calculate the normal vector from the vertex to the centroid
    vertex_normal = np.array(G.nodes[vertex]['coords']) - centroid

    return vertex_normal

def normals(G):
    for vertex in G.nodes:
        G.nodes[vertex]['normal']=vertex_normal(G, vertex)

def divisions(G):
    tot_div=0
    for i in G.nodes:
        tot_div+=G.nodes[i]["divisions"]
    return tot_div//2

def ref_area(step, scale):
    return const.base_area*(1+const.growth_rate*step)*(1-(1-scale)/2)

def ref_perim(step, scale):
    return const.base_perimeter*np.sqrt(1+const.growth_rate*step)*(1-(1-scale)/2)**2

def ref_volume(step):
    return const.base_volume*(1+const.beating_intensity*np.sin(step*2*3.141592/const.heart_rhythm))
# division

def division(gvertices, gcenters, gaxes, center):
    vertices = gcenters.nodes[center]["vertices"]
    num_vert = len(vertices)
    if num_vert>3:
        coords_z=[[vertex, gvertices.nodes[vertex]["coords"][2]] for vertex in vertices]
        for vertex in vertices:
            coords_z.append([vertex, gvertices.nodes[vertex]["coords"][2]])
        top_cell = randint(0,num_vert)
        vertices = vertices[top_cell:]+vertices[:top_cell]
        other_vert = vertices[num_vert//2] # take the opposite one in the cell
        newcenter = gcenters.__len__()
        axis_point = gcenters.nodes[center]["axis"]
        alpha = gvertices.__len__()
        beta = alpha+1 
        centers_alpha = [center for center in gvertices.nodes[vertices[0]]['centers'] if center in gvertices.nodes[vertices[-1]]['centers']]
        centers_alpha.append(newcenter)
        centers_beta = [center for center in gvertices.nodes[vertices[num_vert//2]]['centers'] if center in gvertices.nodes[vertices[(num_vert//2)-1]]['centers']]
        centers_beta.append(newcenter)
        # proper division : 1/ setting up centers, axes, and vertices
        gcenters.nodes[center]['divisions']+=1
        gcenters.add_node(newcenter, **{"axis":axis_point, "divisions":gcenters.nodes[center]['divisions'], "scale":gcenters.nodes[center]['scale'], "area":gcenters.nodes[center]['area'], "perimeter":gcenters.nodes[center]['perimeter']})
        gvertices.add_node(alpha, **{"coords" : ((gvertices.nodes[vertices[0]]['coords'] + gvertices.nodes[vertices[-1]]['coords']))/2, 'centers': centers_alpha, 'fixed':gvertices.nodes[vertices[0]]['fixed'], 'force':[0,0,0]})
        gvertices.add_node(beta, **{"coords" : ((gvertices.nodes[vertices[num_vert//2]]['coords'] + gvertices.nodes[vertices[(num_vert//2)-1]]['coords']))/2, 'centers': centers_beta, 'fixed':gvertices.nodes[vertices[num_vert//2]]['fixed'], 'force':[0,0,0]})
        gaxes.nodes[axis_point]["centers"].append(newcenter)
        gcenters.nodes[center]["vertices"]=[alpha]+vertices[0:num_vert//2]+[beta]
        gcenters.nodes[newcenter]["vertices"]=[alpha]+[beta]+vertices[num_vert//2:]
        if len(centers_alpha)==3:
            kappa = [centeri for centeri in centers_alpha if (centeri != center and centeri != newcenter)][0]
            for idx_loc in range(len(gcenters.nodes[kappa]['vertices'])):
                if gcenters.nodes[kappa]['vertices'][idx_loc]==vertices[0]:
                    gcenters.nodes[kappa]['vertices']= gcenters.nodes[kappa]['vertices'][:idx_loc+1]+[alpha]+gcenters.nodes[kappa]['vertices'][idx_loc+1:]
                    break
        if len(centers_beta)==3:
            lambd = [centeri for centeri in centers_beta if (centeri != center and centeri != newcenter)][0]
            for idx_loc in range(len(gcenters.nodes[lambd]['vertices'])):
                if gcenters.nodes[lambd]['vertices'][idx_loc]==vertices[num_vert//2]:
                    gcenters.nodes[lambd]['vertices']= gcenters.nodes[lambd]['vertices'][:idx_loc+1]+[beta]+gcenters.nodes[lambd]['vertices'][idx_loc+1:]
                    break
        for vertex in vertices[num_vert//2:]:
            for idx_center in range(len(gvertices.nodes[vertex]["centers"])):
                if gvertices.nodes[vertex]["centers"][idx_center]==center:
                    gvertices.nodes[vertex]["centers"].pop(idx_center)
                    break
            gvertices.nodes[vertex]["centers"].append(newcenter)
        # proper division : 2/ setting up edges of vertices
        gvertices.remove_edge(vertices[0], vertices[-1])
        gvertices.remove_edge(vertices[num_vert//2], vertices[(num_vert//2)-1])
        gvertices.add_edge(vertices[-1], alpha)
        gvertices.add_edge(vertices[0], alpha)
        gvertices.add_edge(beta, alpha)
        gvertices.add_edge(vertices[num_vert//2], beta)
        gvertices.add_edge(vertices[(num_vert//2)-1], beta)


# forces

def scale_size(z):
    #if 0.25<z/const.length<0.75:
    #    return 0.75+abs(z/const.length-0.5)
    return 1

def force_dir(point1, point2):
    # unity vector from point 1 to point 2
    force_dir_loc = np.array(point2)-np.array(point1)
    force_dir_loc/= np.linalg.norm(force_dir_loc)
    return force_dir_loc

def vertex_forces(gvertices, vertex):
    edges = gvertices.edges(vertex)
    coords = gvertices.nodes[vertex]['coords']
    for edge in edges:
        gvertices.nodes[vertex]['force']-=const.sigma*(np.array(gvertices.nodes[edge[1]]['coords'])-np.array(coords))

def area_forces(gvertices, gcenters, vertex, center, area, good_edges, coords, step, scale):
    local_area_width = np.linalg.norm(np.array(gvertices.nodes[good_edges[1]]['coords'])-np.array(gvertices.nodes[good_edges[0]]['coords']))
    loc_force_dir = force_dir(gcenters.nodes[center]['coords'], coords)
    gvertices.nodes[vertex]['force']-= const.a_alpha*(area-ref_area(step,scale))*loc_force_dir*local_area_width/1.5**gcenters.nodes[center]["divisions"]

def volume_forces(gvertices, gcenters, gaxes, center, vertex, coords, volume, divs, step):
    loc_force_dir = force_dir(gaxes.nodes[gcenters.nodes[center]["axis"]]['coords']+np.array([0.05,0,0]), coords) # +np.array([0.15,0,0])
    gvertices.nodes[vertex]['force']-= const.vol_alpha*loc_force_dir/3*(volume-ref_volume(step)*(const.n_cells_tot+divs)/const.n_cells_tot)

def innervolume(gcenters, gaxes, step):
    volume=0
    for axis_point in gaxes.nodes:
     # the area and perimeter of a given point / region
        cent_idx = gaxes.nodes[axis_point]['centers']
        coords_axis = gaxes.nodes[axis_point]['coords']
        centers = []
        num_centers = len(cent_idx)
        for i in cent_idx:
            centers.append(gcenters.nodes[i]['coords'])
        # perim = 0
        loc_area = 0
        for i in range(len(centers)):
            #area 
            loc_area += triangle_area(coords_axis, centers[i], centers[(i+1)%num_centers])
            #perim
            # perim += dist(centers[i],centers[(i+1)%centers])
        volume += loc_area*const.length/const.n_cells_vert
    print('volume : ',volume, ' ref_volume : ', ref_volume(step))
    return volume

def bending_forces(gvertices, gcenters, gaxes, vertex, center):
    norm_bend_loc = np.linalg.norm(gvertices.nodes[vertex]['normal']*const.c_alpha)
    loc_base_norm_dir = -force_dir(gaxes.nodes[gcenters.nodes[center]["axis"]]['coords'], gvertices.nodes[vertex]['coords'])
    loc_base_norm_dir /= np.linalg.norm(loc_base_norm_dir)
    return -((gvertices.nodes[vertex]['normal']+ loc_base_norm_dir * 0.0856)*const.c_alpha)*norm_bend_loc*const.nu/0.05

def compute_force(gvertices, gcenters, gaxes, step):
    volume = innervolume(gcenters, gaxes, step)
    for vertex in gvertices.nodes():
        #initialisation and bending, vertex forces
        coords = gvertices.nodes[vertex]['coords']
        centers = gvertices.nodes[vertex]['centers']
        gvertices.nodes[vertex]['force']=bending_forces(gvertices, gcenters, gaxes, vertex, gvertices.nodes[vertex]['centers'][0])
        if const.sigma != 0:
            vertex_forces(gvertices, vertex)
        for center in centers:
            scale = gcenters.nodes[center]['scale']
            A_alpha = gcenters.nodes[center]['area']
            L_alpha = gcenters.nodes[center]['perimeter']
            edges = gvertices.edges(vertex)
            coords = gvertices.nodes[vertex]['coords']
            # volume force
            if const.vol_alpha !=0 :
                divs=divisions(gcenters)
                volume_forces(gvertices, gcenters, gaxes, center, vertex, coords, volume, divs, step)
            # good edges : edges to the vertex, on the cell
            good_edges=[]
            for edge in edges:
                if edge[1] in gcenters.nodes[center]['vertices']:
                    good_edges.append(edge[1])
            # Perimeter force
            if len(good_edges)==2:
                gvertices.nodes[vertex]['force']-=2*const.b_alpha*(L_alpha-ref_perim(step, scale))*(force_dir(gvertices.nodes[good_edges[0]]['coords'], coords)+force_dir(gvertices.nodes[good_edges[1]]['coords'], coords))
                """ if gvertices.nodes[vertex]['coords'][0]<0:
                    gvertices.nodes[vertex]['force']+=2*const.b_alpha*L_alpha*(force_dir(gvertices.nodes[good_edges[0]]['coords'], coords)+force_dir(gvertices.nodes[good_edges[1]]['coords'], coords))
                else:
                    gvertices.nodes[vertex]['force']+=0.2*2*const.b_alpha*L_alpha*(force_dir(gvertices.nodes[good_edges[0]]['coords'], coords)+force_dir(gvertices.nodes[good_edges[1]]['coords'], coords))
 """
            else:
                raise(ValueError('Wrong number of neighbours'))
            # Area force
            area_forces(gvertices, gcenters, vertex, center, A_alpha, good_edges, coords, step, scale)
    bend_avg=0
    for vertex in gvertices.nodes():
        bend_avg+=np.linalg.norm(bending_forces(gvertices, gcenters, gaxes, vertex, gvertices.nodes[vertex]['centers'][0]))
    print('average of bending forces : ', bend_avg/len(gvertices.nodes()))

def compute_length_ij(gvertices, vertex1, vertex2):
    coords1 = np.array(gvertices.nodes[vertex1]['coords'])
    coords2 = np.array(gvertices.nodes[vertex2]['coords'])
    length = np.linalg.norm(coords1 - coords2)
    if np.isnan(length):
        print(f'NaN value found in length between vertices {vertex1} and {vertex2}')
    return length

def reorder_vertices(vertices, coords):
    center = np.mean(coords, axis=0)
    angles = np.arctan2(coords[:,1] - center[1], coords[:,0] - center[0])
    return [vertices[i] for i in np.argsort(angles)]

def compute_A_hat(gvertices, gcenters, center):
    vertices = gcenters.nodes[center]['vertices']
    coords = np.array([gvertices.nodes[vertex]['coords'] for vertex in vertices])
    vertices = reorder_vertices(vertices, coords)
    A_hat = np.zeros(3)
    for i in range(len(vertices)):
        cross_product = np.cross(gvertices.nodes[vertices[i]]['coords'], gvertices.nodes[vertices[(i+1)%len(vertices)]]['coords'])
        if np.isnan(cross_product).any():
            print(f'NaN value found in cross product for vertices {vertices[i]} and {vertices[(i+1)%len(vertices)]}')
        A_hat += cross_product
    A_hat /= 2
    if np.isnan(A_hat).any():
        print(f'NaN value found in A_hat for center {center}')
    return A_hat

def t1_transition(gcenters, gvertices, vertex1, vertex2, dist):
    if gvertices.nodes[vertex1]['coords'][2]>gvertices.nodes[vertex2]['coords'][2]:
    # making sure vertex1 is at the bottom to fit the schema
        _ = vertex1
        vertex1=vertex2
        vertex2=_

    centers1=gvertices.nodes[vertex1]['centers']
    centers2=gvertices.nodes[vertex2]['centers']
    #shared centers : alpha is linked only to v1, beta to v2, gamma and delta are shared.
    tmp = []
    for center in centers1:
        if center in centers2:
            tmp.append(center)
        else:
            c_alpha = center
    x1, y1 = gcenters.nodes[tmp[0]]['coords'][:2]
    x2, y2 = gcenters.nodes[tmp[1]]['coords'][:2]
    if (atan2(y1, x1)<atan2(y2, x2)) != (y1>0 and y2<0 and x1<0):
    # if the first center is more to the "left" (anglewise)
    # chirality !!!
        c_delta, c_gamma = tmp[0], tmp[1]
    else :
        c_delta, c_gamma = tmp[1], tmp[0]
    
    if len(gcenters.nodes[c_delta]['vertices'])==3 or len(gcenters.nodes[c_gamma]['vertices'])==3:
        print('3 vertices !')
        return True
    
    for center in centers2:
        if center not in tmp:
            c_beta = center
    
    # identifying vertices:

    for vertex in gcenters.nodes[c_alpha]['vertices']:
        if vertex in gcenters.nodes[c_delta]['vertices'] and vertex!=vertex1 and vertex!=vertex2:
            vertex1delta = vertex
        elif vertex in gcenters.nodes[c_gamma]['vertices'] and vertex!=vertex1 and vertex!=vertex2:
            vertex1gamma = vertex

    for vertex in gcenters.nodes[c_beta]['vertices']:
        if vertex in gcenters.nodes[c_delta]['vertices'] and vertex!=vertex1 and vertex!=vertex2:
            vertex2delta = vertex
        elif vertex in gcenters.nodes[c_gamma]['vertices'] and vertex!=vertex1 and vertex!=vertex2:
            vertex2gamma = vertex

    # operate transition
    for i in range(len(gcenters.nodes[c_delta]['vertices'])):
        if gcenters.nodes[c_delta]['vertices'][i]==vertex2:
            gcenters.nodes[c_delta]['vertices']=np.delete(gcenters.nodes[c_delta]['vertices'],i)
            break
    for i in range(len(gvertices.nodes[vertex2]['centers'])):
        if gvertices.nodes[vertex2]['centers'][i]==c_delta:
            gvertices.nodes[vertex2]['centers']=np.delete(gvertices.nodes[vertex2]['centers'],i)
            break
    for i in range(len(gcenters.nodes[c_gamma]['vertices'])):
        if gcenters.nodes[c_gamma]['vertices'][i]==vertex1:
            gcenters.nodes[c_gamma]['vertices']=np.delete(gcenters.nodes[c_gamma]['vertices'],i)
            break
    for i in range(len(gvertices.nodes[vertex1]['centers'])):
        if gvertices.nodes[vertex1]['centers'][i]==c_gamma:
            gvertices.nodes[vertex1]['centers']=np.delete(gvertices.nodes[vertex1]['centers'],i)
            break
    gcenters.nodes[c_beta]['vertices']=np.append(gcenters.nodes[c_beta]['vertices'],vertex1)
    gvertices.nodes[vertex1]['centers']=np.append(gvertices.nodes[vertex1]['centers'],c_beta)
    gcenters.nodes[c_alpha]['vertices']=np.append(gcenters.nodes[c_alpha]['vertices'],vertex2)
    gvertices.nodes[vertex2]['centers']=np.append(gvertices.nodes[vertex2]['centers'],c_alpha)
    gvertices.add_edge(vertex1, vertex2delta)
    gvertices.add_edge(vertex2, vertex1gamma)
    gvertices.remove_edge(vertex1, vertex1gamma)
    gvertices.remove_edge(vertex2, vertex2delta)

    #new coordinates, inducing chirality : slightly more to the left, impact is handled by const.t1coeff
    gvertices.nodes[vertex2]['coords']=(gvertices.nodes[vertex2]['coords']+gvertices.nodes[vertex1]['coords'])/2
    direction  = gcenters.nodes[c_delta]['coords']-gvertices.nodes[vertex2]['coords']
    direction/=np.linalg.norm(direction)
    direction*=dist
    gvertices.nodes[vertex1]['coords']=gvertices.nodes[vertex2]['coords']+const.t1_coeff*direction
    gvertices.nodes[vertex1]['force']=np.zeros(3)
    gvertices.nodes[vertex2]['force']=np.zeros(3)
    print('T1 transition between vertices ', vertex1, ' and ', vertex2)
    return True

def rotation(gvertices, vertex, angle):
    coords = gvertices.nodes[vertex]['coords']
    x,y = coords[0], coords[1]
    phi = atan2(y,x)
    r=sqrt(x**2+y**2)
    gvertices.nodes[vertex]['coords'][0]=r*cos(phi+angle)
    gvertices.nodes[vertex]['coords'][1]=r*sin(phi+angle)

def update_positions(gvertices, gcenters, gaxes, step, steps):
    compute_force(gvertices, gcenters, gaxes, step)
    norm_avg = 0
    for vertex in gvertices.nodes():
        norm_avg+=np.linalg.norm(gvertices.nodes[vertex]['force'])
    print('average of forces : ', norm_avg/len(gvertices.nodes()))
    for edge in gvertices.edges:
        loc_dist = dist(gvertices.nodes[edge[0]]['coords'], gvertices.nodes[edge[1]]['coords'])
        if loc_dist< const.t1_min_dist: # T1 transition only if vertices are too close
            t1_transition(gcenters, gvertices, edge[0], edge[1], loc_dist)
    for vertex in gvertices.nodes():
        if gvertices.nodes[vertex]['fixed']==False:
            scale=gcenters.nodes[gvertices.nodes[vertex]["centers"][0]]['scale']
            gvertices.nodes[vertex]['coords'] += scale**2*const.nu * gvertices.nodes[vertex]['force']
            #rotation(gvertices, vertex, 10*const.angle*(1-scale))
        else:
            # rotation(gvertices, vertex, const.angle)
            gvertices.nodes[vertex]['coords'] += const.displacement 
    if const.tau_div !=0:   
        tmp_centers = gcenters.copy()
        for center in tmp_centers.nodes: 
            if random()<1/const.tau_div/steps:
                division(gvertices, gcenters, gaxes, center)
    for center in gcenters.nodes:
        att=gcenters.nodes[center]
        gcenters.nodes[center]['coords']=coords_center(att['vertices'], gvertices)
    for axis_point in gaxes.nodes:
        att = gaxes.nodes[axis_point]
        gaxes.nodes[axis_point]['coords']=coords_center(att['centers'], gcenters)
    normals(gvertices)
    add_area_perim(gcenters, gvertices)

def simulate(gvertices, gcenters, gaxes, steps):
    for step in range(steps):
        print ("step ",step,' of ', steps)
        print('total divisions : ', divisions(gcenters))
        update_positions(gvertices, gcenters, gaxes, step, steps)