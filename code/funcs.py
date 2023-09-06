import numpy as np
import globals as const
import networkx as nx
from random import random, randint
from math import atan2, sqrt, cos, sin
from copy import deepcopy


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
    return np.array([sum([graph.nodes[loc_vertices[i]]['coords'][0] for i in range(l)])/l, sum([graph.nodes[loc_vertices[i]]['coords'][1] for i in range(l)])/l, sum([graph.nodes[loc_vertices[i]]['coords'][2] for i in range(l)])/l])

def def_graphs(scale_artif):
    """
    Generate the graph of vertices, with connecting edges.
    """
    gvertices = nx.Graph()
    gcenters = nx.Graph()
    gaxes = nx.Graph()
    for axis_point in range(const.n_cells_vert):
        gaxes.add_node(axis_point, **{"coords" : [], 'centers': []})
        for cell in range(const.n_cells_hor):
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
        gvertices.nodes[vertex]['t1']=0
        att=gvertices.nodes[vertex]
        theta = 2 * np.pi * att['coords'][0] / const.width
        z = att['coords'][1]
        x = cylinder_radius * np.cos(theta)
        y = cylinder_radius * np.sin(theta)
        att['coords']=[x+const.random_noise*random()+z*const.tilting/const.length,y+const.random_noise*random(),z+const.random_noise*random()]
        #boundary vertices
        if (z < const.length/2/const.n_cells_vert) or (const.length-z < const.length/2/const.n_cells_vert*2):
            att['fixed']=True
        else:
            att['fixed']=False

    for center in gcenters.nodes:
        scale_size(gcenters, gvertices, center, scale_artif)

    for axis_point in gaxes.nodes:
        att = gaxes.nodes[axis_point]
        gaxes.nodes[axis_point]['coords']=coords_center(att['centers'], gcenters)

    add_area_perim(gaxes, gcenters, gvertices)

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

def center_area_perim_normal(center, gaxes, gcenters, gvertices):

    # the area and perimeter of a given point / region
    vert_idx = gcenters.nodes[center]['vertices']
    center_coords = gcenters.nodes[center]['coords']
    vertices = []
    num_vertices = len(vert_idx)
    for i in vert_idx:
        vertices.append(gvertices.nodes[i]['coords'])
    perim = 0
    total_area = 0
    normal=np.array([0.,0.,0.])
    for i in range(len(vertices)):
        #area 
        total_area += triangle_area(center_coords, vertices[i], vertices[(i+1)%num_vertices])
        #perim
        perim += dist(vertices[i],vertices[(i+1)%num_vertices])
        #normal
        normal += np.cross(np.array(center_coords)-np.array(vertices[i]), np.array(center_coords)-np.array(vertices[(i+1)%num_vertices]))
    if np.linalg.norm(normal)!= 0:
        normal/np.linalg.norm(normal)
    if np.dot(normal, np.array(center_coords)-gaxes.nodes[gcenters.nodes[center]["axis"]]['coords'])<0:
        normal= -normal
    return total_area, perim, normal

def add_area_perim(gaxes,gcenters, gvertices):
    for center in gcenters.nodes():
        c_area, c_perim, c_normal = center_area_perim_normal(center, gaxes, gcenters, gvertices)
        gcenters.nodes[center]['area']=c_area
        gcenters.nodes[center]['perimeter']=c_perim
        gcenters.nodes[center]['normal']=c_normal

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
    return const.base_area*(1+const.growth_rate*step)*scale

def ref_perim(step, scale):
    return const.base_perimeter*np.sqrt(1+const.growth_rate*step)*np.sqrt(scale)

def ref_volume(step, rhythm):
    if const.beating_volume_change !=0:
        return const.base_volume*(1+2*step/const.steps+const.beating_volume_change*np.sin(step*2*3.141592/rhythm))
    # return const.base_volume*(1+2*step/const.steps)
    return const.base_volume


def division(gvertices, gcenters, gaxes, center):
    vertices = np.array(gcenters.nodes[center]["vertices"]).tolist()
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
                    gvertices.nodes[vertex]["centers"]=np.delete(gvertices.nodes[vertex]["centers"], idx_center)
                    break
            gvertices.nodes[vertex]["centers"]=np.insert(gvertices.nodes[vertex]["centers"],0,newcenter)
        # proper division : 2/ setting up edges of vertices
        gvertices.remove_edge(vertices[0], vertices[-1])
        gvertices.remove_edge(vertices[num_vert//2], vertices[(num_vert//2)-1])
        gvertices.add_edge(vertices[-1], alpha)
        gvertices.add_edge(vertices[0], alpha)
        gvertices.add_edge(beta, alpha)
        gvertices.add_edge(vertices[num_vert//2], beta)
        gvertices.add_edge(vertices[(num_vert//2)-1], beta)


# forces

def scale_size(gcenters, gvertices, center, scale_artif):
    # using the same cell shapes as in michaela's paper
    att=gcenters.nodes[center]
    gcenters.nodes[center]['coords']=coords_center(att['vertices'], gvertices)
    x=att['coords'][0]
    z=att['coords'][2]
    if (x-z*const.tilting/const.length<0.1 and att['coords'][2]>const.length*0.61) or (x-z*const.tilting/const.length>-0.1 and z<const.length*0.45):
        gcenters.nodes[center]['scale']=(1-1.8*scale_artif)
    elif (x-z*const.tilting/const.length>0.1 and const.length*0.9>z>const.length*0.61) or (x-z*const.tilting/const.length<-0.1 and const.length*0.45>z>const.length*0.1):
        gcenters.nodes[center]['scale']=1+1.2*scale_artif
    elif (const.length*0.61>z>const.length*0.45):
        gcenters.nodes[center]['scale']=1-2*scale_artif
    else:
        gcenters.nodes[center]['scale']=1

def force_dir(point1, point2):
    # unit vector from point 1 to point 2
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
    gvertices.nodes[vertex]['force']-= const.a_alpha*(area-ref_area(step,scale))*loc_force_dir*local_area_width

def volume_forces(gvertices, vertex, volume, divs, step, scale, c_normal, vol_alpha, rhythm):
    # gvertices.nodes[vertex]['force']-= scale*vol_alpha*c_normal/3*(volume-scale*ref_volume(step, rhythm)*(const.n_cells_tot+divs)/const.n_cells_tot)
    gvertices.nodes[vertex]['force']-= vol_alpha*c_normal/3*(volume-ref_volume(step, rhythm)*(const.n_cells_tot+divs)/const.n_cells_tot)
    
def innervolume(gcenters, gaxes, step, rhythm):
    volume=0
    if step==0:
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
        # print('volume : ',volume, ' ref_volume : ', ref_volume(step, rhythm))
        return volume
    else:
        volume=0
        for axispoint in gaxes.nodes:
            if axispoint!=const.n_cells_vert-1:
                volume+= gaxes.nodes[axispoint]['area']*dist(gaxes.nodes[axispoint]['coords'], gaxes.nodes[axispoint+1]['coords'])
        # print('volume : ',volume, ' ref_volume : ', ref_volume(step, rhythm))
        return volume

def bending_forces(gvertices, gcenters, gaxes, vertex, center):
    norm_bend_loc = np.linalg.norm(gvertices.nodes[vertex]['normal']*const.c_alpha)
    loc_base_norm_dir = -force_dir(gaxes.nodes[gcenters.nodes[center]["axis"]]['coords'], gvertices.nodes[vertex]['coords'])
    loc_base_norm_dir /= np.linalg.norm(loc_base_norm_dir)
    return -((gvertices.nodes[vertex]['normal']+ loc_base_norm_dir * 0.0856)*const.c_alpha)*norm_bend_loc*const.nu/0.05

def compute_force(gvertices, gcenters, gaxes, step, vol_alpha, rhythm):
    volume = innervolume(gcenters, gaxes, step, rhythm)
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
            if vol_alpha !=0 :
                divs=divisions(gcenters)
                c_normal = gcenters.nodes[center]['normal']
                volume_forces(gvertices, vertex, volume, divs, step, scale, c_normal, vol_alpha, rhythm)
            # good edges : edges to the vertex, on the cell
            good_edges=[]
            for edge in edges:
                if edge[1] in gcenters.nodes[center]['vertices']:
                    good_edges.append(edge[1])
            # Perimeter force
            if len(good_edges)==2:
                gvertices.nodes[vertex]['force']-=2*const.b_alpha*(L_alpha-ref_perim(step, scale))*(force_dir(gvertices.nodes[good_edges[0]]['coords'], coords)+force_dir(gvertices.nodes[good_edges[1]]['coords'], coords))
            else:
                print(gvertices.nodes[vertex])
                raise(ValueError('Wrong number of neighbours for vertex ', vertex))
            # Area force
            area_forces(gvertices, gcenters, vertex, center, A_alpha, good_edges, coords, step, scale)
    bend_avg=0
    for vertex in gvertices.nodes():
        bend_avg+=np.linalg.norm(bending_forces(gvertices, gcenters, gaxes, vertex, gvertices.nodes[vertex]['centers'][0]))
    if np.isnan(gvertices.nodes[0]['force'][0]):
        print('NaN force !')
        raise ValueError('NaN found !')
    # print('average of bending forces : ', bend_avg/len(gvertices.nodes()))


def t1_transition(gcenters, gvertices, vertex1, vertex2, dist):
    if gvertices.nodes[vertex1]['coords'][2]>gvertices.nodes[vertex2]['coords'][2]:
    # making sure vertex1 is at the bottom 
        _ = vertex1
        vertex1=vertex2
        vertex2=_

    centers1=gvertices.nodes[vertex1]['centers']
    centers2=gvertices.nodes[vertex2]['centers']
    #shared centers : alpha is linked only to v1, beta to v2, gamma and delta are shared.
    tmp = []
    for center_ in centers1:
        if center_ in centers2:
            tmp.append(center_)
        else:
            c_alpha = center_
    x1, y1 = gcenters.nodes[tmp[0]]['coords'][:2]
    x2, y2 = gcenters.nodes[tmp[1]]['coords'][:2]
    if (atan2(y1, x1)<atan2(y2, x2)) != (y1>0 and y2<0 and x1<0):
    # if the first center is more to the "left" (anglewise)
        c_delta, c_gamma = tmp[0], tmp[1]
    else :
        c_delta, c_gamma = tmp[1], tmp[0]
    
    for center_ in centers2:
        if center_ not in tmp:
            c_beta = center_
    
    # identifying vertices:

    for vertex in gcenters.nodes[c_alpha]['vertices']:
        if vertex in gcenters.nodes[c_gamma]['vertices'] and vertex!=vertex1 and vertex!=vertex2:
            vertex1gamma = vertex

    for vertex in gcenters.nodes[c_beta]['vertices']:
        if vertex in gcenters.nodes[c_delta]['vertices'] and vertex!=vertex1 and vertex!=vertex2:
            vertex2delta = vertex

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

    # rerorder vertices in cells for area and perimeter calculations
    reorder_cell_vertices(gcenters, gvertices, c_alpha)
    reorder_cell_vertices(gcenters, gvertices, c_beta)
    reorder_cell_vertices(gcenters, gvertices, c_gamma)
    reorder_cell_vertices(gcenters, gvertices, c_delta)

    # print('T1 transition between vertices ', vertex1, ' and ', vertex2)
    return True

def rotation(gvertices, vertex, angle):
    coords = gvertices.nodes[vertex]['coords']
    x,y = coords[0], coords[1]
    phi = atan2(y,x)
    r=sqrt(x**2+y**2)
    gvertices.nodes[vertex]['coords'][0]=r*cos(phi+angle)
    gvertices.nodes[vertex]['coords'][1]=r*sin(phi+angle)

def axis_normal(gaxes, axispoint):
    if axispoint!=0 and axispoint!=const.n_cells_vert-1:
        vect = gaxes.nodes[axispoint+1]['coords']-gaxes.nodes[axispoint-1]['coords']
        return vect
    else : return np.array([0,0,1])

def rotation_matrix_from_vectors(vect):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vect: the 3d "source" vector
    :return rotation_matrix: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    z_norm = np.array([0,0,1])
    vec_norm= (vect / np.linalg.norm(vect)).reshape(3)
    v = np.cross(vec_norm,z_norm)
    c = np.dot(vec_norm, z_norm)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def reorder_around_and_area(gaxes, gcenters, axispoint):
    loc_origin = gaxes.nodes[axispoint]['coords']
    loc_coords = []
    loc_centers = gaxes.nodes[axispoint]['centers']
    rot_mat = rotation_matrix_from_vectors(axis_normal(gaxes, axispoint))
    for center_ in gaxes.nodes[axispoint]['centers']:
        loc_coords.append(rot_mat.dot(gcenters.nodes[center_]['coords']-loc_origin))
    angles = [np.arctan2(k[1], k[0]) for k in loc_coords]
    gaxes.nodes[axispoint]['centers'] = [loc_centers[i] for i in np.argsort(angles)]
    loc_centers = gaxes.nodes[axispoint]['centers']
    area = 0
    for i in range(len(loc_centers)):
        v1 = gcenters.nodes[loc_centers[i]]['coords']
        v2 = gcenters.nodes[loc_centers[(i+1) % len(loc_centers)]]['coords']
        area += v1[0]*v2[1] - v1[1]*v2[0]
    return abs(area/2)

def reorder_cell_vertices(gcenters, gvertices, center):
    vertex0= gcenters.nodes[center]['vertices'][0]
    new_vertices = [vertex0]
    edges0 = gvertices.edges(vertex0)
    for edge in edges0:
        if edge[0]==vertex0 and edge[1] in gcenters.nodes[center]['vertices']:
            vertex1=edge[1]
            new_vertices.append(vertex1)
            break
        elif edge[1]==vertex0 and edge[0] in gcenters.nodes[center]['vertices']:
            vertex1=edge[0]
            new_vertices.append(vertex1)
            break
    for k in range(1,len(gcenters.nodes[center]['vertices'])-1):
        edges_ = gvertices.edges(new_vertices[k])
        for edge in edges_:
            if edge[0]==new_vertices[k] and edge[1] in gcenters.nodes[center]['vertices'] and edge[1]!=new_vertices[k-1]:
                vertex_=edge[1]
                new_vertices.append(vertex_)
                break
            elif edge[1]==new_vertices[k] and edge[0] in gcenters.nodes[center]['vertices'] and edge[0]!=new_vertices[k-1]:
                vertex_=edge[1]
                new_vertices.append(vertex_)
                break
    gcenters.nodes[center]['vertices']=new_vertices
    

def reassign_centers(gcenters, gaxes): # if a center is closer to a new axis point, it is reassigned
    for center in gcenters.nodes:
        for axispoint in gaxes.nodes:
            if dist(gcenters.nodes[center]['coords'], gaxes.nodes[gcenters.nodes[center]['axis']]['coords']) > const.reassign_margin*dist(gcenters.nodes[center]['coords'], gaxes.nodes[axispoint]['coords']):
                print("center ", center, " has switched to axis point ", axispoint)
                gaxes.nodes[gcenters.nodes[center]['axis']]['centers']=np.delete(gaxes.nodes[gcenters.nodes[center]['axis']]['centers'], np.where(gaxes.nodes[gcenters.nodes[center]['axis']]['centers']==center))
                gaxes.nodes[axispoint]['centers']= np.append(gaxes.nodes[axispoint]['centers'],center)
                gcenters.nodes[center]['axis']=axispoint
                

def update_positions(gvertices, gcenters, gaxes, step, steps, vol_alpha, rhythm):
    compute_force(gvertices, gcenters, gaxes, step, vol_alpha, rhythm)
    norm_avg = 0
    for vertex in gvertices.nodes():
        norm_avg+=np.linalg.norm(gvertices.nodes[vertex]['force'])
    # print('average of forces : ', norm_avg/len(gvertices.nodes()))
    for edge in gvertices.edges:
        if gvertices.has_edge(edge[0], edge[1]): #edge might have been removed by a previous T1 transition
            loc_dist = dist(gvertices.nodes[edge[0]]['coords'], gvertices.nodes[edge[1]]['coords'])
            ok=1
            if not(loc_dist< const.t1_min_dist and gvertices.nodes[edge[0]]['fixed']==False and gvertices.nodes[edge[1]]['fixed']==False): # T1 transition only if vertices are too close
                ok=0
            else:
                for center_ in gvertices.nodes[edge[0]]['centers']:
                    if len(gcenters.nodes[center_])<=4 or gcenters.nodes[center_]['axis'] in [0,const.n_cells_vert-1]:
                        ok=0
            if ok==1:
                t1_transition(gcenters, gvertices, edge[0], edge[1], loc_dist)
        # else:
            # print("edge has disappeared by previous T1 transition !")
    for vertex in gvertices.nodes():
        if gvertices.nodes[vertex]['fixed']==False:
            gvertices.nodes[vertex]['coords'] += const.nu * gvertices.nodes[vertex]['force']
            #rotation(gvertices, vertex, 10*const.angle*(1-scale))
        else:
            # rotation(gvertices, vertex, const.angle)
            gvertices.nodes[vertex]['coords'] += const.displacement 
    if const.tau_div !=0:   
        tmp_centers = gcenters.copy()
        for center in tmp_centers.nodes: 
            if random()<1/const.tau_div/steps and len(gcenters.nodes[center]['vertices'])>=4:
                # print('Division of center ', center)
                division(gvertices, gcenters, gaxes, center)
    for center in gcenters.nodes:
        att=gcenters.nodes[center]
        gcenters.nodes[center]['coords']=coords_center(att['vertices'], gvertices)
    if const.reassign_margin !=0:
        reassign_centers(gcenters, gaxes)
    for axis_point in gaxes.nodes:
        att = gaxes.nodes[axis_point]
        gaxes.nodes[axis_point]['coords']=coords_center(att['centers'], gcenters)
        gaxes.nodes[axis_point]['area']=reorder_around_and_area(gaxes, gcenters, axis_point)
    
    for center in gcenters.nodes:
        reorder_cell_vertices(gcenters, gvertices, center)
    normals(gvertices)
    add_area_perim(gaxes, gcenters, gvertices)

def axis_length(gaxes):
    length=0
    for idx in range(const.n_cells_vert-1):
        length+=dist(gaxes.nodes[idx]['coords'], gaxes.nodes[idx+1]['coords'])
    return length

def simulate(gvertices, gcenters, gaxes, steps, vol_alpha, rhythm):
    all_data=[[gaxes.copy(), gcenters.copy(), deepcopy(gvertices)]]
    for step in range(steps):
        print ("step ",step,' of ', steps)
        # print('total divisions : ', divisions(gcenters))
        update_positions(gvertices, gcenters, gaxes, step, steps, vol_alpha, rhythm)
        all_data.append([gaxes.copy(), deepcopy(gcenters), deepcopy(gvertices)])
        path = axis_length(gaxes)
    return all_data, path