## Script to draw an embedded graph in an orientable surface represented by a 4g-gon
## By José Frías

import math
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
from sympy import *
import itertools as it

# The fuction computes the distance between two points
def dist2pt(A,P):
    return norm(A-P)

# The fuction 'dist2segment' computes the distance between a given point P and a segment determined
# by points A and B
def dist2segment(A,B,P):
    if all(A==P) or all(B==P):
        return 0
    elif np.arccos(np.dot((P-A)/norm(P-A), (B-A)/norm(B-A))) > np.pi/2:
        return norm(P-A)
    elif np.arccos(np.dot((P-B)/norm(P-B), (A-B)/norm(A-B))) > np.pi/2:
        return norm(P-B)
    return norm(np.cross(B-A, A-P))/norm(B-A)

# The function 'belong' has parameters: B (the vertices of the polygon representing
# the underlying surface),  C (the set of vertices of the graph), and
# P is a point.
# The function returns two values:
# 'ind' encodes to which set the point P is closer, it could be closer to a segment of the polygon
# representing the surface (if ind < 4g-1 ) or to a vertex of the graph.
# The second returned value defines position of point P w.r.t. the edge of the polygon or the vertex,
# to which the point P is closer.
def belong(B,C,P):
    ind=0
    dist = dist2segment(B[0,:],B[1,:],P)
    for i in range(len(B)-1):
        if dist > dist2segment(B[i,:],B[i+1,:],P):
            ind = i
            dist = dist2segment(B[i,:],B[i+1,:],P)
    for i in range(len(C)):
        if dist > dist2pt(C[i],P):
            ind = i + len(B)-1
            dist = dist2pt(C[i],P)
    if ind > len(B)-2:
        angle = np.arccos(np.dot((P-C[ind-len(B)+1])/norm(P-C[ind-len(B)+1]), [2*np.pi/(2*len(C)),0]/norm([2*np.pi/(2*len(C)),0])))
        if (P-C[ind-len(B)+1])[1] < 0:
            return ind, 2*np.pi - angle
        else:
            return ind, angle
    else:
        return ind, dist2pt(B[ind,:],P)

# Two auxiliary functions to check when the array v belongs to A.
def find(A,v):
    for i in range(A.shape[0]):
        if np.array_equal(A[i],v):
            return i
    return -1

def find_arr(A,v):
    for i in range(len(A)):
        if np.array_equal(A[i],v):
            return i
    return -1

# An auxiliary function to count how many times appears value k in A[:,0]
def count(A,k):
    j=0
    for i in range(A.shape[0]):
        if A[i,0] == k:
            j += 1
    return j

# The following function receives as input data the genus of the surface (g) and
# the number of vertices of the embedded graph (v) and it returns  the coordinates
# of the vertices of the polygon that represents the surface and the coordinates of
# the vertices of the embedded graph.
def draw_canvas(g,v):
    # Set the aspect of the canvas to draw the graph
    plt.clf()
    figure, axes = plt.subplots()
    #plt.axis([-10., 10., -10., 10.])
    plt.setp(plt.gca(), autoscale_on=True)
    plt.axes().set_aspect('equal', 'datalim')
    pts = []
    plt.title('Embedded Graph Canvas', fontsize=16)

    # Initialaze the arrays of the vertices of the polygon and the vertices of embedded graph
    polygon = np.zeros((4*g+1,2))
    vertex = np.zeros((v,2))

    # Fill the arrays of vertices of polygon and vertices of the graph. Draw the 'fat' vertices of the graph
    for i in range(4*g+1):
        polygon[i,:] = (10*math.cos(2*i*math.pi/(4*g)),10*math.sin(2*i*math.pi/(4*g)))
    for i in range(v):
        vertex[i,:] = (4*math.cos(2*i*math.pi/v),4*math.sin(2*i*math.pi/v))
        circle= plt.Circle((vertex[i,0],vertex[i,1]), radius = 2*math.pi/(2*v))
        axes=plt.gca()
        axes.add_patch(circle)
    # Draw the polygon that defines the canvas to draw the graph
    for i in range(g):
        plt.plot(polygon[4*i:4*i+5 ,0],polygon[4*i:4*i+5,1], color=(1, 1-i/g ,i/g) , linestyle='solid')
    plt.draw()
    return polygon, vertex

# The function 'Darw_graph()' the canvas (polygon and 'fat' vertices of the graph) to plot the
# embedded graph is drawn. The user can draw the embedded graph in the canvas using the mouse or pointer.
# The output of the function are the endpoints of the arcs provided by user.
def draw_graph():
    # Initialize the array 'segments' that will contain the segments provided by user to construct the embedded graph
    # Initialize the variable 'happy' to stop the input of segments by user
    segments=[]
    happy = False

    # Input the segments of the embedded graph until the stop condition is satisfied
    while not happy:
        pts= np.asarray(plt.ginput(-1, timeout=-1, show_clicks=True, mouse_stop = 3))
        segments.append(pts)
        plt.plot(pts[:,0], pts[:,1], 'r', linestyle='solid')
        plt.draw()
        happy = plt.waitforbuttonpress()
    return segments

# The following auxiliary function uses the array of endpoints of arcs provided by user 'segments' and the arrays
# 'polygon' and 'vertex' containing the vertices of the polygon and the vertices of the graph, respectively,
# to codify each endpoint of an arc indicating the edge or the polygon or the vertex of the embedded graph to
# which it is closer, as well as the order in which it appears in the corresponding set.
def codify_arcs(segments, polygon, vertex):
    # The i-th and (i+seg)-th entries of vector 'segments1' indicate where the endpoints of the i-th arc provided by user intersect,
    # they could intersect a vertex of the graph or an edge of the polygon representing the surface.
    # The i-th and (i+seg)-th entries of vector 'orden' indicate the order of the corresponding endpoints of the i-th
    # arc around the corresponding vertex of the graph or along the edge of the polygon determining the surface.
    seg = len(segments)
    segments1 = np.zeros((2*seg,2), dtype=int)
    orden = np.zeros(2*seg)

    for i in range(seg):
        segments1[i,0],orden[i] = belong(polygon,vertex,segments[i][0])
        segments1[i+seg,0],orden[i+seg] = belong(polygon,vertex,segments[i][len(segments[i])-1])

    for i in range(v+4*g):
        segments1[:,1][segments1[:,0]==i] =np.argsort(np.argsort(orden[segments1[:,0]==i]))
    return segments1, seg

# In the function 'simplify_arcs' the arcs provided by user in the surface canvas, are merged to form the edges of the embedded graph
def simplify_arcs(segments1,seg):
    # The  arrays 'seg_aux' and 'aristas' are used to store the set of arcs provided
    # by user to 'merge' them and get the true edges of the embedded graph. The integer
    # 'edge' has the number of edges of the embedded graph.

    edge = len(segments1[segments1[:,0] > 4*g-1])
    seg_aux = segments1[segments1[:,0] > 4*g-1].tolist()
    aristas = np.zeros((edge,2),dtype=int)

    # In this for loop, all the segments provided by user in the canvas are merged to form the atual
    # edges of the embedded graph.

    for i, w in enumerate(seg_aux):
        aristas[i] = w
        j = find(segments1, w)
        if j < seg:
            j = j + seg
        else:
            j=j - seg
        while segments1[j,0] < 4*g :
            long = count(segments1, segments1[j,0])
            if segments1[j,0] % 4 < 2:
                j = find(segments1, [segments1[j,0]+2, long - 1 - segments1[j,1]] )
                if j < seg:
                    j = j + seg
                else:
                    j=j - seg
            else:
                j = find(segments1, [segments1[j,0]-2, long - 1 - segments1[j,1]] )
                if j < seg:
                    j = j + seg
                else:
                    j=j - seg
        aristas[i+int(edge/2)] =  segments1[j]
        del seg_aux[find_arr(seg_aux,segments1[j])]
    aristas[:,0] = aristas[:,0] - 4*g
    return aristas

# The main function of this module, that uses the functions described above whose input values are
# the parameters provided by user: genus of orientable surface (g) and number of verticees of
# the embeded graph (v).
# The output of function  of the function is an array called 'aristas' of size (2e,2), where 'e' is the number
# of edges of the embedded graph and the values 'aristas[i,:]' and 'aristas[i+e,:]'' represent the endpoints
# of the same arc of the embedded graph. For instance, a possible output of the function could be
# [[1 ,1],[0, 0],[0, 1],[2, 2],[2, 1],[1, 3],[1, 2],[2, 0],[0, 2],[1, 0]]
# the entrances  [0,0] and  [1,2] represent the two endpoint of the same arc, and indicate that the 0-th point
#around the 0-th vertex of the graph is connected to the second point around the vertex with label 1.
def embed_graph_canvas(g,v):
    polygon, vertex = draw_canvas(g,v)
    segments = draw_graph()
    segments1, seg = codify_arcs(segments, polygon, vertex)
    aristas = simplify_arcs(segments1,seg)
    return aristas

# Input the genus of the orientable surface and the number of vertices
print('Provide the genus of surface')
g=int(input())
print('Number of vertices in the graph')
v=int(input())

aristas = embed_graph_canvas(g,v)
print(aristas)
