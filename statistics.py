from json import loads, dumps
from scipy.spatial import ConvexHull
import itertools
import numpy as np
import math
from glob import glob
from natsort import natsorted

def volume_convex_hull(pts):
    """ computes the volume of the convex hull of points in 3-space. """

    def tetrahedron_volume(a, b, c, d):
        return np.abs(np.einsum('ij,ij->i', a - d, np.cross(b - d, c - d))) / 6
    ch = ConvexHull(pts)
    simplices = np.column_stack((np.repeat(ch.vertices[0], ch.nsimplex),
                                 ch.simplices))
    tets = ch.points[simplices]
    return np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1], tets[:, 2], tets[:, 3]))


##################################################################


def unique(lst):
    """ returns unique elements of a list. """

    out = []
    for item in lst:
        out = union(out, item)
    return out


####################################################################


def union(*lsts):
    """ returns the union of multiple lists. """

    def union_utility(lst1, lst2):
        if isinstance(lst1, int) or isinstance(lst1, float):
            return union([lst1], lst2)
        if isinstance(lst2, int) or isinstance(lst2, float):
            return union(lst1, [lst2])
        return sorted(list(set(lst1) | set(lst2)))

    lst = []
    for item in lsts:
        lst = union_utility(lst, item)
    return lst


####################################################################


def Sublist(lst1, lst2):
    """ checks whether the first argument is a sublist of the second. """
    
    lst3 = [value for value in lst1 if value in lst2]
    return len(lst3) == len(lst1)


####################################################################

def Intersection(lst1, lst2):
    """ returns the intersection of two lists. """

    return [value for value in lst1 if value in lst2]


####################################################################


def com_flagness(complex):

    ############################
    ###    test connected    ###
    ############################

    def Test_Connected(ry):
        global flag_complex
        component = loads(dumps(top_faces[0]))
        top_faces.remove(top_faces[0])
        k = 0
        flag_positive = 0
        while k < len(top_faces) and len(component) < ry:
            k = k + 1
            intersection = loads(dumps(top_faces[k - 1]))
            intersection = Intersection(intersection, component)
            if not intersection == []:
                component = list(set(component) | set(top_faces[k - 1]))
                top_faces.remove(top_faces[k - 1])
                k = k - 1
                flag_positive = 1
            if k == len(top_faces) and flag_positive == 1:
                k = 0
                flag_positive = 0
        if not len(component) == ry:
            flag_complex = 0

    ##########################
    ###    test surface    ###
    ##########################

    def Test_Surface():
        global flag_complex
        global faces
        global top_faces
        faces = [[], [], []]
        for facet in test_complex:
            if len(set(facet)) == 3:
                faces[2].append(facet)
            else:
                flag_complex = 0
            # this ensures that the complex is simplicial.
        if flag_complex == 1:
            temp = set()
            for facet in faces[2]:
                for item in itertools.combinations(facet, 2):
                    temp.add(item)
            for item in temp:
                faces[1].append(list(item))
            temp = set()
            for facet in faces[2]:
                for item in itertools.combinations(facet, 1):
                    temp.add(item)
            for item in temp:
                faces[0].append(list(item))
            top_faces = loads(dumps(faces[2]))
            Test_Connected(len(faces[0]))  # this ensures that the complex is connected.
            ### test: complex has pseudomanifold property - an edge lies in precisely two facets. ###
            for edge in faces[1]:
                count = 0
                for facet in faces[2]:
                    if set(edge) < set(facet):
                        count = count + 1
                if not count == 2:
                    flag_complex = 0
            ### test: complex has euler characteristic of sphere ###
            if flag_complex == 1:
                if not len(faces[0]) - len(faces[1]) + len(faces[2]) == 2:
                    flag_complex = 0
            if flag_complex == 1:
                ### processing links ###
                for lnk in faces[0]:
                    if flag_complex == 1:
                        link = [[], []]
                        for facet in faces[2]:
                            if lnk[0] in facet:
                                linkelement = loads(dumps(facet))
                                linkelement.remove(lnk[0])
                                link[1].append(linkelement)  # link[1] is the list of edges of the link of vertex lnk[0].
                                link[1] = sorted(link[1])
                                link[0] = list(set(link[0]) | set(
                                    linkelement))  # link[0] is the list of vertices of the link of vertex lnk[0].
                        top_faces = loads(dumps(link[1]))
                        Test_Connected(len(link[0]))  # tests that the link of each vertex is connected (a circle). If not, flag_complex is set to 0 and the loop breaks.

    ###########################
    ###    test flagness    ###
    ###########################

    def TestFlagness():
        global flag_flag
        triangles = [list(item) for item in itertools.combinations(vertices,
                                                                   3)]  # the list "vertices" is defined within the function "ComputeValence()".
        edges = []
        temp = set()
        for facet in test_complex:
            for item in itertools.combinations(facet, 2):
                temp.add(item)
        for item in temp:
            edges.append(list(item))
        edges.sort()
        for tri in triangles:
            if not tri in test_complex:
                boundary_edges = [list(item) for item in itertools.combinations(tri, 2)]
                if Sublist(boundary_edges, edges):
                    flag_flag = 0  # if the simplicial complex has an empty triangle whose edges are edges of the complex, then the dual complex is not flag*.

    #############################
    ###    compute valence    ###
    #############################

    def ComputeValence():
        global distribution_valence
        global vertices
        global valence
        global list_valence
        vertices = []  # vertices of the complex as a list of integers.
        for facet in test_complex:
            vertices = list(set(vertices) | set(facet))
        valence = []  # a list of integers so that valence[i] is the degree of vertex i.
        for k in range(max_vertex):
            valence.append(0) 
        for facet in test_complex:
            for k in facet:
                valence[k] = valence[k] + 1
        list_valence = []  # a list of integers so that list_valence[i] is the number of vertices having degree i.
        for j in range(len(vertices)):
            list_valence.append(0)
        for k in vertices:
            list_valence[valence[k]] = list_valence[valence[k]] + 1
        distribution_valence = list(item1 + item2 for (item1, item2) in
                                    itertools.zip_longest(distribution_valence, list_valence,
                                                            fillvalue=0))  # lists are added to one another component-wise, even if they have different lengths.


    ##################################################################
    ##################################################################
    ###              Global variables and functions                ###
    ##################################################################
    ##################################################################

    global list_valence
    global distribution_valence
    list_valence, distribution_valence = [], [0 for item in range(35)]

    max_vertex = 100
    count_flag_all = 0

    global flag_flag

    global flag_complex


    ###################################################################
    ###################################################################
    ###                       Main Part                             ###
    ###################################################################
    ###################################################################

    for i in range(len(complex)):
        if not bool(complex[i]):
            continue

        flag_complex = 1
        test_complex = loads(dumps(complex[i]))

        Test_Surface()
        if flag_complex != 0:
            ComputeValence()
            if len(test_complex) > 4:
                if list_valence[3] == 0:  # no vertices of degree 3, hence fundamental.
                    flag_flag = 1
                    TestFlagness()
                    if flag_flag == 1: 
                        count_flag_all += 1

    return count_flag_all * 100 / len(complex) 


################################################################################################


file_out = 'statistics.txt'
log_file = open(file_out, 'w+')
log_file.write('{}, {}, {}, {}, {}, {}\n'.format("step", "flag_perc", "avg_face_degree", "polydispersity", "avg_volume"
                                                 , "time"))

# read all coarsening foam files within the working directory
all_files = []
for fpath in glob('foam*.fe'):
    if '_' in fpath:
        continue
    all_files.append(fpath)
all_files = natsorted(all_files)

for file in all_files:
    print("\r{} out of {}".format(file[:-3], len(all_files)), end="")
    with open(file, "r") as f:
        sentences = [elem for elem in f.read().split('\n') if elem]
    NumVertices = 0
    index1 = sentences.index('vertices')
    index2 = sentences.index('edges')
    index3 = sentences.index(' faces')
    index4 = sentences.index(' bodies')
    index5 = sentences.index('#include "foam.my"')
    index6 = len(sentences) - 1

    # extract the vertices, edges, 2-faces, and 3-faces from the foam files
    vertices_all = []
    for item in sentences[index1 + 1:index2]:
        item = item + ' '
        vertices_all.append([])
        s = 0
        while s < len(item) - 1:
            if item[s].isdigit():
                for r in range(s + 1, len(item)):
                    if not item[r].isdigit() and item[r] != '.':
                        vertices_all[-1].append(float(item[s:r]))
                        s = r
                        break
            else:
                s += 1
        vertices_all[-1] = vertices_all[-1][1:]

    edges_all = []
    for item in sentences[index2 + 1:index3]:
        item = item + ' '
        edges_all.append([])
        s = 0
        while s < len(item) - 1:
            if item[s].isdigit():
                for r in range(s + 1, len(item)):
                    if not item[r].isdigit():
                        edges_all[-1].append(int(item[s:r]))
                        s = r
                        break
            else:
                s += 1
        edges_all[-1] = edges_all[-1][1:]

    faces_all = []
    for item in sentences[index3 + 1:index4]:
        item = item + ' '
        faces_all.append([])
        s = 0
        while s < len(item) - 1:
            if item[s].isdigit():
                for r in range(s + 1, len(item)):
                    if not item[r].isdigit():
                        faces_all[-1].append(int(item[s:r]))
                        s = r
                        break
            else:
                s += 1
        faces_all[-1] = faces_all[-1][1:]

    bodies_all = []
    for item in sentences[index4 + 1:index5]:
        if 'radius' in item:
            end = item.index('radius')
            item = item[:end - 1]
        elif 'old_target' in item:
            end = item.index('old_target')
            item = item[:end - 1]
        item = item + ' '
        bodies_all.append([])
        s = 0
        while s < len(item) - 1:
            if item[s].isdigit():
                for r in range(s + 1, len(item)):
                    if not item[r].isdigit():
                        bodies_all[-1].append(int(item[s:r]))
                        s = r
                        break
            else:
                s += 1
        bodies_all[-1] = bodies_all[-1][1:]

    index = -1
    FacetVertex = []
    complex = []
    for body in bodies_all:
        index += 1
        FacetVertex.append([])
        complex.append([])
        vertices_temp = []
        for facet in body:
            FacetVertex[index].append([])
            for edge in faces_all[facet - 1]:
                FacetVertex[index][-1] = union(FacetVertex[index][-1], edges_all[edge - 1])
                for vertex in edges_all[edge - 1]:
                    vertices_temp = union(vertices_temp, [vertex])
        for vertex in vertices_temp:
            complex[index].append([])
            for i in range(len(FacetVertex[index])):
                if vertex in FacetVertex[index][i]:
                    complex[index][-1].append(i + 1)

    perc_flag = com_flagness(complex)
    AvgR2, AvgR3 = 0, 0
    avg_volume = 0
    for body in bodies_all:
        body_facets = [faces_all[item - 1] for item in body]
        body_edges = [edges_all[item - 1] for item in unique(body_facets)]
        body_vertices = [vertices_all[item - 1] for item in unique(body_edges)]
        if min([v[0] for v in body_vertices]) < 1:
            body_vertices = [[v[0] + 10 * (v[0] < 5), v[1], v[2]] for v in body_vertices]
        if min([v[1] for v in body_vertices]) < 1:
            body_vertices = [[v[0], v[1] + 10 * (v[1] < 5), v[2]] for v in body_vertices]
        if min([v[2] for v in body_vertices]) < 1:
            body_vertices = [[v[0], v[1], v[2] + 10 * (v[2] < 5)] for v in body_vertices]

        volume = volume_convex_hull(body_vertices)
        AvgR2 += (3 * abs(volume)/(4 * math.pi)) ** (2 / 3)
        AvgR3 += 3 * volume/(4 * math.pi)
        avg_volume += volume
    AvgR2 /= len(bodies_all)
    AvgR3 /= len(bodies_all)
    avg_volume /= len(bodies_all)
    Polydispersity = (AvgR3 ** (2/3)) / AvgR2 - 1
    avg = sum([len(body) for body in bodies_all]) / len(bodies_all)
    sentence = sentences[index6]
    Ind = sentence.index("=")
    t = float(sentence[Ind + 1:])
    log_file.write('{}, {}, {}, {}, {}, {}\n'.format(file[:-3], perc_flag, avg, Polydispersity, avg_volume, t))
    log_file.flush()
log_file.close()
