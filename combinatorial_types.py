from json import loads, dumps
import itertools
import math
from glob import glob
from natsort import natsorted


####################################################################


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


def sublist(lst1, lst2):
    """ checks whether the first argument is a sublist of the second. """

    if isinstance(lst1, int) or isinstance(lst1, float):
        return sublist([lst1], lst2)
    if isinstance(lst2, int) or isinstance(lst2, float):
        return sublist(lst1, [lst2])
    lst3 = [value for value in lst1 if value in lst2]
    return len(lst3) == len(lst1)


####################################################################


def setdiff(lst1, lst2):
    """ computes the set difference between the first and second argument. """

    if isinstance(lst1, int) or isinstance(lst1, float):
        return setdiff([lst1], lst2)
    if isinstance(lst2, int) or isinstance(lst2, float):
        return setdiff(lst1, [lst2])
    return sorted(list(set(lst1) - set(lst2)))


##################################################################


def ind(lst, Matrix):
    """ ind(lst, Matrix) returns the list of indices i for which lst is a sublist of Matrix[i]. """

    if not Matrix:
        return []
    if isinstance(Matrix[0], int):
        Matrix = [Matrix]
    if isinstance(lst, int):
        return [index for index in range(len(Matrix)) if (
                (isinstance(Matrix[index], list) and lst in Matrix[index]) or (
                isinstance(Matrix[index], int) and lst == Matrix[index]))]
    elif len(lst) == 1:
        return [index for index in range(len(Matrix)) if (
                (isinstance(Matrix[index], list) and lst[0] in Matrix[index]) or (
                isinstance(Matrix[index], int) and lst[0] == Matrix[index]))]
    else:
        return [index for index in range(len(Matrix)) if set(lst).issubset(set(Matrix[index]))]


##################################################################


def Remove_index(lst, indices):
    """ returns the first argument after removing the indices indicated by the second argument. """

    return [lst[item] for item in setdiff(range(len(lst)), indices)]


##################################################################


def isline(lst, Matrix):
    """ isline(lst, Matrix) is TRUE if lst is a sublist of one of the lists contained 
        in Matrix (thought of as a list of lists) and FALSE otherwise. """

    if not Matrix:
        return False
    if isinstance(Matrix[0], int):
        Matrix = [Matrix]
    if isinstance(lst, int):
        return any(lst in item for item in Matrix)
    else:
        for item in Matrix:
            if set(lst).issubset(set(item)):
                return True
        return False


####################################################################


def sym_2(s):
    """ returns the elements of the symmetric group on two elements:
        sym_2(0) is the identity and sym_2(1) is the flip. """

    def sym_2_1(s):
        if s == 1:
            return 1
        elif s == 2:
            return 2


    def sym_2_2(s):
        if s == 1:
            return 2
        elif s == 2:
            return 1
    
    if s == 0:
        return sym_2_1
    elif s == 1:
        return sym_2_2


####################################################################


def lex_order(facets):
    """ 
        Input: list of facets of a simlicial polytope for some labelling of the vertices.
        Output: list of facets after relabelling so that it's lexicographically minimal.
        Example: for an input lst = [[5, 3, 9], [5, 9, 8], [3, 8, 5], [8, 9, 3]], which 
                 corresponds to the facets of a tetrahedron with labelled vertices 3, 5, 8, 9,
                 the output is [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]], which is 
                 lexicographically minimal.
    """

    vertices = unique(facets)
    edges = [list(item) for item in itertools.combinations(vertices, 2) if isline(item, facets)]
    minimal_facets = loads(dumps(facets))
    for edge in edges:
        for i in range(2):
            for j in range(2):
                new_labels, old_labels = {}, {}
                old_vertices, relabeled_facets = [], []
                for k in range(2):
                    new_labels[edge[k]] = sym_2(i)(k + 1)
                    old_labels[sym_2(i)(k + 1)] = edge[k]
                old_vertices = union(old_vertices, edge)
                start_facets, neighboring_vertices = [], []
                for facet in facets:
                    if sublist(edge, facet):
                        start_facets.append(facet)
                        neighboring_vertices = union(neighboring_vertices, facet)
                start_facets = sorted(start_facets)
                neighboring_vertices = setdiff(neighboring_vertices, edge)
                if j == 0:
                    new_labels[neighboring_vertices[0]] = 3
                    new_labels[neighboring_vertices[1]] = 4
                    old_labels[3] = neighboring_vertices[0]
                    old_labels[4] = neighboring_vertices[1]
                else:
                    new_labels[neighboring_vertices[0]] = 4
                    new_labels[neighboring_vertices[1]] = 3
                    old_labels[3] = neighboring_vertices[1]
                    old_labels[4] = neighboring_vertices[0]
                old_vertices = union(old_vertices, neighboring_vertices)
                for facet in start_facets:
                    new_facet = []
                    for k in facet:
                        new_facet = union(new_facet, new_labels[k])
                    relabeled_facets.append(new_facet)
                relabeled_facets = sorted(relabeled_facets)
                current_facets = loads(dumps(facets))
                boundary = [list(item) for item in itertools.combinations(relabeled_facets[0], 2)]
                for item in itertools.combinations(relabeled_facets[1], 2):
                    if not isline(item, boundary):
                        boundary.append(list(item))
                boundary = sorted(boundary)
                Ind = ind([1, 2], boundary)
                boundary = Remove_index(boundary, Ind)
                for facet in start_facets:
                    current_facets = Remove_index(current_facets, ind(facet, current_facets))
                counter = 4
                facets_to_remove = []
                for facet in current_facets:
                    if sublist(facet, old_vertices):
                        facet_boundary = [list(item) for item in itertools.combinations(facet, 2)]
                        facet_boundary = sorted(facet_boundary)
                        for face in facet_boundary:
                            new_face = []
                            for k in range(2):
                                new_face = union(new_face, new_labels[face[k]])
                            Ind = ind(new_face, boundary)
                            if Ind:
                                boundary = Remove_index(boundary, Ind)
                            else:
                                boundary.append(new_face)
                            boundary = sorted(boundary)
                        facets_to_remove.append(facet)
                        new_facet = []
                        for k in facet:
                            new_facet = union(new_facet, new_labels[k])
                        relabeled_facets.append(new_facet)
                        relabeled_facets = sorted(relabeled_facets)
                for facet in facets_to_remove:
                    current_facets = Remove_index(current_facets, ind(facet, current_facets))
                while current_facets:
                    pivot_vertex = []
                    for facet in current_facets:
                        old_edge = []
                        for k in range(2):
                            old_edge = union(old_edge, old_labels[boundary[0][k]])
                        if sublist(old_edge, facet):
                            pivot_vertex = setdiff(facet, old_edge)[0]
                    counter += 1
                    new_labels[pivot_vertex] = counter
                    old_labels[counter] = pivot_vertex
                    old_vertices = union(old_vertices, pivot_vertex)
                    facets_to_remove = []

                    for facet in current_facets:
                        if sublist(facet, old_vertices):
                            facet_boundary = [list(item) for item in itertools.combinations(facet, 2)]
                            facet_boundary = sorted(facet_boundary)
                            for face in facet_boundary:
                                new_face = []
                                for k in range(2):
                                    new_face = union(new_face, new_labels[face[k]])
                                Ind = ind(new_face, boundary)
                                if Ind:
                                    boundary = Remove_index(boundary, Ind)
                                else:
                                    boundary.append(new_face)
                                boundary = sorted(boundary)
                            facets_to_remove.append(facet)
                            new_facet = []
                            for k in facet:
                                new_facet = union(new_facet, new_labels[k])
                            relabeled_facets.append(new_facet)
                            relabeled_facets = sorted(relabeled_facets)
                    for facet in facets_to_remove:
                        current_facets = Remove_index(current_facets, ind(facet, current_facets))
                if relabeled_facets < minimal_facets:
                    minimal_facets = loads(dumps(relabeled_facets))
    return minimal_facets


####################################################################


def Next(num, step=1, n=math.inf):
    """ returns the next number after a given step. All computations are done modulo n.
        Examples: Next(4) = 5
                  Next(4, 3) = 7
                  Next(4, 3, 5) = 2
                  Next(4, 3, 7) = 0
    """

    return int((((num + step - 1) % n) + 1) % n)


def Previous(j, r=1, n=math.inf):
    """ returns the previous number for a given step. All computations are done modulo n.
        Examples: Previous(4) = 3
                  Previous(4, 3) = 1
                  Previous(8, 3, 3) = 2
                  Previous(8, 3, 5) = 0
    """

    return int((((j - 1 - r) % n) + 1) % n)


####################################################################


# (a) Tetrahedron
cell_tetrahedron = [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]
# (b) Triangular prism
cell_prism_3 = [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 5], [2, 4, 5], [3, 4, 5]]
# (c) Cube
cell_cube = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 6], [3, 5, 6], [4, 5, 6]]
# (d) Pentagonal prism
cell_prism_5 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 6], [3, 5, 7], [3, 6, 7], [4, 5, 7], [4, 6, 7]]
# (e) 4-4-0
cell_config_1 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 6], [3, 5, 7], [3, 6, 7], [4, 5, 8], [4, 6, 8], [5, 7, 8], [6, 7, 8]]
# (f) 3-6-0
cell_config_2 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8], [4, 5, 9], [4, 7, 9], [5, 8, 9], [6, 7, 8], [7, 8, 9]]
# (g) 4-4-1
cell_config_3 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 6], [3, 5, 7], [3, 6, 7], [4, 5, 8], [4, 6, 9], [4, 8, 9], [5, 7, 8], [6, 7, 9], [7, 8, 9]]
# (h) 3-6-1
cell_config_4 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8], [4, 5, 9], [4, 7, 9], [5, 8, 10], [5, 9, 10], [6, 7, 10], [6, 8, 10], [7, 9, 10]]
# (i) 4-4-2
cell_config_5 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 6], [3, 5, 7], [3, 6, 8], [3, 7, 8], [4, 5, 9], [4, 6, 10], [4, 9, 10], [5, 7, 9], [6, 8, 10], [7, 8, 9], [8, 9, 10]]
# (j) 2-8-0
cell_config_6 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8], [4, 5, 9], [4, 7, 9], [5, 8, 9], [6, 7, 10], [6, 8, 10], [7, 9, 10], [8, 9, 10]]
# (k) 2-8-1
cell_config_7 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8], [4, 5, 9], [4, 7, 9], [5, 8, 10], [5, 9, 10], [6, 7, 11], [6, 8, 11], [7, 9, 11], [8, 10, 11], [9, 10, 11]]
# (l) 3-6-2
cell_config_8 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8], [4, 5, 9], [4, 7, 9], [5, 8, 10], [5, 9, 10], [6, 7, 11], [6, 8, 10], [6, 10, 11], [7, 9, 11], [9, 10, 11]]
# (m) 0-12-0 dodecahedron
cell_dodecahedron = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 6], [1, 5, 6], [2, 3, 7], [2, 4, 8], [2, 7, 8], [3, 5, 9], [3, 7, 9], [4, 6, 10], [4, 8, 10], [5, 6, 11], [5, 9, 11], [6, 10, 11], [7, 8, 12], [7, 9, 12], [8, 10, 12], [9, 11, 12], [10, 11, 12]]
# (n) 1-10-2 Matzke
cell_config_9 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 9], [3, 8, 9], [4, 5, 10], [4, 7, 11], [4, 10, 11], [5, 8, 10], [6, 7, 12], [6, 9, 12], [7, 11, 12], [8, 9, 13], [8, 10, 13], [9, 12, 13], [10, 11, 13], [11, 12, 13]]
# (o) 1-10-3
cell_config_10 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8], [4, 5, 9], [4, 7, 10], [4, 9, 10], [5, 8, 11], [5, 9, 11], [6, 7, 12], [6, 8, 13], [6, 12, 13], [7, 10, 12], [8, 11, 13], [9, 10, 14], [9, 11, 14], [10, 12, 14], [11, 13, 14], [12, 13, 14]]
# (p) 2-8-4
cell_config_11 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 9], [3, 8, 9], [4, 5, 10], [4, 7, 11], [4, 10, 11], [5, 8, 12], [5, 10, 12], [6, 7, 13], [6, 9, 13], [7, 11, 13], [8, 9, 12], [9, 12, 14], [9, 13, 14], [10, 11, 14], [10, 12, 14], [11, 13, 14]]
# (q) 2-8-4
cell_config_12 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 9], [3, 8, 9], [4, 5, 10], [4, 7, 11], [4, 10, 11], [5, 8, 10], [6, 7, 12], [6, 9, 12], [7, 11, 13], [7, 12, 13], [8, 9, 14], [8, 10, 14], [9, 12, 14], [10, 11, 13], [10, 13, 14], [12, 13, 14]]
# (r) 1-10-4
cell_config_13 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 9], [3, 8, 9], [4, 5, 10], [4, 7, 11], [4, 10, 11], [5, 8, 12], [5, 10, 12], [6, 7, 13], [6, 9, 13], [7, 11, 13], [8, 9, 14], [8, 12, 14], [9, 13, 14], [10, 11, 15], [10, 12, 15], [11, 13, 15], [12, 14, 15], [13, 14, 15]]
# (s) 3-6-4
cell_config_14 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 9], [3, 8, 9], [4, 5, 10], [4, 7, 11], [4, 10, 11], [5, 8, 12], [5, 10, 12], [6, 7, 13], [6, 9, 13], [7, 11, 13], [8, 9, 12], [9, 12, 13], [10, 11, 12], [11, 12, 13]]
# (t) 3-6-4
cell_config_15 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8], [4, 5, 9], [4, 7, 10], [4, 9, 10], [5, 8, 11], [5, 9, 11], [6, 7, 12], [6, 8, 13], [6, 12, 13], [7, 10, 12], [8, 11, 13], [9, 10, 11], [10, 11, 12], [11, 12, 13]]
# (u) 3-7-4
cell_config_16 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8], [4, 5, 9], [4, 7, 10], [4, 9, 10], [5, 8, 11], [5, 9, 11], [6, 7, 12], [6, 8, 13], [6, 12, 13], [7, 10, 12], [8, 11, 14], [8, 13, 14], [9, 10, 14], [9, 11, 14], [10, 12, 14], [12, 13, 14]]
# (v) 0-12-2
cell_config_17 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 6], [1, 5, 6], [2, 3, 7], [2, 4, 8], [2, 7, 8], [3, 5, 9], [3, 7, 9], [4, 6, 10], [4, 8, 10], [5, 6, 11], [5, 9, 11], [6, 10, 12], [6, 11, 12], [7, 8, 13], [7, 9, 14], [7, 13, 14], [8, 10, 13], [9, 11, 14], [10, 12, 13], [11, 12, 14], [12, 13, 14]]
# (w) 2-8-2
cell_config_18 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8], [4, 5, 9], [4, 7, 10], [4, 9, 10], [5, 8, 11], [5, 9, 11], [6, 7, 12], [6, 8, 12], [7, 10, 12], [8, 11, 12], [9, 10, 11], [10, 11, 12]]
# (x) 2-8-3
cell_config_19 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8], [4, 5, 9], [4, 7, 10], [4, 9, 10], [5, 8, 11], [5, 9, 11], [6, 7, 12], [6, 8, 12], [7, 10, 12], [8, 11, 13], [8, 12, 13], [9, 10, 13], [9, 11, 13], [10, 12, 13]]

# Kelvin cell
cell_kelvin = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 8], [2, 7, 8], [3, 5, 9], [3, 6, 10], [3, 9, 10], [4, 5, 11], [4, 7, 12], [4, 11, 12], [5, 9, 13], [5, 11, 13], [6, 8, 10], [7, 8, 12], [8, 10, 14], [8, 12, 14], [9, 10, 13], [10, 13, 14], [11, 12, 13], [12, 13, 14]]

famous_cells = {"cell_tetrahedron": cell_tetrahedron, "cell_prism_3": cell_prism_3, "cell_cube": cell_cube,
                "cell_prism_5": cell_prism_5, "cell_config_1": cell_config_1, "cell_config_2": cell_config_2,
                "cell_config_3": cell_config_3, "cell_config_4": cell_config_4, "cell_config_5": cell_config_5,
                "cell_config_6": cell_config_6, "cell_config_7": cell_config_7, "cell_config_8": cell_config_8,
                "cell_config_9": cell_config_9, "cell_dodecahedron": cell_dodecahedron,
                "cell_config_10": cell_config_10, "cell_config_11": cell_config_11, "cell_config_12": cell_config_12,
                "cell_config_13": cell_config_13, "cell_config_14": cell_config_14, "cell_config_15": cell_config_15,
                "cell_config_16": cell_config_16, "cell_config_17": cell_config_17, "cell_config_18": cell_config_18,
                "cell_config_19": cell_config_19}


####################################################################


frequent_types = open('combinatorial_types.txt', 'w+')

# read all coarsening foam files within the working directory
all_files = []
for fpath in glob('foam*.fe'):
    if '_' in fpath:
        continue
    all_files.append(fpath)
all_files = natsorted(all_files)

for item in famous_cells:
    frequent_types.write('{}, '.format(item))
frequent_types.write('all_cells, all_combinatorial_types\n')

for file in all_files:
    file_out = open('combinatorial_types_{}.txt'.format(file[:-3]), 'w+')
    print("\r{} out of {}".format(file[:-3], len(all_files)), end="")
    with open(file, "r") as f:
        sentences = [elem for elem in f.read().split('\n') if elem]
    
    index1 = sentences.index('vertices')
    index2 = sentences.index('edges')
    index3 = sentences.index(' faces')
    index4 = sentences.index(' bodies')
    index5 = sentences.index('#include "foam.my"')

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

    # compute the dual (simplicial) polytopes of the foam cells (which are simple polytopes)
    index = -1
    FacetVertex = []
    complex_dual = []
    for body in bodies_all:
        index += 1
        FacetVertex.append([])
        complex_dual.append([])
        vertices_temp = []
        for facet in body:
            FacetVertex[index].append([])
            for edge in faces_all[facet - 1]:
                FacetVertex[index][-1] = union(FacetVertex[index][-1], edges_all[edge - 1])
                for vertex in edges_all[edge - 1]:
                    vertices_temp = union(vertices_temp, [vertex])
        for vertex in vertices_temp:
            complex_dual[index].append([])
            for i in range(len(FacetVertex[index])):
                if vertex in FacetVertex[index][i]:
                    complex_dual[index][-1].append(i + 1)

    # compute the combinatorial types within each coarsening step
    multiplicity, unique_lsts = [], []
    multiplicity_underlying_flag, unique_lsts_underlying_flag = [], []
    print("")
    for u in range(len(complex_dual)):
        print("\r", u + 1, "out of", len(complex_dual), end="")
        cell_ordered = lex_order(complex_dual[u])
        file_out.write('lex_ordered[{}] = {}\n'.format(u + 1, cell_ordered))
        if cell_ordered not in unique_lsts:
            unique_lsts.append(cell_ordered)
            multiplicity.append(1)
        else:
            multiplicity[unique_lsts.index(cell_ordered)] += 1
    file_out.close()
    print("")
    
    for item in famous_cells:
        if famous_cells[item] in unique_lsts:
            Ind = unique_lsts.index(famous_cells[item])
            frequent_types.write('{}, '.format(multiplicity[Ind]))
        else:
            frequent_types.write('{}, '.format(0))
    frequent_types.write('{}, {}\n'.format(len(complex_dual), len(unique_lsts)))
    frequent_types.flush()
frequent_types.close()
