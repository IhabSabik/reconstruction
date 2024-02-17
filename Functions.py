import itertools
from copy import copy, deepcopy
import math
import sys
import numpy as np
from time import time
import networkx as nx
from itertools import combinations
from johnson import simple_cycles
from json import loads, dumps
import random
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from collections import defaultdict


def intersect(*lsts):
    """ returns the intersection of multiple lists. """

    def intersect_utility(lst1, lst2):
        if isinstance(lst1, int) or isinstance(lst1, float):
            lst1 = [lst1]
        if isinstance(lst2, int) or isinstance(lst2, float):
            lst2 = [lst2]
        return [value for value in lst1 if value in lst2]

    lst = lsts[0]
    for item in lsts[1:]:
        lst = intersect_utility(lst, item)
    return lst


####################################################################


def sublist(lst_1, lst_2):
    """ checks whether lst_1 is a sublist of lst_2. """

    if isinstance(lst_1, int) or isinstance(lst_1, float):
        return sublist([lst_1], lst_2)
    if isinstance(lst_2, int) or isinstance(lst_2, float):
        return sublist(lst_1, [lst_2])
    lst_3 = [value for value in lst_1 if value in lst_2]
    return len(lst_3) == len(lst_1)


####################################################################


def setdiff(lst_1, lst_2):
    """ returns the set difference lst_1 - lst_2. """

    if isinstance(lst_1, int) or isinstance(lst_1, float):
        return setdiff([lst_1], lst2)
    if isinstance(lst_2, int) or isinstance(lst_2, float):
        return setdiff(lst_1, [lst_2])
    return sorted(list(set(lst_1) - set(lst_2)))


####################################################################


def union(*lsts):
    """ returns the union of multiple lists. """

    def union_utility(lst_1, lst_2):
        if isinstance(lst_1, int) or isinstance(lst_1, float):
            return union([lst_1], lst_2)
        if isinstance(lst_2, int) or isinstance(lst_2, float):
            return union(lst_1, [lst_2])
        return sorted(list(set(lst_1) | set(lst_2)))

    lst = []
    for item in lsts:
        lst = union_utility(lst, item)
    return lst


####################################################################


def is_tetrahedron(lst, Matrix):
    """ checks whether vertices of lst form a clique of length 4. """

    return len(lst) == 4 and is_clique(lst, Matrix)


####################################################################


def is_clique(lst, Matrix):
    """ checks whether vertices of lst form a clique. """

    for i in lst:
        for j in [item for item in lst if item > i]:
            if i not in Matrix[j]:
                return False
    return True


####################################################################


def components(i, list_temp, Matrix):
    """ returns the components of neighbors of vertex i after removing vertices
        with indices in list_temp. """

    temp = set(Matrix[i]) - set(list_temp)
    temp = temp - {i}
    r = len(temp)
    Components = []
    m = -1
    list_all = []

    while len(list_all) != r:
        m += 1
        temp = set(temp) - set(list_all)
        Components.append(list(temp)[0])
        n = 0
        Components[m] = [Components[m]]
        while n != len(Components[m]):
            n = len(Components[m])
            for j in temp - set(Components[m]):
                for k in Components[m]:
                    if k in Matrix[j]:
                        Components[m] = union(Components[m], [j])
                        break
        for element in Components:
            list_all = set(list_all) | set(element)
    return Components


########################################################


def components_list(lst, Matrix):
    """ returns the connected components of cells with indices in lst. """

    temp = loads(dumps(lst))
    r = len(temp)
    Components = []
    m = -1
    list_all = []
    while len(list_all) != r:
        m += 1
        temp = sorted(setdiff(temp, list_all))
        Components.append(temp[0])
        n = 0
        Components[m] = [Components[m]]
        while n != len(Components[m]):
            n = len(Components[m])
            for j in setdiff(temp, Components[m]):
                for k in Components[m]:
                    if k in Matrix[j]:
                        Components[m] = union(Components[m], [j])
                        break
        for element in Components:
            list_all = union(list_all, element)
    return Components


###############################################################


def components_edges(edges):
    """ returns the connected components of edges. """

    temp = unique(edges)
    r = len(temp)
    Components = []
    m = -1
    list_all = []
    while len(list_all) != r:
        m += 1
        temp = sorted(setdiff(temp, list_all))
        Components.append(temp[0])
        n = 0
        Components[m] = [Components[m]]
        while n != len(Components[m]):
            n = len(Components[m])
            for j in setdiff(temp, Components[m]):
                for k in Components[m]:
                    if isline([j, k], edges):
                        Components[m] = union(Components[m], [j])
                        break
        for element in Components:
            list_all = union(list_all, element)
    return Components


###############################################################


def save_list_of_lists(variable, file_name, default_empty="-1"):
    """ saves variable, which is a list of lists, to a text file named file_name.
        If one of the lists in variable is empty, it gets saved as default_empty. """

    with open(file_name, 'w') as f:
        for item in variable:
            if not item:
                f.writelines(default_empty)
            else:
                for i in range(len(item) - 1):
                    f.writelines("%s, " % item[i])
                f.writelines("%s" % item[-1])
            f.writelines("\n")
    return


###############################################################


def save_list(variable, file_name):
    """ saves the list variable to a text file named file_name. """

    with open(file_name, 'w') as f:
        for item in variable:
            f.writelines("%s\n" % item)
    return


###############################################################


def read_list_of_lists_of_lists(file, is_float=False):
    """ reads nested 3-layered lists such that the final layer consists of homogeneous lists.
        Use to read 'candidate_edges' as well as 'changed_tunnels'. """

    lst_out = []
    with open(file, 'r') as f:
        for line in f:
            if line[0] != "[":
                lst_out.append([])
            else:
                count = 1
                for elem in line:
                    if elem == ",":
                        count += 1
                    elif elem == "]":
                        break
                lst = []
                first_time = True
                n = 0
                for elem in line:
                    if elem.isdigit() or elem in ['-', '.']:
                        if first_time:
                            n = elem
                            first_time = False
                        else:
                            n = n + elem
                    elif n != 0:
                        if is_float:
                            lst.append(float(n))
                        else:
                            lst.append(int(n))
                        first_time = True
                        n = 0
                edges = []
                while lst:
                    edges.append(lst[:count])
                    lst = lst[count:]
                lst_out.append(edges)
    return lst_out


################################################################


def dual_list_in_triangulation(lst, triangulation, file_out, bad_vertices=None):
    """ computes the dual of each of the vertices in lst within the triangulation
        and changes the indices from {x_1, ..., x_n} to {1, ..., n}.
        Duals of vertices in bad_vertices get assigned empty sets. """

    if bad_vertices is None:
        bad_vertices = []
    dual_out = open(file_out, "w")
    complex_dual = []
    for vertex in union(lst, bad_vertices):
        if vertex in bad_vertices:
            dual_out.writelines("[]\n")
            complex_dual.append([])
            continue
        dual = [setdiff(triangulation[item], vertex) for item in ind(vertex, triangulation)]
        lst_new = unique(dual)
        dual = [[lst_new.index(item_sub) + 1 for item_sub in item] for item in dual]
        dual.sort()
        dual_out.writelines("[%s]\n" % ', '.join(str(e) for e in dual))
        complex_dual.append(dual)
    dual_out.close()
    return complex_dual


################################################################


def find_K4s(Matrix, lst=None):
    """ returns all 4-cliques consisting of vertices in lst with respect to the adjacency list Matrix.
        If lst is not given, the function returns all 4-cliques in the adjacency list. """

    if lst is None:
        lst = range(len(Matrix))
    K4s = []
    for i in lst:
        for j in [x for x in Matrix[i] if x > i]:
            list_1 = [x for x in Matrix[j] if x in Matrix[i] and x > j]
            for k in list_1:
                list_2 = [x for x in Matrix[k] if x in list_1 and x > k]
                for l in list_2:
                    K4s.append([i, j, k, l])
    return K4s


################################################################


def find_K5s(Matrix, lst=None):
    """ returns all 5-cliques consisting of vertices in lst with respect to the adjacency list Matrix.
        If lst is not given, the function returns all 5-cliques in the adjacency list. """

    if lst is None:
        lst = range(len(Matrix))
    K5s = []
    for i in lst:
        for j in [x for x in Matrix[i] if x > i]:
            list_1 = [x for x in Matrix[j] if x in Matrix[i] and x > j]
            for k in list_1:
                list_2 = [x for x in Matrix[k] if x in list_1 and x > k]
                for l in list_2:
                    list_3 = [x for x in Matrix[l] if x in list_2 and x > l]
                    for m in list_3:
                        K5s.append([i, j, k, l, m])
    return K5s


##################################################################


def Remove(i, j, Matrix):
    """ removes the connection between vertices i and j from the adjacency list Matrix. """

    N = Matrix[:]
    N[i] = sorted(setdiff(N[i], [j]))
    N[j] = sorted(setdiff(N[j], [i]))
    return N


##################################################################


def Add(i, j, Matrix):
    """ adds a connection between vertices i and j to the adjacency list Matrix. """

    N = Matrix[:]
    N[i] = sorted(union(N[i], [j]))
    N[j] = sorted(union(N[j], [i]))
    return N


##################################################################


def distance_vox(cell_1, cell_2, Centroid, coordinates, voxels, dist_path, sampling=4):
    """ approximates the voxelized distance between two sets of coordinates (cell_1 and cell_2) representing point
        clouds. The result is saved in the file specified by dist_path. The approximation involves randomly selecting a
        number of starting points (equal to the sampling parameter) from the first cloud and computing a zig-zag path to
        the nearest points in the other cloud. Additionally, the centroids of both clouds are included among the
        starting points. """

    if cell_1 == cell_2:
        return 0
    if sampling < 0:
        sampling = 0
    cell_1, cell_2 = min(cell_1, cell_2), max(cell_1, cell_2)
    with open(dist_path) as f:
        computed_distances = [[str_2_int_float(num) for num in line.split(',')] for line in f]
    Ind = ind([cell_1, cell_2], computed_distances)
    Ind_new = []
    for item in Ind:
        if isline([cell_1, cell_2], computed_distances[item][:2]):
            Ind_new.append(item)
            break
    if Ind_new:
        return computed_distances[Ind_new[0]][2]
    for cell_ind in [cell_1, cell_2]:
        if len(coordinates[cell_ind]) == 0:
            coordinates[cell_ind] = np.argwhere(voxels == cell_ind + 1) + 1
    cell, cell_op = coordinates[cell_1], coordinates[cell_2]
    center = np.array(Centroid[cell_2])
    center = center[[2, 1, 0]]
    diff = np.linalg.norm(cell - center, axis=1)
    vertex_1 = cell[np.argmin(diff)]

    center = np.array(Centroid[cell_1])
    center = center[[2, 1, 0]]
    diff = np.linalg.norm(cell_op - center, axis=1)
    vertex_2_op = cell_op[np.argmin(diff)]
    diff = np.linalg.norm(cell - vertex_2_op, axis=1)
    vertex_2 = cell[np.argmin(diff)]
    sample = np.array([vertex_1, vertex_2])

    if sampling > 0:
        pts = random.sample(range(cell.shape[0]), min(cell.shape[0], sampling))
        sample = np.append(sample, [cell[item] for item in pts], axis=0)

    dist_out = math.inf
    vertex_out, vertex_op_out = [], []
    for vertex in sample:
        for count in range(3):
            diff = np.linalg.norm(cell_op - vertex, axis=1)
            Ind = np.argmin(diff)
            vertex_op = cell_op[Ind]

            diff = np.linalg.norm(cell - vertex_op, axis=1)
            Ind = np.argmin(diff)
            dist = diff[Ind]
            vertex = cell[Ind]

            if dist < dist_out:
                dist_out = dist
                vertex_out = vertex
                vertex_op_out = vertex_op

    computed_distances.append([cell_1, cell_2, round(dist_out, 2)] + list(vertex_out) + list(vertex_op_out))
    save_list_of_lists(computed_distances, dist_path)
    return dist_out


##################################################################


def neighbors_within_set(vertex, lst, Matrix):
    """ returns the subset of lst adjacent to the vertex. """

    return intersect(lst, Matrix[vertex])


##################################################################


def neighbors_in_triangulation(lst, triangulation):
    """ returns mutual neighbors of vertices of lst within the triangulation. """

    if isinstance(lst, int):
        lst = [lst]
    tets = [triangulation[item] for item in ind(lst, triangulation)]
    lst_out = setdiff(unique(tets), lst)
    lst_out = [int(item) for item in lst_out]
    return lst_out


##################################################################


def unique(lst):
    """ returns the unique elements of a list. """

    out = []
    for item in lst:
        out = union(out, item)
    return out


##################################################################


def is_edge_in_cycle(edge, cycle):
    """ checks whether the edge lies in the cycle. """

    if not isline(edge, cycle):
        return False
    return cycle[Next(cycle.index(edge[0]), 1, len(cycle))] == edge[1] or cycle[
        Next(cycle.index(edge[1]), 1, len(cycle))] == edge[0]


##################################################################


def istriangulatable(lst, Matrix, list_boundary, Centroid, changed_tunnels=None, try_tunnels=None,
                     with_repair=True, predetermined=None, predetermined_empty=None, post_process=False,
                     post_process_further=False):
    if changed_tunnels is None:
        changed_tunnels = [[], []]
    if try_tunnels is None:
        try_tunnels = []
    if predetermined is None:
        predetermined = []
    if predetermined_empty is None:
        predetermined_empty = []

    def comb(x, A):
        if isinstance(A, bool):
            return [x] + [A]
        return [x] + A

    def comb_lists(x, A):
        if isinstance(A[0], bool):
            return [x] + [A]
        return [x] + A

    if isinstance(lst, list) and len(lst) > 1:
        R_1, triangulation_sub_1, problems_1 = istriangulatable(lst[0], Matrix, list_boundary, Centroid,
                                                                changed_tunnels, try_tunnels, with_repair,
                                                                predetermined, predetermined_empty)
        R_2, triangulation_sub_2, problems_2 = istriangulatable(lst[1:], Matrix, list_boundary, Centroid,
                                                                changed_tunnels, try_tunnels, with_repair,
                                                                predetermined, predetermined_empty)
        if not try_tunnels:
            R = comb(R_1, R_2)
            triangulation_sub = [item for item in triangulation_sub_1]
            for tet in triangulation_sub_2:
                if not isline(tet, triangulation_sub_1):
                    triangulation_sub.append(tet)
            return R, triangulation_sub, union(problems_1, problems_2)
        else:
            R = comb_lists(R_1, R_2)
            return R, predetermined, [union(problems_1[item], problems_2[item]) for item in range(len(problems_1))]
    elif isinstance(lst, list) and len(lst) == 1:
        return istriangulatable(lst[0], Matrix, list_boundary, Centroid, changed_tunnels, try_tunnels,
                                with_repair, predetermined, predetermined_empty)
    elif isinstance(lst, int):
        vertex = lst
        if vertex in list_boundary:
            if not try_tunnels:
                return False, predetermined, [math.inf]
            else:
                n = len(list(combinations(try_tunnels, 2)))
                return [False for item in range(n)], predetermined, [[math.inf] for item in range(n)]
        Matrix_temp = loads(dumps(Matrix))
        link = vertexlink(vertex, Matrix)
        K4_temp = []
        K5_temp = []
        vts = Matrix[vertex]
        final = []
        for i in vts:
            Nei_i = intersect(Matrix[i], vts)
            for j in [item for item in Nei_i if item > i]:
                Nei_i_j = intersect(Nei_i, Matrix[j])
                for k in [item for item in Nei_i_j if item > j]:
                    Nei_i_j_k = intersect(Nei_i_j, Matrix[k])
                    for l in [item for item in Nei_i_j_k if item > k]:
                        if not isline([vertex, i, j, k, l], try_tunnels):
                            K4_temp.append([i, j, k, l])
                        else:
                            final.append([i, j, k, l])
                        Nei_i_j_k_l = intersect(Nei_i_j_k, Matrix[l])
                        for m in [item for item in Nei_i_j_k_l if item > l]:
                            K5_temp.append([i, j, k, l, m])
        K4_temp = K4_temp + final
        # print("final:", final)

        # detect constriction and extra edges
        reduced_link = reducedvertexlink(vertex, Matrix, list_boundary, Centroid)
        reduced_edges = []
        for i in vts:
            Nei_i = Neighbors_link(i, reduced_link)
            for j in [item for item in Nei_i if item > i]:
                reduced_edges.append([i, j])
        open_reduced_edges = []
        # print("reduced_edges:", reduced_edges)
        for edge in reduced_edges:
            if len(ind(edge, reduced_link)) != 2:
                open_reduced_edges.append(edge)
        constricted = [tri for tri in cycles_deg(open_reduced_edges, 3) if not isline(tri, K4_temp)]
        # print("constricted:", constricted)
        Ind = ind(vertex, predetermined)
        triangulation_sub = [predetermined[index] for index in Ind]
        predetermined = Remove_index(predetermined, Ind)
        triangulation_predetermined = loads(dumps(triangulation_sub))
        constricted = [item for item in constricted if not isline(item + [vertex], triangulation_predetermined)]
        # print("constricted:", constricted)
        initial_len = len(triangulation_predetermined)
        remove = []
        Holes = holes(vertex, Matrix, list_boundary, Centroid, True)[0]
        boundary_holes = []
        for hole in Holes:
            if sublist(hole, list_boundary):
                boundary_holes.append(hole)
        # print("boundary_holes:", boundary_holes)
        reduced_vertices = unique(reduced_link)
        K4_reduced = []
        K5_reduced = []
        for i in reduced_vertices:
            Nei_i = intersect(Matrix[i], reduced_vertices)
            for j in [item for item in Nei_i if item > i]:
                Nei_i_j = intersect(Nei_i, Matrix[j])
                for k in [item for item in Nei_i_j if item > j]:
                    Nei_i_j_k = intersect(Nei_i_j, Matrix[k])
                    for l in [item for item in Nei_i_j_k if item > k]:
                        K4_reduced.append([i, j, k, l])
                        Nei_i_j_k_l = intersect(Nei_i_j_k, Matrix[l])
                        for m in [item for item in Nei_i_j_k_l if item > l]:
                            K5_reduced.append([i, j, k, l, m])
        for i in range(len(link)):
            if isline([vertex] + link[i], triangulation_predetermined + predetermined_empty):
                continue
            if isline(link[i], constricted):
                remove.append(i)
                continue
            if isline(link[i], K4_temp):
                continue
            extra_triangle = False
            count_deg_1_edges = 0  # the number of edges of degree 1 must be 1.
            for edge in combinations(link[i], 2):
                if len(ind(edge, link)) < 2 and not isline(edge, Holes):
                    count_deg_1_edges += 1
                    # print("0", link[i])
            if count_deg_1_edges > 1:
                constricted.append(link[i])
                remove.append(i)
                continue
            elif count_deg_1_edges == 0:
                triangulation_sub.append([vertex] + link[i])
                continue
            for v in link[i]:
                triangles = [link[item] for item in ind(v, link)]
                if sum([isline(tri, K4_temp) for tri in triangles]) == len(triangles):
                    # if all triangles surrounding 'v' lie in 4-cliques then the next argument does not work.
                    continue
                if v in reduced_vertices:
                    # print(v)
                    num_edges = len(Neighbors_link(v, reduced_link))
                    Ind = ind(v, K4_reduced)
                    # count = 0
                    # for j in combinations(Ind, 2):
                    #     count += len(intersect(K4_reduced[j[0]], K4_reduced[j[1]])) == 3
                    num_triangles = len(ind(v, reduced_link)) - len(Ind)  # - count
                    if num_triangles > num_edges:
                        extra_triangle = True
                        # print("1", link[i], v, num_triangles, num_edges)
                        break
                if isline(v, K5_temp):
                    continue
                num_edges = len(mutual_neighbors([vertex, v], Matrix))
                Ind = ind(v, K4_temp)
                # count = 0
                # for j in combinations(Ind, 2):
                #     count += len(intersect(K4_temp[j[0]], K4_temp[j[1]])) == 3
                num_triangles = len(ind(v, link)) - len(Ind)  # - count
                # if v == 169:
                #     print(link[i], num_triangles, num_edges)
                if num_triangles > num_edges:
                    extra_triangle = True
                    # print("2", link[i], v, num_triangles, num_edges)
                    break
            if not extra_triangle:
                triangulation_sub.append([vertex] + link[i])
                continue
            constricted.append(link[i])
            remove.append(i)

        # print("constricted:", constricted)
        # print("vts:", vts)
        link = Remove_index(link, remove)
        constricted_vts = []
        if not Holes:
            for i in vts:
                Ind = ind(i, link)
                if not Ind:
                    constricted_vts.append(i)
                    Matrix_temp = Remove(vertex, i, Matrix_temp)
        vts = setdiff(vts, constricted_vts)
        # print("vts:", vts)
        repairable = [[] for item in range(4)]
        Rest = []
        if try_tunnels:
            Rest = K4_temp[-1]
            K4_temp = K4_temp[:-1]
        for i in range(len(K4_temp)):
            repairable[0].append(K4_temp[i])
            repairable[1].append([])
            repairable[2].append([])
            repairable[3].append([])
            K5_current = [vertex] + K4_temp[i]
            # tunnel = tunnel_K5(K5_current, Matrix_temp, Centroid, list_boundary)
            tunnel = tunnel_K5_geometric(K5_current, Centroid)
            Ind = ind(K5_current, changed_tunnels[0])
            if Ind:
                tunnel = changed_tunnels[1][Ind[0]]
            set1 = []
            set2 = []
            if tunnel:
                if vertex in tunnel:
                    apex = setdiff(tunnel, [vertex])[0]
                    base = setdiff(K4_temp[i], apex)
                    set1.append(base)
                    set2.append([apex, base[0], base[1]])
                    set2.append([apex, base[0], base[2]])
                    set2.append([apex, base[1], base[2]])
                    repairable[1][-1] = set1
                    repairable[2][-1] = set2
                else:
                    main_diagonal = tunnel
                    rest = setdiff(K4_temp[i], main_diagonal)
                    repairable[1][-1] = [main_diagonal + [rest[0]], main_diagonal + [rest[1]]]
                    repairable[2][-1] = [rest + [main_diagonal[0]], rest + [main_diagonal[1]]]
                    repairable[3][-1] = [main_diagonal, rest]
            else:
                T, apex = ispyramid(K4_temp[i], vertex, Centroid)
                if T:
                    base = setdiff(K4_temp[i], apex)
                    set1.append(base)
                    set2.append([apex, base[0], base[1]])
                    set2.append([apex, base[0], base[2]])
                    set2.append([apex, base[1], base[2]])
                    repairable[1][-1] = set1
                    repairable[2][-1] = set2
                    continue
                T, diagonals = iscross(K4_temp[i], vertex, Centroid)
                if T:
                    main_diagonal = diagonals[0]
                    rest = setdiff(K4_temp[i], main_diagonal)
                    repairable[1][-1] = [main_diagonal + [rest[0]], main_diagonal + [rest[1]]]
                    repairable[2][-1] = [rest + [main_diagonal[0]], rest + [main_diagonal[1]]]
                    repairable[3][-1] = diagonals
                    continue
                T = is_4_simplex(K4_temp[i], vertex, Centroid)
                if not T:
                    continue
                for tri in combinations(K4_temp[i], 3):
                    set2.append(list(tri))
                repairable[1][-1] = set1
                repairable[2][-1] = set2
        # print("repairable:", repairable)
        if Rest:
            repairable[0].append(Rest)
            repairable[1].append([])
            repairable[2].append([])
            repairable[3].append([])
            K5_current = sorted([vertex] + Rest)
            for j in combinations(K5_current, 2):
                tunnel = list(j)
                # print(tunnel)
                set1 = []
                set2 = []
                if vertex in tunnel:
                    apex = setdiff(tunnel, [vertex])[0]
                    base = setdiff(Rest, apex)
                    set1.append(base)
                    set2.append([apex, base[0], base[1]])
                    set2.append([apex, base[0], base[2]])
                    set2.append([apex, base[1], base[2]])
                    repairable[1][-1].append(set1)
                    repairable[2][-1].append(set2)
                    repairable[3][-1].append([])
                else:
                    main_diagonal = tunnel
                    rest = setdiff(Rest, main_diagonal)
                    repairable[1][-1].append([main_diagonal + [rest[0]], main_diagonal + [rest[1]]])
                    repairable[2][-1].append([rest + [main_diagonal[0]], rest + [main_diagonal[1]]])
                    repairable[3][-1].append([main_diagonal, rest])
            # print(repairable[0])
            # print(repairable[1])
            # print(repairable[2])
            # print(repairable[3])
        affected_facets = [[], []]
        for i in range(len(repairable[0])):
            if len(repairable[1][i]) == 2:
                good = repairable[1][i]
                bad = repairable[2][i]
                T = [True for item in good]
                for j in range(len(good)):
                    if isline([vertex] + good[j], predetermined_empty):
                        T[j] = False
                        continue
                    if tunnel_K5([vertex] + repairable[0][i], Matrix_temp, Centroid, list_boundary) \
                            != tunnel_K5_geometric([vertex] + repairable[0][i], Centroid):
                        continue  # if the tunnel has been changed, we cannot use geometry to decide empty tetrahedra.
                    if vertex not in list_boundary and len(components(vertex, good[j],
                                                                      Matrix_temp)) > 1:  # if a triangle splits the link into two components, then it must be empty.
                        T[j] = False
                    mid = [0, 0, 0]
                    for item in good[j]:
                        mid = add(mid, Centroid[item])
                    mid = multiply_scalar(1 / 3, mid)
                    for k in range(len(link)):
                        if isline(link[k], repairable[0][i]):  # or len(intersect(link[k], good[j])) == 2:
                            continue
                        a, b = barycentric_intersection_in_line(link[k], Centroid[vertex], mid, Centroid)
                        if (a > 0 and b > 0) or (a < 0 and b > 1):
                            I = intersect_triangle_line(link[k], Centroid[vertex], mid, Centroid)
                            a, b, c = barycentric_in_triangle(link[k], I, Centroid)
                            if a > 0 and b > 0 and c > 0:
                                if isline([vertex] + link[k], triangulation_sub):
                                    T[j] = False
                                    break
                affected_facets[0].append(good[0])
                affected_facets[1].append(T[0])
                affected_facets[0].append(good[1])
                affected_facets[1].append(T[1])
                affected_facets[0].append(bad[0])
                affected_facets[1].append(False)
                affected_facets[0].append(bad[1])
                affected_facets[1].append(False)
            elif len(repairable[1][i]) == 1:
                good = repairable[2][i]
                bad = repairable[1][i]
                T = [True for item in good]
                for j in range(len(good)):
                    if isline([vertex] + good[j], predetermined_empty):
                        T[j] = False
                        continue
                    if tunnel_K5([vertex] + repairable[0][i], Matrix_temp, Centroid, list_boundary) != \
                            tunnel_K5_geometric([vertex] + repairable[0][i], Centroid):
                        continue  # if the tunnel has been changed, we cannot use geometry to decide empty tetrahedra.
                    if vertex not in list_boundary and len(components(vertex, good[j],
                                                                      Matrix_temp)) > 1:  # if a triangle splits the link into two components, then it must be empty.
                        T[j] = False
                    mid = [0, 0, 0]
                    for item in good[j]:
                        mid = add(mid, Centroid[item])
                    mid = multiply_scalar(1 / 3, mid)
                    for k in range(len(link)):
                        if isline(link[k], repairable[0][i]):  # or len(intersect(link[k], good[j])) == 2:
                            continue
                        a, b = barycentric_intersection_in_line(link[k], Centroid[vertex], mid, Centroid)
                        if (a > 0 and b > 0) or (a < 0 and b > 1):
                            I = intersect_triangle_line(link[k], Centroid[vertex], mid, Centroid)
                            a, b, c = barycentric_in_triangle(link[k], I, Centroid)
                            if a > 0 and b > 0 and c > 0:
                                if isline([vertex] + link[k], triangulation_sub):
                                    T[j] = False
                                    break
                affected_facets[0].append(good[0])
                affected_facets[1].append(T[0])
                affected_facets[0].append(good[1])
                affected_facets[1].append(T[1])
                affected_facets[0].append(good[2])
                affected_facets[1].append(T[2])
                affected_facets[0].append(bad[0])
                affected_facets[1].append(False)
            elif len(repairable[1][i]) == 0:
                good = repairable[2][i]
                T = [True for item in good]
                for j in range(len(good)):
                    if isline([vertex] + good[j], predetermined_empty):
                        T[j] = False
                        continue
                    if tunnel_K5([vertex] + repairable[0][i], Matrix_temp, Centroid, list_boundary) != \
                            tunnel_K5_geometric([vertex] + repairable[0][i], Centroid):
                        continue  # if the tunnel has been changed, we cannot use geometry to decide empty tetrahedra.
                    if vertex not in list_boundary and len(components(vertex, good[j],
                                                                      Matrix_temp)) > 1:  # if a triangle splits the link into two components, then it must be empty.
                        T[j] = False
                        # print("prob:", good[j])
                        continue
                    mid = [0, 0, 0]
                    for item in good[j]:
                        mid = add(mid, Centroid[item])
                    mid = multiply_scalar(1 / 3, mid)
                    for k in range(len(link)):
                        if isline(link[k], repairable[0][i]):  # or len(intersect(link[k], good[j])) == 2:
                            continue
                        a, b = barycentric_intersection_in_line(link[k], Centroid[vertex], mid, Centroid)
                        if (a > 0 and b > 0) or (a < 0 and b > 1):
                            I = intersect_triangle_line(link[k], Centroid[vertex], mid, Centroid)
                            a, b, c = barycentric_in_triangle(link[k], I, Centroid)
                            if a > 0 and b > 0 and c > 0:
                                if isline([vertex] + link[k], triangulation_sub):
                                    T[j] = False
                                    break
                # print("T:", T)
                affected_facets[0].append(good[0])
                affected_facets[1].append(T[0])
                affected_facets[0].append(good[1])
                affected_facets[1].append(T[1])
                affected_facets[0].append(good[2])
                affected_facets[1].append(T[2])
                affected_facets[0].append(good[3])
                affected_facets[1].append(T[3])
            else:  # here we try to change the tunnel of the given 5-clique 'try_tunnels'
                # print("affected_facets:", affected_facets)
                affected_facets_versions = [loads(dumps(affected_facets)) for item in repairable[1][i]]
                for index in range(len(affected_facets_versions)):
                    if len(repairable[1][i][index]) == 2:
                        good = repairable[1][i][index]
                        bad = repairable[2][i][index]
                        T = [True for item in good]
                        for j in range(len(good)):
                            if isline([vertex] + good[j], predetermined_empty):
                                T[j] = False
                                continue
                            if tunnel_K5([vertex] + repairable[0][i], Matrix_temp, Centroid, list_boundary) != \
                                    tunnel_K5_geometric([vertex] + repairable[0][i], Centroid):
                                continue  # if the tunnel has been changed, we cannot use geometry to decide empty tetrahedra.
                            if vertex not in list_boundary and len(components(vertex, good[j],
                                                                              Matrix_temp)) > 1:  # if a triangle splits the link into two components, then it must be empty.
                                T[j] = False
                            mid = [0, 0, 0]
                            for item in good[j]:
                                mid = add(mid, Centroid[item])
                            mid = multiply_scalar(1 / 3, mid)
                            for k in range(len(link)):
                                if isline(link[k], repairable[0][i]):  # or len(intersect(link[k], good[j])) == 2:
                                    continue
                                a, b = barycentric_intersection_in_line(link[k], Centroid[vertex], mid, Centroid)
                                if (a > 0 and b > 0) or (a < 0 and b > 1):
                                    I = intersect_triangle_line(link[k], Centroid[vertex], mid, Centroid)
                                    a, b, c = barycentric_in_triangle(link[k], I, Centroid)
                                    if a > 0 and b > 0 and c > 0:
                                        if isline([vertex] + link[k], triangulation_sub):
                                            T[j] = False
                                            break
                        affected_facets_versions[index][0].append(good[0])
                        affected_facets_versions[index][1].append(T[0])
                        affected_facets_versions[index][0].append(good[1])
                        affected_facets_versions[index][1].append(T[1])
                        affected_facets_versions[index][0].append(bad[0])
                        affected_facets_versions[index][1].append(False)
                        affected_facets_versions[index][0].append(bad[1])
                        affected_facets_versions[index][1].append(False)
                    elif len(repairable[1][i][index]) == 1:
                        good = repairable[2][i][index]
                        bad = repairable[1][i][index]
                        T = [True for item in good]
                        for j in range(len(good)):
                            if isline([vertex] + good[j], predetermined_empty):
                                T[j] = False
                                continue
                            if tunnel_K5([vertex] + repairable[0][i], Matrix_temp, Centroid, list_boundary) != \
                                    tunnel_K5_geometric([vertex] + repairable[0][i], Centroid):
                                continue  # if the tunnel has been changed, we cannot use geometry to decide empty tetrahedra.
                            if vertex not in list_boundary and len(components(vertex, good[j],
                                                                              Matrix_temp)) > 1:  # if a triangle splits the link into two components, then it must be empty.
                                T[j] = False
                            mid = [0, 0, 0]
                            for item in good[j]:
                                mid = add(mid, Centroid[item])
                            mid = multiply_scalar(1 / 3, mid)
                            for k in range(len(link)):
                                if isline(link[k], repairable[0][i]):  # or len(intersect(link[k], good[j])) == 2:
                                    continue
                                a, b = barycentric_intersection_in_line(link[k], Centroid[vertex], mid, Centroid)
                                if (a > 0 and b > 0) or (a < 0 and b > 1):
                                    I = intersect_triangle_line(link[k], Centroid[vertex], mid, Centroid)
                                    a, b, c = barycentric_in_triangle(link[k], I, Centroid)
                                    if a > 0 and b > 0 and c > 0:
                                        if isline([vertex] + link[k], triangulation_sub):
                                            T[j] = False
                                            break
                        affected_facets_versions[index][0].append(good[0])
                        affected_facets_versions[index][1].append(T[0])
                        affected_facets_versions[index][0].append(good[1])
                        affected_facets_versions[index][1].append(T[1])
                        affected_facets_versions[index][0].append(good[2])
                        affected_facets_versions[index][1].append(T[2])
                        affected_facets_versions[index][0].append(bad[0])
                        affected_facets_versions[index][1].append(False)
                    # print(repairable[1][i][index])
                    # tunnels = [list(item) for item in combinations(K5_current, 2)]
                    # print("tunnel:", tunnels[index], good, bad)
        if not try_tunnels:
            for i in range(len(affected_facets[0])):
                if isline([vertex] + affected_facets[0][i], triangulation_sub + predetermined_empty):
                    continue
                T = [isline(affected_facets[0][i], affected_facets[0][item]) for item in range(len(affected_facets[0]))]
                Ind = [index for index in range(len(T)) if T[index]]
                temp = [affected_facets[1][item] for item in Ind]
                if math.nan in temp:
                    continue
                if sum(temp) == len(temp):
                    triangulation_sub.append([vertex] + affected_facets[0][i])
            # link_vts = unique(link)
            if not vts:
                # print("no vertices left!")
                return False, predetermined + triangulation_sub, [math.inf]
            problems = [i for i in vts if
                        not proper_edge_in_triangulation([vertex, i], triangulation_sub, isline(i, boundary_holes))]
            R = not problems and Euler_chr(vertex, triangulation_sub, vts, boundary_holes != [])
            if problems:
                for tri in constricted:
                    if isline(tri + [vertex], predetermined_empty):
                        continue
                    triangulation_sub.append(tri + [vertex])
                    problems_sub = [i for i in tri if not proper_edge_in_triangulation([vertex, i], triangulation_sub,
                                                                                       isline(i, boundary_holes))]
                    if len(problems_sub) < len(intersect(tri, problems)):
                        problems = setdiff(problems, tri)
                        problems = union(problems, problems_sub)
                    else:
                        triangulation_sub = triangulation_sub[:-1]
                # print("constricted:", constricted)
                problems = [i for i in vts if
                            not proper_edge_in_triangulation([vertex, i], triangulation_sub, isline(i, boundary_holes))]
                # print("problems here:", problems)
                R = not problems and Euler_chr(vertex, triangulation_sub, vts, boundary_holes != [])
                # print("problems:", problems)
            # try to remove and add tetrahedra with keeping the predetermined triangulation
            if not R and with_repair:
                while True:
                    triangles = [tri for tri in link if intersect(tri, problems) and not isline([vertex] + tri,
                                                                                                triangulation_predetermined + predetermined_empty)]
                    if not triangles:
                        break
                    affected = unique(triangles)
                    initial = len([i for i in affected if proper_edge_in_triangulation([vertex, i], triangulation_sub,
                                                                                       isline(i, boundary_holes))])
                    initial_vertices = [v for v in union(affected, [vertex]) if
                                        proper_vertex_in_triangulation(v, triangulation_sub, Matrix[v])]
                    gains, gains_vertices = [], []
                    for tri in triangles:
                        triangulation_test = loads(dumps(triangulation_sub))
                        Ind = ind([vertex] + tri, triangulation_sub)
                        if Ind:
                            triangulation_test = Remove_index(triangulation_test, Ind)
                        else:
                            triangulation_test.append([vertex] + tri)
                        gains.append(len([i for i in affected if proper_edge_in_triangulation([vertex, i],
                                                                                              triangulation_test,
                                                                                              isline(i,
                                                                                                     boundary_holes))]) - initial)
                        gains_vertices.append(len([v for v in union(affected, [vertex]) if
                                                   proper_vertex_in_triangulation(v, triangulation_test,
                                                                                  Matrix[v])]) - len(initial_vertices))
                    m = max(gains)
                    if m >= 0:
                        Ind = [in_sub for in_sub in range(len(gains)) if gains[in_sub] == m]
                        for in_sub in setdiff(range(len(gains_vertices)), Ind):
                            gains_vertices[in_sub] = -math.inf
                        n = max(gains_vertices)
                        if m == 0 and n <= 0:
                            break
                        Ind = gains_vertices.index(n)
                        tri = triangles[Ind]
                        Ind = ind([vertex] + tri, triangulation_sub)
                        if Ind:
                            triangulation_sub = Remove_index(triangulation_sub, Ind)
                        else:
                            triangulation_sub.append([vertex] + tri)
                        problems = [i for i in vts if
                                    not proper_edge_in_triangulation([vertex, i], triangulation_sub,
                                                                     isline(i, boundary_holes))]
                    else:
                        break
                R = not problems and Euler_chr(vertex, triangulation_sub, vts, boundary_holes != [])
                # print(R)
            # print("problems here:", problems)
            # try to remove and add tetrahedra to the predetermined triangulation
            if not R and post_process:
                while True:
                    triangles = [tri for tri in link if intersect(tri, problems)]
                    if not triangles:
                        break
                    lst_intersect = [-len(intersect(tri, problems)) for tri in triangles]
                    indices_sorted = list(np.argsort(np.array(lst_intersect)))  # descending order by length
                    triangles = [triangles[item] for item in indices_sorted]
                    initial_vertices = [v for v in union(unique(triangles), [vertex]) if
                                        proper_vertex_in_triangulation(v, triangulation_sub, Matrix[v])]
                    # print("triangles:", triangles)
                    gains, gains_vertices = [], []
                    for tri in triangles:
                        triangulation_test = loads(dumps(triangulation_sub))
                        Ind = ind([vertex] + tri, triangulation_test)
                        if Ind:
                            triangulation_test = Remove_index(triangulation_test, Ind)
                        else:
                            triangulation_test.append([vertex] + tri)
                        edges_proper_link = [list(item) for item in combinations([vertex] + tri, 2) if
                                             proper_edge_in_triangulation(list(item), triangulation_sub)]
                        initial = len(edges_proper_link)
                        edges_proper_link = [list(item) for item in combinations([vertex] + tri, 2) if
                                             proper_edge_in_triangulation(list(item), triangulation_test)]
                        gains.append(len(edges_proper_link) - initial)
                        gains_vertices.append(len([v for v in [vertex] + tri if
                                                   proper_vertex_in_triangulation(v, triangulation_test,
                                                                                  Matrix[v])]) - len(
                            intersect(initial_vertices, [vertex] + tri)))
                    m = max(gains)
                    if m >= 0:
                        Ind = [in_sub for in_sub in range(len(gains)) if gains[in_sub] == m]
                        for in_sub in setdiff(range(len(gains_vertices)), Ind):
                            gains_vertices[in_sub] = -math.inf
                        n = max(gains_vertices)
                        if m == 0 and n <= 0:
                            break
                        Ind = gains_vertices.index(n)
                        tri = triangles[Ind]
                        Ind = ind([vertex] + tri, triangulation_sub)
                        if Ind:
                            triangulation_sub = Remove_index(triangulation_sub, Ind)
                        else:
                            triangulation_sub.append([vertex] + tri)
                        problems = [i for i in vts if
                                    not proper_edge_in_triangulation([vertex, i], triangulation_sub,
                                                                     isline(i, boundary_holes))]
                    else:
                        break
                R = not problems and Euler_chr(vertex, triangulation_sub, vts, boundary_holes != [])
            # try to remove and add tetrahedra to the predetermined triangulation and then retriangulate
            if not R and post_process_further:
                # start_time = time()
                while True:
                    # t = time() - start_time
                    # if t > 20:
                    #     break
                    # print("problems:", problems)
                    triangles = [tri for tri in link if intersect(tri, problems)]
                    if not triangles:
                        break
                    lst_intersect = [-len(intersect(tri, problems)) for tri in triangles]
                    indices_sorted = list(np.argsort(np.array(lst_intersect)))  # descending order by length
                    triangles = [triangles[item] for item in indices_sorted]
                    gains, gains_vertices, retriangulations = [], [], []
                    triangulation_current = loads(dumps(predetermined + triangulation_sub))
                    lst_edges = edges_of(Matrix[vertex] + [vertex], Matrix)
                    lst_vertices = Matrix[vertex] + [vertex]
                    lst_initial = [edge for edge in lst_edges if
                                   proper_edge_in_triangulation(edge, triangulation_current)]
                    lst_initial_vertices = [v for v in lst_vertices if
                                            proper_vertex_in_triangulation(v, triangulation_current, Matrix[v])]
                    # print(len(triangles))
                    for tri in triangles:
                        # print(time() - start_time)
                        Ind = ind([vertex] + tri, triangulation_current)
                        if Ind:
                            retriangulation = Remove_index(triangulation_current, Ind)
                            for v_sub in [vertex] + tri:
                                retriangulation = istriangulatable(v_sub, Matrix, list_boundary, Centroid,
                                                                   changed_tunnels, predetermined=retriangulation,
                                                                   post_process=True)[1]
                        else:
                            retriangulation = loads(dumps(triangulation_current + [[vertex] + tri]))
                            for v_sub in [vertex] + tri:
                                retriangulation = istriangulatable(v_sub, Matrix, list_boundary, Centroid,
                                                                   changed_tunnels, predetermined=retriangulation,
                                                                   post_process=True)[1]
                        retriangulations.append(retriangulation)
                        lst_final_vertices = [v for v in lst_vertices if
                                              proper_vertex_in_triangulation(v, retriangulation, Matrix[v])]
                        lst_final = [edge for edge in lst_edges if
                                     proper_edge_in_triangulation(edge, retriangulation)]
                        gains.append(len(lst_final) - len(lst_initial))
                        gains_vertices.append(len(lst_final_vertices) - len(lst_initial_vertices))
                    m = max(gains)
                    if m >= 0:
                        Ind = [in_sub for in_sub in range(len(gains)) if gains[in_sub] == m]
                        for in_sub in setdiff(range(len(gains_vertices)), Ind):
                            gains_vertices[in_sub] = -math.inf
                        n = max(gains_vertices)
                        if m == 0 and n <= 0:
                            break
                        Ind = gains_vertices.index(n)
                        Ind_sub = ind(vertex, retriangulations[Ind])
                        triangulation_sub = [retriangulations[Ind][index] for index in Ind_sub]
                        predetermined = Remove_index(retriangulations[Ind], Ind_sub)
                        problems = [i for i in vts if
                                    not proper_edge_in_triangulation([vertex, i], triangulation_sub,
                                                                     isline(i, boundary_holes))]
                    else:
                        break
                R = not problems and Euler_chr(vertex, triangulation_sub, vts, boundary_holes != [])
            # print("Newly added:", triangulation_sub[initial_len:])
            return R, predetermined + triangulation_sub, problems

        else:
            R = [True for item in affected_facets_versions]
            problems = [[] for item in affected_facets_versions]
            triangulation_sub_versions = [loads(dumps(triangulation_sub)) for item in affected_facets_versions]
            # print(triangulation_sub)
            # constricted_versions = [loads(dumps(constricted)) for item in affected_facets_versions]
            for index in range(len(affected_facets_versions)):
                affected_facets_version = affected_facets_versions[index]
                for i in range(len(affected_facets_version[0])):
                    if isline([vertex] + affected_facets_version[0][i],
                              triangulation_sub_versions[index] + predetermined_empty):
                        continue
                    T = [isline(affected_facets_version[0][i], affected_facets_version[0][item]) for item in
                         range(len(affected_facets_version[0]))]
                    Ind = [index_sub for index_sub in range(len(T)) if T[index_sub]]
                    temp = [affected_facets_version[1][item] for item in Ind]
                    if math.nan in temp:
                        continue
                    if sum(temp) == len(temp):
                        triangulation_sub_versions[index].append([vertex] + affected_facets_version[0][i])
                # link_vts = unique(link)
                if not vts:
                    R[index] = False
                    continue
                problems[index] = [i for i in vts if
                                   not proper_edge_in_triangulation([vertex, i], triangulation_sub_versions[index],
                                                                    isline(i, boundary_holes))]
                R[index] = not problems[index] and Euler_chr(vertex, triangulation_sub_versions[index], vts,
                                                             boundary_holes != [])
                if not R[index]:
                    for tri in constricted:
                        triangulation_sub_versions[index].append(tri + [vertex])
                        problems_sub = [i for i in tri if
                                        not proper_edge_in_triangulation([vertex, i], triangulation_sub_versions[index],
                                                                         isline(i, boundary_holes))]
                        if len(problems_sub) < len(intersect(tri, problems[index])):
                            problems[index] = setdiff(problems[index], tri)
                            problems[index] = union(problems[index], problems_sub)
                        else:
                            triangulation_sub_versions[index] = triangulation_sub_versions[index][:-1]
                R[index] = not problems[index] and Euler_chr(vertex, triangulation_sub_versions[index], vts,
                                                             boundary_holes != [])
                if not R[index] and with_repair:
                    while True:
                        triangles = [tri for tri in link if
                                     intersect(tri, intersect(problems[index], Rest)) and not isline([vertex] + tri,
                                                                                                     triangulation_predetermined)]
                        for v in intersect(problems[index], Rest):
                            for tri in [link[item] for item in ind(v, link) if
                                        not isline([vertex] + link[item], triangulation_predetermined)]:
                                triangles.append(tri)
                        if not triangles:
                            break
                        affected = unique(triangles)
                        initial = len(
                            [i for i in affected if
                             proper_edge_in_triangulation([vertex, i], triangulation_sub_versions[index],
                                                          isline(i, boundary_holes))])
                        initial_vertices = [v for v in union(affected, [vertex]) if
                                            proper_vertex_in_triangulation(v, triangulation_sub_versions[index],
                                                                           Matrix[v])]
                        gains, gains_vertices = [], []
                        for tri in triangles:
                            triangulation_test = loads(dumps(triangulation_sub_versions[index]))
                            Ind = ind([vertex] + tri, triangulation_sub_versions[index])
                            if Ind:
                                triangulation_test = Remove_index(triangulation_test, Ind)
                            else:
                                triangulation_test.append([vertex] + tri)
                            gains.append(len([i for i in affected if
                                              proper_edge_in_triangulation([vertex, i], triangulation_test,
                                                                           isline(i, boundary_holes))]) - initial)
                            gains_vertices.append(len([v for v in union(affected, [vertex]) if
                                                       proper_vertex_in_triangulation(v, triangulation_test,
                                                                                      Matrix[v])]) -
                                                  len(initial_vertices))
                        m = max(gains)
                        if m >= 0:
                            Ind = [in_sub for in_sub in range(len(gains)) if gains[in_sub] == m]
                            for in_sub in setdiff(range(len(gains_vertices)), Ind):
                                gains_vertices[in_sub] = -math.inf
                            n = max(gains_vertices)
                            if m == 0 and n <= 0:
                                break
                            Ind = gains_vertices.index(n)
                            tri = triangles[Ind]
                            Ind = ind([vertex] + tri, triangulation_sub_versions[index])
                            if Ind:
                                triangulation_sub_versions[index] = Remove_index(triangulation_sub_versions[index], Ind)
                            else:
                                triangulation_sub_versions[index].append([vertex] + tri)
                            problems[index] = [i for i in vts if
                                               not proper_edge_in_triangulation([vertex, i],
                                                                                triangulation_sub_versions[index],
                                                                                isline(i, boundary_holes))]
                        else:
                            break
                    R[index] = not problems[index] and Euler_chr(vertex, triangulation_sub_versions[index], vts,
                                                                 boundary_holes != [])

            return R, [predetermined + triangulation_sub_versions[index] for index in
                       range(len(triangulation_sub_versions))], problems


##################################################################


def constricted_triangles_in_link(vertex, Matrix, list_boundary, Centroid):
    """ detects constrictions and extra edges in the link of the vertex in the clique complex. """

    vts = Matrix[vertex]
    link = vertexlink(vertex, Matrix)
    K4_temp = find_K4s(Matrix, vts)
    K5_temp = find_K5s(Matrix, vts)
    reduced_link = reducedvertexlink(vertex, Matrix, list_boundary, Centroid, with_constriction=True)
    reduced_edges = []
    for i in vts:
        Nei_i = Neighbors_link(i, reduced_link)
        for j in [item for item in Nei_i if item > i]:
            reduced_edges.append([i, j])
    open_reduced_edges = []
    for edge in reduced_edges:
        if len(ind(edge, reduced_link)) != 2:
            open_reduced_edges.append(edge)
    constricted = [tri for tri in cycles_deg(open_reduced_edges, 3) if not isline(tri, K4_temp)]
    Holes = holes(vertex, Matrix, list_boundary, Centroid, True)[0]
    reduced_vertices = unique(reduced_link)
    K4_reduced = [item for item in K4_temp if sublist(item, reduced_vertices)]
    for i in range(len(link)):
        if isline(link[i], K4_temp) or isline(link[i], constricted):
            continue
        if min([len(Neighbors_link(item, link)) for item in link[i]]) < 3:
            continue  # all vertices of a constricted triangle should have degree >= 3
        extra_triangle = False
        count_deg_1_edges = 0  # the number of edges of degree 1 must be 1.
        for edge in combinations(link[i], 2):
            if len(ind(edge, link)) < 2 and not isline(edge, Holes):
                count_deg_1_edges += 1
        if count_deg_1_edges >= 1:
            constricted.append(link[i])
            continue
        elif count_deg_1_edges == 0:
            continue
        for v in link[i]:
            triangles = [link[item] for item in ind(v, link)]
            if sum([isline(tri, K4_temp) for tri in triangles]) == len(triangles):
                # if all triangles surrounding 'v' lie in 4-cliques then the next argument does not work.
                continue
            if v in reduced_vertices:
                num_edges = len(Neighbors_link(v, reduced_link))
                Ind = ind(v, K4_reduced)
                num_triangles = len(ind(v, reduced_link)) - len(Ind)
                if num_triangles > num_edges:
                    extra_triangle = True
                    break
            if isline(v, K5_temp):
                continue
            num_edges = len(mutual_neighbors([vertex, v], Matrix))
            Ind = ind(v, K4_temp)
            num_triangles = len(ind(v, link)) - len(Ind)
            if num_triangles > num_edges:
                extra_triangle = True
                break
        if not extra_triangle:
            continue
        constricted.append(link[i])
    return constricted


##################################################################


def proper_edge_in_triangulation(edge, triangulation, is_boundary_edge=False):
    """ checks whether the edge has a proper link (i.e. a cycle) in the triangulation. """

    edges = edge_link_in_triangulation(edge, triangulation)
    if not edges:
        return False
    for j in range(len(edges)):
        if len(edges[j]) != 2:
            return False
    C = cycles(edges)
    if len(C) > 1:
        return False
    elif len(C) == 0:
        if not is_boundary_edge:
            return False
    else:
        if len(C[0]) != len(edges):
            return False
    return True


##################################################################


def edge_link_in_triangulation(edge, triangulation):
    """ returns the link of edge in the triangulation. """

    return [setdiff(triangulation[item], edge) for item in ind(edge, triangulation)]


##################################################################


def proper_vertex_in_triangulation(vertex, triangulation, predetermined_neighbors=None, return_problems=False):
    """ checks whether the vertex has a proper link (i.e. a 2-sphere) in the triangulation. """

    if predetermined_neighbors is None:
        predetermined_neighbors = []
    T, problems = True, []
    vts = neighbors_in_triangulation(vertex, triangulation)
    if predetermined_neighbors:
        vts = predetermined_neighbors
    if not vts:
        if not return_problems:
            return False
        else:
            return False, []
    for v_sub in vts:
        if not proper_edge_in_triangulation([vertex, v_sub], triangulation):
            if not return_problems:
                return False
            else:
                T = False
                problems.append(v_sub)
    if not problems:
        if not Euler_chr(vertex, triangulation, predetermined_neighbors):
            if not return_problems:
                return False
            else:
                T = False
                problems.append(-1)
    if return_problems:
        return T, problems
    return True


##################################################################


def Euler_chr(vertex, triangulation, predetermined_neighbors=None, punctured_sphere=False):
    """ checks whether the link of a vertex in a triangulation has the Euler characteristic of a 2-sphere. """

    if predetermined_neighbors is None:
        predetermined_neighbors = []
    vts = neighbors_in_triangulation(vertex, triangulation)
    if predetermined_neighbors:
        vts = predetermined_neighbors
    link = [setdiff(triangulation[item], vertex) for item in ind(vertex, triangulation)]
    edges = []
    for tri in link:
        for edge in combinations(tri, 2):
            if not isline(edge, edges):
                edges.append(edge)
    if punctured_sphere:
        return len(vts) - len(edges) + len(link) < 2
    return len(vts) - len(edges) + len(link) == 2


##################################################################


def vertex_link_in_triangulation(vertex, triangulation, forbidden_tetrahedra=None):
    """ returns the link of a vertex in a triangulation. """
    
    if forbidden_tetrahedra is None:
        return [setdiff(triangulation[item], vertex) for item in ind(vertex, triangulation)]
    return [setdiff(triangulation[item], vertex) for item in ind(vertex, triangulation) if
            not isline(triangulation[item], forbidden_tetrahedra)]


##################################################################


def candidate_cycles(cycles, edges_1, edges_2):
    """ computes which cycles contain all the edges in edges_1 and none of the edges in edges_2. """
    
    T = [True for item in cycles]
    for j in range(len(cycles)):
        for k in range(len(edges_1)):
            T[j] = T[j] and isline(edges_1[k], cycles[j])
    for j in range(len(cycles)):
        if not T[j]:
            continue
        for k in range(len(edges_1)):
            T[j] = T[j] and is_edge_in_cycle(edges_1[k], cycles[j])
        for k in range(len(edges_2)):
            if not isline(edges_2[k], cycles[j]):
                continue
            T[j] = T[j] and not is_edge_in_cycle(edges_2[k], cycles[j])
    return T


##################################################################


def Recognize_combinatorial(Matrix, list_boundary):
    """ reconstructs a partial triangulation on the adjacency list Matrix using only combinatorics. """
    
    Matrix_out = loads(dumps(Matrix))
    tetrahedra_certified = []
    tetrahedra_empty = []
    numbPoints = len(Matrix_out)
    K2 = []
    for i in range(numbPoints):
        Nei_i = Matrix_out[i]
        for j in [item for item in Nei_i if item > i]:
            if not sublist([i, j], list_boundary):
                K2.append([i, j])
    K4 = find_K4s(Matrix_out)
    K5 = find_K5s(Matrix_out)

    # Links of edges and their cycles...
    Links_of_edges = [loads(dumps(K2)), [], []]
    for i in range(len(K2)):
        Links_of_edges[1].append([])
        for j in ind(K2[i], K4):
            Links_of_edges[1][i].append(setdiff(K4[j], Links_of_edges[0][i]))
        Links_of_edges[2].append([])
        if not Links_of_edges[2][i]:
            Links_of_edges[2][i] = cycles(Links_of_edges[1][i])

    # Identifying tetrahedra...
    deg_4_vertices = []
    for vertex in range(numbPoints):
        if len(Matrix_out[vertex]) == 4:
            deg_4_vertices.append(vertex)

    Certified = [[] for item in range(numbPoints)]
    Decertified = [[] for item in range(numbPoints)]
    Undecided = [[K4[i] for i in ind([vertex], K4)] for vertex in range(numbPoints)]
    Links = [[setdiff(K4[i], vertex) for i in ind([vertex], K4)] for vertex in range(numbPoints)]

    first_time = True
    affected = []
    Repeat_indices = range(len(Links_of_edges[0]))
    while True:
        Out = True
        iteration = 0
        while True:
            out = True
            iteration += 1
            if first_time:
                affected_prev = loads(dumps(K4))
            else:
                affected_prev = loads(dumps(affected))
            affected = []
            first_time = False
            i = -1
            while i < len(K2) - 1:
                i += 1
                print("\r", i, "out of", len(K2) - 1, "iteration", iteration, end="")
                if i not in Repeat_indices or not isline(K2[i], affected_prev):
                    continue
                if not Links_of_edges[2][i]:
                    continue
                edge = K2[i]
                if intersect(edge, list_boundary):
                    continue
                edges_certified = [setdiff(item, edge) for item in intersect(Certified[edge[0]], Certified[edge[1]])]
                edges_decertified = [setdiff(item, edge) for item in
                                     intersect(Decertified[edge[0]], Decertified[edge[1]])]
                edges_undecided = [setdiff(item, edge) for item in intersect(Undecided[edge[0]], Undecided[edge[1]])]
                if not edges_undecided:
                    Repeat_indices = setdiff(Repeat_indices, [i])
                    continue
                repeat = False

                T = candidate_cycles(Links_of_edges[2][i], edges_certified, edges_decertified)
                if sum(T) == 1:
                    Ind = T.index(True)
                    for k in range(len(Links_of_edges[2][i][Ind])):
                        l = Next(k, 1, len(Links_of_edges[2][i][Ind]))
                        tet = Links_of_edges[0][i] + [Links_of_edges[2][i][Ind][k], Links_of_edges[2][i][Ind][l]]
                        tet.sort()
                        if not isline(tet, tetrahedra_certified):
                            tetrahedra_certified.append(tet)
                            length_added = 1
                            retry = True
                            for index in range(- length_added, 0):
                                for triangle in combinations(tetrahedra_certified[index], 3):
                                    Ind_triangle = ind(triangle, tetrahedra_certified)
                                    if len(Ind_triangle) > 2 and retry:
                                        retry = False
                                        tetrahedra_certified = tetrahedra_certified[:-length_added]
                            if retry:
                                affected.append(tet)
                                edges_certified.append(Links_of_edges[2][i][Ind][k])
                                edges_decertified = [item for item in edges_decertified if not isline(item, [
                                    Links_of_edges[2][i][Ind][k], Links_of_edges[2][i][Ind][l]])]
                                edges_undecided = [item for item in edges_undecided if not isline(item, [
                                    Links_of_edges[2][i][Ind][k], Links_of_edges[2][i][Ind][l]])]
                                for vertex_sub in tet:
                                    Undecided[vertex_sub].remove(tet)
                                    Certified[vertex_sub].append(tet)
                                out = False
                                repeat = True
                                Out = False
                    for k in range(len(Links_of_edges[1][i])):
                        if isedgeincycle(Links_of_edges[1][i][k], Links_of_edges[2][i][Ind]):
                            continue
                        tet = sorted(Links_of_edges[0][i] + Links_of_edges[1][i][k])
                        if not isline(tet, tetrahedra_empty):
                            tetrahedra_empty.append(tet)
                            for vertex_sub in tet:
                                Undecided[vertex_sub].remove(tet)
                                Decertified[vertex_sub].append(tet)
                                Links[vertex_sub].remove(sorted(setdiff(tet, vertex_sub)))
                            affected.append(tet)
                            edges_certified = [item for item in edges_certified if item != Links_of_edges[1][i][k]]
                            edges_decertified.append(Links_of_edges[1][i][k])
                            edges_undecided = [item for item in edges_undecided if item != Links_of_edges[1][i][k]]
                            out = False
                            repeat = True
                            Out = False
                    if repeat:
                        i -= 1
                        continue
                else:
                    # tag an unknown edge in the link as lying in a certified
                    # tetrahedron. If this leads to inconsistency, the edge must
                    # lie in an empty tetrahedron.
                    T = [[] for item in edges_undecided]
                    for j in range(len(edges_undecided)):
                        edges_temp = edges_certified + [edges_undecided[j]]
                        T[j] = candidate_cycles(Links_of_edges[2][i], edges_temp, edges_decertified)
                    Ind = [index for index in range(len(T)) if sum(T[index]) == 0]
                    for j in Ind:
                        tet = sorted(Links_of_edges[0][i] + edges_undecided[j])
                        if not isline(tet, tetrahedra_empty):
                            tetrahedra_empty.append(tet)
                            for vertex_sub in tet:
                                Undecided[vertex_sub].remove(tet)
                                Decertified[vertex_sub].append(tet)
                                Links[vertex_sub].remove(sorted(setdiff(tet, vertex_sub)))
                            affected.append(tet)
                            edges_certified = [item for item in edges_certified if item != edges_undecided[j]]
                            edges_decertified.append(edges_undecided[j])
                            out = False
                            repeat = True
                            Out = False
                    if repeat:
                        i -= 1
                        continue
                    # tag an unknown edge in the link as lying in an empty
                    # tetrahedron. If this leads to inconsistency, the edge must
                    # lie in a certified tetrahedron.
                    T = [[] for item in edges_undecided]
                    for j in range(len(edges_undecided)):
                        edges_temp = edges_decertified + [edges_undecided[j]]
                        T[j] = candidate_cycles(Links_of_edges[2][i], edges_certified, edges_temp)
                    Ind = [index for index in range(len(T)) if sum(T[index]) == 0]
                    for j in Ind:
                        tet = sorted(Links_of_edges[0][i] + edges_undecided[j])
                        if not isline(tet, tetrahedra_certified):
                            tetrahedra_certified.append(tet)
                            length_added = 1
                            retry = True
                            for index in range(- length_added, 0):
                                for triangle in combinations(tetrahedra_certified[index], 3):
                                    Ind_triangle = ind(triangle, tetrahedra_certified)
                                    if len(Ind_triangle) > 2 and retry:
                                        retry = False
                                        tetrahedra_certified = tetrahedra_certified[:-length_added]
                            if retry:
                                affected.append(tet)
                                edges_certified.append(edges_undecided[j])
                                edges_decertified = [item for item in edges_decertified if item != edges_undecided[j]]
                                for vertex_sub in tet:
                                    Undecided[vertex_sub].remove(tet)
                                    Certified[vertex_sub].append(tet)
                                out = False
                                repeat = True
                                Out = False
                    if repeat:
                        i -= 1
                        continue
                edges_undecided = [item for item in Links_of_edges[1][i] if
                                   item not in [item for item in Links_of_edges[1][i] if
                                                item in edges_certified or item in edges_decertified]]
                if not edges_undecided:
                    Repeat_indices = setdiff(Repeat_indices, [i])
            if out:
                break
        print("")
        # Each K5 that does not have a degree-4-vertex can have at most 3 certified
        # tetrahedra. This holds only for 5-cliques of actual vertices (without the
        # vertex at infinity)
        for i in range(len(K5)):
            if intersect(K5[i], deg_4_vertices):
                continue
            lst = K5[i]
            lst_temp = []
            for j in combinations(lst, 4):
                j = list(j)
                if not isline(j, tetrahedra_certified):
                    lst_temp.append(j)
            if len(lst_temp) != 2:
                continue
            for s in [0, 1]:
                if not isline(lst_temp[s], tetrahedra_empty):
                    tet = sorted(lst_temp[s])
                    tetrahedra_empty.append(tet)
                    for vertex_sub in tet:
                        Undecided[vertex_sub].remove(tet)
                        Decertified[vertex_sub].append(tet)
                        Links[vertex_sub].remove(sorted(setdiff(tet, vertex_sub)))
                    affected.append(tet)
                    Out = False
        print("")
        if Out:
            break

    edges_defect_link = []
    for i in range(numbPoints):
        Nei = Matrix_out[i]
        for j in [item for item in Nei if item > i]:
            if sublist([i, j], list_boundary):
                continue
            edges = edge_link_in_triangulation([i, j], tetrahedra_certified)
            C = cycles(edges)
            if len(C) != 1:
                edges_defect_link.append([i, j])
                continue
            if len(C[0]) != len(edges):
                edges_defect_link.append([i, j])
                continue
    vts_defect_links = []
    for edge in edges_defect_link:
        vts_defect_links = union(vts_defect_links, edge)
    vts_defect_links = setdiff(vts_defect_links, list_boundary)

    return tetrahedra_certified, tetrahedra_empty, edges_defect_link, vts_defect_links


##################################################################


def tunnel_K5(K5, Matrix, Centroid, list_boundary):
    """ computes the tunnel of a 5-clique. """

    count = 0
    for j in combinations(K5, 2):
        j = list(j)
        if len(mutual_neighbors(j, Matrix)) == 3:
            count += 1
            tunnel = j
    if count == 1 and not sublist(tunnel, list_boundary):
        return tunnel
    tunnel = []
    count = 0
    for edge in combinations(K5, 2):
        edge = list(edge)
        a, b, c = barycentric_intersection_in_triangle(setdiff(K5, edge), edge, Centroid)
        if a > 0 and b > 0 and c > 0:
            count += 1
            tunnel = edge
    if count != 1:
        tunnel = []
    else:
        fine = []
        for j in setdiff(mutual_neighbors(tunnel, Matrix), K5):
            lst = intersect(Matrix[j], setdiff(K5, tunnel))
            if len(lst) >= 2:
                fine = union(fine, lst)
        for j in setdiff(mutual_neighbors(tunnel, Matrix), K5):
            lst = intersect(Matrix[j], setdiff(K5, tunnel))
            if len(lst) == 1:
                lst = lst[0]
                if lst in fine:
                    continue
                eqn = eqn_plane_3_points(tunnel + [lst], Centroid)
                rest = setdiff(K5, tunnel + [lst])
                for k in rest:
                    if evaluate(eqn, Centroid[j]) * evaluate(eqn, Centroid[k]) > 0:
                        tunnel = [lst, k]
                        break
    return tunnel


##################################################################


def tunnel_K5_geometric(K5, Centroid):
    """ computes the geometric tunnel of a 5-clique. """

    tunnel = []
    count = 0
    for edge in combinations(K5, 2):
        edge = list(edge)
        a, b, c = barycentric_intersection_in_triangle(setdiff(K5, edge), edge, Centroid)
        if a > 0 and b > 0 and c > 0:
            count += 1
            tunnel = edge
    if count != 1:
        tunnel = []
    return tunnel


##################################################################


def eqn_plane_3_points(triple, Centroid=None):
    """ computes the equation 'ax + by + cz + d = 0' of the plane passing through 3 points and returns the coefficients
        a, b, c, d. """

    triple_new = []
    for i in triple:
        if isinstance(i, int):
            triple_new.append(np.array(Centroid[i]))
        else:
            triple_new.append(np.array(i))
    p1, p2, p3 = triple_new
    v1 = p3 - p1
    v2 = p2 - p1
    cp = np.cross(v1, v2)
    a, b, c = cp
    d = - np.dot(cp, p3)
    return [a, b, c, d]


##################################################################


def evaluate(eqn, pt, Centroid=None):
    """ evaluates the equation of a plane given by eqn at the point pt. """

    if isinstance(pt, int):
        pt = Centroid[pt]
    return eqn[0] * pt[0] + eqn[1] * pt[1] + eqn[2] * pt[2] + eqn[3]


##################################################################


def plot_edge_link(link_edge, Matrix, Certified=None, Decertified=None):
    """ plots the (1-skeleton of the) link of an edge in the clique complex. """

    if Certified is None:
        Certified = []
    if Decertified is None:
        Decertified = []

    # Compute the nodes
    vts = mutual_neighbors(link_edge, Matrix)

    # Compute the edges
    edges = []
    for i in vts:
        for j in [item for item in intersect(Matrix[i], vts) if item > i]:
            edges.append([i, j])

    # Create an empty graph with no nodes or edges
    G = nx.Graph()

    # Add nodes to the graph, with the node index as the node label
    for i in vts:
        G.add_node(i)

    # Add edges to the graph, with the edge index as the edge label
    for edge in edges:
        G.add_edge(edge[0], edge[1])

    # Use the spring layout to position the nodes in the graph
    pos = nx.spring_layout(G)

    # Figure out the edge labels
    labels = []
    for i in range(len(edges)):
        if isline(edges[i] + link_edge, Certified):
            labels.append("1")
        elif isline(edges[i] + link_edge, Decertified):
            labels.append("0")
        else:
            labels.append("")

    # Create a dictionary of edge labels
    edge_labels = {(edges[item][0], edges[item][1]): labels[item] for item in range(len(edges))}

    # Draw the graph, with the node labels and edge labels
    nx.draw(G, pos, with_labels=True)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

    # Add a title
    plt.title("Link of edge [" + str(link_edge[0]) + "," + str(link_edge[1]) + "]")

    # Show the plot
    plt.show()


##################################################################


def plot_vertex_link(vertex, Matrix, list_boundary, Centroid, Certified=None, Decertified=None, reduced=False,
                     only_certified=False, highlight=None, highlight_edges=None, show_steps=False, view=None):
    """ plots the (2-skeleton of the) link of a vertex in the clique complex. """

    if Certified is None:
        Certified = []
    if Decertified is None:
        Decertified = []
    if highlight is None:
        highlight = []
    if highlight_edges is None:
        highlight_edges = []
    if view is None:
        view = [30, -60]
    if isinstance(highlight, int):
        highlight = [highlight]
    fig = plt.figure()
    ax = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)

    link = vertexlink(vertex, Matrix)
    if reduced:
        link = reducedvertexlink(vertex, Matrix, list_boundary, Centroid, show_steps=show_steps)
    if only_certified:
        link = [setdiff(Certified[item], vertex) for item in ind(vertex, Certified)]
        print("faces:", len(link))
    edges = []
    vts = Matrix[vertex]
    for j in vts:
        Nei_j = intersect(vts, Matrix[j])
        for k in [item for item in Nei_j if item > j]:
            if not isline([j, k], edges):
                edges.append([j, k])
    if only_certified:
        edges = [edge for edge in edges if isline(edge, link)]
        print("edges:", len(edges))
        print("vertices:", len(vts))
    c = "k"
    if vertex in list_boundary:
        c = "r"
    ax.scatter3D(Centroid[vertex][0], Centroid[vertex][1], Centroid[vertex][2], color=c)
    ax.text(Centroid[vertex][0], Centroid[vertex][1], Centroid[vertex][2], str(vertex))
    ax.set_axis_off()
    ax.grid(False)
    lst = setdiff(vts, list_boundary)
    lst = setdiff(lst, highlight)
    x = [Centroid[item][0] for item in lst]
    y = [Centroid[item][1] for item in lst]
    z = [Centroid[item][2] for item in lst]
    ax.scatter3D(x, y, z, color="k")

    lst = intersect(vts, list_boundary)
    lst = setdiff(lst, highlight)
    x = [Centroid[item][0] for item in lst]
    y = [Centroid[item][1] for item in lst]
    z = [Centroid[item][2] for item in lst]
    ax.scatter3D(x, y, z, color="r")

    if highlight:
        x = [Centroid[item][0] for item in highlight]
        y = [Centroid[item][1] for item in highlight]
        z = [Centroid[item][2] for item in highlight]
        ax.scatter3D(x, y, z, color=[0, 1, 1], s=40)
    for vertex_sub in vts:
        ax.text(Centroid[vertex_sub][0], Centroid[vertex_sub][1], Centroid[vertex_sub][2], str(vertex_sub))

    for triangle in link:
        color = divied_scalar(sum(triangle), triangle)
        x = [Centroid[item][0] for item in triangle]
        y = [Centroid[item][1] for item in triangle]
        z = [Centroid[item][2] for item in triangle]
        verts = [list(zip(x, y, z))]
        ax.add_collection3d(Poly3DCollection(verts, alpha=.25, facecolors=color))
        t = ""
        if not only_certified:
            if isline([vertex] + triangle, Certified):
                t = "1"
            elif isline([vertex] + triangle, Decertified):
                t = "0"
        mid = [sum(x) / 3, sum(y) / 3, sum(z) / 3]
        ax.text(mid[0], mid[1], mid[2], t, fontsize=6)

    for edge in edges:
        start = Centroid[edge[0]]
        end = Centroid[edge[1]]
        c, l = "k", 1.5
        if isline(edge, highlight_edges):
            c, l = "r", 3
        ax.plot([start[0], end[0]], [start[1], end[1]], zs=[start[2], end[2]], alpha=.4, color=c, linewidth=l)
    ax.view_init(elev=view[0], azim=view[1])
    plt.show()


##################################################################


def compute_time(start_time, with_return=False):
    """ returns the elapsed time from a given starting time. """

    t = time() - start_time
    if t < 60:
        t = int(t)
        s = "seconds"
        if t == 1:
            s = "second"
        if not with_return:
            print("Elapsed time: {} {}".format(t, s))
        else:
            return "Elapsed time: {} {}".format(t, s)
    elif t < 3600:
        t1 = int(t / 60)
        t2 = int(t % 60) + 1
        s1 = "minutes"
        if t1 == 1:
            s1 = "minute"
        s2 = "seconds"
        if t2 == 1:
            s2 = "second"
        if not with_return:
            print("Elapsed time: {} {} and {} {}".format(t1, s1, t2, s2))
        else:
            return "Elapsed time: {} {} and {} {}".format(t1, s1, t2, s2)
    elif t < 86400:
        t1 = int(t / 3600)
        t2 = int(t % 3600 / 60) + 1
        s1 = "hours"
        if t1 == 1:
            s1 = "hour"
        s2 = "minutes"
        if t2 == 1:
            s2 = "minute"
        if not with_return:
            print("Elapsed time: {} {} and {} {}".format(t1, s1, t2, s2))
        else:
            return "Elapsed time: {} {} and {} {}".format(t1, s1, t2, s2)
    else:
        t1, t = int(t / 86400), int(t % 86400)
        t2 = int(t / 3600)
        t3 = int(t % 3600 / 60) + 1
        s1 = "days"
        if t1 == 1:
            s1 = "day"
        s2 = "hours"
        if t2 == 1:
            s2 = "hour"
        s3 = "minutes"
        if t3 == 1:
            s3 = "minute"
        if not with_return:
            print("Elapsed time: {} {}, {} {} and {} {}".format(t1, s1, t2, s2, t3, s3))
        else:
            return "Elapsed time: {} {}, {} {} and {} {}".format(t1, s1, t2, s2, t3, s3)


##################################################################


def Remove_index(lst, indices):
    """ removes items from lst with given indices. """

    return [lst[item] for item in setdiff(range(len(lst)), indices)]


##################################################################


def mutual_neighbors(lst, Matrix):
    """ returns the list of mutual neighbors of vertices in lst. """

    if isinstance(lst, int):
        lst = [lst]
    output = range(len(Matrix))
    for i in lst:
        output = intersect(output, Matrix[i])
    return setdiff(output, lst)


##################################################################


def holes(vertex, Matrix, list_boundary, Centroid, low_degree_allowed=False, complete_holes_allowed=False,
          holes_degree=None, allow_stacked_tetrahedron=True, reduced=False):
    """ computes the holes in the vertex link in the clique complex. """

    if vertex in list_boundary:
        return [], [], []
    Matrix_original = loads(dumps(Matrix))

    if allow_stacked_tetrahedron:
        # A stacked tetrahedron contains no holes:
        Matrix_temp = loads(dumps(Matrix))
        while True:
            out = True
            lst = Matrix_temp[vertex]
            for i in lst:
                lst_sub = mutual_neighbors([i, vertex], Matrix_temp)
                if len(lst_sub) == 3 and is_clique(lst_sub, Matrix_temp) and len(Matrix_temp[vertex]) > 4:
                    Matrix_temp = Remove(vertex, i, Matrix_temp)
                    out = False
            if out:
                break
        if is_tetrahedron(Matrix_temp[vertex], Matrix_temp):
            return [], [], []

    vts_original = Matrix[vertex]

    # Remove degree-3-vertices forming apexes of pyramids
    deg_3 = []
    for i in vts_original:
        lst_sub = mutual_neighbors([i, vertex], Matrix)
        if len(lst_sub) == 3 and is_clique(lst_sub, Matrix):
            deg_3.append(i)
    for i in deg_3:
        Nei = mutual_neighbors([i, vertex], Matrix)
        if len(Nei) != 3:
            continue
        if 3 in [len(mutual_neighbors([vertex, v], Matrix)) for v in Nei]:
            continue
        tunnel = tunnel_K5_geometric([i, vertex] + Nei, Centroid)
        if tunnel and isline(tunnel, [i, vertex]):
            # print("remove:", i)
            for j in loads(dumps((Matrix[i]))):
                Matrix = Remove(i, j, Matrix)
        else:
            T, apex = ispyramid(Nei + [i], vertex, Centroid)
            if T and apex == i:
                for j in loads(dumps((Matrix[i]))):
                    Matrix = Remove(i, j, Matrix)

    vts = Matrix[vertex]
    link = vertexlink(vertex, Matrix)
    if reduced:
        link = reducedvertexlink(vertex, Matrix, list_boundary, Centroid)
    # print(link)

    low_degree = []
    for j in vts:
        Nei = intersect(Matrix[j], vts)
        if len(Nei) < 3:
            low_degree.append(j)
    if low_degree and not low_degree_allowed:
        return [], [], []

    # compute the edges and the 4-cliques containing each edge
    edges = []
    four_cliques = []
    five_cliques = []
    six_cliques = []
    deg_3_edges_in_cross = []
    for j in vts:
        Nei_j = intersect(vts, Matrix[j])
        for k in [item for item in Nei_j if item > j]:
            if link and (len(intersect(Matrix[j], vts)) < 2 or len(intersect(Matrix[k], vts)) < 2):
                continue
            edges.append([j, k])
            Nei_j_k = intersect(Nei_j, Matrix[k])
            for l in [item for item in Nei_j_k if item > k]:
                Nei_j_k_l = intersect(Nei_j_k, Matrix[l])
                for m in [item for item in Nei_j_k_l if item > l]:
                    four_cliques.append([j, k, l, m])
                    Nei_j_k_l_m = intersect(Nei_j_k_l, Matrix[m])
                    for n in [item for item in Nei_j_k_l_m if item > m]:
                        five_cliques.append([j, k, l, m, n])
                        Nei_j_k_l_m_n = intersect(Nei_j_k_l_m, Matrix[n])
                        for o in [item for item in Nei_j_k_l_m_n if item > n]:
                            six_cliques.append([j, k, l, m, n, o])
    # print("edges:", edges)
    # remove extra degree-one edges
    lst_remove = []
    for index in range(len(edges)):
        edge = edges[index]
        if len(ind(edge, link)) == 1:
            Matrix = Remove(edge[0], edge[1], Matrix)
            lst_remove.append(index)
    if lst_remove:
        link = vertexlink(vertex, Matrix)
        if reduced:
            link = reducedvertexlink(vertex, Matrix, list_boundary, Centroid)
        edges = Remove_index(edges, lst_remove)
    num_4_cliques = [0 for item in edges]
    for clique in four_cliques:
        tunnel = tunnel_K5_geometric([vertex] + clique, Centroid)
        if tunnel:
            if vertex in tunnel:
                apex = setdiff(tunnel, [vertex])[0]
                base = setdiff(clique, [apex])
                for edge in combinations(base, 2):
                    num_4_cliques[ind(list(edge), edges)[0]] += 1
            else:
                diagonals = [tunnel, setdiff(clique, tunnel)]
                for edge in combinations(clique, 2):
                    edge = sorted(list(edge))
                    if isline(edge, diagonals):
                        continue
                    num_4_cliques[ind(edge, edges)[0]] += 1
                    if not isline(edge, deg_3_edges_in_cross):
                        deg_3_edges_in_cross.append(edge)
        else:
            T, apex = ispyramid(clique, vertex, Centroid)
            if T:
                base = setdiff(clique, [apex])
                for edge in combinations(base, 2):
                    num_4_cliques[ind(list(edge), edges)[0]] += 1
                continue
            T, diagonals = iscross(clique, vertex, Centroid)
            if T:
                for edge in combinations(clique, 2):
                    edge = sorted(list(edge))
                    if isline(edge, diagonals):
                        continue
                    num_4_cliques[ind(edge, edges)[0]] += 1
                    if not isline(edge, deg_3_edges_in_cross):
                        deg_3_edges_in_cross.append(edge)
    for clique in five_cliques:
        for edge in combinations(clique, 2):
            num_4_cliques[ind(list(edge), edges)[0]] -= 1
    for clique in six_cliques:
        for edge in combinations(clique, 2):
            num_4_cliques[ind(list(edge), edges)[0]] += 1
    if reduced:
        vts = unique(link)
        edges = []
        for j in vts:
            Nei_j = Neighbors_link(j, link)
            for k in [item for item in Nei_j if item > j]:
                edges.append([j, k])
        num_4_cliques = [0 for item in edges]
    open_edges = []
    for i in range(len(edges)):
        Ind = ind(edges[i], link)
        if len(Ind) < 2 + num_4_cliques[i]:
            open_edges.append(edges[i])
    # print(len(open_edges), open_edges)
    # return open_edges

    if holes_degree is None:
        temp = cycles(open_edges, minimum_allowed=4)  # cycles of open edges (candidate holes)
    else:
        temp = cycles_deg(open_edges, holes_degree)
    # print("candidates:", len(temp))
    Holes = []
    com_original = components_list(vts_original, Matrix_original)
    for i in range(len(temp)):
        priority = []
        if len(temp[i]) <= 3:
            continue
        # if len(temp[i]) == 4 and sum([isline(edge, temp[i]) for edge in edges]) == 5:
        #     continue
        d = sum([isline(item, temp[i]) for item in four_cliques])
        if sum([isline(edge, temp[i]) for edge in edges]) - d >= 2 * len(
                temp[i]) - 3:  # the hole is already triangulated
            continue
        real_hole = True
        count = 0
        for u in temp[i]:
            lst = intersect(temp[i], Matrix[u])
            if len(lst) == 2:
                continue
            for v in [item for item in lst if item > u]:
                if temp[i].index(u) == Next(temp[i].index(v), 1, len(temp[i])) or temp[i].index(u) == Previous(
                        temp[i].index(v), 1, len(temp[i])):
                    continue
                if not isline([u, v], link):  # there is already a chord in the cycle, so it is not a hole.
                    real_hole = False
                    break
                else:
                    mut = mutual_neighbors([u, v], Matrix)
                    if not isline([u, v], four_cliques) or not intersect(mut, temp[i]):
                        count += 1
        if count >= len(temp[i]) - 3 and len(temp[i]) != len(vts):
            # a polygonal hole with n vertices and n - 3 inner edges (diagonals) is a triangulated polygon,
            # so we consider it a hole only if it is the whole link.
            continue
        if not real_hole:
            continue
        c = len(com_original)
        com = components_list(setdiff(vts_original, temp[i]), Matrix_original)
        diff = len(com) - c
        hole = temp[i]
        if diff > 0:
            count = 0
            isolated_vertices = [item for item in com if len(item) == 1]
            for item in isolated_vertices:
                if item not in com_original and len(intersect(Matrix_original[item[0]], vts)) == 1:
                    count += 1
                    priority.append(item[0])
            if count != diff:
                continue
            replace = {}
            for v in priority:
                lst = intersect(Matrix_original[v], temp[i])
                u = lst[0]
                if u in replace:
                    replace[u].extend([v, u])
                else:
                    replace[u] = [u, v, u]
            hole = []
            for u in temp[i]:
                if u in replace:
                    hole.extend(replace[u])
                else:
                    hole.append(u)
        if not is_clique(temp[i], Matrix):
            Holes.append(hole)
        elif complete_holes_allowed:
            Holes.append(hole)

    lst_remove = []
    for i in range(len(Holes)):
        for j in [item for item in setdiff(range(len(Holes)), i)]:
            if isline(Holes[i], Holes[j]):
                lst_remove = union(lst_remove, j)
            if isline(Holes[j], Holes[i]):
                lst_remove = union(lst_remove, i)
    Holes = [Holes[item] for item in setdiff(range(len(Holes)), lst_remove)]

    # order the holes by their size (number of vertices). If two holes have the same size, order by the number of edges.
    Holes_ordered = []
    while Holes:
        m = min([len(item) for item in Holes])
        S = [item for item in Holes if len(item) == m]
        if len(S) == 1:
            Holes_ordered.append(S[0])
            Holes.remove(S[0])
        else:
            while S:
                m_sub = min([sum([isline(edge, item) for edge in edges]) for item in S])
                Ind = [index for index in range(len(S)) if sum([isline(edge, S[index]) for edge in edges]) == m_sub][0]
                Holes_ordered.append(S[Ind])
                Holes.remove(S[Ind])
                S.remove(S[Ind])
    return Holes_ordered, deg_3_edges_in_cross, open_edges


##################################################################


def exposed_edges_in_link(vertex, Matrix, Centroid):
    """ returns exposed edges in the link of a vertex in the clique complex. """

    # Remove degree-3-vertices forming apexes of pyramids
    deg_3 = []
    vts_original = Matrix[vertex]
    for i in vts_original:
        lst_sub = mutual_neighbors([i, vertex], Matrix)
        if len(lst_sub) == 3 and is_clique(lst_sub, Matrix):
            deg_3.append(i)
    for i in deg_3:
        Nei = mutual_neighbors([i, vertex], Matrix)
        if len(Nei) != 3:
            continue
        tunnel = tunnel_K5_geometric([i, vertex] + Nei, Centroid)
        if tunnel and isline(tunnel, [i, vertex]):
            for j in loads(dumps((Matrix[i]))):
                Matrix = Remove(i, j, Matrix)
        else:
            T, apex = ispyramid(Nei + [i], vertex, Centroid)
            if T and apex == i:
                for j in loads(dumps((Matrix[i]))):
                    Matrix = Remove(i, j, Matrix)

    vts = Matrix[vertex]
    link = vertexlink(vertex, Matrix)

    # compute the edges and the 4-cliques containing each edge
    edges = []
    four_cliques = []
    five_cliques = []
    six_cliques = []
    deg_3_edges_in_cross = []
    for j in vts:
        Nei_j = intersect(vts, Matrix[j])
        for k in [item for item in Nei_j if item > j]:
            if link and (len(intersect(Matrix[j], vts)) < 2 or len(intersect(Matrix[k], vts)) < 2):
                continue
            edges.append([j, k])
            Nei_j_k = intersect(Nei_j, Matrix[k])
            for l in [item for item in Nei_j_k if item > k]:
                Nei_j_k_l = intersect(Nei_j_k, Matrix[l])
                for m in [item for item in Nei_j_k_l if item > l]:
                    four_cliques.append([j, k, l, m])
                    Nei_j_k_l_m = intersect(Nei_j_k_l, Matrix[m])
                    for n in [item for item in Nei_j_k_l_m if item > m]:
                        five_cliques.append([j, k, l, m, n])
                        Nei_j_k_l_m_n = intersect(Nei_j_k_l_m, Matrix[n])
                        for o in [item for item in Nei_j_k_l_m_n if item > n]:
                            six_cliques.append([j, k, l, m, n, o])
    num_4_cliques = [0 for item in edges]
    for clique in four_cliques:
        tunnel = tunnel_K5_geometric([vertex] + clique, Centroid)
        if tunnel:
            if vertex in tunnel:
                apex = setdiff(tunnel, [vertex])[0]
                base = setdiff(clique, [apex])
                for edge in combinations(base, 2):
                    num_4_cliques[ind(list(edge), edges)[0]] += 1
            else:
                diagonals = [tunnel, setdiff(clique, tunnel)]
                for edge in combinations(clique, 2):
                    edge = sorted(list(edge))
                    if isline(edge, diagonals):
                        continue
                    num_4_cliques[ind(edge, edges)[0]] += 1
                    if not isline(edge, deg_3_edges_in_cross):
                        deg_3_edges_in_cross.append(edge)
        else:
            T, apex = ispyramid(clique, vertex, Centroid)
            if T:
                base = setdiff(clique, [apex])
                for edge in combinations(base, 2):
                    num_4_cliques[ind(list(edge), edges)[0]] += 1
                continue
            T, diagonals = iscross(clique, vertex, Centroid)
            if T:
                for edge in combinations(clique, 2):
                    edge = sorted(list(edge))
                    if isline(edge, diagonals):
                        continue
                    num_4_cliques[ind(edge, edges)[0]] += 1
                    if not isline(edge, deg_3_edges_in_cross):
                        deg_3_edges_in_cross.append(edge)
    for clique in five_cliques:
        for edge in combinations(clique, 2):
            num_4_cliques[ind(list(edge), edges)[0]] -= 1
    for clique in six_cliques:
        for edge in combinations(clique, 2):
            num_4_cliques[ind(list(edge), edges)[0]] += 1

    open_edges = []
    for i in range(len(edges)):
        Ind = ind(edges[i], link)
        if len(Ind) < 2 + num_4_cliques[i]:
            open_edges.append(edges[i])
    return open_edges


##################################################################


def pseudoholes(vertex, Matrix, list_boundary, Centroid):
    """ returns the cycles in the link that are almost holes (having all edges open but one). """

    Matrix_reduced, triangles = reducedvertexlink(vertex, Matrix, list_boundary, Centroid, return_Matrix=True,
                                                  return_triangles=True)[1:]
    tets = [[vertex] + item for item in triangles]
    D, T, _ = isHamiltonian(vertex, Matrix_reduced, list_boundary, Centroid, predetermined_tetrahedra=tets)
    if D and T:
        return []
    open_edges = exposed_edges_in_link(vertex, Matrix, Centroid)
    Holes = []
    Graph = graph(open_edges)
    # print("open_edges:", open_edges)
    for com in components_edges(open_edges):
        # print(com)
        if len(com) < 4:
            continue
        endpoints = [item for item in com if len(ind(item, open_edges)) == 1]
        # print("endpoints:", endpoints)
        for j in combinations(endpoints, 2):
            if j[0] not in Matrix[j[1]]:
                continue
            for k in range(4, len(com) + 1):
                C = [path for path in dfs(Graph, j[0], j[1], [k - 1, k])]
                for path in C:
                    if j[0] not in path:
                        path = [j[0]] + path
                    if j[1] not in path:
                        path = path + [j[1]]
                    if not isline_exact(path, Holes)[0]:
                        Holes.append(path)
    if Holes:
        return Holes
    for com in components_edges(open_edges):
        if len(com) < 4:
            continue
        endpoints = unique(open_edges)
        for j in combinations(endpoints, 2):
            if j[0] not in Matrix[j[1]]:
                continue
            for k in range(4, len(com)):
                C = [path for path in dfs(Graph, j[0], j[1], k)]
                for path in C:
                    if j[0] not in path:
                        path = [j[0]] + path
                    if j[1] not in path:
                        path = path + [j[1]]
                    if not isline(path, Holes):
                        Holes.append(path)
    return Holes


##################################################################


def repairable_link_old(lst, Matrix, list_boundary, Centroid):
    if not isinstance(lst, int):
        return [repairable_link_old(item, Matrix, list_boundary, Centroid) for item in lst]
    vertex = lst
    if vertex in list_boundary:
        return False
    Holes, _, open_edges = holes_old(vertex, Matrix, list_boundary, Centroid, True, True)
    for edge in open_edges:
        included = False
        for Hole in Holes:
            if isline(edge, Hole):
                Ind1 = Hole.index(edge[0])
                Ind2 = Hole.index(edge[1])
                if Ind2 in [Next(Ind1, 1, len(Hole)), Previous(Ind1, 1, len(Hole))]:
                    included = True
                    break
        if not included:
            return False
    return True


##################################################################


def repairable_link(lst, Matrix, list_boundary, Centroid):
    if not isinstance(lst, int):
        return [repairable_link(item, Matrix, list_boundary, Centroid) for item in lst]
    vertex = lst
    # if istriangulatable(vertex, Matrix, list_boundary, Centroid)[0]:
    #     return True
    if vertex in list_boundary:
        return False
    Holes, _, open_edges = holes(vertex, Matrix, list_boundary, Centroid, True, True)
    for edge in open_edges:
        included = False
        for Hole in Holes:
            if isline(edge, Hole):
                Ind1 = [Ind for Ind in range(len(Hole)) if Hole[Ind] == edge[0]]
                Ind2 = [Ind for Ind in range(len(Hole)) if Hole[Ind] == edge[1]]
                for i in Ind1:
                    for j in Ind2:
                        if j in [Next(i, 1, len(Hole)), Previous(i, 1, len(Hole))]:
                            included = True
                            break
        if not included:
            return False
    return True


##################################################################


def vertexlink(vertex, Matrix):  # link of a vertex in the clique complex
    link = []
    Nei_i = Matrix[vertex]
    for j in Nei_i:
        Nei_i_j = intersect(Nei_i, Matrix[j])
        for k in [item for item in Nei_i_j if item > j]:
            Nei_i_j_k = intersect(Nei_i_j, Matrix[k])
            for l in [item for item in Nei_i_j_k if item > k]:
                link.append([j, k, l])
    return link


##################################################################


def edges_vertexlink(vertex, Matrix):  # edges of a link of a vertex in the clique complex
    edges = []
    vts = Matrix[vertex]
    for i in vts:
        for j in [item for item in intersect(Matrix[i], vts) if item > i]:
            edges.append([i, j])
    return edges


##################################################################


def reducedvertexlink(vertex, Matrix, list_boundary, Centroid, show_steps=False, return_Matrix=False,
                      return_triangles=False, with_constriction=False, show_intermediate_plots=False):
    Matrix_out = loads(dumps(Matrix))
    # empty = []
    if with_constriction:
        constricted = []
    else:
        constricted = constricted_triangles_in_link(vertex, Matrix_out, list_boundary, Centroid)
    triangles = []
    while True:
        out = True
        change = False
        link = [item for item in vertexlink(vertex, Matrix_out) if not isline(item, constricted)]
        # link = [item for item in link if not isline(item, empty)]
        vts = unique(link)
        if is_tetrahedron(vts, Matrix_out):
            if return_Matrix and return_triangles:
                return link, Matrix_out, triangles
            if return_Matrix:
                return link, Matrix_out
            return link
        T = False
        for i in vts:
            if T:
                continue
            if len(ind(i, link)) == 3:  # remove degree-3 vertices
                vertices = Neighbors_link(i, link)
                if len(vertices) == 3 and is_clique(vertices, Matrix_out):  # and not isline(vertices, empty):
                    T = True
                    Matrix_out = Remove(vertex, i, Matrix_out)
                    triangles.append(vertices)
                    # empty.append(vertices)
                    # print("0 Removed triangle:", vertices)
                    if show_steps:
                        print("0 Removed vertex:", i, " - neighbors:", vertices)
                        if show_intermediate_plots:
                            print(i, Neighbors_link(i, link), intersect(Matrix_out[i], Matrix_out[vertex]))
                            plot_vertex_link(vertex, Add(vertex, i, Matrix_out), list_boundary, Centroid, highlight=[i])
                    out = False
                    change = True
        if change:
            continue
        K4 = []
        for i in vts:
            Nei_i = intersect(vts, Matrix_out[i])
            for j in [item for item in Nei_i if item > i]:
                Nei_i_j = intersect(Nei_i, Matrix_out[j])
                for k in [item for item in Nei_i_j if item > j]:
                    Nei_i_j_k = intersect(Nei_i_j, Matrix_out[k])
                    for l in [item for item in Nei_i_j_k if item > k]:
                        K4.append([i, j, k, l])
        if not K4:
            if return_Matrix and return_triangles:
                return vertexlink(vertex, Matrix_out), Matrix_out, triangles
            if return_Matrix:
                return vertexlink(vertex, Matrix_out), Matrix_out
            return vertexlink(vertex, Matrix_out)
        T = False
        for i in range(len(K4)):
            if T:
                continue
            lst = K4[i]
            degrees = [[], []]
            for j in lst:
                for k in [item for item in lst if item > j]:
                    degrees[0].append([j, k])
                    degrees[1].append(len(ind([j, k], link)))
            Ind = [index for index in range(6) if degrees[1][index] == 2]
            if len(Ind) == 1:  # if only one edge of a 4-clique has degree 2, remove it
                Ind = Ind[0]
                T = True
                Matrix_out = Remove(degrees[0][Ind][0], degrees[0][Ind][1], Matrix_out)
                if show_steps:
                    print("1 Removed edge:", degrees[0][Ind], "from the clique", lst)
                    if show_intermediate_plots:
                        plot_vertex_link(vertex, Add(degrees[0][Ind][0], degrees[0][Ind][1], Matrix_out), list_boundary,
                                         Centroid, highlight=degrees[0][Ind])
                out = False
                change = True
        if change:
            continue
        T = False
        for i in range(len(K4)):
            if T:
                continue
            lst = K4[i]
            degrees = [[], []]
            for j in lst:
                for k in [item for item in lst if item > j]:
                    degrees[0].append([j, k])
                    degrees[1].append(len(ind([j, k], link)))
            Ind = [index for index in range(6) if degrees[1][index] == 2]
            if len(Ind) == 2:  # if two non-intersecting edges of a 4-clique have degree 2, remove any of them
                if intersect(degrees[0][Ind[0]], degrees[0][Ind[1]]):
                    continue
                T = True
                Matrix_out = Remove(degrees[0][Ind[0]][0], degrees[0][Ind[0]][1], Matrix_out)
                if show_steps:
                    print("2 Removed edge:", degrees[0][Ind[0]], "from the clique", lst)
                    if show_intermediate_plots:
                        plot_vertex_link(vertex, Add(degrees[0][Ind[0]][0], degrees[0][Ind[0]][1], Matrix_out),
                                         list_boundary, Centroid, highlight=degrees[0][Ind[0]])
                out = False
                change = True
        if change:
            continue
        T = False
        for i in range(len(K4)):
            if T:
                continue
            lst = K4[i]
            degrees = [[], []]
            for j in lst:
                for k in [item for item in lst if item > j]:
                    degrees[0].append([j, k])
                    degrees[1].append(len(ind([j, k], link)))
            Ind = [index for index in range(6) if degrees[1][index] == 2]
            if len(Ind) == 3:  # if three edges of a 4-clique have degree 2 and two of them are disjoint, remove any
                # of those two
                Ind_sub = []
                for j in range(2):
                    for k in range(j + 1, 3):
                        if not intersect(degrees[0][Ind[j]], degrees[0][Ind[k]]):
                            Ind_sub = [Ind[item] for item in [j, k]]
                if Ind_sub:
                    Ind = Ind_sub
                T = True
                Matrix_out = Remove(degrees[0][Ind[0]][0], degrees[0][Ind[0]][1], Matrix_out)
                if show_steps:
                    print("3 Removed edge:", degrees[0][Ind[0]], "from the clique", lst)
                    if show_intermediate_plots:
                        plot_vertex_link(vertex, Add(degrees[0][Ind[0]][0], degrees[0][Ind[0]][1], Matrix_out),
                                         list_boundary, Centroid, highlight=degrees[0][Ind[0]])
                out = False
                change = True
        if change:
            continue

        # T = False
        # for i in vts:
        #     if T:
        #         continue
        #     if len(Neighbors_link(i, link)) == 4:
        #         if isK4(Neighbors_link(i, link), Matrix_out):
        #             T = True
        #             Matrix_out = Remove(vertex, i, Matrix_out)
        #             if show_steps:
        #                 print("4 Removed vertex:", i)
        #             out = False
        #             change = True
        # if change:
        #     continue
        # T = False
        # for edge in combinations(unique(K4), 2):
        #     edge = list(edge)
        #     if not isline(edge, K4):
        #         continue
        #     if T:
        #         continue
        #     if sublist(intersect(mutual_neighbors(edge, Matrix_out), vts), unique(K4)):
        #         T = True
        #         Matrix_out = Remove(edge[0], edge[1], Matrix_out)
        #         if show_steps:
        #             print("5 Removed edge:", edge)
        #         out = False
        #         change = True
        # if change:
        #     continue

        # edges = []
        # for j in vts:
        #     Nei_j = intersect(vts, Neighbors_link(j, link))
        #     for k in [item for item in Nei_j if item > j]:
        #         if not isline([j, k], edges):
        #             edges.append([j, k])
        # for i in range(len(edges)):
        #     if len(ind(edges[i], link)) < 2:
        #         Matrix_out = Remove(edges[i][0], edges[i][1], Matrix_out)
        #         # print("6 Removed edge:", edges[i])
        #         out = False
        #         change = True
        if out:
            if return_Matrix and return_triangles:
                return vertexlink(vertex, Matrix_out), Matrix_out, triangles
            if return_Matrix:
                return vertexlink(vertex, Matrix_out), Matrix_out
            return vertexlink(vertex, Matrix_out)


##################################################################


def Neighbors_link(vertex, link):
    lst = []
    if isinstance(vertex, int):
        Ind = ind(vertex, link)
        for k in Ind:
            lst = union(lst, link[k])
        lst = setdiff(lst, [vertex])
        return lst
    elif isinstance(vertex, list) and len(vertex) == 1:
        return Neighbors_link(vertex[0], link)
    elif isinstance(vertex, list):
        lst = Neighbors_link(vertex[0], link)
        for i in vertex[1:]:
            lst = intersect(lst, Neighbors_link(i, link))
        return lst


##################################################################


def str_2_int_float(s):
    try:
        return int(s)
    except ValueError:
        return float(s)


##################################################################


def str_2_int_float_nan(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return -1


################################################################


def str_2_bool(s):
    if not s or s == "False":
        return False
    return True


##################################################################


def isK4(lst, Matrix):
    if len(lst) != 4:
        return False
    return lst[0] in Matrix[lst[1]] and lst[0] in Matrix[lst[2]] and lst[0] in Matrix[lst[3]] and lst[1] in Matrix[
        lst[2]] and lst[1] in Matrix[lst[3]] and lst[2] in Matrix[lst[3]]


##################################################################


def isline(lst, Matrix):
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


##################################################################


def isline_exact(lst, Matrix):
    # returns True if there exists an element 'item' in 'Matrix' such that 
    # 'set(item) = set(lst)' and 'len(item) = len(lst)'.
    if not Matrix:
        return False, []
    T = False
    Ind = ind(lst, Matrix)
    lst_remove = []
    if isinstance(Matrix[0], int):
        Matrix = [Matrix]
    if isinstance(lst, int):
        lst = [lst]
    for i in Ind:
        if len(Matrix[i]) != len(lst):
            lst_remove.append(i)
    for i in lst_remove:
        Ind.remove(i)
    if Ind:
        T = True
    return T, Ind


##################################################################


def isline_EXACT(lst, Matrix):
    # returns True if there is an element of 'Matrix' that has the same 
    # cyclic order as 'lst'.
    if not Matrix:
        return False
    Ind = isline_exact(lst, Matrix)[1]
    if not Ind:
        return False
    lst_remove = []
    if isinstance(Matrix[0], int):
        Matrix = [Matrix]
    if isinstance(lst, int):
        lst = [lst]
    for i in Ind:
        if not isequal_cyclic(Matrix[i], lst):
            lst_remove.append(i)
    for i in lst_remove:
        Ind.remove(i)
    if Ind:
        return True
    return False


##################################################################


def isline_EXACT_np(v, Matrix):
    """ returns True if there is an element of 'Matrix' that has the same
        cyclic order as 'v' where 'Matrix' is a numpy array. """

    if Matrix.size == 0:
        return False
    return np.any(np.all(v == Matrix, axis=1))


##################################################################


def isequal_cyclic(lst1, lst2):
    if lst1 == lst2:
        return True
    if len(lst1) != len(lst2) or (len(lst1) == 1 and len(lst2) == 1):
        return False
    if lst1[0] not in lst2:
        return False
    a = lst2.index(lst1[0])
    if lst1[1] == lst2[Next(a, 1, len(lst2))]:
        lst = []
        for s in range(len(lst2)):
            lst.append(lst2[Next(a, s, len(lst2))])
        return lst == lst1
    elif lst1[1] == lst2[Previous(a, 1, len(lst2))]:
        lst = []
        for s in range(len(lst2)):
            lst.append(lst2[Previous(a, s, len(lst2))])
        return lst == lst1
    else:
        return False


##################################################################


def ind(lst, Matrix):
    """ returns indices of sublists of Matrix (seen as an adjacency matrix) than contain lst. """

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


def graph(edges):
    for item in edges:
        if len(item) != 2:
            print("Error: the input should consist of edges!")
    vertices = []
    adjacencies = []
    for edge in edges:
        for i in range(2):
            if edge[i] not in vertices:
                vertices.append(edge[i])
                adjacencies.append([edge[Next(i, 1, 2)]])
            else:
                Ind = vertices.index(edge[i])
                adjacencies[Ind] = union(adjacencies[Ind], [edge[Next(i, 1, 2)]])
    Graph = {}
    for index in range(len(vertices)):
        Graph.update({vertices[index]: adjacencies[index]})
    return Graph


##################################################################


def Next(j, r=1, n=math.inf):
    return int((((j + r - 1) % n) + 1) % n)


##################################################################


def Previous(j, r=1, n=math.inf):
    return int((((j - 1 - r) % n) + 1) % n)


##################################################################


def edges_of(vertices, Matrix):
    edges = []
    for i in vertices:
        for j in [item for item in intersect(Matrix[i], vertices) if item > i]:
            edges.append([i, j])
    return edges


##################################################################


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


def sym_2(s):
    if s == 0:
        return sym_2_1
    elif s == 1:
        return sym_2_2


def lex_order(facets):
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


##################################################################


def cycles(edges, minimum_allowed=3, order_matters=False, within_time=False):
    """ Computes the cycles of the given edges within 5 seconds (if 'within_time' is set to 'True'). """

    Graph = graph(edges)
    C = list(simple_cycles(Graph))
    # print(len(C))
    # C_out = []
    # for item in C:
    #     if len(item) == 2:
    #         continue
    #     if isline_EXACT(item, C_out):
    #         continue
    #     C_out.append(item)
    C = [item for item in C if len(item) >= minimum_allowed]
    start_time = time()
    if order_matters:
        lst_remove = []
        for i in range(len(C)):
            if within_time:
                if time() - start_time > 5:
                    break
            if i in lst_remove:
                continue
            for j in setdiff(range(i + 1, len(C)), lst_remove):
                if len(C[i]) == len(C[j]) and set(C[j]) == set(C[i]) and isline_EXACT(C[i], C[j]):
                    lst_remove.append(j)
        return Remove_index(C, lst_remove)
    lst_remove = []
    for i in range(len(C)):
        if within_time:
            if time() - start_time > 5:
                break
        if i in lst_remove:
            continue
        for j in range(i + 1, len(C)):
            if len(C[i]) != len(C[j]):
                continue
            if set(C[j]) == set(C[i]):
                lst_remove.append(j)
    return Remove_index(C, lst_remove)


##################################################################


def dfs(Graph, start, end, lengths=None):
    """ Depth-first search to find paths of length within 'lengths' """

    if lengths is None: lengths = range(3, 100)
    if isinstance(lengths, int): lengths = [lengths]
    fringe = [(start, [])]
    while fringe:
        state, path = fringe.pop()
        if len(path) in lengths and state == end:
            yield path
            continue
        for next_state in Graph[state]:
            if next_state in path:
                continue
            fringe.append((next_state, path + [next_state]))


##################################################################


def dfsHamiltonian(Graph, show_steps, starting_node=None):
    start_time = time()
    nodes = list(Graph.keys())
    nodes = sorted(nodes, reverse=True)
    # print("nodes:", nodes)
    vertices = unique([node[0] for node in nodes])
    if starting_node is None:
        starting_node = nodes[-1]
    # print(len(vertices), " - ", vertices)
    fringe = [(starting_node, [])]
    while fringe:
        t = time() - start_time
        # print(t)
        if t > 20:
            print(" - timed out!")
            yield [], True
            # break
        state, path = fringe.pop()
        if path and path[0][0] != starting_node[0]:
            break
        # print(state, path, fringe)
        if len(path) == len(vertices):
            yield path, False
            continue
        for next_state in nodes:
            if next_state in path:
                continue
            if intersect(path, Graph[next_state]):
                continue
            if path and next_state < path[-1]:
                continue
            if path and vertices.index(next_state[0]) != vertices.index(path[-1][0]) + 1:
                continue
            fringe.append((next_state, path + [next_state]))
        # print(fringe)
        if show_steps:
            print(path)


##################################################################


def aux_graph(vertex, Matrix, list_boundary, Centroid, show_steps=False, predetermined_tetrahedra=None,
              multi_links=False):
    """ Returns a graph whose vertices are of the form '(v, cycle)' where v is a vertex on the link of 'vertex'
        and 'cycle' is a cycle in the link of the edge '[vertex, v]' in the clique complex. Two vertices
        '(v_1, cycle_1)' and '(v_2, cycle_2)' are adjacent iff 'cycle_1' and 'cycle_2' are not compatible
        in a Hamiltonian cycle of the link of 'vertex' in the clique complex. """

    if predetermined_tetrahedra is None:
        predetermined_tetrahedra = []

    K4 = []
    for i in [vertex] + Matrix[vertex]:
        Nei_i = Matrix[i]
        for j in [item for item in Nei_i if item > i]:
            Nei_i_j = intersect(Nei_i, Matrix[j])
            for k in [item for item in Nei_i_j if item > j]:
                Nei_i_j_k = intersect(Nei_i_j, Matrix[k])
                for l in [item for item in Nei_i_j_k if item > k]:
                    K4.append([i, j, k, l])
    Graph = {}
    Candidate_cycles = {}
    constricted = constricted_triangles_in_link(vertex, Matrix, list_boundary, Centroid)
    for i in range(len(Matrix[vertex])):
        u = Matrix[vertex][i]
        edges = edges_of(mutual_neighbors([vertex, u], Matrix), Matrix)
        # print(u, len(edges))
        C = cycles(edges, minimum_allowed=3, order_matters=True, within_time=True)
        # print(u, len(C))
        if len(C) > 1000:
            if show_steps:
                print("a vertex with too many cycles:", u)
            return {}, {}
        tets = [predetermined_tetrahedra[item] for item in ind([u, vertex], predetermined_tetrahedra)]
        must_contain = [setdiff(tet, [vertex, u]) for tet in tets]  # edges that must appear in the cycle around 'u'
        if not multi_links:
            for edge in edges:
                if isline(edge + [u], constricted):
                    continue
                Ind = ind(edge + [u], K4)
                if len(Ind) == 2 and len(
                        intersect(unique([K4[item] for item in Ind]), Matrix[vertex])) == 3 and not sublist(edge + [u],
                                                                                                            list_boundary):
                    if not isline(edge, must_contain):
                        must_contain.append(edge)
        # print(u, must_contain)
        must_contain = []
        Candidate_cycles[u] = [cycle for cycle in C if
                               all([isline(edge, edges_of_cycle(cycle)) for edge in must_contain])]
        # print(u, Candidate_cycles[u])
        if not Candidate_cycles[u]:
            if show_steps:
                print("a vertex with no cycles:", u)
            return {}, {}
        for s in range(len(Candidate_cycles[u])):
            Graph[u, s] = [(u, t) for t in range(len(Candidate_cycles[u])) if t != s]
            for j in range(i):
                v = Matrix[vertex][j]
                if v not in Matrix[u]:
                    continue
                for t in range(len(Candidate_cycles[v])):
                    cycle_1 = Candidate_cycles[u][s]
                    cycle_2 = Candidate_cycles[v][t]
                    cond_1 = (u in cycle_2 and v not in cycle_1) or (u not in cycle_2 and v in cycle_1)
                    cond_2 = u in cycle_2 and v in cycle_1
                    if cond_2:
                        ind_u, ind_v = cycle_2.index(u), cycle_1.index(v)
                        lst_1 = [cycle_2[Next(ind_u, 1, len(cycle_2))], cycle_2[Previous(ind_u, 1, len(cycle_2))]]
                        lst_2 = [cycle_1[Next(ind_v, 1, len(cycle_1))], cycle_1[Previous(ind_v, 1, len(cycle_1))]]
                        cond_2 = not isline(lst_1, lst_2)
                    if cond_1 or cond_2:
                        Graph[u, s].append((v, t))
                        Graph[v, t].append((u, s))
    return Graph, Candidate_cycles


##################################################################


def isHamiltonian(vertex, Matrix, list_boundary, Centroid, show_steps=False, predetermined_tetrahedra=None,
                  multi_links=False):
    """ checks whether a simplicial 2-complex is Hamiltonian. Best to use with reduced links.
        Returns 3 outputs:

        - decidability within the given time (20 seconds),
        - Hamiltonicity which is reliable only if the first output is 'True',
        - a Hamiltonian cycle (if both first and second outputs are 'True'). """

    if predetermined_tetrahedra is None:
        predetermined_tetrahedra = []

    Graph, Candidate_cycles = aux_graph(vertex, Matrix, list_boundary, Centroid, show_steps=show_steps,
                                        predetermined_tetrahedra=predetermined_tetrahedra, multi_links=multi_links)
    if not Graph:
        return True, False, []

    for Ham_cycle, timed_out in dfsHamiltonian(Graph, show_steps):
        if timed_out:
            return False, False, []
        nodes = list(Graph.keys())
        vertices = unique([node[0] for node in nodes])
        check_consistency = True
        surr_tetrahedra = [predetermined_tetrahedra[item] for item in ind(vertex, predetermined_tetrahedra)]
        for tet in surr_tetrahedra:
            tri = setdiff(tet, vertex)
            if not sublist(tri, vertices) or not is_clique(tri, Matrix):
                continue  # might have done multiple vertex removals
            Ind = vertices.index(tri[0])
            cycle = Candidate_cycles[tri[0]][Ham_cycle[Ind][1]]
            edges = edges_of_cycle(cycle)
            if not isline(tri[1:], edges):
                check_consistency = False
                break
        if check_consistency:
            sphere = []
            for pair in Ham_cycle:
                triangles = [[pair[0]] + edge for edge in edges_of_cycle(Candidate_cycles[pair[0]][pair[1]])]
                sphere += [tri for tri in triangles if not isline(tri, sphere)]
            return True, True, sphere
    return True, False, []


##################################################################


def all_Ham_cycles(vertex, Matrix, list_boundary, Centroid, predetermined_tetrahedra=None):
    """ Computes all Hamiltonian cycles of the link that are consistent with 'predetermined_tetrahedra'. """

    if predetermined_tetrahedra is None:
        predetermined_tetrahedra = []

    Graph, Candidate_cycles = aux_graph(vertex, Matrix, list_boundary, Centroid,
                                        predetermined_tetrahedra=predetermined_tetrahedra, multi_links=True)
    if not Graph:
        return []
    Ham_cycles = []
    for Ham_cycle, timed_out in dfsHamiltonian(Graph, show_steps=False):
        if timed_out:
            break
        nodes = list(Graph.keys())
        vertices = unique([node[0] for node in nodes])
        check_consistency = True
        surr_tetrahedra = [predetermined_tetrahedra[item] for item in ind(vertex, predetermined_tetrahedra)]
        for tet in surr_tetrahedra:
            tri = setdiff(tet, vertex)
            if not sublist(tri, vertices) or not is_clique(tri, Matrix):
                continue  # might have done multiple vertex removals
            Ind = vertices.index(tri[0])
            cycle = Candidate_cycles[tri[0]][Ham_cycle[Ind][1]]
            edges = edges_of_cycle(cycle)
            if not isline(tri[1:], edges):
                check_consistency = False
                break
        if check_consistency:
            sphere = []
            for pair in Ham_cycle:
                triangles = [[pair[0]] + edge for edge in edges_of_cycle(Candidate_cycles[pair[0]][pair[1]])]
                sphere += [tri for tri in triangles if not isline(tri, sphere)]
            Ham_cycles.append([[vertex] + tri for tri in sphere])
    return Ham_cycles


##################################################################


def reconstruct_Hamiltonian(lst_vertices, Matrix, list_boundary, Centroid, triangulation_current, start_time, K5,
                            iteration_index=None, initial_length=None, must_contain=None, level_reached=None,
                            already_computed=None):
    """ Tries to find compatible Hamiltonian cycles in the links of the given vertices 'lst_vertices'. """

    if iteration_index is None:
        iteration_index = 1
    if initial_length is None:
        initial_length = len(lst_vertices)
    if must_contain is None:
        must_contain = []
    if level_reached is None:
        level_reached = 0
    if already_computed is None:
        already_computed = {v: [[], []] for v in lst_vertices}

    if time() - start_time > 5400:  # 1,5 hours
        return triangulation_current, False, False, max(iteration_index, level_reached), already_computed
    vertex = lst_vertices[0]
    if vertex == 900:
        save_list_of_lists(triangulation_current, "triangulation_current.txt")
        plot_vertex_link(vertex, Matrix, list_boundary, Centroid, Certified=triangulation_current)
    current_star = [triangulation_current[item] for item in ind(vertex, triangulation_current)]
    computed = False
    for i in range(len(already_computed[vertex][0])):
        star = already_computed[vertex][0][i]
        if len(current_star) != len(star):
            continue
        if not all([isline(item, star) for item in current_star]):
            continue
        Ind = i
        computed = True
        break
    if computed:
        cycles = already_computed[vertex][1][Ind]
    else:
        cycles = all_Ham_cycles(vertex, Matrix, list_boundary, Centroid,
                                predetermined_tetrahedra=triangulation_current)
        already_computed[vertex][0].append(current_star)
        already_computed[vertex][1].append(cycles)
    print("{}  |  {}".format(vertex, len(cycles)))
    # print("\r{}  |  {}".format(vertex, len(cycles)), end="")
    # print("\rProgress: |{}{}|".format("–" * iteration_index, " " * (initial_length - iteration_index)), end="")
    if cycles:
        cycle = cycles[0]
        for tet in cycle:
            if mutual_neighbors(tet, Matrix) and iteration_index != 1:
                continue
            if all([isline(tet, c) for c in cycles[1:]]):
                if not isline(tet, must_contain):
                    must_contain.append(tet)
                if not isline(tet, triangulation_current):
                    triangulation_current.append(tet)
    elif iteration_index > level_reached and len(ind(vertex, K5)) < 15:
        # print("iteration_index:", iteration_index)
        # print("level_reached:", level_reached)
        R, T, _ = isHamiltonian(vertex, Matrix, list_boundary, Centroid, predetermined_tetrahedra=must_contain)
        if R and not T:  # a compatibility issue
            print("a compatibility issue at the vertex {}".format(vertex))
            return triangulation_current, False, True, max(iteration_index, level_reached), already_computed
    for cycle in cycles:
        triangulation_neu = loads(dumps(triangulation_current))
        triangulation_neu += [tet for tet in cycle if not isline(tet, triangulation_neu)]
        # plot_vertex_link(lst_vertices[1], Matrix, list_boundary, Centroid, Certified=triangulation_current)
        if not proper_vertex_in_triangulation(vertex, triangulation_neu, Matrix[vertex]):
            continue
        if len(lst_vertices) == 1:
            return triangulation_neu, True, False, max(iteration_index, level_reached), already_computed
        else:
            triangulation_out, out, break_out, level_reached, already_computed = reconstruct_Hamiltonian(
                lst_vertices[1:], Matrix, list_boundary, Centroid, triangulation_neu, start_time, K5,
                iteration_index + 1, initial_length, must_contain, level_reached, already_computed)
            if break_out or out:
                return triangulation_out, out, break_out, max(iteration_index, level_reached), already_computed
    return triangulation_current, False, False, max(iteration_index, level_reached), already_computed


##################################################################


def cycles_deg(edges, k):
    """ returns cycles of length 'k' in the graph on 'edges'. """

    Graph = graph(edges)
    if k == 4:
        C = []
        for v in Graph:
            for u_1 in [item for item in Graph[v] if item > v]:
                for u_2 in [item for item in Graph[v] if item > max(v, u_1)]:
                    lst = intersect(Graph[u_1], Graph[u_2])
                    for u_3 in [item for item in lst if item > v]:
                        C.append([v, u_1, u_3, u_2])
        return C
    C = [path for node in Graph for path in dfs(Graph, node, node, k)]
    lst_remove = []
    for i in range(len(C)):
        if i in lst_remove:
            continue
        for j in range(i + 1, len(C)):
            if set(C[j]) == set(C[i]):
                lst_remove.append(j)
    return Remove_index(C, lst_remove)


##################################################################


def edge_in_cycle(edge, cycle):
    """ returns 'True' if the edge lies in the cycle. """
    if not sublist(edge, cycle):
        return False
    Ind_1, Ind_2 = cycle.index(edge[0]), cycle.index(edge[1])
    if Ind_2 not in [Previous(Ind_1, 1, len(cycle)), Next(Ind_1, 1, len(cycle))]:
        return False
    return True


##################################################################


def edges_of_cycle(cycle):
    """ returns the edges of the cycle. """
    return [[cycle[item], cycle[Next(item, 1, len(cycle))]] for item in range(len(cycle))]


##################################################################


def ispath(edges):
    """ checks whether the edges form a path and returns the path in that case. """
    endpoints = [item for item in unique(edges) if len(ind(item, edges)) == 1]
    if len(endpoints) != 2:
        return []
    Graph = graph(edges)
    path = [path for path in dfs(Graph, endpoints[0], endpoints[1], len(edges))][0]
    if path[0] != endpoints[0]:
        path = [endpoints[0]] + path
    if path[-1] != endpoints[1]:
        path += [endpoints[1]]
    return path


##################################################################


def ispyramid(tuple_4, vertex, Centroid):
    T = False
    apex = math.nan
    count = 0
    for j in tuple_4:
        [a, b, c] = barycentric_intersection_in_triangle(setdiff(tuple_4, [j]), [vertex, j], Centroid)
        if a > 0 and b > 0 and c > 0:
            apex = j
            count += 1
    if count == 1:
        return True, apex
    return T, math.nan


##################################################################


def iscross(tuple_4, vertex, Centroid):
    T = False
    diagonals = math.nan
    for j in combinations(tuple_4, 2):
        j = list(j)
        rest = setdiff(tuple_4, j)
        [a, b] = barycentric_intersection_in_line(union([vertex], j), Centroid[rest[0]], Centroid[rest[1]], Centroid)
        if a > 0 and b > 0:
            [a, b, c] = barycentric_intersection_in_triangle(union([vertex], j), rest, Centroid)
            if a > 0 and b > 0 and c > 0:
                T = True
                diagonals = [rest, j]  # the first diagonal is the closest one to vertex.
    return T, diagonals


##################################################################


def is_4_simplex(tuple_4, vertex, Centroid):
    count = 0
    for j in tuple_4:
        [a, b, c] = barycentric_intersection_in_triangle(setdiff(tuple_4, [j]), [vertex, j], Centroid)
        if a > 0 and b > 0 and c > 0:
            count += 1
    if count == 4:
        return True
    return False


##################################################################


def barycentric_intersection_in_triangle(X, Y, Centroid):
    # computes the barycentric coordinates of the intersection point of line
    # passing through the points Y(1), Y(2) and the triangle X(1), X(2), X(3).
    # The coordinates are w.r.t. the triangle.
    A, B, C, P1, P2 = Centroid[X[0]], Centroid[X[1]], Centroid[X[2]], Centroid[Y[0]], Centroid[Y[1]]
    n = cross(subtract(B, A), subtract(C, A))
    u = subtract(P2, P1)
    w = subtract(P1, A)
    D = dot(n, u)
    N = - dot(n, w)
    sI = N / D
    I = add(P1, multiply_scalar(sI, u))
    n1 = cross(subtract(C, B), subtract(I, B))
    n2 = cross(subtract(C, A), subtract(C, I))
    n3 = cross(subtract(B, A), subtract(I, A))
    a = dot(n, n1) / (norm(n) ** 2)
    b = dot(n, n2) / (norm(n) ** 2)
    c = dot(n, n3) / (norm(n) ** 2)
    return a, b, c


##################################################################


def barycentric_intersection_in_line(X, P1, P2, Centroid):
    # computes the barycentric coordinates of the intersection point of line
    # passing through the points P1, P2 and the triangle X(1), X(2), X(3).
    # The coordinates are w.r.t. the line.
    if isinstance(P1, int):
        P1 = Centroid[P1]
    if isinstance(P2, int):
        P2 = Centroid[P2]
    A, B, C = Centroid[X[0]], Centroid[X[1]], Centroid[X[2]]
    n = cross(subtract(B, A), subtract(C, A))
    u = subtract(P2, P1)
    w = subtract(P1, A)
    D = dot(n, u)
    N = - dot(n, w)
    sI = N / D
    I = add(P1, multiply_scalar(sI, u))
    a = dot(subtract(I, P2), subtract(P1, P2)) / dot(subtract(P1, P2), subtract(P1, P2))
    b = dot(subtract(I, P1), subtract(P2, P1)) / dot(subtract(P2, P1), subtract(P2, P1))
    return a, b


##################################################################


def barycentric_in_triangle(X, P, Centroid):
    # computes the barycentric coordinates of point P w.r.t. triangle X.
    if isinstance(P, int):
        P = Centroid[P]
    A, B, C = Centroid[X[0]], Centroid[X[1]], Centroid[X[2]]
    n = cross(subtract(B, A), subtract(C, A))
    n1 = cross(subtract(C, B), subtract(P, B))
    n2 = cross(subtract(C, A), subtract(C, P))
    n3 = cross(subtract(B, A), subtract(P, A))
    a = dot(n, n1) / (norm(n) ** 2)
    b = dot(n, n2) / (norm(n) ** 2)
    c = dot(n, n3) / (norm(n) ** 2)
    return a, b, c


##################################################################


def intersect_triangle_line(X, P1, P2, Centroid):
    # find coordinates of the intersection point of segment P1P2 and triangle X = ABC.
    if isinstance(P1, int):
        P1 = Centroid[P1]
    if isinstance(P2, int):
        P2 = Centroid[P2]
    A, B, C = Centroid[X[0]], Centroid[X[1]], Centroid[X[2]]
    n = cross(subtract(B, A), subtract(C, A))
    u = subtract(P2, P1)
    w = subtract(P1, A)
    D = dot(n, u)
    N = - dot(n, w)
    sI = N / D
    return add(P1, multiply_scalar(sI, u))


##################################################################


def cross(u, v):
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]


##################################################################


def dot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


##################################################################


def subtract(u, v, Centroid=None):
    if isinstance(u, int):
        u = Centroid[u]
    if isinstance(v, int):
        v = Centroid[v]
    output = []
    for i in range(len(u)):
        output.append(u[i] - v[i])
    return output


##################################################################


def add(u, v):
    output = []
    for i in range(len(u)):
        output.append(u[i] + v[i])
    return output


##################################################################


def multiply_vectors(u, v):
    output = []
    for i in range(len(u)):
        output.append(u[i] * v[i])
    return output


##################################################################


def multiply_scalar(a, u):
    output = []
    for i in range(len(u)):
        output.append(a * u[i])
    return output


##################################################################


def divied_scalar(a, u):
    output = []
    for i in range(len(u)):
        output.append(u[i] / a)
    return output


##################################################################


def divide(u, v):
    output = []
    for i in range(len(u)):
        output.append(u[i] / v[i])
    return output


##################################################################


def norm(u):
    return (u[0] ** 2 + u[1] ** 2 + u[2] ** 2) ** (1 / 2)


##################################################################


def normalized(u):
    return divied_scalar(norm(u), u)


##################################################################


def distance_Hausdorff(pt1, pt2, Centroid, Diameters):
    o1 = Centroid[pt1]
    o2 = Centroid[pt2]
    r1 = Diameters[pt1]
    r2 = Diameters[pt2]
    x1 = add(o1, multiply_scalar(r1, normalized(subtract(o1, o2))))
    x2 = add(o1, multiply_scalar(r1, normalized(subtract(o2, o1))))
    y1 = add(o2, multiply_scalar(r2, normalized(subtract(o1, o2))))
    y2 = add(o2, multiply_scalar(r2, normalized(subtract(o2, o1))))
    return max(norm(subtract(x1, y1)), norm(subtract(x2, y2)))

##################################################################
