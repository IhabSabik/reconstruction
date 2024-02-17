import sys

from skimage import io
from Functions import *

start_time = time()
print(start_time)

# suffix = '_meanDiameter_ASM_water.txt'
suffix = '_meanDiameter_tristan_4_8.txt'
# suffix = '_waterbord.txt'
# path = suffix[14:-3] + 'tif'
path = 'labeledpores' + suffix[-16:-3] + 'tif'
# path = suffix[1:-3] + 'tif'
voxels = io.imread(path)
# dist_path = 'computed_distances' + suffix[13:]
dist_path = 'computed_distances' + suffix[-16:]
# dist_path = 'computed_distances' + suffix
f = open(dist_path, 'a')
f.close()
# log_file = open('log_file' + suffix[13:], 'a+')
log_file = open('log_file' + suffix[-16:], 'a+')
# log_file = open('log_file' + suffix, 'a+')
# log_file_detailed = open('log_file_detailed' + suffix[13:], 'a+')
log_file_detailed = open('log_file_detailed' + suffix[-16:], 'a+')
# log_file_detailed = open('log_file_detailed' + suffix, 'a+')

with open('Centroid' + suffix, 'r') as f:
    Centroid = [[float(num) for num in line.split(',')] for line in f]
numbPoints = len(Centroid)

coordinates = [[] for item in range(numbPoints)]

# with open('list_boundary_alpha_shape' + suffix[13:], 'r') as f:
#     list_boundary = [int(elem) - 1 for elem in f.read().split(',') if elem]
with open('list_boundary_alpha_shape' + suffix[-8:], 'r') as f:
    list_boundary = [int(elem) - 1 for elem in f.read().split(',') if elem]
# with open('list_boundary_alpha_shape' + suffix, 'r') as f:
#     list_boundary = [int(elem) - 1 for elem in f.read().split(',') if elem]

with open('list_boundary' + suffix, 'r') as f:
    list_boundary_vox = [int(elem) - 1 for elem in f.read().split(',') if elem]

with open('Diameters' + suffix, 'r') as f:
    Radii = [float(elem) / 2 for elem in f.read().split('\n') if elem]
ignore = [vertex for vertex in range(numbPoints) if Radii[vertex] < (3 * 450 / (4 * math.pi)) ** (1/3)]

list_interior = setdiff(setdiff(range(numbPoints), ignore), list_boundary)

center = sum(np.array(Centroid)) / len(Centroid)
dists = [np.linalg.norm(temp - center) for temp in Centroid]
indices_sorted = list(np.argsort(np.array(dists)))
lst_1 = [int(item) for item in indices_sorted if item in list_interior]
lst_2 = [int(item) for item in indices_sorted if item in list_boundary]
order = lst_1 + lst_2

max_num_open_edges = 25
index_inf = math.inf
threshold = 1

Register_Added = []
Register_Removed = []

Matrix_out = [[] for item in Centroid]
avg_deg = 0

log_file_detailed.write("Compute initial adjacencies...")
log_file_detailed.flush()
log_file.write("Compute initial adjacencies...")
log_file.flush()
while avg_deg < 6:
    threshold += 1
    for vertex in setdiff(range(numbPoints), ignore):
        for vertex_sub in setdiff(range(vertex + 1, numbPoints), ignore):
            if vertex_sub in Matrix_out[vertex]:
                continue
            print("\rCompute initial adjacencies (tolerance = {}):".format(threshold), vertex, "out of", numbPoints - 2,
                  " | ", vertex_sub, "out of", numbPoints - 1, end="")
            if norm(np.array(Centroid[vertex]) - np.array((Centroid[vertex_sub]))) > Radii[vertex] + Radii[vertex_sub] + 10 * threshold:
                continue
            dist = distance_vox(vertex, vertex_sub, Centroid, coordinates, voxels, dist_path)
            if dist < threshold:
                Matrix_out = Add(vertex, vertex_sub, Matrix_out)
    avg_deg = 0
    for vertex in list_interior:
        avg_deg += len(Matrix_out[vertex])
    avg_deg /= len(list_interior)
    print("\nThreshold = {}  |  average face degree = {}".format(threshold, round(avg_deg, 3)))
    log_file.write("\nThreshold = {}  |  average face degree = {}".format(threshold, round(avg_deg, 3)))
    log_file.flush()
    log_file_detailed.write("\nThreshold = {}  |  average face degree = {}".format(threshold, round(avg_deg, 3)))
    log_file_detailed.flush()
print("")

save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_-1.txt")

# with open("Matrix_out" + suffix[:-4] + "_step_-1.txt", 'r') as f:
#     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# for i in range(len(Matrix_out)):
#     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed

log_file_detailed.write("\n\nCheck for isolated vertices...")
log_file_detailed.flush()
log_file.write("\n\nCheck for isolated vertices...")
log_file.flush()
isolated = [item for item in setdiff(range(numbPoints), ignore) if not Matrix_out[item]]
print("Initially isolated:", len(isolated))
log_file_detailed.write("\nInitially isolated: {}".format(len(isolated)))
log_file_detailed.flush()
log_file.write("\nInitially isolated: {}".format(len(isolated)))
log_file.flush()
count = 0
for vertex in isolated:
    count += 1
    print("\rAdd connections to isolated vertices: processing {} out of {}".format(count, len(isolated)), end="")
    distances, vts = [], []
    for vertex_sub in setdiff(range(numbPoints), ignore):
        if vertex_sub == vertex:
            continue
        if norm(np.array(Centroid[vertex]) - np.array((Centroid[vertex_sub]))) > Radii[vertex] + Radii[vertex_sub] + 15 * threshold:
            continue
        dist = distance_vox(vertex, vertex_sub, Centroid, coordinates, voxels, dist_path)
        distances.append(dist)
        vts.append(vertex_sub)
    m = min(distances)
    vertex_sub = vts[distances.index(m)]
    Matrix_out = Add(vertex, vertex_sub, Matrix_out)
print("")

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
avg_deg /= len(list_interior)
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg, 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg, 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_0.txt")

# with open("Matrix_out" + suffix[:-4] + "_step_0.txt", 'r') as f:
#     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# for i in range(len(Matrix_out)):
#     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed

# Add connections to endpoints of edges with low degree
start_time = time()
print("\nAdd connections to endpoints of edges with low degree...")
log_file_detailed.write("\n\nAdd connections to endpoints of edges with low degree...")
log_file_detailed.flush()
log_file.write("\n\nAdd connections to endpoints of edges with low degree...")
log_file.flush()
for u in range(numbPoints):
    for v in [item for item in Matrix_out[u] if item > u]:
        if sublist([u, v], list_boundary):
            continue
        while len(mutual_neighbors([u, v], Matrix_out)) < 3:
            pairs = []
            for w in setdiff(Matrix_out[u], [v] + Matrix_out[v]):
                pairs.append([v, w])
            for w in setdiff(Matrix_out[v], [u] + Matrix_out[u]):
                pairs.append([u, w])
            if not pairs:
                break
            distances = [distance_vox(item[0], item[1], Centroid, coordinates, voxels, dist_path) for item in pairs]
            Ind = distances.index(min(distances))
            if distances[Ind] > 2 * threshold:
                break
            Register_Added.append(pairs[Ind])
            Matrix_out = Add(pairs[Ind][0], pairs[Ind][1], Matrix_out)
            print("Added:", pairs[Ind], distances[Ind])
            log_file_detailed.write("\nAdded: {} with distance {}".format(pairs[Ind], distances[Ind]))
            log_file_detailed.flush()

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_1.txt")
save_list_of_lists(Register_Added, "Register_Added" + suffix[:-4] + "_step_1.txt")

start_time = time()
print("\n\nCompute the square holes in vertex links...")
log_file_detailed.write("\n\nCompute the square holes in vertex links...")
log_file_detailed.flush()
log_file.write("\n\nCompute the square holes in vertex links...")
log_file.flush()
computed_holes = [[] for item in range(numbPoints)]
deg_3_edges_in_cross = [[] for item in range(numbPoints)]
computed_open_edges = [[] for item in range(numbPoints)]
for vertex in range(numbPoints):
    print("\rCompute the square holes in vertex links:", vertex, "out of", numbPoints - 1,
            end="")
    Holes, edges_in_cross, open_edges = holes(vertex, Matrix_out, list_boundary, Centroid, True, holes_degree=4)
    computed_holes[vertex] = Holes
    deg_3_edges_in_cross[vertex] = edges_in_cross
    computed_open_edges[vertex] = open_edges
print("")

log_file_detailed.write("\n\nClose the square holes in vertex links...")
log_file_detailed.flush()
log_file.write("\n\nClose the square holes in vertex links...")
log_file.flush()
iteration = 0
while True:
    out = True
    iteration += 1
    index = -1
    while index < len(order) - 1:
        index += 1
        # print("\r", "Close the square holes:", index, "out of", numbPoints - 1, "iteration", iteration, end="")
        vertex = order[index]
        if len(computed_open_edges[vertex]) < max_num_open_edges:  # only close squared holes in complex links
            continue
        Holes = loads(dumps(computed_holes[vertex]))
        for i in range(len(Holes)):
            if len(Holes[i]) != 4:
                continue
            priority = []
            for j in range(len(Holes[i])):
                if isline([Holes[i][item] for item in [j, Next(j, 1, len(Holes[i]))]],
                            deg_3_edges_in_cross[vertex]) and isline(
                    [Holes[i][item] for item in [j, Previous(j, 1, len(Holes[i]))]], deg_3_edges_in_cross[vertex]):
                    priority.append(Holes[i][j])
            for j in range(len(Holes[i]) - 2):  # extract priority vertices
                if Holes[i][j] == Holes[i][j + 2]:
                    priority = union(priority, Holes[i][j + 1])
            indices = []
            for j in combinations(range(len(Holes[i])), 2):
                j = list(j)
                edge = [Holes[i][item] for item in j]
                if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                    continue
                if priority and not intersect(edge, priority) and len(Holes[i]) != 4:
                    continue
                indices.append(j)
            if not indices:
                continue

            distances, outer_vertices, other_edges = [], [], []
            for j in indices:
                edge = [Holes[i][j[0]], Holes[i][j[1]]]
                distances.append(distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path))
                for k in mutual_neighbors(edge, Matrix_out):
                    if k not in [vertex] + Matrix_out[vertex]:
                        outer_vertices = union(outer_vertices, k)
                    else:
                        for l in Holes[i]:
                            if k not in Matrix_out[l] and not isline([k, l], other_edges) and not isline([k, l],
                                                                                                         Holes[i]):
                                other_edges.append([k, l])
            for v in outer_vertices:
                distances.append(distance_vox(vertex, v, Centroid, coordinates, voxels, dist_path))
            for edge in other_edges:
                distances.append(distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path))
            Ind = distances.index(min(distances))
            if Ind < len(indices):
                edge = [Holes[i][item] for item in indices[Ind]]
            elif Ind < len(indices) + len(outer_vertices):
                edge = [vertex, outer_vertices[Ind - len(indices)]]
            else:
                edge = other_edges[Ind - len(indices) - len(outer_vertices)]
            if distances[Ind] > 3 * threshold:
                continue
            Register_Added.append(edge)
            Matrix_out = Add(edge[0], edge[1], Matrix_out)
            print("Added:", vertex, Holes[i], edge, distances[Ind])
            log_file_detailed.write("\nAdded: linkvertex {}, hole {}, edge {} with distance {}".format(vertex, Holes[i], edge, distances[Ind]))
            log_file_detailed.flush()
            out = False
            All = mutual_neighbors(edge, Matrix_out)
            for vertex_sub in All:
                computed_holes[vertex_sub], deg_3_edges_in_cross[vertex_sub], computed_open_edges[vertex_sub] = holes(vertex_sub, Matrix_out,
                                                                                        list_boundary, Centroid,
                                                                                        True, holes_degree=4)
            index -= 1
            break
    if out:
        break
print("")

save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_2.txt")
save_list_of_lists(Register_Added, "Register_Added" + suffix[:-4] + "_step_2.txt")

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

# with open("Matrix_out" + suffix[:-4] + "_step_2.txt", 'r') as f:
#     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# for i in range(len(Matrix_out)):
#     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed
#
# with open("Register_Added" + suffix[:-4] + "_step_2.txt", 'r') as f:
#     Register_Added = [[int(num) for num in line.split(',')] for line in f]

start_time = time()
print("\n\nCompute candidate edges to be added to close the holes...")
log_file_detailed.write("\n\nCompute candidate edges to be added to close the holes...")
log_file_detailed.flush()
log_file.write("\n\nCompute candidate edges to be added to close the holes...")
log_file.flush()
candidates_edges = [[] for item in range(numbPoints)]
candidates_distances = [[] for item in range(numbPoints)]
lst_compute = range(numbPoints)
while lst_compute:
    vertex = lst_compute[0]
    candidates_edges[vertex] = []
    candidates_distances[vertex] = []
    print("\rCompute the candidates:", vertex, "out of", numbPoints - 1, end="")
    if vertex in ignore:
        lst_compute = lst_compute[1:]
        continue
    Holes, edges_in_cross = [], []
    open_edges = open_edges_in_link(vertex, Matrix_out, Centroid)
    if len(open_edges) <= 75:
        Holes, edges_in_cross, _ = holes(vertex, Matrix_out, list_boundary, Centroid, True)
    elif vertex not in list_boundary:
        Holes = [unique(open_edges)]
    outer_vertices, other_edges = [], []
    for i in range(len(Holes)):
        priority = []
        for j in range(len(Holes[i])):
            if edges_in_cross and isline([Holes[i][item] for item in [j, Next(j, 1, len(Holes[i]))]],
                      edges_in_cross) and isline(
                [Holes[i][item] for item in [j, Previous(j, 1, len(Holes[i]))]], edges_in_cross):
                priority.append(Holes[i][j])
        for j in range(len(Holes[i]) - 2):  # extract priority vertices
            if Holes[i][j] == Holes[i][j + 2]:
                priority = union(priority, Holes[i][j + 1])
        indices = []
        for j in combinations(range(len(Holes[i])), 2):
            j = list(j)
            edge = [Holes[i][item] for item in j]
            if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                continue
            if priority and not intersect(edge, priority) and len(Holes[i]) != 4:
                continue
            indices.append(j)
        for j in range(len(indices)):
            edge = [Holes[i][indices[j][0]], Holes[i][indices[j][1]]]
            if all([v in list_boundary_vox for v in Holes[i]]):
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                if dist > 2 * threshold:
                    continue
            if not isline(edge, candidates_edges[vertex]):
                candidates_edges[vertex].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                candidates_distances[vertex].append(dist)
                for k in mutual_neighbors(edge, Matrix_out):
                    if k not in [vertex] + Matrix_out[vertex]:
                        outer_vertices = union(outer_vertices, k)
                    else:
                        for l in Holes[i]:
                            if k not in Matrix_out[l] and not isline([k, l], other_edges) and not isline([k, l], Holes[i]):
                                dist = distance_vox(k, l, Centroid, coordinates, voxels, dist_path)
                                if dist > 1.5 * threshold:
                                    continue
                                other_edges.append([k, l])
    for v in outer_vertices:
        edge = [vertex, v]
        if not isline(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
            candidates_distances[vertex].append(dist)
    for edge in other_edges:
        if not isline(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
            candidates_distances[vertex].append(dist)
    # lst_remove = [index for index in range(len(candidates_edges[vertex])) if
    #               isline(candidates_edges[vertex][index], list_boundary_vox) and candidates_distances[vertex][
    #                   index] > 2 * threshold]
    # candidates_edges[vertex] = Remove_index(candidates_edges[vertex], lst_remove)
    # candidates_distances[vertex] = Remove_index(candidates_distances[vertex], lst_remove)
    lst_compute = lst_compute[1:]
print("")

save_list_of_lists(candidates_edges, "candidates_edges" + suffix[:-4] + "_step_3.txt")
save_list_of_lists(candidates_distances, "candidates_distances" + suffix[:-4] + "_step_3.txt")

log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

# with open("Matrix_out" + suffix[:-4] + "_step_2.txt", 'r') as f:
#     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# for i in range(len(Matrix_out)):
#     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed
#
# with open("Register_Added" + suffix[:-4] + "_step_2.txt", 'r') as f:
#     Register_Added = [[int(num) for num in line.split(',')] for line in f]
#
# with open("candidates_distances" + suffix[:-4] + "_step_3.txt", 'r') as f:
#     candidates_distances = [[str_2_int_float(num) for num in line.split(',')] for line in f]
# for i in range(len(candidates_distances)):
#     candidates_distances[i] = [x for x in candidates_distances[i] if x != -1]  # -1's get removed
#
# candidates_edges = read_list_of_lists_of_lists("candidates_edges" + suffix[:-4] + "_step_3.txt")

start_time = time()
print("\n\nAdd short candidate edges to close the holes...")
log_file_detailed.write("\n\nAdd short candidate edges to close the holes...")
log_file_detailed.flush()
log_file.write("\n\nAdd short candidate edges to close the holes...")
log_file.flush()
count = 0
while True:
    minima = [min(candidates_distances[vertex], default=math.inf) for vertex in range(numbPoints)]
    if min(minima) > 6 * threshold:  # == math.inf:
        break
    vertex = minima.index(min(minima))
    Ind = candidates_distances[vertex].index(min(candidates_distances[vertex]))
    edge = candidates_edges[vertex][Ind]
    Register_Added.append(edge)
    Matrix_out = Add(edge[0], edge[1], Matrix_out)
    print("Added:", vertex, edge, candidates_distances[vertex][Ind])
    log_file_detailed.write("\nAdded: linkvertex {}, edge {} with distance {}".format(vertex, edge, candidates_distances[vertex][Ind]))
    log_file_detailed.flush()
    All = union(edge, mutual_neighbors(edge, Matrix_out))
    for v_sub in All:
        candidates_edges[v_sub] = []
        candidates_distances[v_sub] = []
        Holes, edges_in_cross = [], []
        open_edges = open_edges_in_link(v_sub, Matrix_out, Centroid)
        if len(open_edges) <= 75:
            Holes, edges_in_cross, _ = holes(v_sub, Matrix_out, list_boundary, Centroid, True)
        elif v_sub not in list_boundary:
            Holes = [unique(open_edges)]
        outer_vertices, other_edges = [], []
        for i in range(len(Holes)):
            priority = []
            for j in range(len(Holes[i])):
                if edges_in_cross and isline([Holes[i][item] for item in [j, Next(j, 1, len(Holes[i]))]],
                        edges_in_cross) and isline(
                    [Holes[i][item] for item in [j, Previous(j, 1, len(Holes[i]))]], edges_in_cross):
                    priority.append(Holes[i][j])
            for j in range(len(Holes[i]) - 2):  # extract priority vertices
                if Holes[i][j] == Holes[i][j + 2]:
                    priority = union(priority, Holes[i][j + 1])
            indices = []
            for j in combinations(range(len(Holes[i])), 2):
                j = list(j)
                edge = [Holes[i][item] for item in j]
                if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                    continue
                if priority and not intersect(edge, priority) and len(Holes[i]) != 4:
                    continue
                indices.append(j)
            for j in range(len(indices)):
                edge = [Holes[i][indices[j][0]], Holes[i][indices[j][1]]]
                if all([v in list_boundary_vox for v in Holes[i]]):
                    dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                    if dist > 2 * threshold:
                        continue
                if not isline(edge, candidates_edges[v_sub]):
                    candidates_edges[v_sub].append(edge)
                    dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                    candidates_distances[v_sub].append(dist)
                    for k in mutual_neighbors(edge, Matrix_out):
                        if k not in [v_sub] + Matrix_out[v_sub]:
                            outer_vertices = union(outer_vertices, k)
                        else:
                            for l in Holes[i]:
                                if k not in Matrix_out[l] and not isline([k, l], other_edges) and not isline([k, l], Holes[i]):
                                    dist = distance_vox(k, l, Centroid, coordinates, voxels, dist_path)
                                    if dist > 1.5 * threshold:
                                        continue
                                    other_edges.append([k, l])
        for v in outer_vertices:
            edge = [v_sub, v]
            if not isline(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                candidates_distances[v_sub].append(dist)
        for edge in other_edges:
            if not isline(edge, candidates_edges[vertex]):
                candidates_edges[v_sub].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                candidates_distances[v_sub].append(dist)
        # lst_remove = [index for index in range(len(candidates_edges[v_sub])) if
        #               isline(candidates_edges[v_sub][index], list_boundary_vox) and candidates_distances[v_sub][
        #                   index] > 2 * threshold]
        # candidates_edges[v_sub] = Remove_index(candidates_edges[v_sub], lst_remove)
        # candidates_distances[v_sub] = Remove_index(candidates_distances[v_sub], lst_remove)

save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_4.txt")
save_list_of_lists(Register_Added, "Register_Added" + suffix[:-4] + "_step_4.txt")

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

# with open("Matrix_out" + suffix[:-4] + "_step_4.txt", 'r') as f:
#     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# for i in range(len(Matrix_out)):
#     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed
#
# with open("Register_Added" + suffix[:-4] + "_step_4.txt", 'r') as f:
#     Register_Added = [[int(num) for num in line.split(',')] for line in f]
#
# Register_Removed = []

start_time = time()
print("\n\nAdd edges to avoid digons...")
log_file_detailed.write("\n\nAdd edges to avoid digons...")
log_file_detailed.flush()
log_file.write("\n\nAdd edges to avoid digons...")
log_file.flush()
# Add edges to avoid digons
print("Compute candidate edges to avoid digons...")
K2 = []
for i in range(numbPoints):
    Nei_i = Matrix_out[i]
    for j in [item for item in Nei_i if item > i]:
        K2.append([i, j])
candidates_edges = [[] for item in range(len(K2))]
candidates_distances = [[] for item in range(len(K2))]
lst_compute = range(len(K2))
count = 0
for index in lst_compute:
    candidate_edge = K2[index]
    candidates_edges[index] = []
    candidates_distances[index] = []
    count += 1
    print("\rComputing", count, "out of", len(K2), end="")
    if sublist(candidate_edge, list_boundary):
        continue
    lst = mutual_neighbors(candidate_edge, Matrix_out)
    if len(lst) > 6:
        continue
    if not cycles(edges_of(mutual_neighbors(candidate_edge, Matrix_out), Matrix_out), minimum_allowed=3):
        for w in setdiff(Matrix_out[candidate_edge[0]], [candidate_edge[1]] + Matrix_out[candidate_edge[1]]):
            candidates_edges[index].append([candidate_edge[1], w])
            dist = distance_vox(candidate_edge[1], w, Centroid, coordinates, voxels, dist_path)
            candidates_distances[index].append(dist)
        for w in setdiff(Matrix_out[candidate_edge[1]], [candidate_edge[0]] + Matrix_out[candidate_edge[0]]):
            candidates_edges[index].append([candidate_edge[0], w])
            dist = distance_vox(candidate_edge[0], w, Centroid, coordinates, voxels, dist_path)
            candidates_distances[index].append(dist)
        for edge in combinations(lst, 2):
            if edge[0] not in Matrix_out[edge[1]]:
                candidates_edges[index].append(list(edge))
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                candidates_distances[index].append(dist)

save_list_of_lists(candidates_edges, "candidates_edges" + suffix[:-4] + "_step_5.txt")
save_list_of_lists(candidates_distances, "candidates_distances" + suffix[:-4] + "_step_5.txt")

# with open("Matrix_out" + suffix[:-4] + "_step_4.txt", 'r') as f:
#     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# for i in range(len(Matrix_out)):
#     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed
#
# with open("Register_Added" + suffix[:-4] + "_step_4.txt", 'r') as f:
#     Register_Added = [[int(num) for num in line.split(',')] for line in f]
#
# K2 = []
# for i in range(numbPoints):
#     Nei_i = Matrix_out[i]
#     for j in [item for item in Nei_i if item > i]:
#         K2.append([i, j])
#
# with open("candidates_distances" + suffix[:-4] + "_step_5.txt", 'r') as f:
#     candidates_distances = [[float(num) for num in line.split(',')] for line in f]
# for i in range(len(candidates_distances)):
#     candidates_distances[i] = [x for x in candidates_distances[i] if x != -1]  # -1's get removed
#
# candidates_edges = read_list_of_lists_of_lists("candidates_edges" + suffix[:-4] + "_step_5.txt")

# Add short candidate edges...
print("\nAdd short candidate edges...")
while True:
    minima = [min(candidates_distances[index], default=math.inf) for index in range(len(K2))]
    if min(minima) > 6 * threshold:  # == math.inf:
        break
    index_edge = minima.index(min(minima))
    Ind = candidates_distances[index_edge].index(min(candidates_distances[index_edge]))
    edge = candidates_edges[index_edge][Ind]
    Register_Added.append(edge)
    Matrix_out = Add(edge[0], edge[1], Matrix_out)
    print("Added:", K2[index_edge], edge, candidates_distances[index_edge][Ind])
    log_file_detailed.write("\nAdded to avoid a digon between {} and {}: edge {} with distance {}".format(K2[index_edge][0], K2[index_edge][1], edge, candidates_distances[index_edge][Ind]))
    log_file_detailed.flush()
    K2.append(edge)
    candidates_edges.append([])
    candidates_distances.append([])
    All = union(edge, mutual_neighbors(edge, Matrix_out))
    All_ind = unique([ind(item, K2) for item in All])
    for ind_sub in All_ind:
        candidates_edges[ind_sub] = []
        candidates_distances[ind_sub] = []
        candidate_edge = K2[ind_sub]
        if sublist(candidate_edge, list_boundary):
            continue
        lst = mutual_neighbors(candidate_edge, Matrix_out)
        if len(lst) > 6:
            continue
        if not cycles(edges_of(mutual_neighbors(candidate_edge, Matrix_out), Matrix_out), minimum_allowed=3):
            for w in setdiff(Matrix_out[candidate_edge[0]], [candidate_edge[1]] + Matrix_out[candidate_edge[1]]):
                candidates_edges[ind_sub].append([candidate_edge[1], w])
                dist = distance_vox(candidate_edge[1], w, Centroid, coordinates, voxels, dist_path)
                candidates_distances[ind_sub].append(dist)
            for w in setdiff(Matrix_out[candidate_edge[1]], [candidate_edge[0]] + Matrix_out[candidate_edge[0]]):
                candidates_edges[ind_sub].append([candidate_edge[0], w])
                dist = distance_vox(candidate_edge[0], w, Centroid, coordinates, voxels, dist_path)
                candidates_distances[ind_sub].append(dist)
            for edge in combinations(lst, 2):
                if edge[0] not in Matrix_out[edge[1]]:
                    candidates_edges[ind_sub].append(list(edge))
                    dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                    candidates_distances[ind_sub].append(dist)

save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_6.txt")
save_list_of_lists(Register_Added, "Register_Added" + suffix[:-4] + "_step_6.txt")

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

# with open("Matrix_out" + suffix[:-4] + "_step_6.txt", 'r') as f:
#     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# for i in range(len(Matrix_out)):
#     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed
#
# with open("Register_Added" + suffix[:-4] + "_step_6.txt", 'r') as f:
#     Register_Added = [[int(num) for num in line.split(',')] for line in f]

start_time = time()
print("\n\nAdd edges of degree 3...")
log_file_detailed.write("\n\nAdd edges of degree 3...")
log_file_detailed.flush()
log_file.write("\n\nAdd edges of degree 3...")
log_file.flush()
# Add edges of degree 3
for vertex in range(numbPoints):
    link = vertexlink(vertex, Matrix_out)
    for tri in link:
        lst = setdiff(mutual_neighbors(tri, Matrix_out), [vertex] + Matrix_out[vertex])
        for outer_vertex in [item for item in lst if item > vertex]:
            dist = distance_vox(vertex, outer_vertex, Centroid, coordinates, voxels, dist_path)
            if dist < 3 * threshold:
                Matrix_out = Add(vertex, outer_vertex, Matrix_out)
                Register_Added.append([vertex, outer_vertex])
                print("Added around a triangle:", vertex, tri, [vertex, outer_vertex], dist)
                log_file_detailed.write("\nAdded around the triangle {}, edge {} with distance {}".format(tri, [vertex, outer_vertex], dist))
                log_file_detailed.flush()

save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_7.txt")
save_list_of_lists(Register_Added, "Register_Added" + suffix[:-4] + "_step_7.txt")

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

start_time = time()
print("\n\nCompute candidate edges within pseudoholes to be added...")
log_file_detailed.write("\n\nCompute candidate edges within pseudoholes to be added...")
log_file_detailed.flush()
log_file.write("\n\nCompute candidate edges within pseudoholes to be added...")
log_file.flush()
candidates_edges = [[] for item in range(numbPoints)]
candidates_distances = [[] for item in range(numbPoints)]
lst_compute = range(numbPoints)
while lst_compute:
    vertex = lst_compute[0]
    candidates_edges[vertex] = []
    candidates_distances[vertex] = []
    print("\rCompute the candidates from pseudoholes:", vertex, "out of", numbPoints - 1, end="")
    if vertex in ignore:
        lst_compute = lst_compute[1:]
        continue
    Holes = []
    open_edges = open_edges_in_link(vertex, Matrix_out, Centroid)
    if len(open_edges) <= 75:
        Holes = pseudoholes(vertex, Matrix_out, list_boundary, Centroid) + holes(vertex, Matrix_out, list_boundary, Centroid, True)[0]
    elif vertex not in list_boundary:
        Holes = [unique(open_edges)]
    outer_vertices, other_edges = [], []
    for i in range(len(Holes)):
        indices = []
        for j in combinations(range(len(Holes[i])), 2):
            j = list(j)
            edge = [Holes[i][item] for item in j]
            if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                continue
            indices.append(j)
        for j in range(len(indices)):
            edge = [Holes[i][indices[j][0]], Holes[i][indices[j][1]]]
            if all([v in list_boundary_vox for v in Holes[i]]):
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                if dist > 2 * threshold:
                    continue
            if not isline(edge, candidates_edges[vertex]):
                candidates_edges[vertex].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                candidates_distances[vertex].append(dist)
                for k in mutual_neighbors(edge, Matrix_out):
                    if k not in [vertex] + Matrix_out[vertex]:
                        outer_vertices = union(outer_vertices, k)
                    else:
                        for l in Holes[i]:
                            if k not in Matrix_out[l] and not isline([k, l], other_edges) and not isline([k, l], Holes[i]):
                                dist = distance_vox(k, l, Centroid, coordinates, voxels, dist_path)
                                if dist > 1.5 * threshold:
                                    continue
                                other_edges.append([k, l])
    for v in outer_vertices:
        edge = [vertex, v]
        if not isline(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
            candidates_distances[vertex].append(dist)
    for edge in other_edges:
        if not isline(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
            candidates_distances[vertex].append(dist)
    lst_compute = lst_compute[1:]
print("")

save_list_of_lists(candidates_edges, "candidates_edges" + suffix[:-4] + "_step_8.txt")
save_list_of_lists(candidates_distances, "candidates_distances" + suffix[:-4] + "_step_8.txt")

log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

# with open("Matrix_out" + suffix[:-4] + "_step_5.txt", 'r') as f:
#     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# for i in range(len(Matrix_out)):
#     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed
# for i in range(len(Matrix_out)):
#     Matrix_out[i] = setdiff(Matrix_out[i], [i])
#
# with open("Register_Added" + suffix[:-4] + "_step_5.txt", 'r') as f:
#     Register_Added = [[int(num) for num in line.split(',')] for line in f]
#
# with open("Register_Removed" + suffix[:-4] + "_step_5.txt", 'r') as f:
#     Register_Removed = [[int(num) for num in line.split(',')] for line in f]
#
# with open("candidates_distances" + suffix[:-4] + "_step_5.txt", 'r') as f:
#     candidates_distances = [[float(num) for num in line.split(',')] for line in f]
# for i in range(len(candidates_distances)):
#     candidates_distances[i] = [x for x in candidates_distances[i] if x != -1]  # -1's get removed
#
# candidates_edges = read_list_of_lists_of_lists("candidates_edges" + suffix[:-4] + "_step_5.txt")

start_time = time()
print("\n\nAdd short candidate edges to close the pseudoholes...")
log_file_detailed.write("\n\nAdd short candidate edges to close the pseudoholes...")
log_file_detailed.flush()
log_file.write("\n\nAdd short candidate edges to close the pseudoholes...")
log_file.flush()
count = 0
while True:
    minima = [min(candidates_distances[vertex], default=math.inf) for vertex in range(numbPoints)]
    if min(minima) > 5 * threshold:  # == math.inf:
        break
    vertex = minima.index(min(minima))
    Ind = candidates_distances[vertex].index(min(candidates_distances[vertex]))
    edge = candidates_edges[vertex][Ind]
    Register_Added.append(edge)
    Matrix_out = Add(edge[0], edge[1], Matrix_out)
    print("Added:", vertex, edge, candidates_distances[vertex][Ind])
    log_file_detailed.write("\nAdded: linkvertex {}, edge {} with distance {}".format(vertex, edge, candidates_distances[vertex][Ind]))
    log_file_detailed.flush()
    All = union(edge, mutual_neighbors(edge, Matrix_out))
    for v_sub in All:
        candidates_edges[v_sub] = []
        candidates_distances[v_sub] = []
        Holes = []
        open_edges = open_edges_in_link(v_sub, Matrix_out, Centroid)
        if len(open_edges) <= 75:
            Holes = pseudoholes(v_sub, Matrix_out, list_boundary, Centroid) + holes(v_sub, Matrix_out, list_boundary, Centroid, True)[0]
        elif v_sub not in list_boundary:
            Holes = [unique(open_edges)]
        outer_vertices, other_edges = [], []
        for i in range(len(Holes)):
            indices = []
            for j in combinations(range(len(Holes[i])), 2):
                j = list(j)
                edge = [Holes[i][item] for item in j]
                if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                    continue
                indices.append(j)
            for j in range(len(indices)):
                edge = [Holes[i][indices[j][0]], Holes[i][indices[j][1]]]
                if all([v in list_boundary_vox for v in Holes[i]]):
                    dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                    if dist > 2 * threshold:
                        continue
                if not isline(edge, candidates_edges[v_sub]):
                    candidates_edges[v_sub].append(edge)
                    dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                    candidates_distances[v_sub].append(dist)
                    for k in mutual_neighbors(edge, Matrix_out):
                        if k not in [v_sub] + Matrix_out[v_sub]:
                            outer_vertices = union(outer_vertices, k)
                        else:
                            for l in Holes[i]:
                                if k not in Matrix_out[l] and not isline([k, l], other_edges) and not isline([k, l], Holes[i]):
                                    dist = distance_vox(k, l, Centroid, coordinates, voxels, dist_path)
                                    if dist > 1.5 * threshold:
                                        continue
                                    other_edges.append([k, l])
        for v in outer_vertices:
            edge = [v_sub, v]
            if not isline(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                candidates_distances[v_sub].append(dist)
        for edge in other_edges:
            if not isline(edge, candidates_edges[vertex]):
                candidates_edges[v_sub].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                candidates_distances[v_sub].append(dist)

save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_9.txt")
save_list_of_lists(Register_Added, "Register_Added" + suffix[:-4] + "_step_9.txt")

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

# with open("Matrix_out" + suffix[:-4] + "_step_8.txt", 'r') as f:
#     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# for i in range(len(Matrix_out)):
#     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed
#
# Register_Removed = []
#
# with open("Register_Added" + suffix[:-4] + "_step_8.txt", 'r') as f:
#     Register_Added = [[int(num) for num in line.split(',')] for line in f]

log_file_detailed.write("\n\nRemove excessive edges...")
log_file_detailed.flush()
log_file.write("\nRemove excessive edges...")
log_file.flush()
print("\n\nRemove excessive edges...")

for vertex in range(numbPoints):
    # print("\rProgress:", vertex, "out of", numbPoints - 1, end="")
    for triangle in vertexlink(vertex, Matrix_out):
        for edge_test in edges_vertexlink(vertex, Matrix_out):
            if intersect(edge_test, triangle):
                continue
            A, B = barycentric_intersection_in_line(triangle, edge_test[0], edge_test[1], Centroid)
            if min(A, B) >= 0.3:
                a, b, c = barycentric_intersection_in_triangle(triangle, edge_test, Centroid)
                if min(a, b, c) > 0:
                    dist = distance_vox(edge_test[0], edge_test[1], Centroid, coordinates, voxels, dist_path)
                    dists = [distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path) for edge in combinations(triangle, 2)]
                    if dist > max([item + 0.5 * threshold for item in dists] + [2 * threshold]):
                        Register_Added = Remove_index(Register_Added, ind(edge_test, Register_Added))
                        Register_Removed.append(edge_test)
                        Matrix_out = Remove(edge_test[0], edge_test[1], Matrix_out)
                        print("Removed:", vertex, triangle, edge_test, dist, dists)
                        log_file_detailed.write("\nRemoved excessive edge {} blocked by {} in the link of {}".format(edge_test, triangle, vertex))
                        log_file_detailed.flush()

save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_10.txt")
save_list_of_lists(Register_Added, "Register_Added" + suffix[:-4] + "_step_10.txt")
save_list_of_lists(Register_Removed, "Register_Removed" + suffix[:-4] + "_step_10.txt")

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

# with open("Matrix_out" + suffix[:-4] + "_step_10.txt", 'r') as f:
#     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# for i in range(len(Matrix_out)):
#     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed
# with open("Register_Added" + suffix[:-4] + "_step_10.txt", 'r') as f:
#     Register_Added = [[int(num) for num in line.split(',')] for line in f]
#
# with open("Register_Removed" + suffix[:-4] + "_step_10.txt", 'r') as f:
#     Register_Removed = [[int(num) for num in line.split(',')] for line in f]
#
# with open('list_boundary_alpha_shape' + suffix[-8:], 'r') as f:
#     list_boundary = [int(elem) - 1 for elem in f.read().split(',') if elem]
#
# list_interior = setdiff(range(numbPoints), list_boundary)

# Collect candidates to close all remaining holes...
start_time = time()
print("\n\nCollect candidates to close all remaining holes...")
log_file_detailed.write("\n\nCollect candidates to close all remaining holes...")
log_file_detailed.flush()
log_file.write("\n\nCollect candidates to close all remaining holes...")
log_file.flush()
candidates_edges = [[] for item in range(numbPoints)]
candidates_distances = [[] for item in range(numbPoints)]
lst_compute = range(numbPoints)
while lst_compute:
    vertex = lst_compute[0]
    candidates_edges[vertex] = []
    candidates_distances[vertex] = []
    print("\rCompute the candidates:", vertex, "out of", numbPoints - 1, end="")
    if vertex in union(list_boundary, list_boundary_vox, ignore):
        lst_compute = lst_compute[1:]
        continue
    Holes = holes(vertex, Matrix_out, list_boundary, Centroid, True)[0]
    Matrix_reduced = reducedvertexlink(vertex, Matrix_out, list_boundary, Centroid, return_Matrix=True)[1]
    Holes += [hole for hole in holes(vertex, Matrix_reduced, list_boundary, Centroid, True)[0] if
              not isline_exact(hole, Holes)[0]]
    scale = 1
    if not Holes:
        Holes = pseudoholes(vertex, Matrix_out, list_boundary, Centroid)
        scale = 0.5
    outer_vertices, other_edges = [], []
    for i in range(len(Holes)):
        boundary_hole = len(intersect(Holes[i], union(list_boundary, list_boundary_vox))) >= len(Holes[i]) - 1
        for j in combinations(range(len(Holes[i])), 2):
            j = list(j)
            edge = [Holes[i][item] for item in j]
            if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                continue
            dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
            if dist > scale * 10 * threshold:
                continue
            if all([v in list_boundary_vox for v in Holes[i]]):
                if dist > 3 * threshold:
                    continue
            if boundary_hole:
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                if dist > 6 * threshold:
                    continue
            if not isline(edge, candidates_edges[vertex]):
                candidates_edges[vertex].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                candidates_distances[vertex].append(dist)
                for k in mutual_neighbors(edge, Matrix_out):
                    if k not in [vertex] + Matrix_out[vertex]:
                        outer_vertices = union(outer_vertices, k)
                    else:
                        for l in Holes[i]:
                            if k not in Matrix_out[l] and not isline([k, l], other_edges) and not isline([k, l], Holes[i]):
                                dist = distance_vox(k, l, Centroid, coordinates, voxels, dist_path)
                                if dist > 1.5 * threshold:
                                    continue
                                other_edges.append([k, l])
    for v in outer_vertices:
        edge = [vertex, v]
        dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
        if dist > scale * 10 * threshold:
            continue
        if not isline(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            candidates_distances[vertex].append(dist)
    for edge in other_edges:
        dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
        if dist > scale * 10 * threshold:
            continue
        if not isline(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            candidates_distances[vertex].append(dist)
    lst_compute = lst_compute[1:]
print("")

save_list_of_lists(candidates_edges, "candidates_edges" + suffix[:-4] + "_step_11.txt")
save_list_of_lists(candidates_distances, "candidates_distances" + suffix[:-4] + "_step_11.txt")

log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()


with open("Matrix_out" + suffix[:-4] + "_step_10.txt", 'r') as f:
    Matrix_out = [[int(num) for num in line.split(',')] for line in f]
for i in range(len(Matrix_out)):
    Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed

with open("Register_Added" + suffix[:-4] + "_step_10.txt", 'r') as f:
    Register_Added = [[int(num) for num in line.split(',')] for line in f]

with open("Register_Removed" + suffix[:-4] + "_step_10.txt", 'r') as f:
    Register_Removed = [[int(num) for num in line.split(',')] for line in f]

with open("candidates_distances" + suffix[:-4] + "_step_11.txt", 'r') as f:
    candidates_distances = [[float(num) for num in line.split(',')] for line in f]
for i in range(len(candidates_distances)):
    candidates_distances[i] = [x for x in candidates_distances[i] if x != -1]  # -1's get removed

candidates_edges = read_list_of_lists_of_lists("candidates_edges" + suffix[:-4] + "_step_11.txt")

start_time = time()
print("\n\nAdd short candidate edges to close remaining holes...")
log_file_detailed.write("\n\nAdd short candidate edges to close remaining holes...")
log_file_detailed.flush()
log_file.write("\n\nAdd short candidate edges to close remaining holes...")
log_file.flush()
while True:
    minima = [min(candidates_distances[vertex], default=math.inf) for vertex in range(numbPoints)]
    if min(minima) > 10 * threshold:  # == math.inf:
        break
    vertex = minima.index(min(minima))
    Ind = candidates_distances[vertex].index(min(candidates_distances[vertex]))
    edge = candidates_edges[vertex][Ind]
    Register_Added.append(edge)
    Matrix_out = Add(edge[0], edge[1], Matrix_out)
    print("Added:", vertex, edge, candidates_distances[vertex][Ind])
    log_file_detailed.write("\nAdded: linkvertex {}, edge {} with distance {}".format(vertex, edge, candidates_distances[vertex][Ind]))
    log_file_detailed.flush()
    All = union(edge, mutual_neighbors(edge, Matrix_out))
    for v_sub in All:
        if v_sub in union(list_boundary, list_boundary_vox) or v_sub in ignore:
            continue
        candidates_edges[v_sub] = []
        candidates_distances[v_sub] = []
        Holes = holes(v_sub, Matrix_out, list_boundary, Centroid, True)[0]
        Matrix_temp = reducedvertexlink(v_sub, Matrix_out, list_boundary, Centroid, return_Matrix=True)[1]
        Holes += [hole for hole in holes(v_sub, Matrix_temp, list_boundary, Centroid, True)[0] if
                  not isline_exact(hole, Holes)[0]]
        scale = 1
        if not Holes:
            Holes = pseudoholes(v_sub, Matrix_out, list_boundary, Centroid)
            scale = 0.5
        outer_vertices, other_edges = [], []
        for i in range(len(Holes)):
            boundary_hole = len(intersect(Holes[i], union(list_boundary, list_boundary_vox))) >= len(Holes[i]) - 1
            for j in combinations(range(len(Holes[i])), 2):
                j = list(j)
                edge = [Holes[i][item] for item in j]
                if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                    continue
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                if dist > scale * 10 * threshold:
                    continue
                if all([v in list_boundary_vox for v in Holes[i]]):
                    if dist > 3 * threshold:
                        continue
                if boundary_hole:
                    dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                    if dist > 6 * threshold:
                        continue
                if not isline(edge, candidates_edges[v_sub]):
                    candidates_edges[v_sub].append(edge)
                    dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                    candidates_distances[v_sub].append(dist)
                    for k in mutual_neighbors(edge, Matrix_out):
                        if k not in [v_sub] + Matrix_out[v_sub]:
                            outer_vertices = union(outer_vertices, k)
                        else:
                            for l in Holes[i]:
                                if k not in Matrix_out[l] and not isline([k, l], other_edges) and not isline([k, l], Holes[i]):
                                    dist = distance_vox(k, l, Centroid, coordinates, voxels, dist_path)
                                    if dist > 1.5 * threshold:
                                        continue
                                    other_edges.append([k, l])
        for v in outer_vertices:
            edge = [v_sub, v]
            dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
            if dist > scale * 10 * threshold:
                continue
            if not isline(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                candidates_distances[v_sub].append(dist)
        for edge in other_edges:
            dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
            if dist > scale * 10 * threshold:
                continue
            if not isline(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                candidates_distances[v_sub].append(dist)

save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_12.txt")
save_list_of_lists(Register_Added, "Register_Added" + suffix[:-4] + "_step_12.txt")

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

# Collect candidates to close all remaining boundary holes...
start_time = time()
print("\n\nCollect candidates to close all remaining boundary holes...")
log_file_detailed.write("\n\nCollect candidates to close all remaining boundary holes...")
log_file_detailed.flush()
log_file.write("\n\nCollect candidates to close all remaining boundary holes...")
log_file.flush()
candidates_edges = [[] for item in range(numbPoints)]
candidates_distances = [[] for item in range(numbPoints)]
lst_compute = range(numbPoints)
while lst_compute:
    vertex = lst_compute[0]
    candidates_edges[vertex] = []
    candidates_distances[vertex] = []
    print("\rCompute the candidates:", vertex, "out of", numbPoints - 1, end="")
    if vertex in union(list_boundary, list_boundary_vox, ignore):
        lst_compute = lst_compute[1:]
        continue
    Holes = holes(vertex, Matrix_out, list_boundary, Centroid, True)[0]
    Matrix_reduced = reducedvertexlink(vertex, Matrix_out, list_boundary, Centroid, return_Matrix=True)[1]
    Holes += [hole for hole in holes(vertex, Matrix_reduced, list_boundary, Centroid, True)[0] if
              not isline_exact(hole, Holes)[0]]
    outer_vertices, other_edges = [], []
    for i in range(len(Holes)):
        for j in combinations(range(len(Holes[i])), 2):
            j = list(j)
            edge = [Holes[i][item] for item in j]
            if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                continue
            dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
            if dist > 10 * threshold:
                continue
            if not isline(edge, candidates_edges[vertex]):
                candidates_edges[vertex].append(edge)
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                candidates_distances[vertex].append(dist)
                for k in mutual_neighbors(edge, Matrix_out):
                    if k not in [vertex] + Matrix_out[vertex]:
                        outer_vertices = union(outer_vertices, k)
                    else:
                        for l in Holes[i]:
                            if k not in Matrix_out[l] and not isline([k, l], other_edges) and not isline([k, l], Holes[i]):
                                dist = distance_vox(k, l, Centroid, coordinates, voxels, dist_path)
                                if dist > 1.5 * threshold:
                                    continue
                                other_edges.append([k, l])
    for v in outer_vertices:
        edge = [vertex, v]
        dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
        if dist > 10 * threshold:
            continue
        if not isline(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            candidates_distances[vertex].append(dist)
    for edge in other_edges:
        dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
        if dist > 10 * threshold:
            continue
        if not isline(edge, candidates_edges[vertex]):
            candidates_edges[vertex].append(edge)
            candidates_distances[vertex].append(dist)
    lst_compute = lst_compute[1:]
print("")

save_list_of_lists(candidates_edges, "candidates_edges" + suffix[:-4] + "_step_13.txt")
save_list_of_lists(candidates_distances, "candidates_distances" + suffix[:-4] + "_step_13.txt")

log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

# with open("Matrix_out" + suffix[:-4] + "_step_10.txt", 'r') as f:
#     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# for i in range(len(Matrix_out)):
#     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed
#
# with open("Register_Added" + suffix[:-4] + "_step_10.txt", 'r') as f:
#     Register_Added = [[int(num) for num in line.split(',')] for line in f]
#
# with open("Register_Removed" + suffix[:-4] + "_step_10.txt", 'r') as f:
#     Register_Removed = [[int(num) for num in line.split(',')] for line in f]
#
# with open("candidates_distances" + suffix[:-4] + "_step_11.txt", 'r') as f:
#     candidates_distances = [[float(num) for num in line.split(',')] for line in f]
# for i in range(len(candidates_distances)):
#     candidates_distances[i] = [x for x in candidates_distances[i] if x != -1]  # -1's get removed
#
# candidates_edges = read_list_of_lists_of_lists("candidates_edges" + suffix[:-4] + "_step_11.txt")

start_time = time()
print("\n\nAdd short candidate edges to close remaining boundary holes...")
log_file_detailed.write("\n\nAdd short candidate edges to close remaining boundary holes...")
log_file_detailed.flush()
log_file.write("\n\nAdd short candidate edges to close remaining boundary holes...")
log_file.flush()
while True:
    minima = [min(candidates_distances[vertex], default=math.inf) for vertex in range(numbPoints)]
    if min(minima) > 10 * threshold:  # == math.inf:
        break
    vertex = minima.index(min(minima))
    Ind = candidates_distances[vertex].index(min(candidates_distances[vertex]))
    edge = candidates_edges[vertex][Ind]
    Register_Added.append(edge)
    Matrix_out = Add(edge[0], edge[1], Matrix_out)
    print("Added:", vertex, edge, candidates_distances[vertex][Ind])
    log_file_detailed.write("\nAdded: linkvertex {}, edge {} with distance {}".format(vertex, edge, candidates_distances[vertex][Ind]))
    log_file_detailed.flush()
    All = union(edge, mutual_neighbors(edge, Matrix_out))
    for v_sub in All:
        if v_sub in union(list_boundary, list_boundary_vox) or v_sub in ignore:
            continue
        candidates_edges[v_sub] = []
        candidates_distances[v_sub] = []
        Holes = holes(v_sub, Matrix_out, list_boundary, Centroid, True)[0]
        Matrix_temp = reducedvertexlink(v_sub, Matrix_out, list_boundary, Centroid, return_Matrix=True)[1]
        Holes += [hole for hole in holes(v_sub, Matrix_temp, list_boundary, Centroid, True)[0] if
                  not isline_exact(hole, Holes)[0]]
        outer_vertices, other_edges = [], []
        for i in range(len(Holes)):
            for j in combinations(range(len(Holes[i])), 2):
                j = list(j)
                edge = [Holes[i][item] for item in j]
                if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
                    continue
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                if dist > 10 * threshold:
                    continue
                if not isline(edge, candidates_edges[v_sub]):
                    candidates_edges[v_sub].append(edge)
                    dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                    candidates_distances[v_sub].append(dist)
                    for k in mutual_neighbors(edge, Matrix_out):
                        if k not in [v_sub] + Matrix_out[v_sub]:
                            outer_vertices = union(outer_vertices, k)
                        else:
                            for l in Holes[i]:
                                if k not in Matrix_out[l] and not isline([k, l], other_edges) and not isline([k, l], Holes[i]):
                                    dist = distance_vox(k, l, Centroid, coordinates, voxels, dist_path)
                                    if dist > 1.5 * threshold:
                                        continue
                                    other_edges.append([k, l])
        for v in outer_vertices:
            edge = [v_sub, v]
            dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
            if dist > 10 * threshold:
                continue
            if not isline(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                candidates_distances[v_sub].append(dist)
        for edge in other_edges:
            dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
            if dist > 10 * threshold:
                continue
            if not isline(edge, candidates_edges[v_sub]):
                candidates_edges[v_sub].append(edge)
                candidates_distances[v_sub].append(dist)

save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_14.txt")
save_list_of_lists(Register_Added, "Register_Added" + suffix[:-4] + "_step_14.txt")

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
log_file_detailed.flush()
log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n{}".format(compute_time(start_time, True)))
log_file.flush()

with open("Matrix_out" + suffix[:-4] + "_step_14.txt", 'r') as f:
    Matrix_out = [[int(num) for num in line.split(',')] for line in f]
for i in range(len(Matrix_out)):
    Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed

with open("Register_Added" + suffix[:-4] + "_step_14.txt", 'r') as f:
    Register_Added = [[int(num) for num in line.split(',')] for line in f]
#
# Detect digons and collect candidates to repair non-Hamiltonian links...
start_time = time()
print("\n\nDetect digons and collect candidates to repair non-Hamiltonian links...")
log_file_detailed.write("\n\nDetect digons and collect candidates to repair non-Hamiltonian links...")
log_file_detailed.flush()
log_file.write("\n\nDetect digons and collect candidates to repair non-Hamiltonian links...")
log_file.flush()
candidates_edges = [[] for item in range(numbPoints)]
candidates_distances = [[] for item in range(numbPoints)]
lst_Hamiltonian, lst_not_Hamiltonian = [], []
digons, digons_neighbors = [], []
loop_over = range(numbPoints)
avoid = union(list_boundary, list_boundary_vox, ignore)
while loop_over:
    vertex = loop_over[0]
    print("\rCompute the candidates:", vertex, "out of", numbPoints - 1, end="")
    if vertex in avoid:
        loop_over = loop_over[1:]
        continue
    Matrix_reduced, triangles = reducedvertexlink(vertex, Matrix_out, list_boundary, Centroid, return_Matrix=True,
                                                  return_triangles=True)[1:]
    tets = [[vertex] + item for item in triangles]
    if isHamiltonian(vertex, Matrix_reduced, list_boundary, Centroid, predetermined_tetrahedra=tets)[1]:
        lst_Hamiltonian.append(vertex)
        loop_over = loop_over[1:]
        continue
    else:
        if isHamiltonian(vertex, Matrix_out, list_boundary, Centroid)[1]:
            lst_Hamiltonian.append(vertex)
            loop_over = loop_over[1:]
            continue
        else:
            lst_not_Hamiltonian.append(vertex)
    if len(Matrix_out[vertex]) < 6:
        loop_over = loop_over[1:]
        continue  # this vertex might be a digon
    if len(Matrix_out[vertex]) > 40:
        loop_over = loop_over[1:]
        continue
    lst = [len(mutual_neighbors([vertex, v], Matrix_out)) for v in Matrix_out[vertex]]
    Ind = sorted(range(len(lst)), key=lambda k: lst[k])
    # avoid repairing digons
    neighbors_ordered = [Matrix_out[vertex][item] for item in Ind]
    digon_involved = False
    for v in neighbors_ordered:
        if len(Matrix_out[v]) > 5 or v in list_boundary:
            continue
        if cycles(edges_of(mutual_neighbors([vertex, v], Matrix_out), Matrix_out), minimum_allowed=3):
            continue
        Matrix_temp = Remove(v, vertex, Matrix_out)
        if isHamiltonian(vertex, Matrix_temp, list_boundary, Centroid)[1]:
            digon_involved = True
            digons.append(v)
            ignore.append(v)
            avoid.append(v)
            digon_neighbors = Matrix_out[v]
            digons_neighbors.append(digon_neighbors)
            for u in digon_neighbors:
                Matrix_out = Remove(u, v, Matrix_out)
            print("\ndigon found: {} - neighbors: {}".format(v, digon_neighbors))
            log_file_detailed.write("\ndigon found: {} - neighbors: {}".format(v, digon_neighbors))
            log_file_detailed.flush()
            for v_sub in digon_neighbors:
                if v_sub in union(list_boundary, list_boundary_vox, ignore):
                    continue
                Matrix_reduced, triangles = reducedvertexlink(v_sub, Matrix_out, list_boundary, Centroid,
                                                   return_Matrix=True, return_triangles=True)[1:]
                lst_Hamiltonian = setdiff(lst_Hamiltonian, v_sub)
                lst_not_Hamiltonian = setdiff(lst_not_Hamiltonian, v_sub)
                loop_over = setdiff(loop_over, v_sub)
                tets = [[v_sub] + item for item in triangles]
                if isHamiltonian(v_sub, Matrix_reduced, list_boundary, Centroid, predetermined_tetrahedra=tets)[1]:
                    lst_Hamiltonian.append(v_sub)
                else:
                    if isHamiltonian(v_sub, Matrix_out, list_boundary, Centroid)[1]:
                        lst_Hamiltonian.append(v_sub)
                    else:
                        lst_not_Hamiltonian.append(v_sub)
                        loop_over.append(v_sub)
            break
    if digon_involved:
        continue
    second_neighbors = union(*[Matrix_out[item] for item in Matrix_out[vertex]])
    initial_Hamiltonian = intersect(second_neighbors, lst_Hamiltonian)
    candidates_Hamiltonian, candidates_add, candidates_add_dist = [], [], []
    for i in Matrix_out[vertex]:
        for j in [item for item in setdiff(Matrix_out[vertex], Matrix_out[i]) if item > i]:
            edge = [i, j]
            dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
            if dist > 10 * threshold:
                continue
            Matrix_temp = Add(edge[0], edge[1], Matrix_out)
            Matrix_reduced, triangles = reducedvertexlink(vertex, Matrix_temp, list_boundary, Centroid,
                                                          return_Matrix=True,
                                                          return_triangles=True)[1:]
            tets = [[vertex] + item for item in triangles]
            check = isHamiltonian(vertex, Matrix_reduced, list_boundary, Centroid, predetermined_tetrahedra=tets)[1]
            if not check:
                check = isHamiltonian(vertex, Matrix_temp, list_boundary, Centroid)[1]
            if check:
                candidates_add.append(edge)
                candidates_add_dist.append(dist)
                candidates_Hamiltonian.append([vertex])
                for v in setdiff(second_neighbors, vertex):
                    if v in union(list_boundary, list_boundary_vox, ignore):
                        continue
                    if v in lst_Hamiltonian and v not in edge:
                        candidates_Hamiltonian[-1].append(v)
                        continue
                    if v not in union(edge, mutual_neighbors(edge, Matrix_out)):
                        continue
                    Matrix_reduced, triangles = reducedvertexlink(v, Matrix_temp, list_boundary, Centroid,
                                                                  return_Matrix=True,
                                                                  return_triangles=True)[1:]
                    tets = [[v] + item for item in triangles]
                    if isHamiltonian(v, Matrix_reduced, list_boundary, Centroid, predetermined_tetrahedra=tets)[1]:
                        candidates_Hamiltonian[-1].append(v)
                        continue
                    if isHamiltonian(v, Matrix_temp, list_boundary, Centroid)[1]:
                        candidates_Hamiltonian[-1].append(v)
    diff = [len(item) - len(initial_Hamiltonian) for item in candidates_Hamiltonian]
    m = max(diff, default=0)
    if m > 0:
        Ind = [index for index in range(len(candidates_Hamiltonian)) if diff[index] == m]
        for index in setdiff(range(len(candidates_add_dist)), Ind):
            candidates_add_dist[index] = math.inf
        Ind = candidates_add_dist.index(min(candidates_add_dist))  # add the shortest edge
        edge = candidates_add[Ind]
        dist = candidates_add_dist[Ind]
        candidates_edges[vertex].append(edge)
        candidates_distances[vertex].append(dist)
    loop_over = loop_over[1:]
print("")

# save_list(lst_Hamiltonian, "lst_Hamiltonian.txt")
# save_list(lst_not_Hamiltonian, "lst_not_Hamiltonian.txt")
# save_list(digons, "digons" + suffix)
# save_list_of_lists(digons_neighbors, "digons_neighbors" + suffix)
save_list_of_lists(candidates_edges, "candidates_edges" + suffix[:-4] + "_step_15.txt")
save_list_of_lists(candidates_distances, "candidates_distances" + suffix[:-4] + "_step_15.txt")
#
# log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
# log_file_detailed.flush()
# log_file.write("\n{}".format(compute_time(start_time, True)))
# log_file.flush()
#
# # with open("Matrix_out" + suffix[:-4] + "_step_10.txt", 'r') as f:
# #     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# # for i in range(len(Matrix_out)):
# #     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed
# #
# # with open("Register_Added" + suffix[:-4] + "_step_10.txt", 'r') as f:
# #     Register_Added = [[int(num) for num in line.split(',')] for line in f]
# #
# # with open("Register_Removed" + suffix[:-4] + "_step_10.txt", 'r') as f:
# #     Register_Removed = [[int(num) for num in line.split(',')] for line in f]
# #
# # with open("candidates_distances" + suffix[:-4] + "_step_11.txt", 'r') as f:
# #     candidates_distances = [[float(num) for num in line.split(',')] for line in f]
# # for i in range(len(candidates_distances)):
# #     candidates_distances[i] = [x for x in candidates_distances[i] if x != -1]  # -1's get removed
# #
# # candidates_edges = read_list_of_lists_of_lists("candidates_edges" + suffix[:-4] + "_step_11.txt")
# #
# # with open("digons" + suffix, 'r') as f:
# #     digons = [int(elem) for elem in f.read().split('\n') if elem]
# # ignore = union(ignore, digons)
#
start_time = time()
print("\n\nAdd short candidate edges to repair non-Hamiltonian links...")
#log_file_detailed.write("\n\nAdd short candidate edges to repair non-Hamiltonian links...")
log_file_detailed.flush()
log_file.write("\n\nAdd short candidate edges to repair non-Hamiltonian links...")
log_file.flush()
while True:
    minima = [min(candidates_distances[vertex], default=math.inf) for vertex in range(numbPoints)]
    if min(minima) > 10 * threshold:  # == math.inf:
        break
    vertex = minima.index(min(minima))
    Ind = candidates_distances[vertex].index(min(candidates_distances[vertex]))
    edge = candidates_edges[vertex][Ind]
    Register_Added.append(edge)
    Matrix_out = Add(edge[0], edge[1], Matrix_out)
    Register_Removed = Remove_index(Register_Removed, ind(edge, Register_Removed))
    Register_Added.append(edge)
    dist = candidates_distances[vertex][Ind]
    print("Repair the link of {} by adding edge {} with length {}".format(vertex, edge, dist))
    # log_file_detailed.write("\nRepair the link of {} by adding edge {} with length {}".format(vertex, edge, dist))
    # log_file_detailed.flush()
    for v_sub in edge + mutual_neighbors(edge, Matrix_out):
        if v_sub in union(list_boundary, list_boundary_vox, ignore):
            continue
        candidates_edges[v_sub] = []
        candidates_distances[v_sub] = []
        lst_Hamiltonian = setdiff(lst_Hamiltonian, v_sub)
        lst_not_Hamiltonian = setdiff(lst_not_Hamiltonian, v_sub)
        Matrix_reduced, triangles = reducedvertexlink(v_sub, Matrix_out, list_boundary, Centroid, return_Matrix=True,
                                                      return_triangles=True)[1:]
        tets = [[v_sub] + item for item in triangles]
        if isHamiltonian(v_sub, Matrix_reduced, list_boundary, Centroid, predetermined_tetrahedra=tets)[1]:
            lst_Hamiltonian = union(lst_Hamiltonian, v_sub)
            continue
        else:
            if isHamiltonian(v_sub, Matrix_out, list_boundary, Centroid)[1]:
                lst_Hamiltonian = union(lst_Hamiltonian, v_sub)
                continue
            else:
                lst_not_Hamiltonian = union(lst_not_Hamiltonian, v_sub)
        if len(Matrix_out[v_sub]) < 6:
            loop_over = loop_over[1:]
            continue  # this vertex might be a digon
        if len(Matrix_out[v_sub]) > 40:
            loop_over = loop_over[1:]
            continue
        second_neighbors = union(*[Matrix_out[item] for item in Matrix_out[v_sub]])
        initial_Hamiltonian = intersect(second_neighbors, lst_Hamiltonian)
        candidates_Hamiltonian, candidates_add, candidates_add_dist = [], [], []
        for i in Matrix_out[v_sub]:
            for j in [item for item in setdiff(Matrix_out[v_sub], Matrix_out[i]) if item > i]:
                edge = [i, j]
                dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
                if dist > 10 * threshold:
                    continue
                Matrix_temp = Add(edge[0], edge[1], Matrix_out)
                Matrix_reduced, triangles = reducedvertexlink(v_sub, Matrix_temp, list_boundary, Centroid,
                                                              return_Matrix=True,
                                                              return_triangles=True)[1:]
                tets = [[v_sub] + item for item in triangles]
                check = isHamiltonian(v_sub, Matrix_reduced, list_boundary, Centroid,
                                      predetermined_tetrahedra=tets)[1]
                if not check:
                    check = isHamiltonian(v_sub, Matrix_temp, list_boundary, Centroid)[1]
                if check:
                    candidates_add.append(edge)
                    candidates_add_dist.append(dist)
                    candidates_Hamiltonian.append([v_sub])
                    for v in setdiff(second_neighbors, v_sub):
                        if v in union(list_boundary, list_boundary_vox, ignore):
                            continue
                        if v in lst_Hamiltonian and v not in edge:
                            candidates_Hamiltonian[-1].append(v)
                            continue
                        if v not in union(edge, mutual_neighbors(edge, Matrix_out)):
                            continue
                        Matrix_reduced, triangles = reducedvertexlink(v, Matrix_temp, list_boundary, Centroid,
                                                                      return_Matrix=True,
                                                                      return_triangles=True)[1:]
                        tets = [[v] + item for item in triangles]
                        if isHamiltonian(v, Matrix_reduced, list_boundary, Centroid, predetermined_tetrahedra=tets)[1]:
                            candidates_Hamiltonian[-1].append(v)
                            continue
                        if isHamiltonian(v, Matrix_temp, list_boundary, Centroid)[1]:
                            candidates_Hamiltonian[-1].append(v)
        diff = [len(item) - len(initial_Hamiltonian) for item in candidates_Hamiltonian]
        m = max(diff, default=0)
        if m > 0:
            Ind = [index for index in range(len(candidates_Hamiltonian)) if diff[index] == m]
            for index in setdiff(range(len(candidates_add_dist)), Ind):
                candidates_add_dist[index] = math.inf
            Ind = candidates_add_dist.index(min(candidates_add_dist))  # add the shortest edge
            edge = candidates_add[Ind]
            dist = candidates_add_dist[Ind]
            candidates_edges[v_sub].append(edge)
            candidates_distances[v_sub].append(dist)

save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_16.txt")
save_list_of_lists(Register_Added, "Register_Added" + suffix[:-4] + "_step_16.txt")
sys.exit()
# avg_deg = 0
# for vertex in list_interior:
#     avg_deg += len(Matrix_out[vertex])
# avg_deg += sum([len(item) for item in digons_neighbors])
# print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
# log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
# log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
# log_file_detailed.flush()
# log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
# log_file.write("\n{}".format(compute_time(start_time, True)))
# log_file.flush()
#
# # Add exceptionally long edges...
# start_time = time()
# print("\n\nAdd exceptionally long edges...")
# log_file_detailed.write("\n\nAdd exceptionally long edges...")
# log_file_detailed.flush()
# log_file.write("\n\nAdd exceptionally long edges...")
# log_file.flush()
# print("Compute the candidate edges...")
# log_file_detailed.write("\nCompute the candidate edges...")
# log_file_detailed.flush()
# log_file.write("\nCompute the candidate edges...")
# log_file.flush()
#
# candidates_edges = [[] for item in range(numbPoints)]
# candidates_distances = [[] for item in range(numbPoints)]
# lst_compute = range(numbPoints)
# while lst_compute:
#     vertex = lst_compute[0]
#     candidates_edges[vertex] = []
#     candidates_distances[vertex] = []
#     print("\rCompute the candidates:", vertex, "out of", numbPoints - 1, end="")
#     if vertex in union(list_boundary, list_boundary_vox, ignore):
#         lst_compute = lst_compute[1:]
#         continue
#     Holes = holes(vertex, Matrix_out, list_boundary, Centroid, True)[0]
#     Matrix_reduced = reducedvertexlink(vertex, Matrix_out, list_boundary, Centroid, return_Matrix=True)[1]
#     Holes += [hole for hole in holes(vertex, Matrix_reduced, list_boundary, Centroid, True)[0] if
#               not isline_exact(hole, Holes)[0]]
#     scale = 1
#     if not Holes:
#         Holes = pseudoholes(vertex, Matrix_out, list_boundary, Centroid)
#         scale = 0.25
#     outer_vertices, other_edges = [], []
#     for i in range(len(Holes)):
#         # boundary_hole = len(intersect(Holes[i], union(list_boundary, list_boundary_vox))) >= len(Holes[i]) - 1
#         for j in combinations(range(len(Holes[i])), 2):
#             j = list(j)
#             edge = [Holes[i][item] for item in j]
#             if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
#                 continue
#             dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
#             if dist > scale * 20 * threshold:
#                 continue
#             if not isline(edge, candidates_edges[vertex]):
#                 candidates_edges[vertex].append(edge)
#                 candidates_distances[vertex].append(dist)
#                 for k in mutual_neighbors(edge, Matrix_out):
#                     if k not in [vertex] + Matrix_out[vertex]:
#                         outer_vertices = union(outer_vertices, k)
#                     else:
#                         for l in Holes[i]:
#                             if k not in Matrix_out[l] and not isline([k, l], other_edges) and not isline([k, l], Holes[i]):
#                                 dist = distance_vox(k, l, Centroid, coordinates, voxels, dist_path)
#                                 if dist > 1.5 * threshold:
#                                     continue
#                                 other_edges.append([k, l])
#     for v in outer_vertices:
#         edge = [vertex, v]
#         dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
#         if dist > scale * 20 * threshold:
#             continue
#         if not isline(edge, candidates_edges[vertex]):
#             candidates_edges[vertex].append(edge)
#             candidates_distances[vertex].append(dist)
#     for edge in other_edges:
#         dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
#         if dist > scale * 20 * threshold:
#             continue
#         if not isline(edge, candidates_edges[vertex]):
#             candidates_edges[vertex].append(edge)
#             candidates_distances[vertex].append(dist)
#     lst_compute = lst_compute[1:]
# print("")
#
# save_list_of_lists(candidates_edges, "candidates_edges" + suffix[:-4] + "_step_17.txt")
# save_list_of_lists(candidates_distances, "candidates_distances" + suffix[:-4] + "_step_17.txt")
#
# start_time = time()
# print("Add (probably long) candidate edges to close remaining holes...")
# log_file_detailed.write("\nAdd (probably long) candidate edges to close remaining holes...")
# log_file_detailed.flush()
# log_file.write("\nAdd (probably long) candidate edges to close remaining holes...")
# log_file.flush()
# count = 0
# while True:
#     minima = [min(candidates_distances[vertex], default=math.inf) for vertex in range(numbPoints)]
#     if min(minima) > 20 * threshold:  # == math.inf:
#         break
#     vertex = minima.index(min(minima))
#     Ind = candidates_distances[vertex].index(min(candidates_distances[vertex]))
#     edge = candidates_edges[vertex][Ind]
#     Register_Removed = Remove_index(Register_Removed, ind(edge, Register_Removed))
#     Register_Added.append(edge)
#     Matrix_out = Add(edge[0], edge[1], Matrix_out)
#     print("Added:", vertex, edge, candidates_distances[vertex][Ind])
#     log_file_detailed.write("\nAdded: linkvertex {}, edge {} with distance {}".format(vertex, edge, candidates_distances[vertex][Ind]))
#     log_file_detailed.flush()
#     All = union(edge, mutual_neighbors(edge, Matrix_out))
#     for v_sub in All:
#         if v_sub in union(list_boundary, list_boundary_vox) or vertex in ignore:
#             continue
#         candidates_edges[v_sub] = []
#         candidates_distances[v_sub] = []
#         Holes = holes(v_sub, Matrix_out, list_boundary, Centroid, True)[0]
#         Matrix_reduced = reducedvertexlink(v_sub, Matrix_out, list_boundary, Centroid, return_Matrix=True)[1]
#         Holes += [hole for hole in holes(v_sub, Matrix_reduced, list_boundary, Centroid, True)[0] if
#                   not isline_exact(hole, Holes)[0]]
#         scale = 1
#         if not Holes:
#             Holes = pseudoholes(v_sub, Matrix_out, list_boundary, Centroid)
#             scale = 0.25
#         outer_vertices, other_edges = [], []
#         for i in range(len(Holes)):
#             # boundary_hole = len(intersect(Holes[i], union(list_boundary, list_boundary_vox))) >= len(Holes[i]) - 1
#             for j in combinations(range(len(Holes[i])), 2):
#                 j = list(j)
#                 edge = [Holes[i][item] for item in j]
#                 if edge[0] in Matrix_out[edge[1]] or edge[0] == edge[1]:
#                     continue
#                 dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
#                 if dist > scale * 20 * threshold:
#                     continue
#                 if not isline(edge, candidates_edges[v_sub]):
#                     candidates_edges[v_sub].append(edge)
#                     candidates_distances[v_sub].append(dist)
#                     for k in mutual_neighbors(edge, Matrix_out):
#                         if k not in [v_sub] + Matrix_out[v_sub]:
#                             outer_vertices = union(outer_vertices, k)
#                         else:
#                             for l in Holes[i]:
#                                 if k not in Matrix_out[l] and not isline([k, l], other_edges) and not isline([k, l], Holes[i]):
#                                     dist = distance_vox(k, l, Centroid, coordinates, voxels, dist_path)
#                                     if dist > 1.5 * threshold:
#                                         continue
#                                     other_edges.append([k, l])
#         for v in outer_vertices:
#             edge = [v_sub, v]
#             dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
#             if dist > scale * 20 * threshold:
#                 continue
#             if not isline(edge, candidates_edges[v_sub]):
#                 candidates_edges[v_sub].append(edge)
#                 candidates_distances[v_sub].append(dist)
#         for edge in other_edges:
#             dist = distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path)
#             if dist > scale * 20 * threshold:
#                 continue
#             if not isline(edge, candidates_edges[v_sub]):
#                 candidates_edges[v_sub].append(edge)
#                 candidates_distances[v_sub].append(dist)
#
# save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_18.txt")
# save_list_of_lists(Register_Added, "Register_Added" + suffix[:-4] + "_step_18.txt")
# save_list_of_lists(Register_Removed, "Register_Removed" + suffix[:-4] + "_step_18.txt")
#
# avg_deg = 0
# for vertex in list_interior:
#     avg_deg += len(Matrix_out[vertex])
# avg_deg += sum([len(item) for item in digons_neighbors])
# print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
# log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
# log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
# log_file_detailed.flush()
# log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
# log_file.write("\n{}".format(compute_time(start_time, True)))
# log_file.flush()
#
# remaining_interior_holes = []
# for vertex in range(numbPoints):
#     if vertex in union(list_boundary, list_boundary_vox, ignore):
#         continue
#     Holes = holes(vertex, Matrix_out, list_boundary, Centroid, True)[0]
#     if Holes:
#         remaining_interior_holes.append(vertex)
#         continue
#     Holes = pseudoholes(vertex, Matrix_out, list_boundary, Centroid)
#     if Holes:
#         remaining_interior_holes.append(vertex)
#         continue
#     Matrix_reduced = reducedvertexlink(vertex, Matrix_out, list_boundary, Centroid, return_Matrix=True)[1]
#     Holes = holes(vertex, Matrix_reduced, list_boundary, Centroid, True)[0]
#     if Holes:
#         remaining_interior_holes.append(vertex)
# print("{} remaining holes".format(len(remaining_interior_holes)))
# log_file_detailed.write("\n\n{} interior vertices might still have holes in their links".format(len(remaining_interior_holes)))
# log_file_detailed.flush()
# log_file.write("\n{} interior vertices might still have holes in their links".format(len(remaining_interior_holes)))
# log_file.flush()
#
# save_list(remaining_interior_holes, "remaining_interior_holes" + suffix[:-4] + "_step_19.txt")
#
# log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
# log_file_detailed.flush()
# log_file.write("\n{}".format(compute_time(start_time, True)))
# log_file.flush()
#
# K6 = []
# for i in range(numbPoints):
#     Nei_i = Matrix_out[i]
#     for j in [item for item in Nei_i if item > i]:
#         Nei_i_j = intersect(Nei_i, Matrix_out[j])
#         for k in [item for item in Nei_i_j if item > j]:
#             Nei_i_j_k = intersect(Nei_i_j, Matrix_out[k])
#             for l in [item for item in Nei_i_j_k if item > k]:
#                 Nei_i_j_k_l = intersect(Nei_i_j_k, Matrix_out[l])
#                 for m in [item for item in Nei_i_j_k_l if item > l]:
#                     Nei_i_j_k_l_m = intersect(Nei_i_j_k_l, Matrix_out[m])
#                     for n in [item for item in Nei_i_j_k_l_m if item > m]:
#                         K6.append([i, j, k, l, m, n])
#
# # Remove excessive edges from 6-cliques...
# start_time = time()
# print("\n\nRemove excessive edges from 6-cliques...")
# log_file_detailed.write("\n\nRemove excessive edges from 6-cliques...")
# log_file_detailed.flush()
# log_file.write("\n\nRemove excessive edges from 6-cliques...")
# log_file.flush()
# print("Compute the candidate edges...")
# log_file_detailed.write("\nCompute the candidate edges...")
# log_file_detailed.flush()
# log_file.write("\nCompute the candidate edges...")
# log_file.flush()
#
# candidates_edges = [[] for item in K6]
# candidates_distances = [[] for item in K6]
# for i in range(len(K6)):
#     print("\rCompute the candidates:", i + 1, "out of", len(K6), end="")
#     lst = K6[i]
#     Holes_all = []
#     for v in lst:
#         Holes = holes(v, Matrix_out, list_boundary, Centroid, True)[0]
#         Matrix_reduced = reducedvertexlink(v, Matrix_out, list_boundary, Centroid, return_Matrix=True)[1]
#         Holes += [hole for hole in holes(v, Matrix_reduced, list_boundary, Centroid, True)[0] if not isline_exact(hole, Holes)[0]]
#         Holes_all.append(Holes)
#     edges, distances = [], []
#     for edge in combinations(lst, 2):
#         lst_sub = mutual_neighbors(edge, Matrix_out)
#         lst_sub = setdiff(lst_sub, union(list_boundary, list_boundary_vox))
#         Matrix_temp = Remove(edge[0], edge[1], Matrix_out)
#         removable = True
#         for v in lst_sub:
#             if isHamiltonian(v, Matrix_out, list_boundary, Centroid)[1] and not isHamiltonian(v, Matrix_temp,
#                                                                                            list_boundary, Centroid)[1]:
#                 removable = False
#                 break
#             for j in range(len(lst)):
#                 v = lst[j]
#                 if isline(v, union(list_boundary, list_boundary_vox, ignore)):
#                     continue
#                 Holes = holes(v, Matrix_temp, list_boundary, Centroid, True)[0]
#                 Matrix_reduced = reducedvertexlink(v, Matrix_temp, list_boundary, Centroid, return_Matrix=True)[1]
#                 Holes += [hole for hole in holes(v, Matrix_reduced, list_boundary, Centroid, True)[0] if
#                           not isline_exact(hole, Holes)[0]]
#                 if any([not isline_exact(hole, Holes_all[j])[0] for hole in Holes]):
#                     removable = False
#                     break
#             if not removable:
#                 break
#         if removable:
#             edges.append(edge)
#             distances.append(distance_vox(edge[0], edge[1], Centroid, coordinates, voxels, dist_path))
#     if distances:
#         Ind = distances.index(max(distances))
#         if distances[Ind] < threshold:
#             candidates_distances[i].append(-math.inf)
#             continue
#         candidates_edges[i].append(list(edges[Ind]))
#         candidates_distances[i].append(distances[Ind])
#     else:
#         candidates_distances[i].append(-math.inf)
# print("")
#
# save_list_of_lists(candidates_edges, "candidates_edges" + suffix[:-4] + "_step_20.txt")
# save_list_of_lists(candidates_distances, "candidates_distances" + suffix[:-4] + "_step_20.txt")
#
# print("Remove excessive edges...")
# log_file_detailed.write("\nRemove excessive edges...")
# log_file_detailed.flush()
# log_file.write("\nRemove excessive edges...")
# log_file.flush()
# count = 0
# while True:
#     lst = [candidates_distances[item][0] for item in range(len(candidates_distances))]
#     m = max(lst)
#     if m == -math.inf:
#         break
#     Ind = lst.index(m)
#     edge = candidates_edges[Ind][0]
#     Register_Added = Remove_index(Register_Added, ind(edge, Register_Added))
#     Register_Removed.append(edge)
#     Matrix_out = Remove(edge[0], edge[1], Matrix_out)
#     print("Removed from 6-clique:", K6[Ind], edge, candidates_distances[Ind][0])
#     log_file_detailed.write("\nRemoved from 6-clique {} - edge {} with distance {}".format(K6[Ind], edge, candidates_distances[Ind][0]))
#     log_file_detailed.flush()
#     Ind = ind(edge, K6)
#     for i in Ind:
#         candidates_edges[i] = []
#         candidates_distances[i] = [-math.inf]
#
# save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_21.txt")
# save_list_of_lists(Register_Added, "Register_Added" + suffix[:-4] + "_step_21.txt")
# save_list_of_lists(Register_Removed, "Register_Removed" + suffix[:-4] + "_step_21.txt")
#
# avg_deg = 0
# for vertex in list_interior:
#     avg_deg += len(Matrix_out[vertex])
# print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
# avg_deg += sum([len(item) for item in digons_neighbors])
# log_file_detailed.write("\n\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
# log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
# log_file_detailed.flush()
# log_file.write("\nAverage face degree = {}".format(round(avg_deg / len(list_interior), 3)))
# log_file.write("\n{}".format(compute_time(start_time, True)))
# log_file.flush()
#
# start_time = time()
# print("\n\nCompute the list of interior vertices with repairable or proper links...")
# log_file_detailed.write("\n\nCompute the list of interior vertices with repairable or proper links...")
# log_file_detailed.flush()
# log_file.write("\n\nCompute the list of interior vertices with repairable or proper links...")
# log_file.flush()
# lst_repairable = [False for item in range(numbPoints)]
# lst_proper = [False for item in range(numbPoints)]
# for vertex in range(numbPoints):
#     print("\rCompute the list of vertices with repairable or proper links:", vertex, "out of", numbPoints - 1,
#           end="")
#     if repairable_link(vertex, Matrix_out, list_boundary, Centroid):
#         lst_repairable[vertex] = True
#     if propervertexlink(vertex, Matrix_out, list_boundary, Centroid):
#         lst_proper[vertex] = True
# print("")
# print("# repairable:", len([item for item in list_interior if lst_repairable[item]]), "out of", len(list_interior))
# print("# proper:", len([item for item in list_interior if lst_proper[item]]), "out of", len(list_interior))
# log_file_detailed.write("\n# repairable: {} out of {}".format(len([item for item in list_interior if lst_repairable[item]]), len(list_interior)))
# log_file_detailed.write("\n# proper: {} out of {}".format(len([item for item in list_interior if lst_proper[item]]), len(list_interior)))
# log_file_detailed.flush()
# log_file.write("\n# repairable: {} out of {}".format(len([item for item in list_interior if lst_repairable[item]]), len(list_interior)))
# log_file.write("\n# proper: {} out of {}".format(len([item for item in list_interior if lst_proper[item]]), len(list_interior)))
# log_file.flush()
#
# save_list(lst_repairable, "lst_repairable" + suffix[:-4] + "_step_22.txt")
# save_list(lst_proper, "lst_proper" + suffix[:-4] + "_step_22.txt")
#
# # with open("Matrix_out" + suffix[:-4] + "_step_21.txt", 'r') as f:
# #     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# # for i in range(len(Matrix_out)):
# #     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed
# #
# # with open("Register_Added" + suffix[:-4] + "_step_21.txt", 'r') as f:
# #     Register_Added = [[int(num) for num in line.split(',')] for line in f]
# #
# # with open("Register_Removed" + suffix[:-4] + "_step_21.txt", 'r') as f:
# #     Register_Removed = [[int(num) for num in line.split(',')] for line in f]
# #
# # with open("lst_repairable" + suffix[:-4] + "_step_22.txt", 'r') as f:
# #     lst_repairable = [[str_2_bool(num) for num in line.split('\n')][0] for line in f]
# #
# # with open("lst_proper" + suffix[:-4] + "_step_22.txt", 'r') as f:
# #     lst_proper = [[str_2_bool(num) for num in line.split('\n')][0] for line in f]
# #
# # with open("digons" + suffix, 'r') as f:
# #     digons = [int(elem) for elem in f.read().split('\n') if elem]
# # ignore = union(ignore, digons)
#
# K5 = find_K5s(Matrix_out)
#
# remaining = []
# for vertex in [item for item in list_interior if not lst_repairable[item] and not lst_proper[item]]:
#     if vertex not in union(list_boundary_vox, list_boundary):
#         remaining.append(vertex)
#         continue
#     _, _, problems = istriangulatable(vertex, Matrix_out, list_boundary, Centroid)
#     if setdiff(problems, union(list_boundary_vox, list_boundary)):
#         remaining.append(vertex)
#
# log_file_detailed.write("\n\nChange tunnels of 5-cliques...")
# log_file_detailed.flush()
# log_file.write("\n\nChange tunnels of 5-cliques...")
# log_file.flush()
# print("changing tunnels...")
# # changed_tunnels = [[], []]
# # for vertex in remaining:
# #     K5_temp = [K5[item] for item in ind(vertex, K5)]
# #     T, _, problems = istriangulatable(vertex, Matrix_out, list_boundary, Centroid, changed_tunnels)
# #     if T:
# #         continue
# #     K5_temp = [item for item in K5_temp if intersect(item, problems)]
# #     for K5_current in K5_temp:
# #         tunnels = [list(item) for item in combinations(K5_current, 2)]
# #         K5_current_sub = setdiff(K5_current, index_inf)
# #         diff = []
# #         tunnel_initial = tunnel_K5_geometric(K5_current, Centroid)
# #         Ind = ind(K5_current, changed_tunnels[0])
# #         if Ind:
# #             tunnel_initial = changed_tunnels[1][Ind[0]]
# #         T_repaired = istriangulatable(K5_current_sub, Matrix_out, list_boundary, Centroid, changed_tunnels, try_tunnels=K5_current)[0]
# #         if tunnel_initial in tunnels:
# #             Ind = tunnels.index(tunnel_initial)
# #             initial_repaired = sum([T_repaired[item][Ind] for item in range(len(K5_current_sub))])
# #         else:
# #             initial_repaired = sum(istriangulatable(K5_current_sub, Matrix_out, list_boundary, Centroid, changed_tunnels)[0])
# #         T_repaired = np.array(T_repaired)
# #         diff = np.sum(T_repaired, axis=0) - initial_repaired
# #         m = max(diff)
# #         print(vertex, list(diff), m)
# #         if m <= 0:
# #             continue
# #         print(("Tunnel changed: vertex {}, 5-clique: {}, gains: {}  |  {}".format(vertex, K5_current, list(diff), m)))
# #         log_file_detailed.write("\nTunnel changed: vertex {}, 5-clique: {}, gains: {}  |  {}".format(vertex, K5_current, list(diff), m))
# #         log_file_detailed.flush()
# #         Ind_repaired = [item for item in range(len(diff)) if diff[item] == m]
# #         if len(Ind_repaired) == 1:
# #             tunnel = tunnels[Ind_repaired[0]]
# #         else:
# #             T, _, problems = istriangulatable(K5_current_sub, Matrix_out, list_boundary, Centroid, changed_tunnels, try_tunnels=K5_current, with_repair=False)
# #             num_of_problems = [len(item) for item in problems]
# #             m = min(num_of_problems)
# #             Ind_not_repaired = [item for item in range(len(num_of_problems)) if num_of_problems[item] == m]
# #             Ind = intersect(Ind_repaired, Ind_not_repaired)
# #             if Ind:
# #                 Ind = Ind[0]
# #             else:
# #                 Ind = Ind_repaired[0]
# #             tunnel = tunnels[Ind]
# #
# #         Ind = ind(K5_current, changed_tunnels[0])
# #         if Ind:
# #             Ind = Ind[0]
# #             changed_tunnels[1][Ind] = tunnel
# #         else:
# #             changed_tunnels[0].append(K5_current)
# #             changed_tunnels[1].append(tunnel)
# #         if istriangulatable(vertex, Matrix_out, list_boundary, Centroid, changed_tunnels)[0]:
# #             break
# #
# # save_list_of_lists(changed_tunnels, "changed_tunnels" + suffix[:-4] + "_step_23.txt")
#
# changed_tunnels = [[], []]
# for vertex in remaining:
#     K5_temp = [K5[item] for item in ind(vertex, K5)]
#     T, _, problems = istriangulatable(vertex, Matrix_out, list_boundary, Centroid, changed_tunnels)
#     if T:
#         continue
#     K5_temp = [item for item in K5_temp if intersect(item, problems)]
#     for K5_current in K5_temp:
#         tunnels = [list(item) for item in combinations(K5_current, 2)]
#         K5_current_sub = setdiff(K5_current, index_inf)
#         diff = []
#         tunnel_initial = tunnel_K5_geometric(K5_current, Centroid)
#         Ind = ind(K5_current, changed_tunnels[0])
#         if Ind:
#             tunnel_initial = changed_tunnels[1][Ind[0]]
#         T_repaired = istriangulatable(K5_current_sub, Matrix_out, list_boundary, Centroid, changed_tunnels, try_tunnels=K5_current)[0]
#         if tunnel_initial in tunnels:
#             Ind = tunnels.index(tunnel_initial)
#             initial_repaired = sum([T_repaired[item][Ind] for item in range(len(K5_current_sub))])
#         else:
#             initial_repaired = sum(istriangulatable(K5_current_sub, Matrix_out, list_boundary, Centroid, changed_tunnels)[0])
#         T_repaired = np.array(T_repaired)
#         diff = np.sum(T_repaired, axis=0) - initial_repaired
#         m = max(diff)
#         print(vertex, list(diff), m)
#         if m <= 0:
#             continue
#         print(("Tunnel changed: vertex {}, 5-clique: {}, gains: {}  |  {}".format(vertex, K5_current, list(diff), m)))
#         log_file_detailed.write("\nTunnel changed: vertex {}, 5-clique: {}, gains: {}  |  {}".format(vertex, K5_current, list(diff), m))
#         log_file_detailed.flush()
#         Ind_repaired = [item for item in range(len(diff)) if diff[item] == m]
#         if len(Ind_repaired) == 1:
#             tunnel = tunnels[Ind_repaired[0]]
#         else:
#             T, _, problems = istriangulatable(K5_current_sub, Matrix_out, list_boundary, Centroid, changed_tunnels, try_tunnels=K5_current, with_repair=False)
#             num_of_problems = [len(item) for item in problems]
#             m = min(num_of_problems)
#             Ind_not_repaired = [item for item in range(len(num_of_problems)) if num_of_problems[item] == m]
#             Ind = intersect(Ind_repaired, Ind_not_repaired)
#             if Ind:
#                 Ind = Ind[0]
#             else:
#                 Ind = Ind_repaired[0]
#             tunnel = tunnels[Ind]
#
#         Ind = ind(K5_current, changed_tunnels[0])
#         if Ind:
#             Ind = Ind[0]
#             changed_tunnels[1][Ind] = tunnel
#         else:
#             changed_tunnels[0].append(K5_current)
#             changed_tunnels[1].append(tunnel)
#         if istriangulatable(vertex, Matrix_out, list_boundary, Centroid, changed_tunnels)[0]:
#             break
#
# save_list_of_lists(changed_tunnels, "changed_tunnels" + suffix[:-4] + "_step_23_temp.txt")
#
# # changed_tunnels = read_list_of_lists_of_lists("changed_tunnels" + suffix[:-4] + "_step_23.txt")
#
# log_file_detailed.write("\n\nCompute the list of interior vertices with repairable or proper links...")
# log_file_detailed.flush()
# log_file.write("\n\nCompute the list of interior vertices with repairable or proper links...")
# log_file.flush()
# lst_proper = [False for item in range(numbPoints)]
# for vertex in range(numbPoints):
#     print("\rCompute the list of vertices with proper links:", vertex, "out of", numbPoints - 1, end="")
#     if propervertexlink(vertex, Matrix_out, list_boundary, Centroid, changed_tunnels):
#         lst_proper[vertex] = True
# print("")
# print("# repairable:", len([item for item in list_interior if lst_repairable[item]]), "out of", len(list_interior))
# print("# proper:", len([item for item in list_interior if lst_proper[item]]), "out of", len(list_interior))
# log_file_detailed.write("\n# repairable: {} out of {}".format(len([item for item in list_interior if lst_repairable[item]]), len(list_interior)))
# log_file_detailed.write("\n# proper: {} out of {}".format(len([item for item in list_interior if lst_proper[item]]), len(list_interior)))
# log_file_detailed.write("\n{}".format(compute_time(start_time, True)))
# log_file_detailed.flush()
# log_file.write("\n# repairable: {} out of {}".format(len([item for item in list_interior if lst_repairable[item]]), len(list_interior)))
# log_file.write("\n# proper: {} out of {}".format(len([item for item in list_interior if lst_proper[item]]), len(list_interior)))
# log_file.write("\n{}".format(compute_time(start_time, True)))
# log_file.flush()
#
# save_list(lst_repairable, "lst_repairable" + suffix[:-4] + "_step_24.txt")
# save_list(lst_proper, "lst_proper" + suffix[:-4] + "_step_24.txt")
#
# # Reconstruct the triangulation
# start_time = time()
# print("Reconstruct the triangulation at a first try...")
# log_file_detailed.write("\n\nReconstruct the triangulation at a first try...")
# log_file_detailed.flush()
# log_file.write("\n\nReconstruct the triangulation at a first try...")
# log_file.flush()
#
# triangulation_reconstructed = []
# lst_repeat = []
# for vertex in order:
#     initial_len = len(triangulation_reconstructed)
#     T, triangulation_reconstructed, _ = istriangulatable(vertex, Matrix_out, list_boundary, Centroid, changed_tunnels, with_repair=False, predetermined=triangulation_reconstructed)
#     if not T:
#         triangulation_reconstructed = triangulation_reconstructed[:initial_len]
#         lst_repeat.append(vertex)
# for vertex in lst_repeat:
#     _, triangulation_reconstructed, _ = istriangulatable(vertex, Matrix_out, list_boundary, Centroid, changed_tunnels, predetermined=triangulation_reconstructed)
#
# save_list_of_lists(triangulation_reconstructed, "triangulation_test_step_1" + suffix)
#
# # with open("Matrix_out" + suffix[:-4] + "_temp.txt", 'r') as f:
# #     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# # for i in range(len(Matrix_out)):
# #     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed
# #
# # with open("triangulation_test_1.txt") as f:
# #     triangulation_reconstructed = [[int(num) for num in line.split(',')] for line in f]
# #
# # with open("digons" + suffix, 'r') as f:
# #     digons = [int(elem) for elem in f.read().split('\n') if elem]
# # ignore = union(ignore, digons)
#
# edges_defect_link, edges_defect_link_interior = [], []
# count, count_interior = 0, 0
# for i in range(numbPoints):
#     Nei = Matrix_out[i]
#     for j in [item for item in Nei if item > i]:
#         if sublist([i, j], union(list_boundary, list_boundary_vox)):
#             continue
#         if not intersect([i, j], union(list_boundary, list_boundary_vox)):
#             count += 1
#             count_interior += 1
#             if not proper_edge_in_triangulation([i, j], triangulation_reconstructed):
#                 edges_defect_link.append([i, j])
#                 edges_defect_link_interior.append([i, j])
#         else:
#             count += 1
#             if not proper_edge_in_triangulation([i, j], triangulation_reconstructed):
#                 edges_defect_link.append([i, j])
# vts_defect_link = setdiff(unique(edges_defect_link), list_boundary)
# count, count_interior = max(1, count), max(1, count_interior)
#
# print("\nAll non-boundary edges: {}   -   reconstructed at a first try: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# print("All interior edges: {}   -   reconstructed at a first try in the interior: {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
#
# log_file_detailed.write("\nNumber of non-boundary edges: {}, properly reconstructed at a first try: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# log_file_detailed.write("\nNumber of interior edges: {}, properly reconstructed at a first try: {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
# log_file_detailed.flush()
# log_file.write("\nNumber of non-boundary edges: {}, properly reconstructed at a first try: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# log_file.write("\nNumber of interior edges: {}, properly reconstructed at a first try: {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
# log_file.flush()
# log_file_detailed.write("\n\n{}".format(compute_time(start_time, True)))
# log_file_detailed.flush()
# log_file.write("\n{}".format(compute_time(start_time, True)))
# log_file.flush()
#
# start_time = time()
# print("\n\nApply a post process step to try to repair the triangulation...")
# log_file_detailed.write("\n\nApply a post process step to try to repair the triangulation...")
# log_file_detailed.flush()
# log_file.write("\n\nApply a post process step to try to repair the triangulation...")
# log_file.flush()
#
# for vertex in range(numbPoints):
#     if vertex in ignore:
#         continue
#     if proper_vertex_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex]):
#         continue
#     _, triangulation_reconstructed, _ = istriangulatable(vertex, Matrix_out, list_boundary, Centroid, changed_tunnels, predetermined=triangulation_reconstructed, post_process=True)
#
# save_list_of_lists(triangulation_reconstructed, "triangulation_test_step_2" + suffix)
#
# edges_defect_link, edges_defect_link_interior = [], []
# count, count_interior = 0, 0
# for i in range(numbPoints):
#     Nei = Matrix_out[i]
#     for j in [item for item in Nei if item > i]:
#         if sublist([i, j], union(list_boundary, list_boundary_vox)):
#             continue
#         if not intersect([i, j], union(list_boundary, list_boundary_vox)):
#             count += 1
#             count_interior += 1
#             if not proper_edge_in_triangulation([i, j], triangulation_reconstructed):
#                 edges_defect_link.append([i, j])
#                 edges_defect_link_interior.append([i, j])
#         else:
#             count += 1
#             if not proper_edge_in_triangulation([i, j], triangulation_reconstructed):
#                 edges_defect_link.append([i, j])
# vts_defect_link = setdiff(unique(edges_defect_link), list_boundary)
# count, count_interior = max(1, count), max(1, count_interior)
#
# print("\nAll non-boundary edges: {}   -   reconstructed with a post process: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# print("All interior edges: {}   -   reconstructed with a post process in the interior: {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
#
# log_file_detailed.write("\nNumber of non-boundary edges: {}, properly reconstructed with a post process: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# log_file_detailed.write("\nNumber of interior edges: {}, properly reconstructed with a post process: {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
# log_file_detailed.flush()
# log_file.write("\nNumber of non-boundary edges: {}, properly reconstructed with a post process: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# log_file.write("\nNumber of interior edges: {}, properly reconstructed with a post process: {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
# log_file.flush()
# log_file_detailed.write("\n\n{}".format(compute_time(start_time, True)))
# log_file_detailed.flush()
# log_file.write("\n{}".format(compute_time(start_time, True)))
# log_file.flush()
#
# # with open("Matrix_out" + suffix[:-4] + "_step_21.txt", 'r') as f:
# #     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# # for i in range(len(Matrix_out)):
# #     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed
# #
# # with open("Register_Added" + suffix[:-4] + "_step_21.txt", 'r') as f:
# #     Register_Added = [[int(num) for num in line.split(',')] for line in f]
# #
# # with open("Register_Removed" + suffix[:-4] + "_step_21.txt", 'r') as f:
# #     Register_Removed = [[int(num) for num in line.split(',')] for line in f]
# #
# # changed_tunnels = read_list_of_lists_of_lists("changed_tunnels" + suffix[:-4] + "_step_23.txt")
# #
# # with open("triangulation_test_step_2" + suffix, 'r') as f:
# #     triangulation_reconstructed = [[int(num) for num in line.split(',')] for line in f]
# #
# # with open("digons" + suffix, 'r') as f:
# #     digons = [int(elem) for elem in f.read().split('\n') if elem]
# # ignore = union(ignore, digons)
#
# start_time = time()
# print("\n\nChange links of non-proper edges and detect digons and multi-connected cells...")
# log_file_detailed.write("\n\nChange links of non-proper edges and detect digons and multi-connected cells...")
# log_file_detailed.flush()
# log_file.write("\n\nChange links of non-proper edges and detect digons and multi-connected cells...")
# log_file.flush()

# multi_edges = []
# caused_by_digons = []
#
# K2, K4 = [], []
# for i in range(numbPoints):
#     Nei_i = Matrix_out[i]
#     for j in [item for item in Nei_i if item > i]:
#         if not sublist([i, j], union(list_boundary, list_boundary_vox)):
#             K2.append([i, j])
#         Nei_i_j = intersect(Nei_i, Matrix_out[j])
#         for k in [item for item in Nei_i_j if item > j]:
#             Nei_i_j_k = intersect(Nei_i_j, Matrix_out[k])
#             for l in [item for item in Nei_i_j_k if item > k]:
#                 K4.append([i, j, k, l])
#
# empty = []
# processed = []
# count = 0
# loop_over = K2
#
# while loop_over:
#     edge = loop_over[0]
#     count += 1
#     multi = False
#     # print("\rProcessing: {} out of {}".format(count, len(edges_defect_link)), end="")
#     if proper_edge_in_triangulation(edge, triangulation_reconstructed):
#         # if all([proper_vertex_in_triangulation(v, triangulation_reconstructed) for v in setdiff(edge + mutual_neighbors(edge, Matrix_out), union(list_boundary, list_boundary_vox))]):
#         loop_over = loop_over[1:]
#         continue
#     lst = mutual_neighbors(edge, Matrix_out)
#     if len(lst) < 3 or (len(lst) < 6 and not cycles(edges_of(lst, Matrix_out), minimum_allowed=3)):
#         caused_by_digons.append(edge)
#         loop_over = loop_over[1:]
#         continue
#     vts = mutual_neighbors(edge, Matrix_out)
#     edges = edges_of(vts, Matrix_out)
#     Matrix_temp = Remove(edge[0], edge[1], Matrix_out)
#     Holes_1 = holes(edge[0], Matrix_temp, list_boundary, Centroid, True)[0]
#     Holes_1 = [hole for hole in Holes_1 if all([isline(e, edges) for e in edges_of_cycle(hole)])]
#     Matrix_reduced = reducedvertexlink(edge[0], Matrix_temp, list_boundary, Centroid, return_Matrix=True)[1]
#     Holes_1_red = holes(edge[0], Matrix_reduced, list_boundary, Centroid, True)[0]
#     Holes_1_red = [hole for hole in Holes_1_red if all([isline(e, edges) for e in edges_of_cycle(hole)])]
#     Holes_1 += [hole for hole in Holes_1_red if not isline_exact(hole, Holes_1)[0]]
#
#     Holes_2 = holes(edge[1], Matrix_temp, list_boundary, Centroid, True)[0]
#     Holes_2 = [hole for hole in Holes_2 if all([isline(e, edges) for e in edges_of_cycle(hole)])]
#     Matrix_reduced = reducedvertexlink(edge[1], Matrix_temp, list_boundary, Centroid, return_Matrix=True)[1]
#     Holes_2_red = holes(edge[1], Matrix_reduced, list_boundary, Centroid, True)[0]
#     Holes_2_red = [hole for hole in Holes_2_red if all([isline(e, edges) for e in edges_of_cycle(hole)])]
#     Holes_2 += [hole for hole in Holes_2_red if not isline_exact(hole, Holes_2)[0]]
#     if len(Holes_1) > 1 and len(Holes_2) > 1:
#         # consistent_holes = all([isline_EXACT(hole, Holes_2) for hole in Holes_1]) and len(Holes_1) == len(Holes_2)
#         multi = True  # a probable multi edge
#         if len(edges) > 25:  # too large to process
#             if not isline(edge, multi_edges):
#                 multi_edges.append(edge)
#             loop_over = loop_over[1:]
#             continue
#         # try to fix the triangulation
#         edges = edges_of(union(*(Holes_1 + Holes_2)), Matrix_out)
#         C = cycles(edges, minimum_allowed=4, order_matters=True)
#         edges = edges_of(union(*(Holes_1 + Holes_2), edge), Matrix_out)
#         lst_remove = []
#         for i in range(len(C)):
#             if i in lst_remove:
#                 continue
#             for j in range(i + 1, len(C)):
#                 if len(C[i]) != len(C[j]):
#                     continue
#                 if set(C[j]) == set(C[i]) and isline_EXACT(C[i], C[j]):
#                     lst_remove.append(j)
#         C = Remove_index(C, lst_remove)
#         lengths = unique([len(item) for item in C])
#         C_ordered = []
#         for i in lengths:
#             C_ordered.extend([item for item in C if len(item) == i])
#         candidates = []
#         for cycle in C_ordered:
#             Com_1 = components(edge[0], cycle, Matrix_temp)
#             Com_2 = components(edge[1], cycle, Matrix_temp)
#             if len(Com_1) != 1 or len(Com_2) != 1:
#                 continue
#             Matrix_temp_sub = loads(dumps(Matrix_out))
#             for i in Matrix_temp_sub[edge[0]]:
#                 if i in cycle:
#                     continue
#                 Matrix_temp_sub = Remove(i, edge[1], Matrix_temp_sub)
#             Holes = holes(edge[0], Matrix_temp_sub, list_boundary, Centroid, True)[0]
#             # Holes += pseudoholes(edge[0], Matrix_temp_sub, list_boundary, Centroid)
#             # Matrix_reduced = reducedvertexlink(edge[0], Matrix_temp_sub, list_boundary, Centroid, return_Matrix=True)[1]
#             # Holes += holes(edge[0], Matrix_reduced, list_boundary, Centroid, True)[0]
#             Holes = [hole for hole in Holes if all([isline(e, edges) for e in edges_of_cycle(hole)])]
#             if Holes:
#                 continue
#             Matrix_temp_sub = loads(dumps(Matrix_out))
#             for i in Matrix_temp_sub[edge[1]]:
#                 if i in cycle:
#                     continue
#                 Matrix_temp_sub = Remove(i, edge[0], Matrix_temp_sub)
#             Holes = holes(edge[1], Matrix_temp_sub, list_boundary, Centroid, True)[0]
#             # Holes += pseudoholes(edge[1], Matrix_temp_sub, list_boundary, Centroid)
#             # Matrix_reduced = reducedvertexlink(edge[1], Matrix_temp_sub, list_boundary, Centroid, return_Matrix=True)[1]
#             # Holes += holes(edge[1], Matrix_reduced, list_boundary, Centroid, True)[0]
#             Holes = [hole for hole in Holes if all([isline(e, edges) for e in edges_of_cycle(hole)])]
#             if Holes:
#                 continue
#             candidates.append(cycle)
#         lst_remove = []
#         for index in range(len(candidates)):
#             if len(setdiff(range(len(candidates)), lst_remove)) <= 10:
#                 break
#             cycle = candidates[index]
#             d = sum([isline(item, cycle) for item in K4])
#             if sum([isline(edge, cycle) for edge in edges]) - d >= 2 * len(cycle) - 3:  # the hole is already triangulated
#                 lst_remove.append(index)
#                 continue
#             Matrix_temp_sub = loads(dumps(Matrix_out))
#             for i in Matrix_temp_sub[edge[0]]:
#                 if i in cycle:
#                     continue
#                 Matrix_temp_sub = Remove(i, edge[1], Matrix_temp_sub)
#             Centroid_temp = loads(dumps(Centroid))
#             center = [0, 0, 0]
#             for i in cycle:
#                 center += np.array(Centroid[i])
#             center /= len(cycle)
#             center = 2 * center - np.array(Centroid[edge[0]])
#             Centroid_temp[edge[1]] = list(center)
#             Holes = holes(edge[0], Matrix_temp_sub, list_boundary, Centroid_temp, True)[0]
#             # Holes += pseudoholes(edge[0], Matrix_temp_sub, list_boundary, Centroid_temp)
#             # Matrix_reduced = reducedvertexlink(edge[0], Matrix_temp_sub, list_boundary, Centroid, return_Matrix=True)[1]
#             # Holes += holes(edge[0], Matrix_reduced, list_boundary, Centroid, True)[0]
#             Holes = [hole for hole in Holes if all([isline(e, edges) for e in edges_of_cycle(hole)])]
#             if Holes:
#                 lst_remove.append(index)
#                 continue
#             Matrix_temp_sub = loads(dumps(Matrix_out))
#             for i in Matrix_temp_sub[edge[1]]:
#                 if i in cycle:
#                     continue
#                 Matrix_temp_sub = Remove(i, edge[0], Matrix_temp_sub)
#             Centroid_temp = loads(dumps(Centroid))
#             center = [0, 0, 0]
#             for i in cycle:
#                 center += np.array(Centroid[i])
#             center /= len(cycle)
#             center = 2 * center - np.array(Centroid[edge[1]])
#             Centroid_temp[edge[0]] = list(center)
#             Holes = holes(edge[1], Matrix_temp_sub, list_boundary, Centroid_temp, True)[0]
#             # Holes += pseudoholes(edge[1], Matrix_temp_sub, list_boundary, Centroid_temp)
#             # Matrix_reduced = reducedvertexlink(edge[1], Matrix_temp_sub, list_boundary, Centroid, return_Matrix=True)[1]
#             # Holes += holes(edge[1], Matrix_reduced, list_boundary, Centroid, True)[0]
#             Holes = [hole for hole in Holes if all([isline(e, edges) for e in edges_of_cycle(hole)])]
#             if Holes:
#                 lst_remove.append(index)
#                 continue
#         candidates = Remove_index(candidates, lst_remove)
#         print("\nnumber of candidates for {} = {}".format(edge, len(candidates)))
#         if not candidates or len(candidates) > 20:
#             if multi and not isline(edge, multi_edges):
#                 multi_edges.append(edge)
#             loop_over = loop_over[1:]
#             continue
#     else:
#         if len(edges) > 25:  # too large to process
#             loop_over = loop_over[1:]
#             continue
#         if Holes_1 + Holes_2:
#             C = cycles(edges, minimum_allowed=4, order_matters=True)
#         else:
#             C = cycles(edges, minimum_allowed=3, order_matters=True)
#         edges = edges_of(union(*(Holes_1 + Holes_2), edge), Matrix_out)
#         lengths = unique([len(item) for item in C])
#         C_ordered = []
#         for i in lengths:
#             C_ordered.extend([item for item in C if len(item) == i])
#         candidates = []
#         for cycle in C_ordered:
#             Com_1 = components(edge[0], cycle, Matrix_temp)
#             Com_2 = components(edge[1], cycle, Matrix_temp)
#             if len(Com_1) != 1 or len(Com_2) != 1:
#                 continue
#             Matrix_temp_sub = loads(dumps(Matrix_out))
#             for i in Matrix_temp_sub[edge[0]]:
#                 if i in cycle:
#                     continue
#                 Matrix_temp_sub = Remove(i, edge[1], Matrix_temp_sub)
#             Holes = holes(edge[0], Matrix_temp_sub, list_boundary, Centroid, True)[0]
#             # Holes += pseudoholes(edge[0], Matrix_temp_sub, list_boundary, Centroid)
#             # Matrix_reduced = reducedvertexlink(edge[0], Matrix_temp_sub, list_boundary, Centroid, return_Matrix=True)[1]
#             # Holes += holes(edge[0], Matrix_reduced, list_boundary, Centroid, True)[0]
#             Holes = [hole for hole in Holes if all([isline(e, edges) for e in edges_of_cycle(hole)])]
#             if Holes:
#                 continue
#             Matrix_temp_sub = loads(dumps(Matrix_out))
#             for i in Matrix_temp_sub[edge[1]]:
#                 if i in cycle:
#                     continue
#                 Matrix_temp_sub = Remove(i, edge[0], Matrix_temp_sub)
#             Holes = holes(edge[1], Matrix_temp_sub, list_boundary, Centroid, True)[0]
#             # Holes += pseudoholes(edge[1], Matrix_temp_sub, list_boundary, Centroid)
#             # Matrix_reduced = reducedvertexlink(edge[1], Matrix_temp_sub, list_boundary, Centroid, return_Matrix=True)[1]
#             # Holes += holes(edge[1], Matrix_reduced, list_boundary, Centroid, True)[0]
#             Holes = [hole for hole in Holes if all([isline(e, edges) for e in edges_of_cycle(hole)])]
#             if Holes:
#                 continue
#             candidates.append(cycle)
#         # print("\n", candidates)
#         lst_remove = []
#         for index in range(len(candidates)):
#             if len(setdiff(range(len(candidates)), lst_remove)) <= 10:
#                 break
#             cycle = candidates[index]
#             d = sum([isline(item, cycle) for item in K4])
#             if len(cycle) > 3 and sum([isline(edge, cycle) for edge in edges]) - d >= 2 * len(cycle) - 3:  # the hole is already triangulated
#                 lst_remove.append(index)
#             Matrix_temp_sub = loads(dumps(Matrix_out))
#             for i in Matrix_temp_sub[edge[0]]:
#                 if i in cycle:
#                     continue
#                 Matrix_temp_sub = Remove(i, edge[1], Matrix_temp_sub)
#             Centroid_temp = loads(dumps(Centroid))
#             center = [0, 0, 0]
#             for i in cycle:
#                 center += np.array(Centroid[i])
#             center /= len(cycle)
#             center = 2 * center - np.array(Centroid[edge[0]])
#             Centroid_temp[edge[1]] = list(center)
#             Holes = holes(edge[0], Matrix_temp_sub, list_boundary, Centroid_temp, True)[0]
#             # Holes += pseudoholes(edge[0], Matrix_temp_sub, list_boundary, Centroid_temp)
#             # Matrix_reduced = reducedvertexlink(edge[0], Matrix_temp_sub, list_boundary, Centroid, return_Matrix=True)[1]
#             # Holes += holes(edge[0], Matrix_reduced, list_boundary, Centroid, True)[0]
#             Holes = [hole for hole in Holes if all([isline(e, edges) for e in edges_of_cycle(hole)])]
#             if Holes:
#                 lst_remove.append(index)
#                 continue
#             Matrix_temp_sub = loads(dumps(Matrix_out))
#             for i in Matrix_temp_sub[edge[1]]:
#                 if i in cycle:
#                     continue
#                 Matrix_temp_sub = Remove(i, edge[0], Matrix_temp_sub)
#             Centroid_temp = loads(dumps(Centroid))
#             center = [0, 0, 0]
#             for i in cycle:
#                 center += np.array(Centroid[i])
#             center /= len(cycle)
#             center = 2 * center - np.array(Centroid[edge[1]])
#             Centroid_temp[edge[0]] = list(center)
#             Holes = holes(edge[1], Matrix_temp_sub, list_boundary, Centroid_temp, True)[0]
#             # Holes += pseudoholes(edge[1], Matrix_temp_sub, list_boundary, Centroid_temp)
#             # Matrix_reduced = reducedvertexlink(edge[1], Matrix_temp_sub, list_boundary, Centroid, return_Matrix=True)[1]
#             # Holes += holes(edge[1], Matrix_reduced, list_boundary, Centroid, True)[0]
#             Holes = [hole for hole in Holes if all([isline(e, edges) for e in edges_of_cycle(hole)])]
#             if Holes:
#                 lst_remove.append(index)
#                 continue
#         candidates = Remove_index(candidates, lst_remove)
#         print("\nnumber of candidates for {} = {}".format(edge, len(candidates)))
#         if not candidates or len(candidates) > 20:
#             candidates = Holes_1 + Holes_2
#             if candidates:
#                 print("number of candidates reduced to", len(Holes_1 + Holes_2))
#         if not candidates:
#             loop_over = loop_over[1:]
#             continue
#     # we consider the second neighbors
#     lst_vertices = union(Matrix_out[edge[0]], Matrix_out[edge[1]])
#     lst_vertices = unique([Matrix_out[item] for item in lst_vertices])
#     lst_edges = edges_of(lst_vertices, Matrix_out)
#     # print("all:", len(lst))
#     proper_edges_initial = [edge for edge in lst_edges if proper_edge_in_triangulation(edge, triangulation_reconstructed)]
#     proper_vertices_initial = [v for v in lst_vertices if proper_vertex_in_triangulation(v, triangulation_reconstructed, Matrix_out[v])]
#     proper_vertices_initial_int = setdiff(proper_vertices_initial, list_boundary)
#     # print("initial:", len(lst_initial))
#     gains_edges, gains_vertices, gains_vertices_int, retriangulations, lsts, candidates_unfold, emptys, losses = [], [], [], [], [], [], [], []
#     for cycle in candidates:
#         triangulation_test = loads(dumps(triangulation_reconstructed))
#         Ind = []
#         for u in mutual_neighbors(edge, Matrix_out):
#             Ind = union(Ind, ind([edge[0], u], triangulation_test))
#             Ind = union(Ind, ind([edge[1], u], triangulation_test))
#         for edge_sub in combinations(cycle, 2):
#             Ind = union(Ind, ind(edge_sub, triangulation_test))
#         triangulation_test = Remove_index(triangulation_test, Ind)
#         Centroid_temp_0, Centroid_temp_1 = loads(dumps(Centroid)), loads(dumps(Centroid))
#         center = [0, 0, 0]
#         for i in cycle:
#             center += np.array(Centroid[i])
#         center /= len(cycle)
#         center_0 = 2 * center - np.array(Centroid[edge[1]])
#         center_1 = 2 * center - np.array(Centroid[edge[0]])
#         Centroid_temp_0[edge[1]] = list(center_1)
#         Centroid_temp_1[edge[0]] = list(center_0)
#
#         edges = edges_of_cycle(cycle)
#         star = [edge + item for item in edges]
#         empty_sub = [K4[item] for item in ind(edge, K4) if not isline(K4[item], star)]
#         empty_test = empty + [item for item in empty_sub if not isline(item, empty)]
#         triangles = []
#         # if len(cycle) > 3:
#         #     triangles = cycles_deg(edges_of(cycle, Matrix_out), 3)
#         if triangles and len(triangles) < 4:
#             for num in range(2 ** len(triangles)):
#                 candidates_unfold.append(cycle)
#                 pos = bin(num)[2:].zfill(len(triangles))
#                 empty_test_sub = loads(dumps(empty_test))
#                 retriangulation = triangulation_test + star
#                 # emp = []
#                 for index_sub in range(len(pos)):
#                     if pos[index_sub] == "1":
#                         empty_test_sub += [item for item in [[edge[0]] + triangles[index_sub], [edge[1]] + triangles[index_sub]] if not isline(item, empty_test_sub)]
#                         # emp.append(triangles[index_sub])
#                 emptys.append(empty_test_sub)
#                 for v_sub in union(unique(candidates), edge):
#                     # if proper_vertex_in_triangulation(v_sub, retriangulation, Matrix_out[v_sub]):
#                     #     continue
#                     retriangulation = istriangulatable(v_sub, Matrix_out, list_boundary, Centroid,
#                                                        changed_tunnels, predetermined=retriangulation, predetermined_empty=empty_test_sub)[1]
#                 retriangulations.append(retriangulation)
#
#                 proper_vertices_final = [v for v in lst_vertices if proper_vertex_in_triangulation(v, retriangulation, Matrix_out[v])]
#                 proper_edges_final = [edge for edge in lst_edges if proper_edge_in_triangulation(edge, retriangulation)]
#                 proper_vertices_final_int = setdiff(proper_vertices_final, list_boundary)
#                 lsts.append(proper_edges_final)
#                 losses.append([item for item in proper_edges_initial if item not in proper_edges_final])
#                 gains_edges.append(len(proper_edges_final) - len(proper_edges_initial))
#                 gains_vertices.append(len(proper_vertices_final) - len(proper_vertices_initial))
#                 gains_vertices_int.append(len(proper_vertices_final_int) - len(proper_vertices_initial_int))
#                 # if sublist(proper_vertices_initial, proper_vertices_final):
#                 #     proper_edges_final = [edge for edge in lst_edges if proper_edge_in_triangulation(edge, retriangulation)]
#                 #     proper_vertices_final_int = setdiff(proper_vertices_final, list_boundary)
#                 #     lsts.append(proper_edges_final)
#                 #     losses.append([item for item in proper_edges_initial if item not in proper_edges_final])
#                 #     gains_edges.append(len(proper_edges_final) - len(proper_edges_initial))
#                 #     gains_vertices.append(len(proper_vertices_final) - len(proper_vertices_initial))
#                 #     gains_vertices_int.append(len(proper_vertices_final_int) - len(proper_vertices_initial_int))
#                 #     # if gains[-1] > 0:
#                 #     #     print(gains[-1], "- cycle:", cycle, "- empty tris:", emp)
#                 # else:
#                 #     lsts.append([])
#                 #     losses.append([])
#                 #     gains_edges.append(-1)
#                 #     gains_vertices.append(-1)
#                 #     gains_vertices_int.append(-1)
#         else:
#             candidates_unfold.append(cycle)
#             emptys.append(empty_test)
#             retriangulation = triangulation_test + star
#             for v_sub in union(unique(candidates), edge):
#                 # if proper_vertex_in_triangulation(v_sub, retriangulation, Matrix_out[v_sub]):
#                 #     continue
#                 retriangulation = istriangulatable(v_sub, Matrix_out, list_boundary, Centroid,
#                                                    changed_tunnels, predetermined=retriangulation, predetermined_empty=empty_test)[1]
#             retriangulations.append(retriangulation)
#
#             proper_vertices_final = [v for v in lst_vertices if proper_vertex_in_triangulation(v, retriangulation, Matrix_out[v])]
#             proper_edges_final = [edge for edge in lst_edges if proper_edge_in_triangulation(edge, retriangulation)]
#             proper_vertices_final_int = setdiff(proper_vertices_final, list_boundary)
#             lsts.append(proper_edges_final)
#             losses.append([item for item in proper_edges_initial if item not in proper_edges_final])
#             gains_edges.append(len(proper_edges_final) - len(proper_edges_initial))
#             gains_vertices.append(len(proper_vertices_final) - len(proper_vertices_initial))
#             gains_vertices_int.append(len(proper_vertices_final_int) - len(proper_vertices_initial_int))
#             # if sublist(proper_vertices_initial, proper_vertices_final):
#             #     proper_edges_final = [edge for edge in lst_edges if proper_edge_in_triangulation(edge, retriangulation)]
#             #     proper_vertices_final_int = setdiff(proper_vertices_final, list_boundary)
#             #     lsts.append(proper_edges_final)
#             #     losses.append([item for item in proper_edges_initial if item not in proper_edges_final])
#             #     gains_edges.append(len(proper_edges_final) - len(proper_edges_initial))
#             #     gains_vertices.append(len(proper_vertices_final) - len(proper_vertices_initial))
#             #     gains_vertices_int.append(len(proper_vertices_final_int) - len(proper_vertices_initial_int))
#             #     # print("\n", candidates_unfold[-1])
#             #     # print(gains[-1])
#             #     # print("gained:", setdiff(lst_final, lst_initial))
#             #     # print("lost:", setdiff(lst_initial, lst_final))
#             #     # if gains[-1] > 0:
#             #     #     print(gains[-1], "- cycle:", cycle)
#             # else:
#             #     lsts.append([])
#             #     losses.append([])
#             #     gains_edges.append(-1)
#             #     gains_vertices.append(-1)
#             #     gains_vertices_int.append(-1)
#         candidates_unfold.append(cycle)
#         emptys.append(empty_test)
#         retriangulation = triangulation_test + star
#         for v_sub in union(unique(candidates), edge):
#             # if proper_vertex_in_triangulation(v_sub, retriangulation, Matrix_out[v_sub]):
#             #     continue
#             Centroid_current = Centroid
#             if v_sub == edge[0]:
#                 Centroid_current = Centroid_temp_0
#             elif v_sub == edge[1]:
#                 Centroid_current = Centroid_temp_1
#             retriangulation = istriangulatable(v_sub, Matrix_out, list_boundary, Centroid_current,
#                                                changed_tunnels, predetermined=retriangulation, predetermined_empty=empty_test)[1]
#         retriangulations.append(retriangulation)
#
#         proper_vertices_final = [v for v in lst_vertices if proper_vertex_in_triangulation(v, retriangulation, Matrix_out[v])]
#         proper_edges_final = [edge for edge in lst_edges if proper_edge_in_triangulation(edge, retriangulation)]
#         proper_vertices_final_int = setdiff(proper_vertices_final, list_boundary)
#         lsts.append(proper_edges_final)
#         losses.append([item for item in proper_edges_initial if item not in proper_edges_final])
#         gains_edges.append(len(proper_edges_final) - len(proper_edges_initial))
#         gains_vertices.append(len(proper_vertices_final) - len(proper_vertices_initial))
#         gains_vertices_int.append(len(proper_vertices_final_int) - len(proper_vertices_initial_int))
#         # if sublist(proper_vertices_initial, proper_vertices_final):
#         #     proper_edges_final = [edge for edge in lst_edges if proper_edge_in_triangulation(edge, retriangulation)]
#         #     proper_vertices_final_int = setdiff(proper_vertices_final, list_boundary)
#         #     lsts.append(proper_edges_final)
#         #     losses.append([item for item in proper_edges_initial if item not in proper_edges_final])
#         #     gains_edges.append(len(proper_edges_final) - len(proper_edges_initial))
#         #     gains_vertices.append(len(proper_vertices_final) - len(proper_vertices_initial))
#         #     gains_vertices_int.append(len(proper_vertices_final_int) - len(proper_vertices_initial_int))
#         # else:
#         #     lsts.append([])
#         #     losses.append([])
#         #     gains_edges.append(-1)
#         #     gains_vertices.append(-1)
#         #     gains_vertices_int.append(-1)
#
#     m = max(gains_edges)
#     if m > 0:
#         Ind = [index for index in range(len(gains_edges)) if gains_edges[index] == m]
#         print("gains_edges: {}  -  {}".format(m, Ind))
#         for index in setdiff(range(len(gains_vertices)), Ind):
#             gains_vertices[index] = -math.inf
#         m = max(gains_vertices)
#         Ind = [index for index in range(len(gains_vertices)) if gains_vertices[index] == m]
#         print("gains_vertices: {}  -  {}".format(m, Ind))
#         for index in setdiff(range(len(gains_vertices_int)), Ind):
#             gains_vertices_int[index] = -math.inf
#         Ind = gains_vertices_int.index(max(gains_vertices_int))
#         print("gains_vertices_interior: {}  -  {}".format(m, [Ind]))
#         triangulation_reconstructed = loads(dumps(retriangulations[Ind]))
#         empty = loads(dumps(emptys[Ind]))
#         loop_over.extend(losses[Ind])
#         print("best found:", candidates_unfold[Ind])
#         log_file_detailed.write("\nConsidered edge: {}, best link found: {}, number of edges gained: {}".format(edge, candidates_unfold[Ind], max(gains_edges)))
#         log_file_detailed.flush()
#     else:
#         if multi and not isline(edge, multi_edges):
#             multi_edges.append(edge)
#     loop_over = loop_over[1:]
#
# save_list_of_lists(triangulation_reconstructed, "triangulation_test_step_3" + suffix)
#
# edges_defect_link, edges_defect_link_interior = [], []
# count, count_interior, count_reconstructable, count_reconstructable_interior = 0, 0, 0, 0
# for i in range(numbPoints):
#     Nei = Matrix_out[i]
#     for j in [item for item in Nei if item > i]:
#         if sublist([i, j], union(list_boundary, list_boundary_vox)):
#             continue
#         if not intersect([i, j], union(list_boundary, list_boundary_vox)):
#             count += 1
#             count_interior += 1
#             if not isline([i, j], multi_edges + caused_by_digons):
#                 count_reconstructable += 1
#                 count_reconstructable_interior += 1
#             if not proper_edge_in_triangulation([i, j], triangulation_reconstructed):
#                 edges_defect_link.append([i, j])
#                 edges_defect_link_interior.append([i, j])
#         else:
#             count += 1
#             if not isline([i, j], multi_edges + caused_by_digons):
#                 count_reconstructable += 1
#             if not proper_edge_in_triangulation([i, j], triangulation_reconstructed):
#                 edges_defect_link.append([i, j])
# count, count_interior = max(1, count), max(1, count_interior)
# count_reconstructable, count_reconstructable_interior = max(1, count_reconstructable), max(1, count_reconstructable_interior)
#
# save_list_of_lists(caused_by_digons, "caused_by_digons" + suffix)
# save_list_of_lists(multi_edges, "multi_edges" + suffix)
# save_list_of_lists(triangulation_reconstructed, "triangulation_test_step_4" + suffix)
#
# print("\nAll non-boundary edges: {}   -   reconstructed after changing links of non-proper edges: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# print("All interior edges: {}   -   reconstructed after changing links of non-proper edges (interior): {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
# print("Number of multi-edges:", len(multi_edges))
# print("Number of non-properly reconstructed edges around digons:", len(caused_by_digons))
# leftovers = [item for item in edges_defect_link if not isline(item, multi_edges + caused_by_digons)]
# leftovers_interior = [item for item in edges_defect_link_interior if not isline(item, multi_edges + caused_by_digons)]
# print("All reconstructable non-boundary edges: {}   -   reconstructed among reconstructable edges: {} %".format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
# print("All reconstructable interior edges: {}   -   reconstructed among reconstructable edges (interior): {} %".format(count_reconstructable_interior, round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))
#
#
# log_file_detailed.write("\n\nNumber of non-boundary edges: {}, properly reconstructed after changing links of non-proper edges: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# log_file_detailed.write("\nNumber of interior edges: {}, properly reconstructed after changing links of non-proper edges: {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
# log_file.write("\nNumber of non-boundary edges: {}, properly reconstructed after changing links of non-proper edges: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# log_file.write("\nNumber of interior edges: {}, properly reconstructed after changing links of non-proper edges: {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
# s = "pairs"
# if len(multi_edges) == 1:
#     s = "pair"
# log_file_detailed.write("\n\nNumber of multi-connected cells: {} {}".format(len(multi_edges), s))
# log_file_detailed.write("\nNumber of non-properly reconstructed edges around digons: {}".format(len(caused_by_digons)))
# log_file_detailed.write("\n\nNumber of reconstructable non-boundary edges: {}, properly reconstructed among reconstructable edges: {} %".format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
# log_file_detailed.write("\nNumber of reconstructable interior edges: {}, properly reconstructed after changing links of non-proper edges: {} %".format(count_reconstructable_interior, round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))
# log_file_detailed.write("\n\n{}".format(compute_time(start_time, True)))
# log_file_detailed.flush()
# log_file.write("\nNumber of multi-connected cells: {} {}".format(len(multi_edges), s))
# log_file.write("\nNumber of non-properly reconstructed edges around digons: {}".format(len(caused_by_digons)))
# log_file.write("\nNumber of reconstructable non-boundary edges: {}, properly reconstructed among reconstructable edges: {} %".format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
# log_file.write("\nNumber of reconstructable interior edges: {}, properly reconstructed after changing links of non-proper edges: {} %".format(count_reconstructable_interior, round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))
# log_file.write("\n{}".format(compute_time(start_time, True)))
# log_file.flush()
#
# start_time = time()
# print("\n\nChange tunnels in the triangulation...")
# log_file_detailed.write("\n\nChange tunnels in the triangulation...")
# log_file_detailed.flush()
# log_file.write("\n\nChange tunnels in the triangulation...")
# log_file.flush()
#
# K5 = find_K5s(Matrix_out)
#
# loop_over = [v for v in order if isline(v, leftovers)]
# while loop_over:
#     vertex = loop_over[0]
#     T, problems = proper_vertex_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex], return_problems=True)
#     if T:
#         loop_over = loop_over[1:]
#         continue
#     problems = [item for item in problems if not isline([vertex, item], multi_edges + caused_by_digons)]
#     if not problems:
#         loop_over = loop_over[1:]
#         continue
#     updated = False
#     K5_all = [K5[item] for item in ind(vertex, K5)]
#     lst_vertices = unique(K5_all)
#     proper_vertices_initial = [v for v in lst_vertices if proper_vertex_in_triangulation(v, triangulation_reconstructed, Matrix_out[v])]
#     proper_vertices_initial_int = setdiff(proper_vertices_initial, list_boundary)
#     K5_targeted = [item for item in K5_all if intersect(item, problems)]
#     for K5_current in K5_targeted:
#         if proper_vertex_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex]):
#             continue
#         lst_edges, tets = [], []
#         for clique in K5_all:
#             if len(intersect(clique, K5_current)) >= 4:
#                 lst_edges += [list(item) for item in combinations(clique, 2) if not isline(item, lst_edges)]
#                 tets += [list(item) for item in combinations(clique, 4) if not isline(item, tets)]
#         proper_edges_initial = [edge for edge in lst_edges if proper_edge_in_triangulation(edge, triangulation_reconstructed)]
#         gains_edges, gains_vertices, gains_vertices_int, retriangulations, lsts, emptys, losses = [], [], [], [], [], [], []
#         triangulation_test = loads(dumps(triangulation_reconstructed))
#         tets = [list(item) for item in combinations(K5_current, 4)]
#         for tet in tets:
#             triangulation_test = Remove_index(triangulation_test, ind(tet, triangulation_test))
#         for num in range(2 ** 5 - 1):  # K5 with 'all-ones' is not possible, so we disregard 31 = 11111
#             pos = bin(num)[2:].zfill(5)
#             retriangulation = loads(dumps(triangulation_test))
#             retriangulation += [tets[index] for index in range(len(pos)) if pos[index] == "1"]
#             emp = [tets[index] for index in range(len(pos)) if pos[index] == "0"]
#             empty_test = empty + [item for item in emp if not isline(item, empty)]
#             emptys.append(loads(dumps(empty_test)))
#             retriangulations.append(retriangulation)
#
#             proper_vertices_final = [v for v in lst_vertices if proper_vertex_in_triangulation(v, retriangulation, Matrix_out[v])]
#             proper_edges_final = [edge for edge in lst_edges if proper_edge_in_triangulation(edge, retriangulation)]
#             proper_vertices_final_int = setdiff(proper_vertices_final, list_boundary)
#             lsts.append(proper_edges_final)
#             losses.append([item for item in proper_edges_initial if item not in proper_edges_final])
#             gains_edges.append(len(proper_edges_final) - len(proper_edges_initial))
#             gains_vertices.append(len(proper_vertices_final) - len(proper_vertices_initial))
#             gains_vertices_int.append(len(proper_vertices_final_int) - len(proper_vertices_initial_int))
#             # if sublist(proper_vertices_initial, proper_vertices_final):
#             #     proper_edges_final = [edge for edge in lst_edges if proper_edge_in_triangulation(edge, retriangulation)]
#             #     proper_vertices_final_int = setdiff(proper_vertices_final, list_boundary)
#             #     lsts.append(proper_edges_final)
#             #     losses.append([item for item in proper_edges_initial if item not in proper_edges_final])
#             #     gains_edges.append(len(proper_edges_final) - len(proper_edges_initial))
#             #     gains_vertices.append(len(proper_vertices_final) - len(proper_vertices_initial))
#             #     gains_vertices_int.append(len(proper_vertices_final_int) - len(proper_vertices_initial_int))
#             # else:
#             #     lsts.append([])
#             #     losses.append([])
#             #     gains_edges.append(-1)
#             #     gains_vertices.append(-1)
#             #     gains_vertices_int.append(-1)
#
#         m = max(gains_edges)
#         if m > 0:
#             print("vertex: {}, 5-clique: {}".format(vertex, K5_current))
#             Ind = [index for index in range(len(gains_edges)) if gains_edges[index] == m]
#             print("gains_edges: {}  -  {}".format(m, Ind))
#             for index in setdiff(range(len(gains_vertices)), Ind):
#                 gains_vertices[index] = -math.inf
#             m = max(gains_vertices)
#             Ind = [index for index in range(len(gains_vertices)) if gains_vertices[index] == m]
#             print("gains_vertices: {}  -  {}".format(m, Ind))
#             for index in setdiff(range(len(gains_vertices_int)), Ind):
#                 gains_vertices_int[index] = -math.inf
#             Ind = gains_vertices_int.index(max(gains_vertices_int))
#             print("gains_vertices_interior: {}  -  {}".format(m, [Ind]))
#             triangulation_reconstructed = loads(dumps(retriangulations[Ind]))
#             empty = loads(dumps(emptys[Ind]))
#             loop_over.extend(unique(losses[Ind]))
#             updated = True
#             log_file_detailed.write("\nConsidered 5-clique: {}, number of edges gained after a tunnel change: {}".format(K5_current, max(gains_edges)))
#             log_file_detailed.flush()
#     if not updated:
#         loop_over = loop_over[1:]
#
# edges_defect_link, edges_defect_link_interior = [], []
# count, count_interior, count_reconstructable, count_reconstructable_interior = 0, 0, 0, 0
# for i in range(numbPoints):
#     Nei = Matrix_out[i]
#     for j in [item for item in Nei if item > i]:
#         if sublist([i, j], union(list_boundary, list_boundary_vox)):
#             continue
#         if not intersect([i, j], union(list_boundary, list_boundary_vox)):
#             count += 1
#             count_interior += 1
#             if not isline([i, j], multi_edges + caused_by_digons):
#                 count_reconstructable += 1
#                 count_reconstructable_interior += 1
#             if not proper_edge_in_triangulation([i, j], triangulation_reconstructed):
#                 edges_defect_link.append([i, j])
#                 edges_defect_link_interior.append([i, j])
#         else:
#             count += 1
#             if not isline([i, j], multi_edges + caused_by_digons):
#                 count_reconstructable += 1
#             if not proper_edge_in_triangulation([i, j], triangulation_reconstructed):
#                 edges_defect_link.append([i, j])
#
# count, count_interior = max(1, count), max(1, count_interior)
# count_reconstructable, count_reconstructable_interior = max(1, count_reconstructable), max(1, count_reconstructable_interior)
#
# print("\nAll non-boundary edges: {}   -   reconstructed after changing tunnels: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# print("All interior edges: {}   -   reconstructed after changing tunnels (interior): {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
# print("Number of multi-edges:", len(multi_edges))
# print("Number of non-properly reconstructed edges around digons:", len(caused_by_digons))
# leftovers = [item for item in edges_defect_link if not isline(item, multi_edges + caused_by_digons)]
# leftovers_interior = [item for item in edges_defect_link_interior if not isline(item, multi_edges + caused_by_digons)]
# print("All reconstructable non-boundary edges: {}   -   reconstructed among reconstructable edges: {} %".format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
# print("All reconstructable interior edges: {}   -   reconstructed among reconstructable edges (interior): {} %".format(count_reconstructable_interior, round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))
#
#
# log_file_detailed.write("\n\nNumber of non-boundary edges: {}, properly reconstructed after changing tunnels: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# log_file_detailed.write("\nNumber of interior edges: {}, properly reconstructed after changing tunnels: {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
# log_file.write("\nNumber of non-boundary edges: {}, properly reconstructed after changing tunnels: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# log_file.write("\nNumber of interior edges: {}, properly reconstructed after changing tunnels: {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
#
# log_file_detailed.write("\n\nNumber of reconstructable non-boundary edges: {}, properly reconstructed among reconstructable edges: {} %".format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
# log_file_detailed.write("\nNumber of reconstructable interior edges: {}, properly reconstructed among reconstructable edges: {} %".format(count_reconstructable_interior, round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))
#
# log_file.write("\nNumber of reconstructable non-boundary edges: {}, properly reconstructed among reconstructable edges: {} %".format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
# log_file.write("\nNumber of reconstructable interior edges: {}, properly reconstructed among reconstructable edges: {} %".format(count_reconstructable_interior, round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))
#
# save_list_of_lists(triangulation_reconstructed, "triangulation_test_step_5" + suffix)
#
"""
This step is skipped -- Date: 23.01.2024
# vertices_defect_link = []
# count = 0
# for vertex in range(numbPoints):
#     if vertex in union(list_boundary, list_boundary_vox, ignore):
#         continue
#     if vertex in union(unique(caused_by_digons), unique(multi_edges)):
#         continue
#     count += 1
#     if not proper_vertex_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex]):
#         vertices_defect_link.append(vertex)
#
# start_time = time()
# print("\n\nTry a further post-process...")
# log_file_detailed.write("\n\nTry a further post-process...")
# log_file_detailed.flush()
# log_file.write("\n\nTry a further post-process...")
# log_file.flush()
#
# # Try a further post-process:
# loop_over = vertices_defect_link
# while loop_over:
#     vertex = loop_over[0]
#     print("\rcurrent vertex: {}, number of remaining vertices {}".format(vertex, len(loop_over)), end="")
#     if len(vertexlink(vertex, Matrix_out)) > 50:  # too large to process
#         loop_over = loop_over[1:]
#         continue
#     if proper_vertex_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex]):
#         loop_over = loop_over[1:]
#         continue
#     triangulation_temp = istriangulatable(vertex, Matrix_out, list_boundary, Centroid, changed_tunnels, predetermined=triangulation_reconstructed, post_process=True, post_process_further=True)[1]
#
#     lst = Matrix_out[vertex]
#     for v in Matrix_out[vertex]:
#         lst = union(lst, Matrix_out[v])
#     lst_initial = [v for v in lst if proper_vertex_in_triangulation(v, triangulation_reconstructed, Matrix_out[v])]
#     lst_final = [v for v in lst if proper_vertex_in_triangulation(v, triangulation_temp, Matrix_out[v])]
#     lost = setdiff(setdiff(lst_initial, lst_final), union(list_boundary, list_boundary_vox))
#     gained = setdiff(setdiff(lst_final, lst_initial), union(list_boundary, list_boundary_vox))
#     print("\nlost:", lost)
#     print("gained:", gained)
#     if len(gained) > len(lost):
#         triangulation_reconstructed = loads(dumps(triangulation_temp))
#         loop_over = union(loop_over, setdiff(lost, union(unique(caused_by_digons), unique(multi_edges))))
#     else:
#         loop_over = loop_over[1:]
#
# edges_defect_link, edges_defect_link_interior = [], []
# count, count_interior, count_reconstructable, count_reconstructable_interior = 0, 0, 0, 0
# for i in range(numbPoints):
#     Nei = Matrix_out[i]
#     for j in [item for item in Nei if item > i]:
#         if sublist([i, j], union(list_boundary, list_boundary_vox)):
#             continue
#         if not intersect([i, j], union(list_boundary, list_boundary_vox)):
#             count += 1
#             count_interior += 1
#             if not isline([i, j], multi_edges + caused_by_digons):
#                 count_reconstructable += 1
#                 count_reconstructable_interior += 1
#             if not proper_edge_in_triangulation([i, j], triangulation_reconstructed):
#                 edges_defect_link.append([i, j])
#                 edges_defect_link_interior.append([i, j])
#         else:
#             count += 1
#             if not isline([i, j], multi_edges + caused_by_digons):
#                 count_reconstructable += 1
#             if not proper_edge_in_triangulation([i, j], triangulation_reconstructed):
#                 edges_defect_link.append([i, j])
#
# count, count_interior = max(1, count), max(1, count_interior)
# count_reconstructable, count_reconstructable_interior = max(1, count_reconstructable), max(1, count_reconstructable_interior)
#
# print("\nAll non-boundary edges: {}   -   reconstructed after a further post-process: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# print("All interior edges: {}   -   reconstructed after a further post-process (interior): {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
# print("Number of multi-edges:", len(multi_edges))
# print("Number of non-properly reconstructed edges around digons:", len(caused_by_digons))
# leftovers = [item for item in edges_defect_link if not isline(item, multi_edges + caused_by_digons)]
# leftovers_interior = [item for item in edges_defect_link_interior if not isline(item, multi_edges + caused_by_digons)]
# print("All reconstructable non-boundary edges: {}   -   reconstructed among reconstructable edges: {} %".format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
# print("All reconstructable interior edges: {}   -   reconstructed among reconstructable edges (interior): {} %".format(count_reconstructable_interior, round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))
#
#
# log_file_detailed.write("\n\nNumber of non-boundary edges: {}, properly reconstructed after a further post-process: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# log_file_detailed.write("\nNumber of interior edges: {}, properly reconstructed after a further post-process: {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
# log_file.write("\nNumber of non-boundary edges: {}, properly reconstructed after a further post-process: {} %".format(count, round(100 * (1 - len(edges_defect_link) / count), 2)))
# log_file.write("\nNumber of interior edges: {}, properly reconstructed after a further post-process: {} %".format(count_interior, round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
#
# log_file_detailed.write("\n\nNumber of reconstructable non-boundary edges: {}, properly reconstructed among reconstructable edges: {} %".format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
# log_file_detailed.write("\nNumber of reconstructable interior edges: {}, properly reconstructed among reconstructable edges: {} %".format(count_reconstructable_interior, round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))
#
# log_file.write("\nNumber of reconstructable non-boundary edges: {}, properly reconstructed among reconstructable edges: {} %".format(count_reconstructable, round(100 * (1 - len(leftovers) / count_reconstructable), 2)))
# log_file.write("\nNumber of reconstructable interior edges: {}, properly reconstructed among reconstructable edges: {} %".format(count_reconstructable_interior, round(100 * (1 - len(leftovers_interior) / count_reconstructable_interior), 2)))
#
# save_list_of_lists(triangulation_reconstructed, "triangulation_test_step_6" + suffix)
#
# log_file_detailed.write("\n\n{}\n".format(compute_time(start_time, True)))
# log_file_detailed.flush()
# log_file.write("\n{}\n".format(compute_time(start_time, True)))
# log_file.flush()
"""

with open("digons" + suffix, 'r') as f:
    digons = [int(elem) for elem in f.read().split('\n') if elem]
ignore = union(ignore, digons)

with open("Matrix_out" + suffix[:-4] + "_step_18.txt", 'r') as f:
    Matrix_out = [[int(num) for num in line.split(',')] for line in f]
for i in range(len(Matrix_out)):
    Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed

with open("Register_Removed" + suffix[:-4] + "_step_21.txt", 'r') as f:
    Register_Removed = [[int(num) for num in line.split(',')] for line in f]

with open("caused_by_digons" + suffix, 'r') as f:
    caused_by_digons = [[int(num) for num in line.split(',')] for line in f]

with open("multi_edges" + suffix, 'r') as f:
    multi_edges = [[int(num) for num in line.split(',')] for line in f]

with open("triangulation_test_step_5_meanDiameter_tristan_4_8.txt", 'r') as f:
    triangulation_reconstructed = [[int(num) for num in line.split(',')] for line in f]



with open("triangulation_reconstructed" + suffix, 'r') as f:
    triangulation_reconstructed = [[int(num) for num in line.split(',')] for line in f]
K2, K4 = [], []
for i in range(numbPoints):
    Nei_i = Matrix_out[i]
    for j in [item for item in Nei_i if item > i]:
        K2.append([i, j])
print("Edges:", len(K2))
sys.exit()



with open("digons_neighbors" + suffix, 'r') as f:
    digons_neighbors = [[int(num) for num in line.split(',')] for line in f]

start_time = time()
print("\n\nTry to remove excessive edges to repair defect links...")
log_file_detailed.write("\n\nTry to remove excessive edges to repair defect links...")
log_file_detailed.flush()
log_file.write("\n\nTry to remove excessive edges to repair defect links...")
log_file.flush()

vertices_defect_link = []
for vertex in range(numbPoints):
    if vertex in union(list_boundary, list_boundary_vox, ignore):
        continue
    if vertex in union(unique(caused_by_digons), unique(multi_edges)):
        continue
    if not proper_vertex_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex]):
        vertices_defect_link.append(vertex)

for vertex in vertices_defect_link:
    if proper_vertex_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex]):
        continue
    candidates, improvements, distances = [], [], []
    for u in Matrix_out[vertex]:
        Matrix_temp = Remove(vertex, u, Matrix_out)
        triangulation_temp = Remove_index(triangulation_reconstructed, ind([vertex, u], triangulation_reconstructed))
        if not proper_vertex_in_triangulation(vertex, triangulation_temp, Matrix_temp[vertex]):
            continue
        lst = union(mutual_neighbors([vertex, u], Matrix_out), u)
        T1 = [proper_vertex_in_triangulation(v, triangulation_reconstructed, Matrix_out[v]) for v in lst]
        T2 = [proper_vertex_in_triangulation(v, triangulation_temp, Matrix_temp[v]) for v in lst]
        if not all([not T1[item] or T2[item] for item in range(len(T1))]):
            continue
        candidates.append(u)
        improvements.append(sum(T2) - sum(T1))
        distances.append(distance_vox(vertex, u, Centroid, coordinates, voxels, dist_path))
    if not candidates:
        continue
    Ind = [item for item in range(len(candidates)) if improvements[item] == max(improvements)]
    Ind = distances.index(max([distances[item] for item in Ind]))
    u, dist = candidates[Ind], distances[Ind]
    Matrix_out = Remove(vertex, u, Matrix_out)
    Register_Removed.append([vertex, u])
    # log_file_detailed.write(
    #     "\nRemoved excessive edge {} of length {} to repair the link of {}".format([vertex, u], dist, vertex))
    # log_file_detailed.flush()
    print("Removed excessive edge {} of length {} to repair the link of {}".format([vertex, u], dist, vertex))
    triangulation_reconstructed = Remove_index(triangulation_reconstructed, ind([vertex, u], triangulation_reconstructed))

vertices_defect_link = []
lst_good = []
count_vts = 0
for vertex in range(numbPoints):
    if vertex in union(list_boundary, list_boundary_vox, ignore):
        continue
    if vertex in union(unique(caused_by_digons), unique(multi_edges)):
        continue
    count_vts += 1
    if proper_vertex_in_triangulation(vertex, triangulation_reconstructed, Matrix_out[vertex]):
        lst_good.append(vertex)
    else:
        vertices_defect_link.append(vertex)

save_list_of_lists(triangulation_reconstructed, "triangulation_reconstructed" + suffix)
save_list(vertices_defect_link, "vertices_defect_link" + suffix)
save_list_of_lists(Matrix_out, "Matrix_out" + suffix[:-4] + "_step_22.txt")
save_list_of_lists(Register_Removed, "Register_Removed" + suffix[:-4] + "_step_22.txt")

avg_deg = 0
for vertex in list_interior:
    avg_deg += len(Matrix_out[vertex])
print("\nAverage face degree =", round(avg_deg / len(list_interior), 3))
avg_deg += sum([len(item) for item in digons_neighbors])

# with open("Matrix_out" + suffix[:-4] + "_step_22.txt", 'r') as f:
#     Matrix_out = [[int(num) for num in line.split(',')] for line in f]
# for i in range(len(Matrix_out)):
#     Matrix_out[i] = [x for x in Matrix_out[i] if x != -1]  # -1's get removed

edges_defect_link, edges_defect_link_interior = [], []
count, count_interior = 0, 0
for i in range(numbPoints):
    Nei = Matrix_out[i]
    for j in [item for item in Nei if item > i]:
        if sublist([i, j], union(list_boundary, list_boundary_vox)):
            continue
        if intersect([i, j], union(unique(caused_by_digons), unique(multi_edges))):
            continue
        if not intersect([i, j], union(list_boundary, list_boundary_vox)):
            count += 1
            count_interior += 1
            if not proper_edge_in_triangulation([i, j], triangulation_reconstructed):
                edges_defect_link.append([i, j])
                edges_defect_link_interior.append([i, j])
        else:
            count += 1
            if not proper_edge_in_triangulation([i, j], triangulation_reconstructed):
                edges_defect_link.append([i, j])
vts_defect_link = setdiff(unique(edges_defect_link), list_boundary)
count, count_interior = max(1, count), max(1, count_interior)
##
print("Num tetrahedra:", len(triangulation_reconstructed))
##
print("Multi-connected pairs of cells: {}".format(len(multi_edges)))
print("Number of non-properly reconstructed edges around digons: {}".format(len(caused_by_digons)))
print("Non-boundary edges: {}".format(count))
print("Edges in the strict interior: {}".format(count_interior))
print("Percentage of properly reconstructed edges (non-boundary) {} %".format(round(100 * (1 - len(edges_defect_link) / count), 2)))
print("Percentage of properly reconstructed edges (strict interior) {} %".format(round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
print("\nAll cells:", numbPoints)
print("All interior vertices:", len(setdiff(range(numbPoints), union(list_boundary, list_boundary_vox))))
print("Reconstructable interior vertices (among reconstructables): {} %".format(
                                                        round(100 * (1 - len(vertices_defect_link) / count_vts), 2)))
print("Average face degree of interior cells = {}".format(round(avg_deg / len(list_interior), 3)))
sys.exit()
log_file_detailed.write("\n\nMulti-connected pairs of cells: {}".format(len(multi_edges)))
log_file_detailed.write("\nNumber of non-properly reconstructed edges around digons: {}".format(len(caused_by_digons)))
log_file_detailed.write("\nNon-boundary edges: {}".format(count))
log_file_detailed.write("\nEdges in the strict interior: {}".format(count_interior))
log_file_detailed.write("\nPercentage of properly reconstructed edges (non-boundary) {} %".format(round(100 * (1 - len(edges_defect_link) / count), 2)))
log_file_detailed.write("\nPercentage of properly reconstructed edges (strict interior) {} %".format(round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
log_file_detailed.write("\n\nAll cells: {}".format(numbPoints))
log_file_detailed.write("\nAll interior vertices: {}".format(len(setdiff(range(numbPoints), union(list_boundary,
                                                                                                  list_boundary_vox)))))
log_file_detailed.write("\nReconstructable interior vertices (among reconstructables): {} %".format(
                                                        round(100 * (1 - len(vertices_defect_link) / count_vts), 2)))
log_file_detailed.write("\nAverage face degree of interior cells = {}".format(round(avg_deg / len(list_interior), 3)))
log_file_detailed.write("\n\n{}\n".format(compute_time(start_time, True)))
log_file_detailed.flush()

log_file.write("\n\nMulti-connected pairs of cells: {}".format(len(multi_edges)))
log_file.write("\nNumber of non-properly reconstructed edges around digons: {}".format(len(caused_by_digons)))
log_file.write("\nNon-boundary edges: {}".format(count))
log_file.write("\nEdges in the strict interior: {}".format(count_interior))
log_file.write("\nPercentage of properly reconstructed edges (non-boundary) {} %".format(round(100 * (1 - len(edges_defect_link) / count), 2)))
log_file.write("\nPercentage of properly reconstructed edges (strict interior) {} %".format(round(100 * (1 - len(edges_defect_link_interior) / count_interior), 2)))
log_file.write("\n\nAll cells: {}".format(numbPoints))
log_file.write("\nAll interior vertices: {}".format(len(setdiff(range(numbPoints),
                                                                union(list_boundary, list_boundary_vox)))))
log_file.write("\nReconstructable interior vertices (among reconstructables): {} %".format(
                                                        round(100 * (1 - len(vertices_defect_link) / count_vts), 2)))
log_file.write("\nAverage face degree of interior cells = {}".format(round(avg_deg / len(list_interior), 3)))
log_file.write("\n\n{}\n".format(compute_time(start_time, True)))
log_file.flush()

degrees = {}
lst = union(list_boundary, list_boundary_vox)
for vertex in range(len(Matrix_out)):
    if vertex in lst:
        continue
    if len(Matrix_out[vertex]) in degrees:
        degrees[len(Matrix_out[vertex])] += 1
    else:
        degrees[len(Matrix_out[vertex])] = 1
count = 0
for i in range(4, 400):
    if i in degrees:
        count += degrees[i]
for i in range(4, 400):
    if i in degrees:
        print("{} vertices: {}%".format(i, round(100 * degrees[i] / count, 2)))
        log_file_detailed.write("\n{} vertices: {}%".format(i, round(100 * degrees[i] / count, 2)))
        log_file_detailed.flush()
        log_file.write("\n{} vertices: {}%".format(i, round(100 * degrees[i] / count, 2)))
        log_file.flush()

# cell_tetrahedron = [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]
# cell_prism_3 = [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 5], [2, 4, 5], [3, 4, 5]]
# cell_cube = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 6], [3, 5, 6], [4, 5, 6]]
# cell_prism_5 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 6], [3, 5, 7], [3, 6, 7], [4, 5, 7],
#                 [4, 6, 7]]
# cell_config_1 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 6], [3, 5, 7], [3, 6, 7], [4, 5, 8],
#                  [4, 6, 8], [5, 7, 8], [6, 7, 8]]
# cell_config_2 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8],
#                  [4, 5, 9], [4, 7, 9], [5, 8, 9], [6, 7, 8], [7, 8, 9]]
# cell_config_3 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 6], [3, 5, 7], [3, 6, 7], [4, 5, 8],
#                  [4, 6, 9], [4, 8, 9], [5, 7, 8], [6, 7, 9], [7, 8, 9]]
# cell_config_4 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8],
#                  [4, 5, 9], [4, 7, 9], [5, 8, 10], [5, 9, 10], [6, 7, 10], [6, 8, 10], [7, 9, 10]]
# cell_config_5 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 6], [3, 5, 7], [3, 6, 8], [3, 7, 8],
#                  [4, 5, 9], [4, 6, 10], [4, 9, 10], [5, 7, 9], [6, 8, 10], [7, 8, 9], [8, 9, 10]]
# cell_config_6 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8],
#                  [4, 5, 9], [4, 7, 9], [5, 8, 9], [6, 7, 10], [6, 8, 10], [7, 9, 10], [8, 9, 10]]
# cell_config_7 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8],
#                  [4, 5, 9], [4, 7, 9], [5, 8, 10], [5, 9, 10], [6, 7, 11], [6, 8, 11], [7, 9, 11], [8, 10, 11],
#                  [9, 10, 11]]
# cell_config_8 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8],
#                  [4, 5, 9], [4, 7, 9], [5, 8, 10], [5, 9, 10], [6, 7, 11], [6, 8, 10], [6, 10, 11], [7, 9, 11],
#                  [9, 10, 11]]
# cell_dodecahedron = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 6], [1, 5, 6], [2, 3, 7], [2, 4, 8], [2, 7, 8], [3, 5, 9],
#                      [3, 7, 9], [4, 6, 10], [4, 8, 10], [5, 6, 11], [5, 9, 11], [6, 10, 11], [7, 8, 12], [7, 9, 12],
#                      [8, 10, 12], [9, 11, 12], [10, 11, 12]]
# cell_config_9 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 9],
#                  [3, 8, 9], [4, 5, 10], [4, 7, 11], [4, 10, 11], [5, 8, 10], [6, 7, 12], [6, 9, 12], [7, 11, 12],
#                  [8, 9, 13], [8, 10, 13], [9, 12, 13], [10, 11, 13], [11, 12, 13]]
# cell_config_10 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8],
#                   [4, 5, 9], [4, 7, 10], [4, 9, 10], [5, 8, 11], [5, 9, 11], [6, 7, 12], [6, 8, 13], [6, 12, 13],
#                   [7, 10, 12], [8, 11, 13], [9, 10, 14], [9, 11, 14], [10, 12, 14], [11, 13, 14], [12, 13, 14]]
# cell_config_11 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 9],
#                   [3, 8, 9], [4, 5, 10], [4, 7, 11], [4, 10, 11], [5, 8, 12], [5, 10, 12], [6, 7, 13], [6, 9, 13],
#                   [7, 11, 13], [8, 9, 12], [9, 12, 14], [9, 13, 14], [10, 11, 14], [10, 12, 14], [11, 13, 14]]
# cell_config_12 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 9],
#                   [3, 8, 9], [4, 5, 10], [4, 7, 11], [4, 10, 11], [5, 8, 10], [6, 7, 12], [6, 9, 12], [7, 11, 13],
#                   [7, 12, 13], [8, 9, 14], [8, 10, 14], [9, 12, 14], [10, 11, 13], [10, 13, 14], [12, 13, 14]]
# cell_config_13 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 9],
#                   [3, 8, 9], [4, 5, 10], [4, 7, 11], [4, 10, 11], [5, 8, 12], [5, 10, 12], [6, 7, 13], [6, 9, 13],
#                   [7, 11, 13], [8, 9, 14], [8, 12, 14], [9, 13, 14], [10, 11, 15], [10, 12, 15], [11, 13, 15],
#                   [12, 14, 15], [13, 14, 15]]
# cell_config_14 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 9],
#                   [3, 8, 9], [4, 5, 10], [4, 7, 11], [4, 10, 11], [5, 8, 12], [5, 10, 12], [6, 7, 13], [6, 9, 13],
#                   [7, 11, 13], [8, 9, 12], [9, 12, 13], [10, 11, 12], [11, 12, 13]]
# cell_config_15 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8],
#                   [4, 5, 9], [4, 7, 10], [4, 9, 10], [5, 8, 11], [5, 9, 11], [6, 7, 12], [6, 8, 13], [6, 12, 13],
#                   [7, 10, 12], [8, 11, 13], [9, 10, 11], [10, 11, 12], [11, 12, 13]]
# cell_config_16 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8],
#                   [4, 5, 9], [4, 7, 10], [4, 9, 10], [5, 8, 11], [5, 9, 11], [6, 7, 12], [6, 8, 13], [6, 12, 13],
#                   [7, 10, 12], [8, 11, 14], [8, 13, 14], [9, 10, 14], [9, 11, 14], [10, 12, 14], [12, 13, 14]]
# cell_config_17 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 6], [1, 5, 6], [2, 3, 7], [2, 4, 8], [2, 7, 8], [3, 5, 9],
#                   [3, 7, 9], [4, 6, 10], [4, 8, 10], [5, 6, 11], [5, 9, 11], [6, 10, 12], [6, 11, 12], [7, 8, 13],
#                   [7, 9, 14], [7, 13, 14], [8, 10, 13], [9, 11, 14], [10, 12, 13], [11, 12, 14], [12, 13, 14]]
# cell_config_18 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8],
#                   [4, 5, 9], [4, 7, 10], [4, 9, 10], [5, 8, 11], [5, 9, 11], [6, 7, 12], [6, 8, 12], [7, 10, 12],
#                   [8, 11, 12], [9, 10, 11], [10, 11, 12]]
# cell_config_19 = [[1, 2, 3], [1, 2, 4], [1, 3, 5], [1, 4, 5], [2, 3, 6], [2, 4, 7], [2, 6, 7], [3, 5, 8], [3, 6, 8],
#                   [4, 5, 9], [4, 7, 10], [4, 9, 10], [5, 8, 11], [5, 9, 11], [6, 7, 12], [6, 8, 12], [7, 10, 12],
#                   [8, 11, 13], [8, 12, 13], [9, 10, 13], [9, 11, 13], [10, 12, 13]]
#
# famous_cells = {"cell_tetrahedron": cell_tetrahedron, "cell_prism_3": cell_prism_3, "cell_cube": cell_cube,
#                 "cell_prism_5": cell_prism_5, "cell_config_1": cell_config_1, "cell_config_2": cell_config_2,
#                 "cell_config_3": cell_config_3, "cell_config_4": cell_config_4, "cell_config_5": cell_config_5,
#                 "cell_config_6": cell_config_6, "cell_config_7": cell_config_7, "cell_config_8": cell_config_8,
#                 "cell_config_9": cell_config_9, "cell_dodecahedron": cell_dodecahedron,
#                 "cell_config_10": cell_config_10, "cell_config_11": cell_config_11, "cell_config_12": cell_config_12,
#                 "cell_config_13": cell_config_13, "cell_config_14": cell_config_14, "cell_config_15": cell_config_15,
#                 "cell_config_16": cell_config_16, "cell_config_17": cell_config_17, "cell_config_18": cell_config_18,
#                 "cell_config_19": cell_config_19}
#
# complex_dual = dual_list_in_triangulation(lst_good, triangulation_reconstructed, "dual_cells" + suffix,
#                                           setdiff(range(numbPoints), lst_good))
# frequent_types = open('combinatorial_types_frequent{}'.format(suffix), 'w+')
# frequent_types_underlying_flag = open('combinatorial_types_underlying_flag_frequent{}'.format(suffix), 'w+')
#
# for item in famous_cells:
#     frequent_types.write('{}, '.format(item))
#     frequent_types_underlying_flag.write('{}, '.format(item))
# frequent_types.write('all_reconstructed_interior_cells, all_combinatorial_types\n')
# frequent_types_underlying_flag.write('all_reconstructed_interior_cells, all_combinatorial_types\n')
# frequent_types.flush()
# frequent_types_underlying_flag.flush()
#
# file_out = open('combinatorial_types{}'.format(suffix), 'w+')
# file_out_underlying_flag = open('combinatorial_types_underlying_flag{}'.format(suffix), 'w+')
#
# multiplicity, unique_lsts = [], []
# multiplicity_underlying_flag, unique_lsts_underlying_flag = [], []
# print("\nCompute combinatorial types of interior cells...")
# log_file_detailed.write("\n\nCompute combinatorial types of interior cells...\n")
# log_file.write("\n\nCompute combinatorial types of interior cells...\n")
# for u in range(len(complex_dual)):
#     print("\r", u + 1, "out of", len(complex_dual), end="")
#     if not complex_dual[u]:
#         file_out.write('lex_ordered[{}] = {}\n'.format(u + 1, []))
#         file_out.flush()
#         file_out_underlying_flag.write('lex_ordered[{}] = {}\n'.format(u + 1, []))
#         file_out_underlying_flag.flush()
#         continue
#     cell_underlying_flag = loads(dumps(complex_dual[u]))
#     while True:
#         degrees = [sum([index in item for item in cell_underlying_flag]) for index in unique(cell_underlying_flag)]
#         if 3 not in degrees or len(degrees) == 4:
#             break
#         v = degrees.index(3) + 1
#         surr_v = [item for item in cell_underlying_flag if v in item]
#         cell_underlying_flag = [item for item in cell_underlying_flag if v not in item]
#         cell_underlying_flag.append(setdiff(unique(surr_v), v))
#         # relabel
#         cell_underlying_flag = [[x - 1 if x > v else x for x in item] for item in cell_underlying_flag]
#     cell_ordered = lex_order(complex_dual[u])
#     cell_underlying_flag_ordered = lex_order(cell_underlying_flag)
#     file_out.write('lex_ordered[{}] = {}\n'.format(u + 1, cell_ordered))
#     file_out.flush()
#     file_out_underlying_flag.write('lex_ordered[{}] = {}\n'.format(u + 1, cell_underlying_flag_ordered))
#     file_out_underlying_flag.flush()
#     if cell_ordered not in unique_lsts:
#         unique_lsts.append(cell_ordered)
#         multiplicity.append(1)
#     else:
#         multiplicity[unique_lsts.index(cell_ordered)] += 1
#     if cell_underlying_flag_ordered not in unique_lsts_underlying_flag:
#         unique_lsts_underlying_flag.append(cell_underlying_flag_ordered)
#         multiplicity_underlying_flag.append(1)
#     else:
#         multiplicity_underlying_flag[unique_lsts_underlying_flag.index(cell_underlying_flag_ordered)] += 1
# print("\n")
# file_out.close()
# file_out_underlying_flag.close()
#
# print("Number of interior reconstructed cells:", len(lst_good))
# print("Number of combinatorial types:", len(unique_lsts))
# print("----")
# print("Number of combinatorial types of underlying flag types:", len(unique_lsts_underlying_flag))
# log_file_detailed.write("Number of interior reconstructed cells: {}\nNumber of combinatorial types: {}\n"
#                         "Number of classes of underlying flag types: {}\n".format(
#                                                 len(lst_good), len(unique_lsts), len(unique_lsts_underlying_flag)))
# log_file.write("Number of interior reconstructed cells: {}\nNumber of combinatorial types: {}\n"
#                "Number of classes of underlying flag types: {}\n".format(
#                                                 len(lst_good), len(unique_lsts), len(unique_lsts_underlying_flag)))
# print("\n")
# for i in range(len(multiplicity_underlying_flag)):
#     if multiplicity_underlying_flag[i] > 10:
#         out = unique_lsts_underlying_flag[i]
#         if not out:
#             continue
#         for key in famous_cells:
#             if out == famous_cells[key]:
#                 out = key
#                 break
#         print("{} - {}".format(multiplicity_underlying_flag[i], out))
# print("")
# for item in famous_cells:
#     if famous_cells[item] in unique_lsts:
#         Ind = unique_lsts.index(famous_cells[item])
#         frequent_types.write('{}, '.format(multiplicity[Ind]))
#     else:
#         frequent_types.write('{}, '.format(0))
#     if famous_cells[item] in unique_lsts_underlying_flag:
#         Ind = unique_lsts_underlying_flag.index(famous_cells[item])
#         frequent_types_underlying_flag.write('{}, '.format(multiplicity_underlying_flag[Ind]))
#     else:
#         frequent_types_underlying_flag.write('{}, '.format(0))
# frequent_types.write('{}, {}\n'.format(len(lst_good), len(unique_lsts)))
# frequent_types.flush()
# frequent_types_underlying_flag.write('{}, {}\n'.format(len(lst_good), len(unique_lsts_underlying_flag)))
# frequent_types_underlying_flag.flush()
#
# frequent_types.close()
# frequent_types_underlying_flag.close()
log_file_detailed.close()
log_file.close()
