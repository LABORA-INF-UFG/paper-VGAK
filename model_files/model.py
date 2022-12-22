import time
import json
from docplex.mp.model import Model
# import networkx as nx
# import matplotlib.pyplot as plt
# from matplotlib.lines import Line2D


class Path:
    def __init__(self, id, source, target, seq, p1, p2, p3, delay_p1, delay_p2, delay_p3):
        self.id = id
        self.source = source
        self.target = target
        self.seq = seq
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.delay_p1 = delay_p1
        self.delay_p2 = delay_p2
        self.delay_p3 = delay_p3

    def __str__(self):
        return "ID: {}\tSEQ: {}\t P1: {}\t P2: {}\t P3: {}\t dP1: {}\t dP2: {}\t dP3: {}".format(self.id, self.seq, self.p1, self.p2, self.p3, self.delay_p1, self.delay_p2, self.delay_p3)


class CR:
    def __init__(self, id, cpu, num_BS):
        self.id = id
        self.cpu = cpu
        self.num_BS = num_BS

    def __str__(self):
        return "ID: {}\tCPU: {}".format(self.id, self.cpu)


class DRC:
    def __init__(self, id, cpu_CU, cpu_DU, cpu_RU, ram_CU, ram_DU, ram_RU, Fs_CU, Fs_DU, Fs_RU, delay_BH, delay_MH,
                 delay_FH, bw_BH, bw_MH, bw_FH):
        self.id = id

        self.cpu_CU = cpu_CU
        self.ram_CU = ram_CU
        self.Fs_CU = Fs_CU

        self.cpu_DU = cpu_DU
        self.ram_DU = ram_DU
        self.Fs_DU = Fs_DU

        self.cpu_RU = cpu_RU
        self.ram_RU = ram_RU
        self.Fs_RU = Fs_RU

        self.delay_BH = delay_BH
        self.delay_MH = delay_MH
        self.delay_FH = delay_FH

        self.bw_BH = bw_BH
        self.bw_MH = bw_MH
        self.bw_FH = bw_FH


class FS:
    def __init__(self, id, f_cpu, f_ram):
        self.id = id
        self.f_cpu = f_cpu
        self.f_ram = f_ram


class RU:
    def __init__(self, id, CR):
        self.id = id
        self.CR = CR

    def __str__(self):
        return "RU: {}\tCR: {}".format(self.id, self.CR)


links = []
capacity = {}
delay = {}
crs = {}
paths = {}
conj_Fs = {}


def read_topology():
    with open('model_files/T2_files/25_CRs_links_HC.json') as json_file:
        data = json.load(json_file)
        json_links = data["links"]
        for item in json_links:
            link = item
            source_node = link["fromNode"]
            destination_node = link["toNode"]
            if source_node < destination_node:
                capacity[(source_node, destination_node)] = link["capacity"]
                delay[(source_node, destination_node)] = link["delay"]
                links.append((source_node, destination_node))
            else:
                capacity[(destination_node, source_node)] = link["capacity"]
                delay[(destination_node, source_node)] = link["delay"]
                links.append((destination_node, source_node))
        with open('model_files/T2_files/25_CRs_nodes_HC.json') as json_file:
            data = json.load(json_file)
            json_nodes = data["nodes"]
            for item in json_nodes:
                node = item
                CR_id = node["nodeNumber"]
                node_CPU = node["cpu"]
                cr = CR(CR_id, node_CPU, 0)
                crs[CR_id] = cr
        crs[0] = CR(0, 0, 0)
        with open('paths.json') as json_paths_file:
            json_paths_f = json.load(json_paths_file)
            json_paths = json_paths_f["paths"]
            for item in json_paths:
                path = json_paths[item]
                path_id = path["id"]
                path_source = path["source"]
                if path_source == "CN":
                    path_source = 0
                path_target = path["target"]
                path_seq = path["seq"]
                paths_p = [path["p1"], path["p2"], path["p3"]]
                list_p1 = []
                list_p2 = []
                list_p3 = []
                for path_p in paths_p:
                    aux = ""
                    sum_delay = 0
                    for tup in path_p:
                        aux += tup
                        tup_aux = tup
                        tup_aux = tup_aux.replace('(', '')
                        tup_aux = tup_aux.replace(')', '')
                        tup_aux = tuple(map(int, tup_aux.split(', ')))
                        if path_p == path["p1"]:
                            list_p1.append(tup_aux)
                        elif path_p == path["p2"]:
                            list_p2.append(tup_aux)
                        elif path_p == path["p3"]:
                            list_p3.append(tup_aux)
                        sum_delay += delay[tup_aux]
                    if path_p == path["p1"]:
                        delay_p1 = sum_delay
                    elif path_p == path["p2"]:
                        delay_p2 = sum_delay
                    elif path_p == path["p3"]:
                        delay_p3 = sum_delay
                    if path_seq[0] == 0:
                        delay_p1 = 0
                    if path_seq[1] == 0:
                        delay_p2 = 0
                p = Path(path_id, path_source, path_target, path_seq, list_p1, list_p2, list_p3, delay_p1, delay_p2, delay_p3)
                paths[path_id] = p


def DRC_structure():

    DRC1 = DRC(1, 0.1, 0.42, 0.48, 0.01, 0.01, 0.01,
               ['f8'], ['f7', 'f6', 'f5', 'f4', 'f3', 'f2'], ['f1', 'f0'], 10, 10, 0.25, 9.9, 13.2, 42.6)
    DRC2 = DRC(2, 0.2, 0.32, 0.48, 0.01, 0.01, 0.01,
               ['f8', 'f7'], ['f6', 'f5', 'f4', 'f3', 'f2'], ['f1', 'f0'], 10, 10, 0.25, 9.9, 13.2, 42.6)
    DRC4 = DRC(4, 0.1, 0.25, 0.65, 0.01, 0.01, 0.01,
               ['f8'], ['f7', 'f6', 'f5', 'f4', 'f3'], ['f2', 'f1', 'f0'], 10, 10, 0.25, 9.9, 13.2, 13.6)
    DRC5 = DRC(5, 0.2, 0.15, 0.65, 0.01, 0.01, 0.01,
               ['f8', 'f7'], ['f6', 'f5', 'f4', 'f3'], ['f2', 'f1', 'f0'], 10, 10, 0.25, 9.9, 13.2, 13.6)
    DRC9 = DRC(9, 0, 0.52, 0.48, 0, 0.01, 0.01,
               [0], ['f8', 'f7', 'f6', 'f5', 'f4', 'f3', 'f2'], ['f1', 'f0'], 0, 10, 0.25, 0, 9.9, 42.6)
    DRC10 = DRC(10, 0, 0.35, 0.65, 0, 0.01, 0.01,
                [0], ['f8', 'f7', 'f6', 'f5', 'f4', 'f3'], ['f2', 'f1', 'f0'], 0, 10, 0.25, 0, 3, 13.6)
    DRC8 = DRC(8, 0, 0, 1, 0, 0, 0.01,
               [0], [0], ['f8', 'f7', 'f6', 'f5', 'f4', 'f3', 'f2', 'f1', 'f0'], 0, 0, 10, 0, 0, 9.9)

    DRCs = {1: DRC1, 2: DRC2, 4: DRC4, 5: DRC5, 8: DRC8, 9: DRC9, 10: DRC10}

    return DRCs


def make_links(link_list):
    for i, n in enumerate(link_list):
        for j, k in enumerate(n):
            if k:
                yield (i, j)


def geo_neighbors(geo_links, rus):
    neighbors = {}

    for ru in rus:
        neighbors[ru] = []
        for l in geo_links:
            if l[0] == ru:
                neighbors[ru].append(l[1])
            if l[1] == ru:
                neighbors[ru].append(l[0])

    print(neighbors)
    return neighbors


def RU_location():
    rus = {}
    count = 1
    with open('model_files/T2_files/25_CRs_nodes_HC.json') as json_file:
        data = json.load(json_file)
        json_crs = data["nodes"]
        for item in json_crs:
            node = item
            num_rus = node["RU"]
            num_cr = node["nodeNumber"]
            for i in range(0, num_rus):
                rus[count] = RU(count, int(num_cr))
                count += 1
    return rus


DRC_f1 = 0
f1_vars = []
f2_vars = []


def run_PlaceRAN():
    print("Running Stage - 1")
    print("-----------------------------------------------------------------------------------------------------------")
    alocation_time_start = time.time()
    read_topology()
    DRCs = DRC_structure()
    rus = RU_location()

    RUs_demand = {1 : 610.08, 2 : 637.44, 3 : 637.44, 4 : 637.44, 5 : 590.29, 6 : 656.00, 7 : 620.13, 8 : 620.13, 9 : 620.13, 10 : 620.13, 11 : 637.44, 12 : 637.44, 13 : 620.13, 14 : 670.84, 15 : 656.00, 16 : 620.13, 17 : 598.95, 18 : 637.44, 19 : 610.08, 20 : 620.13, 21 : 610.08 }

    geo_links = [(1, 9), (1, 10), (1, 4), (1, 11), (1, 14), (1, 2), (9, 10), (9, 13), (9, 4), (10, 11), (12, 13),
                 (12, 4), (12, 19), (12, 3), (12, 14), (13, 4), (4, 14), (18, 19), (18, 3), (18, 8), (19, 3), (3, 7),
                 (3, 14), (3, 8), (5, 11), (5, 15), (5, 2), (5, 16), (5, 17), (11, 15), (11, 2), (15, 17), (7, 14),
                 (7, 2), (7, 8), (7, 21), (7, 16), (14, 2), (2, 16), (20, 8), (20, 21), (20, 6), (8, 21), (21, 16),
                 (6, 21), (6, 16), (6, 17), (16, 17)]#, (18, 20)] # ADICIONEI 18 E 20

    neighbors = geo_neighbors(geo_links, rus)

    F1 = FS('f8', 2, 2)
    F2 = FS('f7', 2, 2)
    F3 = FS('f6', 2, 2)
    F4 = FS('f5', 2, 2)
    F5 = FS('f4', 2, 2)
    F6 = FS('f3', 2, 2)
    F7 = FS('f2', 2, 2)
    F8 = FS('f1', 2, 2)
    F9 = FS('f0', 2, 2)

    conj_Fs = {'f8': F1, 'f7': F2, 'f6': F3, 'f5': F4, 'f4': F5, 'f3': F6, 'f2': F7}


    mdl = Model(name='NGRAN Problem', log_output=True)
    mdl.parameters.mip.tolerances.mipgap = 0

    i = [(p, d, b) for p in paths for d in DRCs for b in rus if paths[p].seq[2] == rus[b].CR]

    mdl.x = mdl.binary_var_dict(i, name='x')

    count_clusters = mdl.sum(mdl.min(1, mdl.sum(mdl.x[it] for it in i if paths[it[0]].seq[1] == c and it[1] != 8)) for c in crs)

    count_clusters_D_RAN = mdl.sum(mdl.min(1, mdl.sum(mdl.x[it] for it in i if paths[it[0]].seq[2] == c and it[1] == 8)) for c in crs)

    count_agg = mdl.sum(
        mdl.sum(mdl.max(0, (mdl.sum(mdl.x[it] for it in i if ((o in DRCs[it[1]].Fs_CU and paths[it[0]].seq[0] == crs[c].id) or (
                    o in DRCs[it[1]].Fs_DU and paths[it[0]].seq[1] == crs[c].id) or (
                                                                          o in DRCs[it[1]].Fs_RU and paths[it[0]].seq[
                                                                      2] == crs[c].id))) - 1)) for o in conj_Fs) for c in crs)

    mdl.minimize(count_clusters + count_clusters_D_RAN)

    mdl.add_constraint(mdl.sum(mdl.x[it] for it in i if it[2] == 18 and paths[it[0]].seq[1] == 4) == 0)

    for b in rus:
        mdl.add_constraint(mdl.sum(mdl.x[it] for it in i if it[2] == b) == 1, 'unicity')
    mdl.add_constraint(mdl.sum(mdl.x[it] for it in i if paths[it[0]].target != rus[it[2]].CR) == 0, 'path')
    mdl.add_constraint(mdl.sum(mdl.x[it] for it in i if paths[it[0]].seq[0] != 0 and (it[1] == 6 or it[1] == 7 or it[1] == 8 or it[1] == 9 or it[1] == 10)) == 0, 'DRCs_path_pick')
    mdl.add_constraint(mdl.sum(mdl.x[it] for it in i if paths[it[0]].seq[0] == 0 and it[1] != 6 and it[1] != 7 and it[1] != 8 and it[1] != 9 and it[1] != 10) == 0, 'DRCs_path_pick2')
    mdl.add_constraint(mdl.sum(mdl.x[it] for it in i if paths[it[0]].seq[0] == 0 and paths[it[0]].seq[1] == 0 and it[1] != 8) == 0, 'DRCs_path_pick3')
    mdl.add_constraint(mdl.sum(mdl.x[it] for it in i if paths[it[0]].seq[0] == 0 and paths[it[0]].seq[1] != 0 and it[1] == 8) == 0, 'DRCs_path_pick4')
    for ru in rus:
        mdl.add_constraint(mdl.sum(mdl.x[it] for it in i if paths[it[0]].seq[2] != rus[ru].CR and it[2] == rus[ru].id) == 0)
    for l in links:
        k = (l[1], l[0])
        mdl.add_constraint(
            mdl.sum(mdl.x[it] * DRCs[it[1]].bw_BH for it in i if l in paths[it[0]].p1 or k in paths[it[0]].p1)
            + mdl.sum(mdl.x[it] * DRCs[it[1]].bw_MH for it in i if l in paths[it[0]].p2 or k in paths[it[0]].p2)
            + mdl.sum(mdl.x[it] * DRCs[it[1]].bw_FH for it in i if l in paths[it[0]].p3 or k in paths[it[0]].p3)
            <= capacity[l], 'links_bw')
    for it in i:
        mdl.add_constraint((mdl.x[it] * paths[it[0]].delay_p1) <= DRCs[it[1]].delay_BH, 'delay_req_p1')
    for it in i:
        mdl.add_constraint((mdl.x[it] * paths[it[0]].delay_p2) <= DRCs[it[1]].delay_MH, 'delay_req_p2')
    for it in i:
        mdl.add_constraint((mdl.x[it] * paths[it[0]].delay_p3 <= DRCs[it[1]].delay_FH), 'delay_req_p3')
    for c in crs:
        mdl.add_constraint(mdl.sum(mdl.x[it] * DRCs[it[1]].cpu_CU * RUs_demand[it[2]] for it in i
                                   if c == paths[it[0]].seq[0])
                           + mdl.sum(mdl.x[it] * DRCs[it[1]].cpu_DU * RUs_demand[it[2]] for it in i
                                   if c == paths[it[0]].seq[1])
                           + mdl.sum(mdl.x[it] * DRCs[it[1]].cpu_RU * RUs_demand[it[2]] for it in i
                                     if c == paths[it[0]].seq[2]) <= crs[c].cpu, 'crs_cpu_usage')

    # NEW CONSTRAINT GEO NEIGHBORS
    for it in i:
        if it[1] != 8:
            mdl.add_constraint(mdl.x[it] <= mdl.sum(mdl.x[it1] for it1 in i if it1[2] in neighbors[it[2]]
                                                    and paths[it[0]].seq[1] == paths[it1[0]].seq[1]))

    alocation_time_end = time.time()
    start_time = time.time()
    mdl.solve()
    end_time = time.time()
    print("Stage 1 - Alocation Time: {}".format(alocation_time_end - alocation_time_start))
    print("Stage 1 - Enlapsed Time: {}".format(end_time - start_time))
    print("Stage 1 -Optimal Solution: {}".format(mdl.solution.get_objective_value()))
    print("Stage 1 - Alocation Time: {}".format(alocation_time_end - alocation_time_start))
    print("Stage 1 - Enlapsed Time: {}".format(end_time - start_time))
    print("Stage 1 -Optimal Solution: {}".format(mdl.solution.get_objective_value()))
    global f1_vars
    for it in i:
        if mdl.x[it].solution_value > 0:
            f1_vars.append(it)

    DUs_map = {}
    for it in i:
        if mdl.x[it].solution_value > 0.8:
            if it[1] != 8:
                DUs_map[it[2]] = paths[it[0]].seq[1]
            else:
                DUs_map[it[2]] = paths[it[0]].seq[2]
            print("x{} -> {}".format(it, mdl.x[it].solution_value))
            print(paths[it[0]].seq)

    disp_Fs = {}

    for cr in crs:
        disp_Fs[cr] = {'f8': 0, 'f7': 0, 'f6': 0, 'f5': 0, 'f4': 0, 'f3': 0, 'f2': 0, 'f1': 0, 'f0': 0}

    for it in i:
        for cr in crs:
            if mdl.x[it].solution_value > 0.8:
                if cr in paths[it[0]].seq:
                    seq = paths[it[0]].seq
                    if cr == seq[0]:
                        Fs = DRCs[it[1]].Fs_CU
                        for o in Fs:
                            if o != 0:
                                dct = disp_Fs[cr]
                                dct["{}".format(o)] += 1
                                disp_Fs[cr] = dct

                    if cr == seq[1]:
                        Fs = DRCs[it[1]].Fs_DU
                        for o in Fs:
                            if o != 0:
                                dct = disp_Fs[cr]
                                dct["{}".format(o)] += 1
                                disp_Fs[cr] = dct

                    if cr == seq[2]:
                        Fs = DRCs[it[1]].Fs_RU
                        for o in Fs:
                            if o != 0:
                                dct = disp_Fs[cr]
                                dct["{}".format(o)] += 1
                                disp_Fs[cr] = dct

    print("FO: {}".format(mdl.solution.get_objective_value()))

    with open("stage_1_solution.json", "w") as stage_3_result:
        result_list = {"Solution": []}
        for it in i:
            if mdl.x[it].solution_value > 0.1:
                result = {"RU_id": 0, "RU_DRC": 0, "CU_loc": 0, "DU_loc": 0, "RU_loc": 0, "path": []}
                path_sol = []
                sol_dsg = it[1]
                ru_id = it[2]
                if paths[it[0]].p1:
                    for item in paths[it[0]].p1:
                        path_sol.append(item)
                if paths[it[0]].p2:
                    for item in paths[it[0]].p2:
                        path_sol.append(item)
                if paths[it[0]].p3:
                    for item in paths[it[0]].p3:
                        path_sol.append(item)
                result["path"] = path_sol
                cu_loc = paths[it[0]].seq[0]
                du_loc = paths[it[0]].seq[1]
                ru_loc = paths[it[0]].seq[2]
                result["RU_id"] = ru_id
                if du_loc == 0:
                    du_loc = ru_loc
                    cu_loc = du_loc
                elif cu_loc == 0 and it[1] > 8:
                    cu_loc = du_loc
                elif cu_loc == 0 and it[1] < 9:
                    cu_loc = du_loc
                    du_loc = ru_loc
                result["RU_id"] = ru_id
                result["RU_DRC"] = sol_dsg
                result["CU_loc"] = cu_loc
                result["DU_loc"] = du_loc
                result["RU_loc"] = ru_loc
                result["path"] = path_sol

                result_list["Solution"].append(result)
        json.dump(result_list, stage_3_result)

        geo_links = [(1, 9), (1, 10), (1, 4), (1, 11), (1, 14), (1, 2), (9, 10), (9, 13), (9, 4), (10, 11), (12, 13),
                     (12, 4), (12, 19), (12, 3), (12, 14), (13, 4), (4, 14), (18, 19), (18, 3), (18, 8), (19, 3),
                     (3, 7),
                     (3, 14), (3, 8), (5, 11), (5, 15), (5, 2), (5, 16), (5, 17), (11, 15), (11, 2), (15, 17), (7, 14),
                     (7, 2), (7, 8), (7, 21), (7, 16), (14, 2), (2, 16), (20, 8), (20, 21), (20, 6), (8, 21), (21, 16),
                     (6, 21), (6, 16), (6, 17), (16, 17)]

        G = nx.Graph()
        G.add_edges_from(geo_links)
        color_map = []

        print(DUs_map)

        for node in G:
            if DUs_map[node] == 1:
                color_map.append("blue")
            elif DUs_map[node] == 2:
                color_map.append("green")
            elif DUs_map[node] == 3:
                color_map.append("gray")
            elif DUs_map[node] == 4:
                color_map.append("m")
            elif DUs_map[node] == 5:
                color_map.append("darkcyan")
            elif DUs_map[node] == 6:
                color_map.append("lightseagreen")
            elif DUs_map[node] == 7:
                color_map.append("plum")
            elif DUs_map[node] == 8:
                color_map.append("violet")
            elif DUs_map[node] == 9:
                color_map.append("limegreen")
            elif DUs_map[node] == 10:
                color_map.append("lime")
            elif DUs_map[node] == 11:
                color_map.append("peru")
            elif DUs_map[node] == 12:
                color_map.append("chocolate")
            else:
                color_map.append("red")

        fig, ax = plt.subplots(1, figsize=(12, 8), dpi=60)
        nx.draw_networkx(G, node_color=color_map)

        plt.show()

    return mdl.solution.get_objective_value()


if __name__ == '__main__':
    start_all = time.time()
    run_PlaceRAN()
    end_all = time.time()
    print("TOTAL TIME: {}".format(end_all - start_all))