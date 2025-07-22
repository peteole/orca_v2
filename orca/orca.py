#!/usr/bin/env python3
"""
Pure Python implementation of ORCA (Orbit Counting Algorithm)
for counting graphlet orbits in networks.

This is a translation of the original C++ ORCA implementation.
It counts node orbits and edge orbits for graphlets on 4 or 5 nodes.
"""

import sys
import time
import bisect
from typing import Literal, Optional
import numpy as np


class Pair:
    """Represents an undirected edge (pair of nodes)"""

    def __init__(self, a: int, b: int):
        self.a = min(a, b)
        self.b = max(a, b)

    def __lt__(self, other):
        if self.a == other.a:
            return self.b < other.b
        return self.a < other.a

    def __eq__(self, other):
        return self.a == other.a and self.b == other.b

    def __hash__(self):
        return (self.a << 8) ^ (self.b << 0)


class Triple:
    """Represents a sorted triple of nodes"""

    def __init__(self, a: int, b: int, c: int):
        nodes = sorted([a, b, c])
        self.a, self.b, self.c = nodes

    def __lt__(self, other):
        if self.a == other.a:
            if self.b == other.b:
                return self.c < other.c
            return self.b < other.b
        return self.a < other.a

    def __eq__(self, other):
        return self.a == other.a and self.b == other.b and self.c == other.c

    def __hash__(self):
        return (self.a << 16) ^ (self.b << 8) ^ (self.c << 0)


class Orca:
    """ORCA algorithm implementation for counting graphlet orbits"""

    def __init__(self):
        self.n = 0  # number of nodes
        self.m = 0  # number of edges
        self.deg = []  # degrees of individual nodes
        self.edges = []  # list of edges as Pair objects

        # Adjacency representations
        self.adj = []  # adjacency lists
        self.inc = []  # incidence lists: (neighbor, edge_id)
        self.adj_matrix = {}  # adjacency matrix for small graphs
        self.use_adj_matrix = False

        # Orbit counts
        self.orbit = []  # orbit[x][o] - how many times node x participates in orbit o
        self.eorbit = []  # eorbit[x][o] - how many times edge x participates in edge orbit o

        # Common neighbor caches
        self.common2 = {}  # common neighbors for pairs
        self.common3 = {}  # common neighbors for triples

    def adjacent(self, x: int, y: int) -> bool:
        """Check if nodes x and y are adjacent"""
        if self.use_adj_matrix:
            return (x, y) in self.adj_matrix
        else:
            # Binary search in sorted adjacency list
            pos = bisect.bisect_left(self.adj[x], y)
            return pos < len(self.adj[x]) and self.adj[x][pos] == y

    def get_edge_id(self, x: int, y: int) -> int:
        """Get the edge ID for edge (x,y)"""
        pos = bisect.bisect_left(self.adj[x], y)
        if pos < len(self.adj[x]) and self.adj[x][pos] == y:
            return self.inc[x][pos][1]
        return -1

    def common2_get(self, pair: Pair) -> int:
        """Get common neighbor count for a pair"""
        return self.common2.get(pair, 0)

    def common3_get(self, triple: Triple) -> int:
        """Get common neighbor count for a triple"""
        return self.common3.get(triple, 0)

    def setup_graph_from_edges(self, edge_list: np.ndarray, num_nodes: Optional[int] = None) -> bool:
        """Setup graph directly from edge list array"""
        try:
            # Convert edge list to Python ints to avoid numpy int types
            edge_list = edge_list.astype(int)
            
            if num_nodes is None:
                self.n = int(np.max(edge_list)) + 1
            else:
                self.n = num_nodes
            
            self.m = len(edge_list)
            self.deg = [0] * self.n
            self.edges = []

            edge_set = set()
            for i, (a, b) in enumerate(edge_list):
                a, b = int(a), int(b)  # Ensure Python ints
                
                if not (0 <= a < self.n) or not (0 <= b < self.n):
                    raise ValueError("Node ids should be between 0 and n-1.")

                if a == b:
                    raise ValueError("Self loops are not allowed.")

                edge = Pair(a, b)
                if edge in edge_set:
                    raise ValueError("Duplicate edges found.")

                edge_set.add(edge)
                self.edges.append(edge)
                self.deg[a] += 1
                self.deg[b] += 1

            return True

        except Exception as e:
            raise ValueError(f"Error setting up graph: {e}")

    def read_graph(self, filename: str) -> bool:
        """Read graph from input file"""
        try:
            with open(filename, "r") as f:
                line = f.readline().strip()
                self.n, self.m = map(int, line.split())

                self.deg = [0] * self.n
                self.edges = []

                edge_set = set()
                for i in range(self.m):
                    line = f.readline().strip()
                    a, b = map(int, line.split())

                    if not (0 <= a < self.n) or not (0 <= b < self.n):
                        print("Node ids should be between 0 and n-1.")
                        return False

                    if a == b:
                        print("Self loops are not allowed.")
                        return False

                    edge = Pair(a, b)
                    if edge in edge_set:
                        print("Duplicate edges found.")
                        return False

                    edge_set.add(edge)
                    self.edges.append(edge)
                    self.deg[a] += 1
                    self.deg[b] += 1

                d_max = max(self.deg) if self.deg else 0
                print(f"nodes: {self.n}")
                print(f"edges: {self.m}")
                print(f"max degree: {d_max}")

                return True

        except Exception as e:
            print(f"Error reading file: {e}")
            return False

    def setup_adjacency(self):
        """Set up adjacency and incidence lists"""
        # Decide whether to use adjacency matrix
        if self.n * self.n < 100 * 1024 * 1024 * 8:
            self.use_adj_matrix = True
            self.adj_matrix = {}
            for edge in self.edges:
                self.adj_matrix[(edge.a, edge.b)] = True
                self.adj_matrix[(edge.b, edge.a)] = True

        # Initialize adjacency and incidence lists
        self.adj = [[] for _ in range(self.n)]
        self.inc = [[] for _ in range(self.n)]

        # Build adjacency and incidence lists
        for i, edge in enumerate(self.edges):
            a, b = edge.a, edge.b
            self.adj[a].append(b)
            self.adj[b].append(a)
            self.inc[a].append((b, i))
            self.inc[b].append((a, i))

        # Sort adjacency and incidence lists
        for i in range(self.n):
            self.adj[i].sort()
            self.inc[i].sort()

        # Initialize orbit counts
        self.orbit = [[0] * 73 for _ in range(self.n)]
        self.eorbit = [[0] * 68 for _ in range(self.m)]

    def count4(self):
        """Count graphlets on max 4 nodes"""
        start_time = time.time()
        start_time_all = start_time

        print("stage 1 - precomputing common nodes")
        # Precompute triangles that span over edges
        tri = [0] * self.m
        for i in range(self.m):
            if i % (self.m // 100 + 1) == 0:
                progress = 100 * i // self.m
                print(f"{progress}%", end="\r")

            x, y = self.edges[i].a, self.edges[i].b
            xi, yi = 0, 0

            while xi < len(self.adj[x]) and yi < len(self.adj[y]):
                if self.adj[x][xi] == self.adj[y][yi]:
                    tri[i] += 1
                    xi += 1
                    yi += 1
                elif self.adj[x][xi] < self.adj[y][yi]:
                    xi += 1
                else:
                    yi += 1

        end_time = time.time()
        print(f"\n{end_time - start_time:.2f}")
        start_time = end_time

        # Count full graphlets
        print("stage 2 - counting full graphlets")
        C4 = [0] * self.n

        for x in range(self.n):
            if x % (self.n // 100 + 1) == 0:
                progress = 100 * x // self.n
                print(f"{progress}%", end="\r")

            for nx in range(len(self.adj[x])):
                y = self.adj[x][nx]
                if y >= x:
                    break

                neigh = []
                for ny in range(len(self.adj[y])):
                    z = self.adj[y][ny]
                    if z >= y:
                        break
                    if self.adjacent(x, z):
                        neigh.append(z)

                for i in range(len(neigh)):
                    z = neigh[i]
                    for j in range(i + 1, len(neigh)):
                        zz = neigh[j]
                        if self.adjacent(z, zz):
                            C4[x] += 1
                            C4[y] += 1
                            C4[z] += 1
                            C4[zz] += 1

        end_time = time.time()
        print(f"\n{end_time - start_time:.2f}")
        start_time = end_time

        # Build systems of equations
        print("stage 3 - building systems of equations")
        common = [0] * self.n
        common_list = []

        for x in range(self.n):
            if x % (self.n // 100 + 1) == 0:
                progress = 100 * x // self.n
                print(f"{progress}%", end="\r")

            f_12_14 = f_10_13 = 0
            f_13_14 = f_11_13 = 0
            f_7_11 = f_5_8 = 0
            f_6_9 = f_9_12 = f_4_8 = f_8_12 = 0
            f_14 = C4[x]

            # Reset common neighbor counts
            for node in common_list:
                common[node] = 0
            common_list = []

            self.orbit[x][0] = self.deg[x]

            # x - middle node
            for nx1 in range(len(self.inc[x])):
                y, ey = self.inc[x][nx1]

                for ny in range(len(self.inc[y])):
                    z, ez = self.inc[y][ny]
                    if self.adjacent(x, z):  # triangle
                        if z < y:
                            f_12_14 += tri[ez] - 1
                            f_10_13 += (self.deg[y] - 1 - tri[ez]) + (
                                self.deg[z] - 1 - tri[ez]
                            )
                    else:
                        if common[z] == 0:
                            common_list.append(z)
                        common[z] += 1

                for nx2 in range(nx1 + 1, len(self.inc[x])):
                    z, ez = self.inc[x][nx2]
                    if self.adjacent(y, z):  # triangle
                        self.orbit[x][3] += 1
                        f_13_14 += (tri[ey] - 1) + (tri[ez] - 1)
                        f_11_13 += (self.deg[x] - 1 - tri[ey]) + (
                            self.deg[x] - 1 - tri[ez]
                        )
                    else:  # path
                        self.orbit[x][2] += 1
                        f_7_11 += (self.deg[x] - 1 - tri[ey] - 1) + (
                            self.deg[x] - 1 - tri[ez] - 1
                        )
                        f_5_8 += (self.deg[y] - 1 - tri[ey]) + (
                            self.deg[z] - 1 - tri[ez]
                        )

            # x - side node
            for nx1 in range(len(self.inc[x])):
                y, ey = self.inc[x][nx1]

                for ny in range(len(self.inc[y])):
                    z, ez = self.inc[y][ny]
                    if x == z:
                        continue
                    if not self.adjacent(x, z):  # path
                        self.orbit[x][1] += 1
                        f_6_9 += self.deg[y] - 1 - tri[ey] - 1
                        f_9_12 += tri[ez]
                        f_4_8 += self.deg[z] - 1 - tri[ez]
                        f_8_12 += common[z] - 1

            # Solve system of equations
            self.orbit[x][14] = f_14
            self.orbit[x][13] = (f_13_14 - 6 * f_14) // 2
            self.orbit[x][12] = f_12_14 - 3 * f_14
            self.orbit[x][11] = (f_11_13 - f_13_14 + 6 * f_14) // 2
            self.orbit[x][10] = f_10_13 - f_13_14 + 6 * f_14
            self.orbit[x][9] = (f_9_12 - 2 * f_12_14 + 6 * f_14) // 2
            self.orbit[x][8] = (f_8_12 - 2 * f_12_14 + 6 * f_14) // 2
            self.orbit[x][7] = (f_13_14 + f_7_11 - f_11_13 - 6 * f_14) // 6
            self.orbit[x][6] = (2 * f_12_14 + f_6_9 - f_9_12 - 6 * f_14) // 2
            self.orbit[x][5] = 2 * f_12_14 + f_5_8 - f_8_12 - 6 * f_14
            self.orbit[x][4] = 2 * f_12_14 + f_4_8 - f_8_12 - 6 * f_14

        end_time = time.time()
        print(f"\n{end_time - start_time:.2f}")

        end_time_all = time.time()
        print(f"total: {end_time_all - start_time_all:.2f}")

    def ecount4(self):
        """Count edge orbits of graphlets on max 4 nodes"""
        start_time = time.time()
        start_time_all = start_time

        print("stage 1 - precomputing common nodes")
        # Precompute triangles that span over edges
        tri = [0] * self.m
        for i in range(self.m):
            if i % (self.m // 100 + 1) == 0:
                progress = 100 * i // self.m
                print(f"{progress}%", end="\r")

            x, y = self.edges[i].a, self.edges[i].b
            xi, yi = 0, 0

            while xi < len(self.adj[x]) and yi < len(self.adj[y]):
                if self.adj[x][xi] == self.adj[y][yi]:
                    tri[i] += 1
                    xi += 1
                    yi += 1
                elif self.adj[x][xi] < self.adj[y][yi]:
                    xi += 1
                else:
                    yi += 1

        end_time = time.time()
        print(f"\n{end_time - start_time:.2f}")
        start_time = end_time

        # Count full graphlets
        print("stage 2 - counting full graphlets")
        C4 = [0] * self.m
        neighx = [-1] * self.n  # lookup table - edges to neighbors of x

        for x in range(self.n):
            if x % (self.n // 100 + 1) == 0:
                progress = 100 * x // self.n
                print(f"{progress}%", end="\r")

            # Set up neighbor lookup
            for nx in range(len(self.inc[x])):
                y, xy = self.inc[x][nx]
                neighx[y] = xy

            for nx in range(len(self.inc[x])):
                y, xy = self.inc[x][nx]
                if y >= x:
                    break

                neigh = []
                neigh_edges = []

                for ny in range(len(self.inc[y])):
                    z, yz = self.inc[y][ny]
                    if z >= y:
                        break
                    if neighx[z] == -1:
                        continue

                    xz = neighx[z]
                    neigh.append(z)
                    neigh_edges.append((xz, yz))

                for i in range(len(neigh)):
                    z = neigh[i]
                    xz, yz = neigh_edges[i]

                    for j in range(i + 1, len(neigh)):
                        w = neigh[j]
                        xw, yw = neigh_edges[j]

                        if self.adjacent(z, w):
                            C4[xy] += 1
                            C4[xz] += 1
                            C4[yz] += 1
                            C4[xw] += 1
                            C4[yw] += 1

            # Reset neighbor lookup
            for nx in range(len(self.inc[x])):
                y, xy = self.inc[x][nx]
                neighx[y] = -1

        # Count full graphlets for the smallest edge
        for x in range(self.n):
            if x % (self.n // 100 + 1) == 0:
                progress = 100 * x // self.n
                print(f"{progress}%", end="\r")

            for nx in range(len(self.inc[x]) - 1, -1, -1):
                y, xy = self.inc[x][nx]
                if y <= x:
                    break

                neigh = []
                for ny in range(len(self.adj[y]) - 1, -1, -1):
                    z = self.adj[y][ny]
                    if z <= y:
                        break
                    if not self.adjacent(x, z):
                        continue
                    neigh.append(z)

                for i in range(len(neigh)):
                    z = neigh[i]
                    for j in range(i + 1, len(neigh)):
                        zz = neigh[j]
                        if self.adjacent(z, zz):
                            C4[xy] += 1

        end_time = time.time()
        print(f"\n{end_time - start_time:.2f}")
        start_time = end_time

        # Build systems of equations
        print("stage 3 - building systems of equations")
        common = [0] * self.n
        common_list = []

        for x in range(self.n):
            if x % (self.n // 100 + 1) == 0:
                progress = 100 * x // self.n
                print(f"{progress}%", end="\r")

            # Reset common neighbor counts
            for node in common_list:
                common[node] = 0
            common_list = []

            # Build common neighbor counts
            for nx in range(len(self.adj[x])):
                y = self.adj[x][nx]
                for ny in range(len(self.adj[y])):
                    z = self.adj[y][ny]
                    if z == x:
                        continue
                    if common[z] == 0:
                        common_list.append(z)
                    common[z] += 1

            for nx in range(len(self.inc[x])):
                y, xy = self.inc[x][nx]
                e = xy

                for n1 in range(len(self.inc[x])):
                    z, xz = self.inc[x][n1]
                    if z == y:
                        continue

                    if self.adjacent(y, z):  # triangle
                        if x < y:
                            self.eorbit[e][1] += 1
                            self.eorbit[e][10] += tri[xy] - 1
                            self.eorbit[e][7] += self.deg[z] - 2
                        self.eorbit[e][9] += tri[xz] - 1
                        self.eorbit[e][8] += self.deg[x] - 2

                for n1 in range(len(self.inc[y])):
                    z, yz = self.inc[y][n1]
                    if z == x:
                        continue

                    if not self.adjacent(x, z):  # path x-y-z
                        self.eorbit[e][0] += 1
                        self.eorbit[e][6] += tri[yz]
                        self.eorbit[e][5] += common[z] - 1
                        self.eorbit[e][4] += self.deg[y] - 2
                        self.eorbit[e][3] += self.deg[x] - 1
                        self.eorbit[e][2] += self.deg[z] - 1

        # Solve system of equations
        for e in range(self.m):
            self.eorbit[e][11] = C4[e]
            self.eorbit[e][10] = (self.eorbit[e][10] - 2 * self.eorbit[e][11]) // 2
            self.eorbit[e][9] = self.eorbit[e][9] - 4 * self.eorbit[e][11]
            self.eorbit[e][8] = (
                self.eorbit[e][8]
                - self.eorbit[e][9]
                - 4 * self.eorbit[e][10]
                - 4 * self.eorbit[e][11]
            )
            self.eorbit[e][7] = (
                self.eorbit[e][7] - self.eorbit[e][9] - 2 * self.eorbit[e][11]
            )
            self.eorbit[e][6] = (self.eorbit[e][6] - self.eorbit[e][9]) // 2
            self.eorbit[e][5] = (self.eorbit[e][5] - self.eorbit[e][9]) // 2
            self.eorbit[e][4] = (
                self.eorbit[e][4]
                - 2 * self.eorbit[e][6]
                - self.eorbit[e][8]
                - self.eorbit[e][9]
            ) // 2
            self.eorbit[e][3] = (
                self.eorbit[e][3]
                - 2 * self.eorbit[e][5]
                - self.eorbit[e][8]
                - self.eorbit[e][9]
            ) // 2
            self.eorbit[e][2] = (
                self.eorbit[e][2]
                - 2 * self.eorbit[e][5]
                - 2 * self.eorbit[e][6]
                - self.eorbit[e][9]
            )

        end_time = time.time()
        print(f"\n{end_time - start_time:.2f}")

        end_time_all = time.time()
        print(f"total: {end_time_all - start_time_all:.2f}")

    def count5(self):
        """Count graphlets on max 5 nodes"""
        import time
        start_time_all = time.time()
        start_time = start_time_all

        print("stage 1 - precomputing common nodes")
        frac_prev = -1

        # Clear existing common neighbor caches
        self.common2.clear()
        self.common3.clear()

        # Precompute common nodes for pairs and triples
        for x in range(self.n):
            frac = int(100 * x / self.n)
            if frac != frac_prev:
                print(f"{frac}%\r", end="")
                frac_prev = frac

            for n1 in range(len(self.adj[x])):
                a = self.adj[x][n1]
                for n2 in range(n1 + 1, len(self.adj[x])):
                    b = self.adj[x][n2]
                    ab = Pair(a, b)
                    self.common2[ab] = self.common2.get(ab, 0) + 1

                    for n3 in range(n2 + 1, len(self.adj[x])):
                        c = self.adj[x][n3]
                        st = (
                            self.adjacent(a, b)
                            + self.adjacent(a, c)
                            + self.adjacent(b, c)
                        )
                        if st < 2:
                            continue
                        abc = Triple(a, b, c)
                        self.common3[abc] = self.common3.get(abc, 0) + 1

        # Precompute triangles
        tri = [0] * self.m
        for i in range(self.m):
            x, y = self.edges[i].a, self.edges[i].b
            xi, yi = 0, 0

            while xi < len(self.adj[x]) and yi < len(self.adj[y]):
                if self.adj[x][xi] == self.adj[y][yi]:
                    tri[i] += 1
                    xi += 1
                    yi += 1
                elif self.adj[x][xi] < self.adj[y][yi]:
                    xi += 1
                else:
                    yi += 1

        end_time = time.time()
        print(f"{end_time - start_time:.2f} sec")
        start_time = end_time

        print("stage 2 - counting full graphlets")
        C5 = [0] * self.n
        neighx = [-1] * self.n
        neigh = []
        neigh_edges = []
        neigh2 = []
        neigh2_edges = []
        frac_prev = -1

        for x in range(self.n):
            frac = int(100 * x / self.n)
            if frac != frac_prev:
                print(f"{frac}%\r", end="")
                frac_prev = frac

            # Set up neighx lookup
            for nx in range(len(self.adj[x])):
                y = self.inc[x][nx][0]
                xy = self.inc[x][nx][1]
                neighx[y] = xy

            for nx in range(len(self.adj[x])):
                y = self.inc[x][nx][0]
                xy = self.inc[x][nx][1]
                if y >= x:
                    break

                neigh.clear()
                neigh_edges.clear()

                for ny in range(len(self.adj[y])):
                    z = self.inc[y][ny][0]
                    yz = self.inc[y][ny][1]
                    if z >= y:
                        break
                    if neighx[z] == -1:
                        continue
                    xz = neighx[z]
                    neigh.append(z)
                    neigh_edges.append((xz, yz))

                for i in range(len(neigh)):
                    z = neigh[i]
                    xz, yz = neigh_edges[i]
                    neigh2.clear()
                    neigh2_edges.clear()

                    for j in range(i + 1, len(neigh)):
                        w = neigh[j]
                        xw, yw = neigh_edges[j]
                        if self.adjacent(z, w):
                            neigh2.append(w)
                            zw = self.get_edge_id(z, w)
                            neigh2_edges.append((xw, yw, zw))

                    for i2 in range(len(neigh2)):
                        z2 = neigh2[i2]
                        z2x, z2y, z2z = neigh2_edges[i2]

                        for j2 in range(i2 + 1, len(neigh2)):
                            z3 = neigh2[j2]
                            z3x, z3y, z3z = neigh2_edges[j2]
                            if self.adjacent(z2, z3):
                                C5[x] += 1
                                C5[y] += 1
                                C5[z] += 1
                                C5[z2] += 1
                                C5[z3] += 1

            # Clear neighx lookup
            for nx in range(len(self.adj[x])):
                y = self.inc[x][nx][0]
                neighx[y] = -1

        end_time = time.time()
        print(f"{end_time - start_time:.2f} sec")
        start_time = end_time

        print("stage 3 - building systems of equations")
        common_x = [0] * self.n
        common_x_list = []
        common_a = [0] * self.n
        common_a_list = []
        frac_prev = -1

        for x in range(self.n):
            frac = int(100 * x / self.n)
            if frac != frac_prev:
                print(f"{frac}%\r", end="")
                frac_prev = frac

            # Reset common neighbor counts
            for node in common_x_list:
                common_x[node] = 0
            common_x_list.clear()

            # Smaller graphlets
            self.orbit[x][0] = len(self.adj[x])

            for nx1 in range(len(self.adj[x])):
                a = self.adj[x][nx1]
                for nx2 in range(nx1 + 1, len(self.adj[x])):
                    b = self.adj[x][nx2]
                    if self.adjacent(a, b):
                        self.orbit[x][3] += 1
                    else:
                        self.orbit[x][2] += 1

                for na in range(len(self.adj[a])):
                    b = self.adj[a][na]
                    if b != x and not self.adjacent(x, b):
                        self.orbit[x][1] += 1
                        if common_x[b] == 0:
                            common_x_list.append(b)
                        common_x[b] += 1

            # Initialize f variables for equations
            f_71 = f_70 = f_67 = f_66 = f_58 = f_57 = 0  # 14
            f_69 = f_68 = f_64 = f_61 = f_60 = f_55 = f_48 = f_42 = f_41 = 0  # 13
            f_65 = f_63 = f_59 = f_54 = f_47 = f_46 = f_40 = 0  # 12
            f_62 = f_53 = f_51 = f_50 = f_49 = f_38 = f_37 = f_36 = 0  # 8
            f_44 = f_33 = f_30 = f_26 = 0  # 11
            f_52 = f_43 = f_32 = f_29 = f_25 = 0  # 10
            f_56 = f_45 = f_39 = f_31 = f_28 = f_24 = 0  # 9
            f_35 = f_34 = f_27 = f_18 = f_16 = f_15 = 0  # 4
            f_17 = 0  # 5
            f_22 = f_20 = f_19 = 0  # 6
            f_23 = f_21 = 0  # 7

            for nx1 in range(len(self.adj[x])):
                a = self.inc[x][nx1][0]
                xa = self.inc[x][nx1][1]

                # Reset common_a
                for node in common_a_list:
                    common_a[node] = 0
                common_a_list.clear()

                for na in range(len(self.adj[a])):
                    b = self.adj[a][na]
                    for nb in range(len(self.adj[b])):
                        c = self.adj[b][nb]
                        if c == a or self.adjacent(a, c):
                            continue
                        if common_a[c] == 0:
                            common_a_list.append(c)
                        common_a[c] += 1

                # x = orbit-14 (tetrahedron)
                for nx2 in range(nx1 + 1, len(self.adj[x])):
                    b = self.inc[x][nx2][0]
                    xb = self.inc[x][nx2][1]
                    if not self.adjacent(a, b):
                        continue
                    for nx3 in range(nx2 + 1, len(self.adj[x])):
                        c = self.inc[x][nx3][0]
                        xc = self.inc[x][nx3][1]
                        if not self.adjacent(a, c) or not self.adjacent(b, c):
                            continue
                        self.orbit[x][14] += 1
                        f_70 += self.common3_get(Triple(a, b, c)) - 1
                        f_71 += (tri[xa] > 2 and tri[xb] > 2) and (self.common3_get(Triple(x, a, b)) - 1) or 0
                        f_71 += (tri[xa] > 2 and tri[xc] > 2) and (self.common3_get(Triple(x, a, c)) - 1) or 0
                        f_71 += (tri[xb] > 2 and tri[xc] > 2) and (self.common3_get(Triple(x, b, c)) - 1) or 0
                        f_67 += tri[xa] - 2 + tri[xb] - 2 + tri[xc] - 2
                        f_66 += self.common2_get(Pair(a, b)) - 2
                        f_66 += self.common2_get(Pair(a, c)) - 2
                        f_66 += self.common2_get(Pair(b, c)) - 2
                        f_58 += len(self.adj[x]) - 3
                        f_57 += len(self.adj[a]) - 3 + len(self.adj[b]) - 3 + len(self.adj[c]) - 3

                # x = orbit-13 (diamond)
                for nx2 in range(len(self.adj[x])):
                    b = self.inc[x][nx2][0]
                    xb = self.inc[x][nx2][1]
                    if not self.adjacent(a, b):
                        continue
                    for nx3 in range(nx2 + 1, len(self.adj[x])):
                        c = self.inc[x][nx3][0]
                        xc = self.inc[x][nx3][1]
                        if not self.adjacent(a, c) or self.adjacent(b, c):
                            continue
                        self.orbit[x][13] += 1
                        f_69 += (tri[xb] > 1 and tri[xc] > 1) and (self.common3_get(Triple(x, b, c)) - 1) or 0
                        f_68 += self.common3_get(Triple(a, b, c)) - 1
                        f_64 += self.common2_get(Pair(b, c)) - 2
                        f_61 += tri[xb] - 1 + tri[xc] - 1
                        f_60 += self.common2_get(Pair(a, b)) - 1
                        f_60 += self.common2_get(Pair(a, c)) - 1
                        f_55 += tri[xa] - 2
                        f_48 += len(self.adj[b]) - 2 + len(self.adj[c]) - 2
                        f_42 += len(self.adj[x]) - 3
                        f_41 += len(self.adj[a]) - 3

                # x = orbit-12 (diamond)
                for nx2 in range(nx1 + 1, len(self.adj[x])):
                    b = self.inc[x][nx2][0]
                    xb = self.inc[x][nx2][1]
                    if not self.adjacent(a, b):
                        continue
                    for na in range(len(self.adj[a])):
                        c = self.inc[a][na][0]
                        ac = self.inc[a][na][1]
                        if c == x or self.adjacent(x, c) or not self.adjacent(b, c):
                            continue
                        self.orbit[x][12] += 1
                        f_65 += (tri[ac] > 1) and self.common3_get(Triple(a, b, c)) or 0
                        f_63 += common_x[c] - 2
                        f_59 += tri[ac] - 1 + self.common2_get(Pair(b, c)) - 1
                        f_54 += self.common2_get(Pair(a, b)) - 2
                        f_47 += len(self.adj[x]) - 2
                        f_46 += len(self.adj[c]) - 2
                        f_40 += len(self.adj[a]) - 3 + len(self.adj[b]) - 3

                # x = orbit-8 (cycle)
                for nx2 in range(nx1 + 1, len(self.adj[x])):
                    b = self.inc[x][nx2][0]
                    xb = self.inc[x][nx2][1]
                    if self.adjacent(a, b):
                        continue
                    for na in range(len(self.adj[a])):
                        c = self.inc[a][na][0]
                        ac = self.inc[a][na][1]
                        if c == x or self.adjacent(x, c) or not self.adjacent(b, c):
                            continue
                        self.orbit[x][8] += 1
                        f_62 += (tri[ac] > 0) and self.common3_get(Triple(a, b, c)) or 0
                        f_53 += tri[xa] + tri[xb]
                        f_51 += tri[ac] + self.common2_get(Pair(c, b))
                        f_50 += common_x[c] - 2
                        f_49 += common_a[b] - 2
                        f_38 += len(self.adj[x]) - 2
                        f_37 += len(self.adj[a]) - 2 + len(self.adj[b]) - 2
                        f_36 += len(self.adj[c]) - 2

                # x = orbit-11 (paw)
                for nx2 in range(nx1 + 1, len(self.adj[x])):
                    b = self.inc[x][nx2][0]
                    xb = self.inc[x][nx2][1]
                    if not self.adjacent(a, b):
                        continue
                    for nx3 in range(len(self.adj[x])):
                        c = self.inc[x][nx3][0]
                        xc = self.inc[x][nx3][1]
                        if c == a or c == b or self.adjacent(a, c) or self.adjacent(b, c):
                            continue
                        self.orbit[x][11] += 1
                        f_44 += tri[xc]
                        f_33 += len(self.adj[x]) - 3
                        f_30 += len(self.adj[c]) - 1
                        f_26 += len(self.adj[a]) - 2 + len(self.adj[b]) - 2

                # x = orbit-10 (paw)
                for nx2 in range(len(self.adj[x])):
                    b = self.inc[x][nx2][0]
                    xb = self.inc[x][nx2][1]
                    if not self.adjacent(a, b):
                        continue
                    for nb in range(len(self.adj[b])):
                        c = self.inc[b][nb][0]
                        bc = self.inc[b][nb][1]
                        if c == x or c == a or self.adjacent(a, c) or self.adjacent(x, c):
                            continue
                        self.orbit[x][10] += 1
                        f_52 += common_a[c] - 1
                        f_43 += tri[bc]
                        f_32 += len(self.adj[b]) - 3
                        f_29 += len(self.adj[c]) - 1
                        f_25 += len(self.adj[a]) - 2

                # x = orbit-9 (paw)
                for na1 in range(len(self.adj[a])):
                    b = self.inc[a][na1][0]
                    ab = self.inc[a][na1][1]
                    if b == x or self.adjacent(x, b):
                        continue
                    for na2 in range(na1 + 1, len(self.adj[a])):
                        c = self.inc[a][na2][0]
                        ac = self.inc[a][na2][1]
                        if c == x or not self.adjacent(b, c) or self.adjacent(x, c):
                            continue
                        self.orbit[x][9] += 1
                        f_56 += (tri[ab] > 1 and tri[ac] > 1) and self.common3_get(Triple(a, b, c)) or 0
                        f_45 += self.common2_get(Pair(b, c)) - 1
                        f_39 += tri[ab] - 1 + tri[ac] - 1
                        f_31 += len(self.adj[a]) - 3
                        f_28 += len(self.adj[x]) - 1
                        f_24 += len(self.adj[b]) - 2 + len(self.adj[c]) - 2

                # x = orbit-4 (path)
                for na in range(len(self.adj[a])):
                    b = self.inc[a][na][0]
                    ab = self.inc[a][na][1]
                    if b == x or self.adjacent(x, b):
                        continue
                    for nb in range(len(self.adj[b])):
                        c = self.inc[b][nb][0]
                        bc = self.inc[b][nb][1]
                        if c == a or self.adjacent(a, c) or self.adjacent(x, c):
                            continue
                        self.orbit[x][4] += 1
                        f_35 += common_a[c] - 1
                        f_34 += common_x[c]
                        f_27 += tri[bc]
                        f_18 += len(self.adj[b]) - 2
                        f_16 += len(self.adj[x]) - 1
                        f_15 += len(self.adj[c]) - 1

                # x = orbit-5 (path)
                for nx2 in range(len(self.adj[x])):
                    b = self.inc[x][nx2][0]
                    xb = self.inc[x][nx2][1]
                    if b == a or self.adjacent(a, b):
                        continue
                    for nb in range(len(self.adj[b])):
                        c = self.inc[b][nb][0]
                        bc = self.inc[b][nb][1]
                        if c == x or self.adjacent(a, c) or self.adjacent(x, c):
                            continue
                        self.orbit[x][5] += 1
                        f_17 += len(self.adj[a]) - 1

                # x = orbit-6 (claw)
                for na1 in range(len(self.adj[a])):
                    b = self.inc[a][na1][0]
                    ab = self.inc[a][na1][1]
                    if b == x or self.adjacent(x, b):
                        continue
                    for na2 in range(na1 + 1, len(self.adj[a])):
                        c = self.inc[a][na2][0]
                        ac = self.inc[a][na2][1]
                        if c == x or self.adjacent(x, c) or self.adjacent(b, c):
                            continue
                        self.orbit[x][6] += 1
                        f_22 += len(self.adj[a]) - 3
                        f_20 += len(self.adj[x]) - 1
                        f_19 += len(self.adj[b]) - 1 + len(self.adj[c]) - 1

                # x = orbit-7 (claw)
                for nx2 in range(nx1 + 1, len(self.adj[x])):
                    b = self.inc[x][nx2][0]
                    xb = self.inc[x][nx2][1]
                    if self.adjacent(a, b):
                        continue
                    for nx3 in range(nx2 + 1, len(self.adj[x])):
                        c = self.inc[x][nx3][0]
                        xc = self.inc[x][nx3][1]
                        if self.adjacent(a, c) or self.adjacent(b, c):
                            continue
                        self.orbit[x][7] += 1
                        f_23 += len(self.adj[x]) - 3
                        f_21 += len(self.adj[a]) - 1 + len(self.adj[b]) - 1 + len(self.adj[c]) - 1

            # Solve equations
            self.orbit[x][72] = C5[x]
            self.orbit[x][71] = (f_71 - 12 * self.orbit[x][72]) // 2
            self.orbit[x][70] = (f_70 - 4 * self.orbit[x][72])
            self.orbit[x][69] = (f_69 - 2 * self.orbit[x][71]) // 4
            self.orbit[x][68] = (f_68 - 2 * self.orbit[x][71])
            self.orbit[x][67] = (f_67 - 12 * self.orbit[x][72] - 4 * self.orbit[x][71])
            self.orbit[x][66] = (f_66 - 12 * self.orbit[x][72] - 2 * self.orbit[x][71] - 3 * self.orbit[x][70])
            self.orbit[x][65] = (f_65 - 3 * self.orbit[x][70]) // 2
            self.orbit[x][64] = (f_64 - 2 * self.orbit[x][71] - 4 * self.orbit[x][69] - 1 * self.orbit[x][68])
            self.orbit[x][63] = (f_63 - 3 * self.orbit[x][70] - 2 * self.orbit[x][68])
            self.orbit[x][62] = (f_62 - 1 * self.orbit[x][68]) // 2
            self.orbit[x][61] = (f_61 - 4 * self.orbit[x][71] - 8 * self.orbit[x][69] - 2 * self.orbit[x][67]) // 2
            self.orbit[x][60] = (f_60 - 4 * self.orbit[x][71] - 2 * self.orbit[x][68] - 2 * self.orbit[x][67])
            self.orbit[x][59] = (f_59 - 6 * self.orbit[x][70] - 2 * self.orbit[x][68] - 4 * self.orbit[x][65])
            self.orbit[x][58] = (f_58 - 4 * self.orbit[x][72] - 2 * self.orbit[x][71] - 1 * self.orbit[x][67])
            self.orbit[x][57] = (f_57 - 12 * self.orbit[x][72] - 4 * self.orbit[x][71] - 3 * self.orbit[x][70] - 1 * self.orbit[x][67] - 2 * self.orbit[x][66])
            self.orbit[x][56] = (f_56 - 2 * self.orbit[x][65]) // 3
            self.orbit[x][55] = (f_55 - 2 * self.orbit[x][71] - 2 * self.orbit[x][67]) // 3
            self.orbit[x][54] = (f_54 - 3 * self.orbit[x][70] - 1 * self.orbit[x][66] - 2 * self.orbit[x][65]) // 2
            self.orbit[x][53] = (f_53 - 2 * self.orbit[x][68] - 2 * self.orbit[x][64] - 2 * self.orbit[x][63])
            self.orbit[x][52] = (f_52 - 2 * self.orbit[x][66] - 2 * self.orbit[x][64] - 1 * self.orbit[x][59]) // 2
            self.orbit[x][51] = (f_51 - 2 * self.orbit[x][68] - 2 * self.orbit[x][63] - 4 * self.orbit[x][62])
            self.orbit[x][50] = (f_50 - 1 * self.orbit[x][68] - 2 * self.orbit[x][63]) // 3
            self.orbit[x][49] = (f_49 - 1 * self.orbit[x][68] - 1 * self.orbit[x][64] - 2 * self.orbit[x][62]) // 2
            self.orbit[x][48] = (f_48 - 4 * self.orbit[x][71] - 8 * self.orbit[x][69] - 2 * self.orbit[x][68] - 2 * self.orbit[x][67] - 2 * self.orbit[x][64] - 2 * self.orbit[x][61] - 1 * self.orbit[x][60])
            self.orbit[x][47] = (f_47 - 3 * self.orbit[x][70] - 2 * self.orbit[x][68] - 1 * self.orbit[x][66] - 1 * self.orbit[x][63] - 1 * self.orbit[x][60])
            self.orbit[x][46] = (f_46 - 3 * self.orbit[x][70] - 2 * self.orbit[x][68] - 2 * self.orbit[x][65] - 1 * self.orbit[x][63] - 1 * self.orbit[x][59])
            self.orbit[x][45] = (f_45 - 2 * self.orbit[x][65] - 2 * self.orbit[x][62] - 3 * self.orbit[x][56])
            self.orbit[x][44] = (f_44 - 1 * self.orbit[x][67] - 2 * self.orbit[x][61]) // 4
            self.orbit[x][43] = (f_43 - 2 * self.orbit[x][66] - 1 * self.orbit[x][60] - 1 * self.orbit[x][59]) // 2
            self.orbit[x][42] = (f_42 - 2 * self.orbit[x][71] - 4 * self.orbit[x][69] - 2 * self.orbit[x][67] - 2 * self.orbit[x][61] - 3 * self.orbit[x][55])
            self.orbit[x][41] = (f_41 - 2 * self.orbit[x][71] - 1 * self.orbit[x][68] - 2 * self.orbit[x][67] - 1 * self.orbit[x][60] - 3 * self.orbit[x][55])
            self.orbit[x][40] = (f_40 - 6 * self.orbit[x][70] - 2 * self.orbit[x][68] - 2 * self.orbit[x][66] - 4 * self.orbit[x][65] - 1 * self.orbit[x][60] - 1 * self.orbit[x][59] - 4 * self.orbit[x][54])
            self.orbit[x][39] = (f_39 - 4 * self.orbit[x][65] - 1 * self.orbit[x][59] - 6 * self.orbit[x][56]) // 2
            self.orbit[x][38] = (f_38 - 1 * self.orbit[x][68] - 1 * self.orbit[x][64] - 2 * self.orbit[x][63] - 1 * self.orbit[x][53] - 3 * self.orbit[x][50])
            self.orbit[x][37] = (f_37 - 2 * self.orbit[x][68] - 2 * self.orbit[x][64] - 2 * self.orbit[x][63] - 4 * self.orbit[x][62] - 1 * self.orbit[x][53] - 1 * self.orbit[x][51] - 4 * self.orbit[x][49])
            self.orbit[x][36] = (f_36 - 1 * self.orbit[x][68] - 2 * self.orbit[x][63] - 2 * self.orbit[x][62] - 1 * self.orbit[x][51] - 3 * self.orbit[x][50])
            self.orbit[x][35] = (f_35 - 1 * self.orbit[x][59] - 2 * self.orbit[x][52] - 2 * self.orbit[x][45]) // 2
            self.orbit[x][34] = (f_34 - 1 * self.orbit[x][59] - 2 * self.orbit[x][52] - 1 * self.orbit[x][51]) // 2
            self.orbit[x][33] = (f_33 - 1 * self.orbit[x][67] - 2 * self.orbit[x][61] - 3 * self.orbit[x][58] - 4 * self.orbit[x][44] - 2 * self.orbit[x][42]) // 2
            self.orbit[x][32] = (f_32 - 2 * self.orbit[x][66] - 1 * self.orbit[x][60] - 1 * self.orbit[x][59] - 2 * self.orbit[x][57] - 2 * self.orbit[x][43] - 2 * self.orbit[x][41] - 1 * self.orbit[x][40]) // 2
            self.orbit[x][31] = (f_31 - 2 * self.orbit[x][65] - 1 * self.orbit[x][59] - 3 * self.orbit[x][56] - 1 * self.orbit[x][43] - 2 * self.orbit[x][39])
            self.orbit[x][30] = (f_30 - 1 * self.orbit[x][67] - 1 * self.orbit[x][63] - 2 * self.orbit[x][61] - 1 * self.orbit[x][53] - 4 * self.orbit[x][44])
            self.orbit[x][29] = (f_29 - 2 * self.orbit[x][66] - 2 * self.orbit[x][64] - 1 * self.orbit[x][60] - 1 * self.orbit[x][59] - 1 * self.orbit[x][53] - 2 * self.orbit[x][52] - 2 * self.orbit[x][43])
            self.orbit[x][28] = (f_28 - 2 * self.orbit[x][65] - 2 * self.orbit[x][62] - 1 * self.orbit[x][59] - 1 * self.orbit[x][51] - 1 * self.orbit[x][43])
            self.orbit[x][27] = (f_27 - 1 * self.orbit[x][59] - 1 * self.orbit[x][51] - 2 * self.orbit[x][45]) // 2
            self.orbit[x][26] = (f_26 - 2 * self.orbit[x][67] - 2 * self.orbit[x][63] - 2 * self.orbit[x][61] - 6 * self.orbit[x][58] - 1 * self.orbit[x][53] - 2 * self.orbit[x][47] - 2 * self.orbit[x][42])
            self.orbit[x][25] = (f_25 - 2 * self.orbit[x][66] - 2 * self.orbit[x][64] - 1 * self.orbit[x][59] - 2 * self.orbit[x][57] - 2 * self.orbit[x][52] - 1 * self.orbit[x][48] - 1 * self.orbit[x][40]) // 2
            self.orbit[x][24] = (f_24 - 4 * self.orbit[x][65] - 4 * self.orbit[x][62] - 1 * self.orbit[x][59] - 6 * self.orbit[x][56] - 1 * self.orbit[x][51] - 2 * self.orbit[x][45] - 2 * self.orbit[x][39])
            self.orbit[x][23] = (f_23 - 1 * self.orbit[x][55] - 1 * self.orbit[x][42] - 2 * self.orbit[x][33]) // 4
            self.orbit[x][22] = (f_22 - 2 * self.orbit[x][54] - 1 * self.orbit[x][40] - 1 * self.orbit[x][39] - 1 * self.orbit[x][32] - 2 * self.orbit[x][31]) // 3
            self.orbit[x][21] = (f_21 - 3 * self.orbit[x][55] - 3 * self.orbit[x][50] - 2 * self.orbit[x][42] - 2 * self.orbit[x][38] - 2 * self.orbit[x][33])
            self.orbit[x][20] = (f_20 - 2 * self.orbit[x][54] - 2 * self.orbit[x][49] - 1 * self.orbit[x][40] - 1 * self.orbit[x][37] - 1 * self.orbit[x][32])
            self.orbit[x][19] = (f_19 - 4 * self.orbit[x][54] - 4 * self.orbit[x][49] - 1 * self.orbit[x][40] - 2 * self.orbit[x][39] - 1 * self.orbit[x][37] - 2 * self.orbit[x][35] - 2 * self.orbit[x][31])
            self.orbit[x][18] = (f_18 - 1 * self.orbit[x][59] - 1 * self.orbit[x][51] - 2 * self.orbit[x][46] - 2 * self.orbit[x][45] - 2 * self.orbit[x][36] - 2 * self.orbit[x][27] - 1 * self.orbit[x][24]) // 2
            self.orbit[x][17] = (f_17 - 1 * self.orbit[x][60] - 1 * self.orbit[x][53] - 1 * self.orbit[x][51] - 1 * self.orbit[x][48] - 1 * self.orbit[x][37] - 2 * self.orbit[x][34] - 2 * self.orbit[x][30]) // 2
            self.orbit[x][16] = (f_16 - 1 * self.orbit[x][59] - 2 * self.orbit[x][52] - 1 * self.orbit[x][51] - 2 * self.orbit[x][46] - 2 * self.orbit[x][36] - 2 * self.orbit[x][34] - 1 * self.orbit[x][29])
            self.orbit[x][15] = (f_15 - 1 * self.orbit[x][59] - 2 * self.orbit[x][52] - 1 * self.orbit[x][51] - 2 * self.orbit[x][45] - 2 * self.orbit[x][35] - 2 * self.orbit[x][34] - 2 * self.orbit[x][27])

        end_time = time.time()
        print(f"{end_time - start_time:.2f} sec")

        end_time_all = time.time()
        print(f"total: {end_time_all - start_time_all:.2f} sec")

    def ecount5(self):
        """Count edge orbits of graphlets on max 5 nodes"""
        import time
        start_time_all = time.time()
        start_time = start_time_all
        
        print("stage 1 - precomputing common nodes")
        frac_prev = -1
        
        # Clear existing common neighbor caches
        self.common2.clear()
        self.common3.clear()
        
        # Precompute common nodes
        for x in range(self.n):
            frac = int(100 * x / self.n)
            if frac != frac_prev:
                print(f"{frac}%\r", end="")
                frac_prev = frac
                
            for n1 in range(len(self.adj[x])):
                a = self.adj[x][n1]
                for n2 in range(n1 + 1, len(self.adj[x])):
                    b = self.adj[x][n2]
                    ab = Pair(a, b)
                    self.common2[ab] = self.common2.get(ab, 0) + 1
                    
                    for n3 in range(n2 + 1, len(self.adj[x])):
                        c = self.adj[x][n3]
                        st = (self.adjacent(a, b) + 
                              self.adjacent(a, c) + 
                              self.adjacent(b, c))
                        if st < 2:
                            continue
                        abc = Triple(a, b, c)
                        self.common3[abc] = self.common3.get(abc, 0) + 1
                        
        # Precompute triangles that span over edges
        tri = [0] * self.m
        for i in range(self.m):
            x, y = self.edges[i].a, self.edges[i].b
            xi, yi = 0, 0
            while xi < len(self.adj[x]) and yi < len(self.adj[y]):
                if self.adj[x][xi] == self.adj[y][yi]:
                    tri[i] += 1
                    xi += 1
                    yi += 1
                elif self.adj[x][xi] < self.adj[y][yi]:
                    xi += 1
                else:
                    yi += 1
                    
        end_time = time.time()
        print(f"{end_time - start_time:.2f} sec")
        start_time = end_time
        
        print("stage 2 - counting full graphlets")
        C5 = [0] * self.m
        neighx = [-1] * self.n  # lookup table - edges to neighbors of x
        neigh = []  # common neighbors of x and y
        neigh_edges = []  # list of common neighbor edge IDs
        neigh2 = []  # second-level neighbors
        neigh2_edges = []  # second-level neighbor edge IDs
        frac_prev = -1
        
        for x in range(self.n):
            frac = int(100 * x / self.n)
            if frac != frac_prev:
                print(f"{frac}%\r", end="")
                frac_prev = frac
                
            # Set up neighx lookup
            for nx in range(len(self.adj[x])):
                y = self.inc[x][nx][0]
                xy = self.inc[x][nx][1]
                neighx[y] = xy
                
            for nx in range(len(self.adj[x])):
                y = self.inc[x][nx][0]
                xy = self.inc[x][nx][1]
                if y >= x:
                    break
                    
                neigh.clear()
                neigh_edges.clear()
                
                for ny in range(len(self.adj[y])):
                    z = self.inc[y][ny][0]
                    yz = self.inc[y][ny][1]
                    if z >= y:
                        break
                    if neighx[z] == -1:
                        continue
                    xz = neighx[z]
                    neigh.append(z)
                    neigh_edges.append((xz, yz))
                    
                for i in range(len(neigh)):
                    z = neigh[i]
                    xz, yz = neigh_edges[i]
                    neigh2.clear()
                    neigh2_edges.clear()
                    
                    for j in range(i + 1, len(neigh)):
                        w = neigh[j]
                        xw, yw = neigh_edges[j]
                        if self.adjacent(z, w):
                            neigh2.append(w)
                            zw = self.get_edge_id(z, w)
                            neigh2_edges.append((xw, yw, zw))
                            
                    for i2 in range(len(neigh2)):
                        z2 = neigh2[i2]
                        z2x, z2y, z2z = neigh2_edges[i2]
                        
                        for j2 in range(i2 + 1, len(neigh2)):
                            z3 = neigh2[j2]
                            z3x, z3y, z3z = neigh2_edges[j2]
                            if self.adjacent(z2, z3):
                                zid = self.get_edge_id(z2, z3)
                                C5[xy] += 1
                                C5[xz] += 1
                                C5[yz] += 1
                                C5[z2x] += 1
                                C5[z2y] += 1
                                C5[z2z] += 1
                                C5[z3x] += 1
                                C5[z3y] += 1
                                C5[z3z] += 1
                                C5[zid] += 1
                                
            # Clear neighx lookup
            for nx in range(len(self.adj[x])):
                y = self.inc[x][nx][0]
                neighx[y] = -1
                
        end_time = time.time()
        print(f"{end_time - start_time:.2f}")
        start_time = end_time
        
        print("stage 3 - building systems of equations")
        common_x = [0] * self.n
        common_x_list = []
        common_y = [0] * self.n
        common_y_list = []
        frac_prev = -1
        
        for x in range(self.n):
            frac = int(100 * x / self.n)
            if frac != frac_prev:
                print(f"{frac}%\r", end="")
                frac_prev = frac
                
            # Clear common_x
            for i in common_x_list:
                common_x[i] = 0
            common_x_list.clear()
            
            # Build common nodes of x
            for nx in range(len(self.adj[x])):
                a = self.adj[x][nx]
                for na in range(len(self.adj[a])):
                    z = self.adj[a][na]
                    if z == x:
                        continue
                    if common_x[z] == 0:
                        common_x_list.append(z)
                    common_x[z] += 1
                    
            for nx in range(len(self.adj[x])):
                y = self.inc[x][nx][0]
                xy = self.inc[x][nx][1]
                e = xy
                if y >= x:
                    break
                    
                # Clear common_y
                for i in common_y_list:
                    common_y[i] = 0
                common_y_list.clear()
                
                # Build common nodes of y
                for ny in range(len(self.adj[y])):
                    a = self.adj[y][ny]
                    for na in range(len(self.adj[a])):
                        z = self.adj[a][na]
                        if z == y:
                            continue
                        if common_y[z] == 0:
                            common_y_list.append(z)
                        common_y[z] += 1
                        
                # Initialize all orbit counts to 0
                for orbit_idx in range(68):
                    self.eorbit[e][orbit_idx] = 0
                        
                # Initialize f variables
                f_66 = f_65 = f_62 = f_61 = f_60 = f_51 = f_50 = 0
                f_64 = f_58 = f_55 = f_48 = f_41 = f_35 = 0
                f_63 = f_59 = f_57 = f_54 = f_53 = f_52 = f_47 = f_40 = f_39 = f_34 = f_33 = 0
                f_45 = f_36 = f_26 = f_23 = f_19 = 0
                f_49 = f_38 = f_37 = f_32 = f_25 = f_22 = f_18 = 0
                f_56 = f_46 = f_44 = f_43 = f_42 = f_31 = f_30 = 0
                f_27 = f_17 = f_15 = 0
                f_20 = f_16 = f_13 = 0
                f_29 = f_28 = f_24 = f_21 = f_14 = f_12 = 0
                
                # Smaller (3-node) graphlets
                self.orbit[x][0] = len(self.adj[x])
                for nx1 in range(len(self.adj[x])):
                    z = self.adj[x][nx1]
                    if z == y:
                        continue
                    if self.adjacent(y, z):
                        self.eorbit[e][1] += 1
                    else:
                        self.eorbit[e][0] += 1
                        
                for ny in range(len(self.adj[y])):
                    z = self.adj[y][ny]
                    if z == x:
                        continue
                    if not self.adjacent(x, z):
                        self.eorbit[e][0] += 1

                # Edge-orbit 11 = (14,14)
                for nx1 in range(len(self.adj[x])):
                    a = self.adj[x][nx1]
                    xa = self.inc[x][nx1][1]
                    if a == y or not self.adjacent(y, a):
                        continue
                    for nx2 in range(nx1 + 1, len(self.adj[x])):
                        b = self.adj[x][nx2]
                        xb = self.inc[x][nx2][1]
                        if b == y or not self.adjacent(y, b) or not self.adjacent(a, b):
                            continue
                        ya = self.get_edge_id(y, a)
                        yb = self.get_edge_id(y, b)
                        ab = self.get_edge_id(a, b)
                        self.eorbit[e][11] += 1
                        f_66 += self.common3_get(Triple(x, y, a)) - 1
                        f_66 += self.common3_get(Triple(x, y, b)) - 1
                        f_65 += self.common3_get(Triple(a, b, x)) - 1
                        f_65 += self.common3_get(Triple(a, b, y)) - 1
                        f_62 += tri[xy] - 2
                        f_61 += (tri[xa] - 2) + (tri[xb] - 2) + (tri[ya] - 2) + (tri[yb] - 2)
                        f_60 += tri[ab] - 2
                        f_51 += (len(self.adj[x]) - 3) + (len(self.adj[y]) - 3)
                        f_50 += (len(self.adj[a]) - 3) + (len(self.adj[b]) - 3)

                # Edge-orbit 10 = (13,13)
                for nx1 in range(len(self.adj[x])):
                    a = self.adj[x][nx1]
                    xa = self.inc[x][nx1][1]
                    if a == y or not self.adjacent(y, a):
                        continue
                    for nx2 in range(nx1 + 1, len(self.adj[x])):
                        b = self.adj[x][nx2]
                        xb = self.inc[x][nx2][1]
                        if b == y or not self.adjacent(y, b) or self.adjacent(a, b):
                            continue
                        ya = self.get_edge_id(y, a)
                        yb = self.get_edge_id(y, b)
                        self.eorbit[e][10] += 1
                        f_64 += self.common3_get(Triple(a, b, x)) - 1
                        f_64 += self.common3_get(Triple(a, b, y)) - 1
                        f_58 += self.common2_get(Pair(a, b)) - 2
                        f_55 += (tri[xa] - 1) + (tri[xb] - 1) + (tri[ya] - 1) + (tri[yb] - 1)
                        f_48 += tri[xy] - 2
                        f_41 += (len(self.adj[a]) - 2) + (len(self.adj[b]) - 2)
                        f_35 += (len(self.adj[x]) - 3) + (len(self.adj[y]) - 3)

                # Edge-orbit 9 = (12,13)
                for nx in range(len(self.adj[x])):
                    a = self.adj[x][nx]
                    xa = self.inc[x][nx][1]
                    if a == y:
                        continue
                    for ny in range(len(self.adj[y])):
                        b = self.adj[y][ny]
                        yb = self.inc[y][ny][1]
                        if b == x or not self.adjacent(a, b):
                            continue
                        adj_ya = self.adjacent(y, a)
                        adj_xb = self.adjacent(x, b)
                        if adj_ya + adj_xb != 1:
                            continue
                        ab = self.get_edge_id(a, b)
                        self.eorbit[e][9] += 1
                        if adj_xb:
                            xb = self.get_edge_id(x, b)
                            f_63 += self.common3_get(Triple(a, b, y)) - 1
                            f_59 += self.common3_get(Triple(a, b, x))
                            f_57 += common_y[a] - 2
                            f_54 += tri[yb] - 1
                            f_53 += tri[xa] - 1
                            f_47 += tri[xb] - 2
                            f_40 += len(self.adj[y]) - 2
                            f_39 += len(self.adj[a]) - 2
                            f_34 += len(self.adj[x]) - 3
                            f_33 += len(self.adj[b]) - 3
                        elif adj_ya:
                            ya = self.get_edge_id(y, a)
                            f_63 += self.common3_get(Triple(a, b, x)) - 1
                            f_59 += self.common3_get(Triple(a, b, y))
                            f_57 += common_x[b] - 2
                            f_54 += tri[xa] - 1
                            f_53 += tri[yb] - 1
                            f_47 += tri[ya] - 2
                            f_40 += len(self.adj[x]) - 2
                            f_39 += len(self.adj[b]) - 2
                            f_34 += len(self.adj[y]) - 3
                            f_33 += len(self.adj[a]) - 3
                        f_52 += tri[ab] - 1

                # Edge-orbit 8 = (10,11)
                for nx in range(len(self.adj[x])):
                    a = self.adj[x][nx]
                    if a == y or not self.adjacent(y, a):
                        continue
                    for nx1 in range(len(self.adj[x])):
                        b = self.adj[x][nx1]
                        if b == y or b == a or self.adjacent(y, b) or self.adjacent(a, b):
                            continue
                        self.eorbit[e][8] += 1
                    for ny1 in range(len(self.adj[y])):
                        b = self.adj[y][ny1]
                        if b == x or b == a or self.adjacent(x, b) or self.adjacent(a, b):
                            continue
                        self.eorbit[e][8] += 1

                # Edge-orbit 7 = (10,10)
                for nx in range(len(self.adj[x])):
                    a = self.adj[x][nx]
                    if a == y or not self.adjacent(y, a):
                        continue
                    for na in range(len(self.adj[a])):
                        b = self.adj[a][na]
                        ab = self.inc[a][na][1]
                        if b == x or b == y or self.adjacent(x, b) or self.adjacent(y, b):
                            continue
                        self.eorbit[e][7] += 1
                        f_45 += common_x[b] - 1
                        f_45 += common_y[b] - 1
                        f_36 += tri[ab]
                        f_26 += len(self.adj[a]) - 3
                        f_23 += len(self.adj[b]) - 1
                        f_19 += (len(self.adj[x]) - 2) + (len(self.adj[y]) - 2)

                # Edge-orbit 6 = (9,11)
                for ny1 in range(len(self.adj[y])):
                    a = self.adj[y][ny1]
                    ya = self.inc[y][ny1][1]
                    if a == x or self.adjacent(x, a):
                        continue
                    for ny2 in range(ny1 + 1, len(self.adj[y])):
                        b = self.adj[y][ny2]
                        yb = self.inc[y][ny2][1]
                        if b == x or self.adjacent(x, b) or not self.adjacent(a, b):
                            continue
                        ab = self.get_edge_id(a, b)
                        self.eorbit[e][6] += 1
                        f_49 += self.common3_get(Triple(y, a, b))
                        f_38 += tri[ab] - 1
                        f_37 += tri[xy]
                        f_32 += (tri[ya] - 1) + (tri[yb] - 1)
                        f_25 += len(self.adj[y]) - 3
                        f_22 += len(self.adj[x]) - 1
                        f_18 += (len(self.adj[a]) - 2) + (len(self.adj[b]) - 2)
                for nx1 in range(len(self.adj[x])):
                    a = self.adj[x][nx1]
                    xa = self.inc[x][nx1][1]
                    if a == y or self.adjacent(y, a):
                        continue
                    for nx2 in range(nx1 + 1, len(self.adj[x])):
                        b = self.adj[x][nx2]
                        xb = self.inc[x][nx2][1]
                        if b == y or self.adjacent(y, b) or not self.adjacent(a, b):
                            continue
                        ab = self.get_edge_id(a, b)
                        self.eorbit[e][6] += 1
                        f_49 += self.common3_get(Triple(x, a, b))
                        f_38 += tri[ab] - 1
                        f_37 += tri[xy]
                        f_32 += (tri[xa] - 1) + (tri[xb] - 1)
                        f_25 += len(self.adj[x]) - 3
                        f_22 += len(self.adj[y]) - 1
                        f_18 += (len(self.adj[a]) - 2) + (len(self.adj[b]) - 2)

                # Edge-orbit 5 = (8,8)
                for nx in range(len(self.adj[x])):
                    a = self.adj[x][nx]
                    xa = self.inc[x][nx][1]
                    if a == y or self.adjacent(y, a):
                        continue
                    for ny in range(len(self.adj[y])):
                        b = self.adj[y][ny]
                        yb = self.inc[y][ny][1]
                        if b == x or self.adjacent(x, b) or not self.adjacent(a, b):
                            continue
                        ab = self.get_edge_id(a, b)
                        self.eorbit[e][5] += 1
                        f_56 += self.common3_get(Triple(x, a, b))
                        f_56 += self.common3_get(Triple(y, a, b))
                        f_46 += tri[xy]
                        f_44 += tri[xa] + tri[yb]
                        f_43 += tri[ab]
                        f_42 += common_x[b] - 2
                        f_42 += common_y[a] - 2
                        f_31 += (len(self.adj[x]) - 2) + (len(self.adj[y]) - 2)
                        f_30 += (len(self.adj[a]) - 2) + (len(self.adj[b]) - 2)

                # Edge-orbit 4 = (6,7)
                for ny1 in range(len(self.adj[y])):
                    a = self.adj[y][ny1]
                    if a == x or self.adjacent(x, a):
                        continue
                    for ny2 in range(ny1 + 1, len(self.adj[y])):
                        b = self.adj[y][ny2]
                        if b == x or self.adjacent(x, b) or self.adjacent(a, b):
                            continue
                        self.eorbit[e][4] += 1
                        f_27 += tri[xy]
                        f_17 += len(self.adj[y]) - 3
                        f_15 += (len(self.adj[a]) - 1) + (len(self.adj[b]) - 1)
                for nx1 in range(len(self.adj[x])):
                    a = self.adj[x][nx1]
                    if a == y or self.adjacent(y, a):
                        continue
                    for nx2 in range(nx1 + 1, len(self.adj[x])):
                        b = self.adj[x][nx2]
                        if b == y or self.adjacent(y, b) or self.adjacent(a, b):
                            continue
                        self.eorbit[e][4] += 1
                        f_27 += tri[xy]
                        f_17 += len(self.adj[x]) - 3
                        f_15 += (len(self.adj[a]) - 1) + (len(self.adj[b]) - 1)

                # Edge-orbit 3 = (5,5)
                for nx in range(len(self.adj[x])):
                    a = self.adj[x][nx]
                    if a == y or self.adjacent(y, a):
                        continue
                    for ny in range(len(self.adj[y])):
                        b = self.adj[y][ny]
                        if b == x or self.adjacent(x, b) or self.adjacent(a, b):
                            continue
                        self.eorbit[e][3] += 1
                        f_20 += tri[xy]
                        f_16 += (len(self.adj[x]) - 2) + (len(self.adj[y]) - 2)
                        f_13 += (len(self.adj[a]) - 1) + (len(self.adj[b]) - 1)

                # Edge-orbit 2 = (4,5)
                for ny in range(len(self.adj[y])):
                    a = self.adj[y][ny]
                    if a == x or self.adjacent(x, a):
                        continue
                    for na in range(len(self.adj[a])):
                        b = self.adj[a][na]
                        ab = self.inc[a][na][1]
                        if b == y or self.adjacent(y, b) or self.adjacent(x, b):
                            continue
                        self.eorbit[e][2] += 1
                        f_29 += common_y[b] - 1
                        f_28 += common_x[b]
                        f_24 += tri[xy]
                        f_21 += tri[ab]
                        f_14 += len(self.adj[a]) - 2
                        f_12 += len(self.adj[b]) - 1
                for nx in range(len(self.adj[x])):
                    a = self.adj[x][nx]
                    if a == y or self.adjacent(y, a):
                        continue
                    for na in range(len(self.adj[a])):
                        b = self.adj[a][na]
                        ab = self.inc[a][na][1]
                        if b == x or self.adjacent(x, b) or self.adjacent(y, b):
                            continue
                        self.eorbit[e][2] += 1
                        f_29 += common_x[b] - 1
                        f_28 += common_y[b]
                        f_24 += tri[xy]
                        f_21 += tri[ab]
                        f_14 += len(self.adj[a]) - 2
                        f_12 += len(self.adj[b]) - 1

                # Solve system of equations (complete from C++ code)
                self.eorbit[e][67] = C5[e]
                self.eorbit[e][66] = (f_66 - 6 * self.eorbit[e][67]) // 2
                self.eorbit[e][65] = (f_65 - 6 * self.eorbit[e][67])
                self.eorbit[e][64] = (f_64 - 2 * self.eorbit[e][66])
                self.eorbit[e][63] = (f_63 - 2 * self.eorbit[e][65]) // 2
                self.eorbit[e][62] = (f_62 - 2 * self.eorbit[e][66] - 3 * self.eorbit[e][67])
                self.eorbit[e][61] = (f_61 - 2 * self.eorbit[e][65] - 4 * self.eorbit[e][66] - 12 * self.eorbit[e][67])
                self.eorbit[e][60] = (f_60 - 1 * self.eorbit[e][65] - 3 * self.eorbit[e][67])
                self.eorbit[e][59] = (f_59 - 2 * self.eorbit[e][65]) // 2
                self.eorbit[e][58] = (f_58 - 1 * self.eorbit[e][64] - 1 * self.eorbit[e][66])
                self.eorbit[e][57] = (f_57 - 2 * self.eorbit[e][63] - 2 * self.eorbit[e][64] - 2 * self.eorbit[e][65])
                self.eorbit[e][56] = (f_56 - 2 * self.eorbit[e][63]) // 2
                self.eorbit[e][55] = (f_55 - 4 * self.eorbit[e][62] - 2 * self.eorbit[e][64] - 4 * self.eorbit[e][66])
                self.eorbit[e][54] = (f_54 - 1 * self.eorbit[e][61] - 2 * self.eorbit[e][63] - 2 * self.eorbit[e][65]) // 2
                self.eorbit[e][53] = (f_53 - 2 * self.eorbit[e][59] - 2 * self.eorbit[e][64] - 2 * self.eorbit[e][65])
                self.eorbit[e][52] = (f_52 - 2 * self.eorbit[e][59] - 2 * self.eorbit[e][63] - 2 * self.eorbit[e][65])
                self.eorbit[e][51] = (f_51 - 1 * self.eorbit[e][61] - 2 * self.eorbit[e][62] - 1 * self.eorbit[e][65] - 4 * self.eorbit[e][66] - 6 * self.eorbit[e][67])
                self.eorbit[e][50] = (f_50 - 2 * self.eorbit[e][60] - 1 * self.eorbit[e][61] - 2 * self.eorbit[e][65] - 2 * self.eorbit[e][66] - 6 * self.eorbit[e][67])
                self.eorbit[e][49] = (f_49 - 1 * self.eorbit[e][59]) // 3
                self.eorbit[e][48] = (f_48 - 2 * self.eorbit[e][62] - 1 * self.eorbit[e][66]) // 3
                self.eorbit[e][47] = (f_47 - 2 * self.eorbit[e][59] - 1 * self.eorbit[e][61] - 2 * self.eorbit[e][65]) // 2
                self.eorbit[e][46] = (f_46 - 1 * self.eorbit[e][57] - 1 * self.eorbit[e][63])
                self.eorbit[e][45] = (f_45 - 1 * self.eorbit[e][52] - 4 * self.eorbit[e][58] - 4 * self.eorbit[e][60])
                self.eorbit[e][44] = (f_44 - 2 * self.eorbit[e][56] - 1 * self.eorbit[e][57] - 2 * self.eorbit[e][63])
                self.eorbit[e][43] = (f_43 - 2 * self.eorbit[e][56] - 1 * self.eorbit[e][63])
                self.eorbit[e][42] = (f_42 - 2 * self.eorbit[e][56] - 1 * self.eorbit[e][57] - 2 * self.eorbit[e][63]) // 2
                self.eorbit[e][41] = (f_41 - 1 * self.eorbit[e][55] - 2 * self.eorbit[e][58] - 2 * self.eorbit[e][62] - 2 * self.eorbit[e][64] - 2 * self.eorbit[e][66])
                self.eorbit[e][40] = (f_40 - 2 * self.eorbit[e][54] - 1 * self.eorbit[e][55] - 1 * self.eorbit[e][57] - 1 * self.eorbit[e][61] - 2 * self.eorbit[e][63] - 2 * self.eorbit[e][64] - 2 * self.eorbit[e][65])
                self.eorbit[e][39] = (f_39 - 1 * self.eorbit[e][52] - 1 * self.eorbit[e][53] - 1 * self.eorbit[e][57] - 2 * self.eorbit[e][59] - 2 * self.eorbit[e][63] - 2 * self.eorbit[e][64] - 2 * self.eorbit[e][65])
                self.eorbit[e][38] = (f_38 - 3 * self.eorbit[e][49] - 1 * self.eorbit[e][56] - 1 * self.eorbit[e][59])
                self.eorbit[e][37] = (f_37 - 1 * self.eorbit[e][53] - 1 * self.eorbit[e][59])
                self.eorbit[e][36] = (f_36 - 1 * self.eorbit[e][52] - 2 * self.eorbit[e][60]) // 2
                self.eorbit[e][35] = (f_35 - 6 * self.eorbit[e][48] - 1 * self.eorbit[e][55] - 4 * self.eorbit[e][62] - 1 * self.eorbit[e][64] - 2 * self.eorbit[e][66])
                self.eorbit[e][34] = (f_34 - 2 * self.eorbit[e][47] - 1 * self.eorbit[e][53] - 1 * self.eorbit[e][55] - 2 * self.eorbit[e][59] - 1 * self.eorbit[e][61] - 2 * self.eorbit[e][64] - 2 * self.eorbit[e][65])
                self.eorbit[e][33] = (f_33 - 2 * self.eorbit[e][47] - 1 * self.eorbit[e][52] - 2 * self.eorbit[e][54] - 2 * self.eorbit[e][59] - 1 * self.eorbit[e][61] - 2 * self.eorbit[e][63] - 2 * self.eorbit[e][65])
                self.eorbit[e][32] = (f_32 - 6 * self.eorbit[e][49] - 1 * self.eorbit[e][53] - 2 * self.eorbit[e][59]) // 2
                self.eorbit[e][31] = (f_31 - 2 * self.eorbit[e][42] - 1 * self.eorbit[e][44] - 2 * self.eorbit[e][46] - 2 * self.eorbit[e][56] - 2 * self.eorbit[e][57] - 2 * self.eorbit[e][63])
                self.eorbit[e][30] = (f_30 - 2 * self.eorbit[e][42] - 2 * self.eorbit[e][43] - 1 * self.eorbit[e][44] - 4 * self.eorbit[e][56] - 1 * self.eorbit[e][57] - 2 * self.eorbit[e][63])
                self.eorbit[e][29] = (f_29 - 2 * self.eorbit[e][38] - 1 * self.eorbit[e][45] - 1 * self.eorbit[e][52]) // 2
                self.eorbit[e][28] = (f_28 - 2 * self.eorbit[e][43] - 1 * self.eorbit[e][45] - 1 * self.eorbit[e][52]) // 2
                self.eorbit[e][27] = (f_27 - 1 * self.eorbit[e][34] - 1 * self.eorbit[e][47])
                self.eorbit[e][26] = (f_26 - 1 * self.eorbit[e][33] - 2 * self.eorbit[e][36] - 1 * self.eorbit[e][50] - 1 * self.eorbit[e][52] - 2 * self.eorbit[e][60]) // 2
                self.eorbit[e][25] = (f_25 - 2 * self.eorbit[e][32] - 1 * self.eorbit[e][37] - 3 * self.eorbit[e][49] - 1 * self.eorbit[e][53] - 1 * self.eorbit[e][59])
                self.eorbit[e][24] = (f_24 - 1 * self.eorbit[e][39] - 1 * self.eorbit[e][45] - 1 * self.eorbit[e][52])
                self.eorbit[e][23] = (f_23 - 2 * self.eorbit[e][36] - 1 * self.eorbit[e][45] - 1 * self.eorbit[e][52] - 2 * self.eorbit[e][58] - 2 * self.eorbit[e][60])
                self.eorbit[e][22] = (f_22 - 1 * self.eorbit[e][37] - 1 * self.eorbit[e][44] - 1 * self.eorbit[e][53] - 1 * self.eorbit[e][56] - 1 * self.eorbit[e][59])
                self.eorbit[e][21] = (f_21 - 2 * self.eorbit[e][38] - 2 * self.eorbit[e][43] - 1 * self.eorbit[e][52]) // 2
                self.eorbit[e][20] = (f_20 - 1 * self.eorbit[e][40] - 1 * self.eorbit[e][54])
                self.eorbit[e][19] = (f_19 - 1 * self.eorbit[e][33] - 2 * self.eorbit[e][41] - 1 * self.eorbit[e][45] - 2 * self.eorbit[e][50] - 1 * self.eorbit[e][52] - 4 * self.eorbit[e][58] - 4 * self.eorbit[e][60])
                self.eorbit[e][18] = (f_18 - 2 * self.eorbit[e][32] - 2 * self.eorbit[e][38] - 1 * self.eorbit[e][44] - 6 * self.eorbit[e][49] - 1 * self.eorbit[e][53] - 2 * self.eorbit[e][56] - 2 * self.eorbit[e][59])
                self.eorbit[e][17] = (f_17 - 2 * self.eorbit[e][25] - 1 * self.eorbit[e][27] - 1 * self.eorbit[e][32] - 1 * self.eorbit[e][34] - 1 * self.eorbit[e][47]) // 3
                self.eorbit[e][16] = (f_16 - 2 * self.eorbit[e][20] - 2 * self.eorbit[e][22] - 1 * self.eorbit[e][31] - 2 * self.eorbit[e][40] - 1 * self.eorbit[e][44] - 2 * self.eorbit[e][54]) // 2
                self.eorbit[e][15] = (f_15 - 2 * self.eorbit[e][25] - 2 * self.eorbit[e][29] - 1 * self.eorbit[e][31] - 2 * self.eorbit[e][32] - 1 * self.eorbit[e][34] - 2 * self.eorbit[e][42] - 2 * self.eorbit[e][47])
                self.eorbit[e][14] = (f_14 - 1 * self.eorbit[e][18] - 2 * self.eorbit[e][21] - 1 * self.eorbit[e][30] - 2 * self.eorbit[e][38] - 1 * self.eorbit[e][39] - 2 * self.eorbit[e][43] - 1 * self.eorbit[e][52]) // 2
                self.eorbit[e][13] = (f_13 - 2 * self.eorbit[e][22] - 2 * self.eorbit[e][28] - 1 * self.eorbit[e][31] - 1 * self.eorbit[e][40] - 2 * self.eorbit[e][44] - 2 * self.eorbit[e][54])
                self.eorbit[e][12] = (f_12 - 2 * self.eorbit[e][21] - 2 * self.eorbit[e][28] - 2 * self.eorbit[e][29] - 2 * self.eorbit[e][38] - 2 * self.eorbit[e][43] - 1 * self.eorbit[e][45] - 1 * self.eorbit[e][52])
                
        end_time = time.time()
        print(f"{end_time - start_time:.2f}")
        
        end_time_all = time.time()
        print(f"total: {end_time_all - start_time_all:.2f}")

    def write_results(self, filename: str, graphlet_size: int = 5):
        """Write node orbit results to file"""
        orbit_counts = {5: 73, 4: 15}
        no = orbit_counts.get(graphlet_size, 73)

        with open(filename, "w") as f:
            for i in range(self.n):
                line_parts = []
                for j in range(no):
                    line_parts.append(str(self.orbit[i][j]))
                f.write(" ".join(line_parts) + "\n")

    def write_edge_results(self, filename: str, graphlet_size: int = 5):
        """Write edge orbit results to file"""
        orbit_counts = {5: 68, 4: 12}
        no = orbit_counts.get(graphlet_size, 68)

        with open(filename, "w") as f:
            for i in range(self.m):
                line_parts = []
                for j in range(no):
                    line_parts.append(str(self.eorbit[i][j]))
                f.write(" ".join(line_parts) + "\n")




def run_orca(
    edge_list: np.ndarray,
    num_nodes: Optional[int] = None,
    mode: Literal["node", "edge"] = "node",
    graphlet_size: Literal[4, 5] = 4,
    debug: bool = False,
) -> np.ndarray:
    """
    Run ORCA orbit counting directly on edge list without file I/O.
    
    Args:
        edge_list: Nx2 numpy array of edges (each row is [node1, node2])
        num_nodes: Number of nodes in graph (if None, inferred from edge_list)
        mode: Either "node" for node orbits or "edge" for edge orbits
        graphlet_size: Size of graphlets to count (4 or 5)
        debug: If True, print debug information
        
    Returns:
        Numpy array of orbit counts. For node mode: shape (num_nodes, num_orbits)
        For edge mode: shape (num_edges, num_orbits)
    """
    # Handle output suppression
    original_stdout = sys.stdout
    suppress_output = False
    
    if not debug:
        try:
            import os
            null_device = os.devnull
            sys.stdout = open(null_device, 'w')
            suppress_output = True
        except Exception:
            # Fallback: just use original stdout
            suppress_output = False
    
    try:
        # Initialize ORCA
        orca = Orca()
        
        # Setup graph from edge list
        orca.setup_graph_from_edges(edge_list, num_nodes)
        
        # Setup adjacency structures
        orca.setup_adjacency()
        
        # Run orbit counting
        if mode == "node":
            if graphlet_size == 4:
                orca.count4()
                num_orbits = 15
            elif graphlet_size == 5:
                orca.count5()
                num_orbits = 73
            else:
                raise ValueError("graphlet_size must be 4 or 5")
            
            # Convert to numpy array
            result = np.zeros((orca.n, num_orbits), dtype=np.int64)
            for i in range(orca.n):
                for j in range(num_orbits):
                    result[i, j] = orca.orbit[i][j]
            
        elif mode == "edge":
            if graphlet_size == 4:
                orca.ecount4()
                num_orbits = 12
            elif graphlet_size == 5:
                orca.ecount5()
                num_orbits = 68
            else:
                raise ValueError("graphlet_size must be 4 or 5")
            
            # Convert to numpy array
            result = np.zeros((orca.m, num_orbits), dtype=np.int64)
            for i in range(orca.m):
                for j in range(num_orbits):
                    result[i, j] = orca.eorbit[i][j]
        else:
            raise ValueError("mode must be 'node' or 'edge'")
            
        return result
        
    finally:
        # Restore stdout
        if suppress_output:
            try:
                sys.stdout.close()
            except Exception:
                pass
        sys.stdout = original_stdout


def orca_nodes(
    edge_list: np.ndarray,
    num_nodes: Optional[int] = None,
    graphlet_size: Literal[4, 5] = 4,
    debug: bool = False,
) -> np.ndarray:
    """
    Convenience function for computing node orbits.
    
    Args:
        edge_list: Nx2 numpy array of edges
        num_nodes: Number of nodes (if None, inferred from edge_list)
        graphlet_size: Size of graphlets to count (4 or 5)
        debug: If True, print debug information
        
    Returns:
        Numpy array of shape (num_nodes, num_orbits) with node orbit counts
    """
    return run_orca(
        edge_list, num_nodes, mode="node", graphlet_size=graphlet_size, debug=debug
    )


def orca_edges(
    edge_list: np.ndarray,
    num_nodes: Optional[int] = None,
    graphlet_size: Literal[4, 5] = 4,
    debug: bool = False,
) -> np.ndarray:
    """
    Convenience function for computing edge orbits.
    
    Args:
        edge_list: Nx2 numpy array of edges
        num_nodes: Number of nodes (if None, inferred from edge_list)
        graphlet_size: Size of graphlets to count (4 or 5)
        debug: If True, print debug information
        
    Returns:
        Numpy array of shape (num_edges, num_orbits) with edge orbit counts
    """
    return run_orca(
        edge_list, num_nodes, mode="edge", graphlet_size=graphlet_size, debug=debug
    )
