% Created with:       FlExtract v1.13
% Element type:       truss
% Number of nodes:    11
% Number of elements: 24

clear all;

% Node coordinates: x, y
X = [
0	0
0	0.333333
0	0.666667
0	1
1	0.333333
1	0.666667
1	1
2	0
2	0.333333
2	0.666667
2	1
];
% Element connectivity: node1_id, node2_id, material_id
IX = [
2	1	1
5	1	1
3	2	1
6	2	1
5	2	1
4	3	1
7	3	1
6	3	1
5	3	1
7	4	1
6	4	1
6	5	1
10	5	1
9	5	1
8	5	1
7	6	1
11	6	1
10	6	1
9	6	1
11	7	1
10	7	1
9	8	1
10	9	1
11	10	1
];
% Element properties: Young's modulus, area
mprop = [
1	1
2	2
];
% Nodal diplacements: node_id, degree of freedom (1 - x, 2 - y), displacement
bound = [
1	1	0
1	2	0
8	2	0
];
% Nodal loads: node_id, degree of freedom (1 - x, 2 - y), load
loads = [
7	2	-0.01
];
% Control parameters
plotdof = 2;
