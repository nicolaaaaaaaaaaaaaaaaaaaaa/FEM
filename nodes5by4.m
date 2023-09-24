% Created with:       FlExtract v1.13
% Element type:       truss
% Number of nodes:    19
% Number of elements: 50

clear all;

% Node coordinates: x, y
X = [
0	0
0	0.333333
0	0.666667
0	1
0.5	0
0.5	0.333333
0.5	0.666667
0.5	1
1	0.333333
1	0.666667
1	1
1.5	0
1.5	0.333333
1.5	0.666667
1.5	1
2	0
2	0.333333
2	0.666667
2	1
];
% Element connectivity: node1_id, node2_id, material_id
IX = [
2	1	1
6	1	1
5	1	1
3	2	1
7	2	1
6	2	1
5	2	1
4	3	1
8	3	1
7	3	1
6	3	1
8	4	1
7	4	1
6	5	1
9	5	1
7	6	1
10	6	1
9	6	1
8	7	1
11	7	1
10	7	1
9	7	1
11	8	1
10	8	1
10	9	1
14	9	1
13	9	1
12	9	1
11	10	1
15	10	1
14	10	1
13	10	1
15	11	1
14	11	1
13	12	1
17	12	1
16	12	1
14	13	1
18	13	1
17	13	1
16	13	1
15	14	1
19	14	1
18	14	1
17	14	1
19	15	1
18	15	1
17	16	1
18	17	1
19	18	1
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
16	2	0
];
% Nodal loads: node_id, degree of freedom (1 - x, 2 - y), load
loads = [
11	2	-0.01
];
% Control parameters
plotdof = 3;
