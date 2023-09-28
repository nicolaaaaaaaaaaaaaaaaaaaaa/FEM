% File assignment1_ex1_input.m
%
% 
% No. of Nodes: 16  
% No. of Elements : 
clear all
clf

% Coordinates of 3 nodes,
X = [  0.00  0.00 
       0.25  0.00 
       0.00  0.25 
       0.25  0.25
       0.00  0.50
       0.25  0.50
       0.00  0.75
       0.25  0.75
       0.00  1.00
       0.25  1.00
       0.50  1.00
       0.75  1.00
       1.00  1.00
       0.50  0.75
       0.75  0.75
       1.00  0.75 ];
       
% Topology matrix IX(node1,node2,propno),
IX = [ 1   2   1 
       1   3   1
       2   3   1 
       2   4   1
       3   4   1
       3   5   1
       3   6   1
       4   6   1
       5   6   1 
       5   7   1
       6   7   1
       6   8   1
       7   8   1
       7   9   1
       9   10  1
       9   8   1
       7   10  1
       8   10  1
       10  11  1
       10  14  1
       8   14  1
       11  14  1
       11  12  1
       11  15  1
       14  12  1
       14  15  1
       12  15  1
       12  13  1
       12  16  1
       15  13  1
       15  16  1
       13  16  1 ];


% Element property matrix mprop = [ E A ],
mprop = [ 70*10^9 0.0002
          2.0 1.0 ];

% Prescribed loads mat(node,dof,force)
loads = [ 16   2  -15*10^3 ];

% Boundary conditions mat(node,dof,disp)   
bound = [ 1  1  0.0
          1  2  0.0
          2  2  0.0 ];

% Control Parameters
%plotdof = 6;
