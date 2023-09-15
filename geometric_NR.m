% File rubber.m, DAY1, modified 3/9 by OS
%
% Example: 2-bar truss
% No. of Nodes: 3  
% No. of Elements : 2
clear all

% Coordinates of 6 nodes,
X = [  0.00  0.40 
       1.50  0.00
       3.00  0.40 ];

% Topology matrix IX(node1,node2,propno),
IX = [ 1  2  1 
       2  3  1 ];
      
% Element property matrix mprop = [ E A],
mprop = [ 1 2 ];

% Prescribed loads mat(node,ldof,fmax,n)
loads = [ 2  2 0.03 200];

% Boundary conditions mat(node,ldof,disp)   
bound = [ 1  1  0.0
          1  2  0.0
          3  1  0.0 
          3  2  0.0];

%Maximum number of iterations

maxn = 100;

%Epsilon for comparing residual
ep = 10^-8;

% Control Parameters
plotdof = 6;

