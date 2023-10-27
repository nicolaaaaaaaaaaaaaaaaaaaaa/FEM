%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Basis truss program                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fea()

close all
clc

%--- Input file ----------------------------------------------------------%
example1_modified               % Input file
%test1                   % Input file

neqn = size(X,1)*size(X,2);         % Number of equations
ne = size(IX,1);                    % Number of elements
disp(['Number of DOF ' sprintf('%d',neqn) ...
    ' Number of elements ' sprintf('%d',ne)]);

%--- Initialize arrays ---------------------------------------------------%
Kmatr=zeros(neqn,neqn);                 % Stiffness matrix
P=zeros(neqn,1);                        % Force vector
D=zeros(neqn,1);                        % Displacement vector
R=zeros(neqn,1);                        % Residual vector
strain=zeros(ne,1);                     % Element strain vector
stress=zeros(ne,1);                     % Element stress vector

%--- Calculate displacements ---------------------------------------------%
[P]=buildload(X,IX,ne,P,loads,mprop);       % Build global load vector

[Kmatr]=buildstiff(X,IX,ne,mprop,Kmatr);    % Build global stiffness matrix

[Kmatr,P]=enforce(Kmatr,P,bound);           % Enforce boundary conditions

D = Kmatr\P;                                            % Solve system of equations

[strain,stress]=recover(mprop,X,IX,D,ne,strain,stress); % Calculate element 
                                                        % stress and strain
                                                        
%--- Plot results --------------------------------------------------------%                                                        
PlotStructure(X,IX,ne,neqn,bound,loads,D,stress)        % Plot structure

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build global load vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P]=buildload(X,IX,ne,P,loads,mprop);
for i=1:size(loads,1)
    node = loads(i,1);
    dir = loads(i, 2);
    mag = loads(i, 3);
    dof = node.*2 - 2 + dir;
    P(dof) = mag;
    end
return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build global stiffness matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K]=buildstiff(X,IX,ne,mprop,K)

% This subroutine builds the global stiffness matrix from
% the local element stiffness matrices

for e=1:ne
    node1 = IX(e,1);
    node2 = IX(e,2);
    node1x = X(node1,1);
    node1y = X(node1,2);
    node2x = X(node2,1);
    node2y = X(node2,2);
    delx = node2x - node1x;
    dely = node2y - node1y;
    [Bo] = (1/(delx^2 + dely^2)).*[-delx -dely delx dely]';
    propno = IX(e, 3);
    [ke] = mprop(propno,1).*mprop(propno,2).*(sqrt(delx^2 + dely^2)).*[Bo]*[Bo]';
    i = node1.*2 - 1;
    j = node2.*2 - 1;
    edof = [i,i+1,j,j+1];
    K(edof,edof) = K(edof,edof) + ke;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Enforce boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Kmatr,P]=enforce(Kmatr,P,bound)

% This subroutine enforces the support boundary conditions

for i=1:size(bound,1)
    dof = bound(i,1).*2 -2 + bound(i,2);
    Kmatr(dof,:) = 0;
    Kmatr(:,dof) = 0;
    Kmatr(dof,dof) = 1;
    P(dof) = 0;
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate element strain and stress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [strain,stress]=recover(mprop,X,IX,D,ne,strain,stress);

% This subroutine recovers the element stress, element strain, 
% and nodal reaction forces
        
for e=1:ne
    node1 = IX(e,1);
    node2 = IX(e,2);
    node1x = X(node1,1);
    node1y = X(node1,2);
    node2x = X(node2,1);
    node2y = X(node2,2);
    delx = node2x - node1x;
    dely = node2y - node1y;
    [Bo] = (1/(delx^2 + dely^2)).*[-delx -dely delx dely]';
    ui = D(node1.*2 - 1);
    vi = D(node1.*2);
    uj = D(node2.*2 - 1);
    vj = D(node2.*2);
    d = [ui vi uj vj]';
    propno = IX(e, 3);
    elstress = mprop(propno,1).*Bo'*d;
    stress(e) = elstress;
    elstrain = Bo'*d;
    strain(e) = elstrain;
    elforce = elstress.*mprop(propno,2);
    force(e) = elforce;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotStructure(X,IX,ne,neqn,bound,loads,D,stress)

% This subroutine plots the undeformed and deformed structure

h1=0;h2=0;
% Plotting Un-Deformed and Deformed Structure
clf
hold on
box on
for e = 1:ne
    xx = X(IX(e,1:2),1);
    yy = X(IX(e,1:2),2);
    h1=plot(xx,yy,'k:','LineWidth',1.);
    edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)];
    xx = xx + D(edof(1:2:4));
    yy = yy + D(edof(2:2:4));
    
    h2=plot(xx,yy,'b','LineWidth',3.5);    
end
plotsupports
plotloads

legend([h1 h2],{'Undeformed state',...
                'Deformed state'})

axis equal;
hold off

return
