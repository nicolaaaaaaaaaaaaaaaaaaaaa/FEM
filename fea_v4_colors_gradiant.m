%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Basis truss program                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function fea()

close all
clc

%--- Input file ----------------------------------------------------------%
%Truss4              % Input file
%test1                   % Input file
try_02

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
force = zeros(ne,1); 

lenght = [0];


%--- Calculate displacements ---------------------------------------------%

[P]=buildload(P,loads);       % Build global load vector

[Kmatr]=buildstiff(X,IX,ne,mprop,Kmatr);    % Build global stiffness matrix

[Kmatr,P]=enforce(Kmatr,P,bound);           % Enforce boundary conditions

D = Kmatr\P;                                            % Solve system of equations

[strain,stress,force]=recover(mprop,X,IX,D,ne,strain,stress,force); % Calculate element 
                                                        % stress, strain
                                                        % and force

lenght = Totlenght(lenght,X,IX,ne);


disp('tot lenght =')
disp(lenght)



disp('stress = ')
disp(stress)

%--- Plot results --------------------------------------------------------%                                                        
PlotStructure(X,IX,ne,neqn,bound,loads,D,stress)        % Plot structure

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% tot lenght %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [l] = Totlenght(l,X,IX,ne)
for e=1:ne

    % find corresponding nodes
    node1 = IX(e,1);
    node2 = IX(e,2);

    % evaluate coordinates of the nodes
    node1x = X(node1,1);
    node1y = X(node1,2);
    node2x = X(node2,1);
    node2y = X(node2,2);

    % compute undeformed lenght
    delx = node2x - node1x;
    dely = node2y - node1y;

    l = l + sqrt(delx^2 + dely^2);
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build global load vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P]=buildload(P,loads)
for i=1:size(loads,1)
    % extract information from load matrix
    node = loads(i,1);
    dir = loads(i, 2);
    mag = loads(i, 3);
    
    % compute corrisponding dof
    dof = node.*2 - 2 + dir;        

    % insert value in matrix
    P(dof) = mag;
end
return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build global stiffness matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K]=buildstiff(X,IX,ne,mprop,K)

% This subroutine builds the global stiffness matrix from
% the local element stiffness matrices

for e=1:ne

    % find corresponding nodes
    node1 = IX(e,1);
    node2 = IX(e,2);

    % evaluate coordinates of the nodes
    node1x = X(node1,1);
    node1y = X(node1,2);
    node2x = X(node2,1);
    node2y = X(node2,2);

    % compute undeformed lenght
    delx = node2x - node1x;
    dely = node2y - node1y;

    % compute linear displacement vector
    [Bo] = (1/(delx^2 + dely^2)).*[-delx -dely delx dely]';

    % compute E and A
    propno = IX(e, 3);
    E = mprop(propno,1);
    A = mprop(propno,2);

    % compute element stiffness matrix
    [ke] = E.*A.*(sqrt(delx^2 + dely^2)).*Bo*Bo';
    
    % compute corrisponding dof
    dofi = node1.*2 - 1;
    dofj = node2.*2 - 1;
    edof = [dofi,dofi+1,dofj,dofj+1];

    % compute global stiffness matrix
    K(edof,edof) = K(edof,edof) + ke;

end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Enforce boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Kmatr,P]=enforce(Kmatr,P,bound)

% This subroutine enforces the support boundary conditions

for i=1:size(bound,1)

    % compute corrisponding dof
    dof = bound(i,1).*2 -2 + bound(i,2);

    % enforce row and column = 0 in stiffness matrix
    Kmatr(dof,:) = 0;
    Kmatr(:,dof) = 0;
    Kmatr(dof,dof) = 1;

    % enforce row = 0 in force vector
    P(dof) = 0;

end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate element strain and stress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [strain,stress,force]=recover(mprop,X,IX,D,ne,strain,stress, force)

% This subroutine recovers the element stress, element strain, 
% and nodal reaction forces
        
for e=1:ne

    % find corresponding nodes
    node1 = IX(e,1);
    node2 = IX(e,2);

    % evaluate coordinates of the nodes
    node1x = X(node1,1);
    node1y = X(node1,2);
    node2x = X(node2,1);
    node2y = X(node2,2);

    % compute undeformed lenght
    delx = node2x - node1x;
    dely = node2y - node1y;

    % compute linear displacement vector
    [Bo] = (1/(delx^2 + dely^2)).*[-delx -dely delx dely]';

    % compute displacement vector
    ui = D(node1.*2 - 1);
    vi = D(node1.*2);
    uj = D(node2.*2 - 1);
    vj = D(node2.*2);
    d = [ui vi uj vj]';

    % compute E and A
    propno = IX(e, 3);
    E = mprop(propno,1);
    A = mprop(propno,2);

    % compute stress matrix
    elstress = E.*Bo'*d;
    stress(e) = elstress;

    % compute strain matrix
    elstrain = Bo'*d;
    strain(e) = elstrain;

    % compute force matrix
    elforce = elstress.*A;
    force(e) = elforce;
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotStructure(X,IX,ne,neqn,bound,loads,D,stress)

% This subroutine plots the undeformed and deformed structure

h1=0;h2=0;
% Plotting Un-Deformed and Deformed Structure
% Color code:
    % green = zero stress
    % blue = tension
    % dark blue = max tension
    % red = compression
    % dark red = max compression

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
    col = [0 0 0];
    elstress = stress(e);
    if elstress == max(stress)
        col = [0, 0, 0.6];  % dark blue
    elseif elstress == min(stress)
        col = [0.6, 0, 0];  %  dark Red
    elseif elstress == 0
        col = [0, 1, 0];  %green
    elseif elstress < 0
        col = [1, 0, 0];  % Red
    elseif elstress > 0
        col = [0, 0, 1];  % blue
    end
    h2=plot(xx,yy, 'Color', col, 'LineWidth',3.5); 
end
plotsupports
plotloads

legend([h1 h2],{'Undeformed state',...
                'Deformed state'})

axis equal;
hold off

return