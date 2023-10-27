%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Nonlinearity Newton Rhapson Method                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fea()

close all
clc

%--- Input file ----------------------------------------------------------%
ex3             % Input file

neqn = size(X,1)*size(X,2);         % Number of equations
ne = size(IX,1);                    % Number of elements
disp(['Number of DOF ' sprintf('%d',neqn) ...
    ' Number of elements ' sprintf('%d',ne)]);

%--- Initialize arrays ---------------------------------------------------%
Kmatr=zeros(neqn,neqn);                 % Stiffness matrix
P=zeros(neqn,1);                        % Force vector
Pfinal = zeros(neqn,1);                 % Final force vector
delP=zeros(neqn,1);                     % Force increment vector
D=zeros(neqn,1);                        % Displacement vector
R=zeros(neqn,1);                        % Residual vector
Rint=zeros(neqn,1);
strain=zeros(ne,1);                     % Element strain vector
stress=zeros(ne,1);                     % Element stress vector
force = zeros(ne,1); 
Darray = [];                             % Store D
Parray = [];                             % Store P  


%--- Calculate displacements ---------------------------------------------%
[delP,Pfinal]=buildload(delP,Pfinal, loads);                           % Build global load vector



for n = 1:loads(1,4)+1
    Parray = [Parray,P];
    P = P + delP;
    Darray = [Darray,D];
    for i = 0:maxn
        [R] = residual(X,IX,ne,mprop,P,D,R,Rint);       % Calculate residual
        [R]=enforceR(R,bound);                          % Enforce boundary conditions on R
        if norm(R) <= ep*norm(Pfinal)                        % Exit iteration loop if residual is small 
            break
        end
        [Kmatr]=buildstiff(X,IX,ne,mprop,Kmatr,D);         % Build tangent stiffness matrix
        [Kmatr,P]=enforce(Kmatr,P,bound);               % Enforce boundary conditions on Kmatr
        delD = -Kmatr\(R);
        D = D + delD;                                      % Solve system of equations
    end
end



%--- Plot results --------------------------------------------------------%                                                        
PlotStructure(X,IX,ne,neqn,bound,loads,D,stress)        % Plot structure

%Plotdisplacement(Parray,Darray)

% plot FEM solution

%plot(Darray(85,:), Parray(85,:))

%hold on

% hold off

P1 = Parray(88,21);
D1 = Darray(88,21);
L = 20;
EI = (P1*L^3)/(3*D1);
disp(EI)
% 
%EI = 0.4977;
% 
Pcrit = (pi^2*EI)/(4*L^2);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build global load vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delP,Pfinal]=buildload(delP,Pfinal,loads)
for i=1:size(loads,1)
    % extract information from load matrix
    node = loads(i,1);
    dir = loads(i, 2);
    del = loads(i, 3)/loads(i, 4);
    Pf = loads(i,3);
    
    % compute corrisponding dof
    dof = node.*2 - 2 + dir;        

    % insert value in matrix
    delP(dof) = del;
    Pfinal(dof) = Pf;
end
return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build global stiffness matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K]=buildstiff(X,IX,ne,mprop,K,D)

% This subroutine builds the global stiffness matrix from
% the local element stiffness matrices

K = K*0;

Bdmult = [1 0 -1 0;
          0 1 0 -1;
          -1 0 1 0;
          0 -1 0 1];

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

    % compute E, A, and Lo

    propno = IX(e, 3);
    A = mprop(propno,2);
    E = mprop(propno,1);
    Lo = sqrt(delx^2 + dely^2);

    ui = D(node1.*2 - 1);
    vi = D(node1.*2);
    uj = D(node2.*2 - 1);
    vj = D(node2.*2);
    d = [ui vi uj vj]';

    % compute linear displacement vector
    [Bo] = (1/Lo^2).*[-delx -dely delx dely]';
    [Bd] = (1/Lo^2).*Bdmult*d;

    
    % compute element stiffness matrix

    sig = E.*(Bo'+0.5.*Bd')*d;
    
    
    [ksige] = (1/Lo^2).*Bdmult.*A.*sig.*Lo;
    [k0e] = E.*A.*Lo.*Bo*Bo';
    [kde] = (A*E*Lo.*Bo*Bd') + (A*E*Lo.*Bd*Bo') + (A*E*Lo.*Bd*Bd');

    [ke] = ksige + k0e + kde;
    
    
    % local to global dof
    for i = 1:2
        edof(2.*i-1) = 2.*IX(e,i)-1;
        edof(2.*i) = 2.*IX(e,i);
    end

    % compute global stiffness matrix
    K(edof,edof) = K(edof,edof) + ke;

end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build residual vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R]=residual(X,IX,ne,mprop,P,D,R,Rint)

% This subroutine builds the residual matrix from
% the local element residual matrices

Rint = Rint*0;

Bdmult = [1 0 -1 0;
          0 1 0 -1;
          -1 0 1 0;
          0 -1 0 1];

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

    Lo = sqrt(delx^2 + dely^2);


    % compute E and A
    propno = IX(e, 3);
    A = mprop(propno,2);
    E = mprop(propno,1);

    ui = D(node1.*2 - 1);
    vi = D(node1.*2);
    uj = D(node2.*2 - 1);
    vj = D(node2.*2);
    d = [ui vi uj vj]';

   
    % compute linear displacement vector
    [Bo] = (1/Lo^2).*[-delx -dely delx dely]';
    [Bd] = (1/Lo^2).*Bdmult*d;
    [Bhat] = Bo + Bd;
    
    sig = E.*(Bo'+0.5.*Bd')*d;

    rint = Bhat*A*sig*Lo;

    for i = 1:2
       edof(2.*i-1) = 2.*IX(e,i)-1;
       edof(2.*i) = 2.*IX(e,i);
    end

    Rint(edof,1) = Rint(edof,1) + rint;

end
R = Rint - P;
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
%%% Enforce boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R]=enforceR(R,bound)

% This subroutine enforces the support boundary conditions on R

for i=1:size(bound,1)

    % compute corrisponding dof
    dof = bound(i,1).*2 -2 + bound(i,2);

    % enforce row = 0 in residual vector
    R(dof) = 0;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot displacement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Plotdisplacement(Parray,Darray)
%figure
hold on
for i = 1:length(Parray(1))
    plot(Darray(i), Parray(i))
end
hold off
return
