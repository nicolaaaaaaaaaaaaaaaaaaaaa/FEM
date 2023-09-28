%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Exercise 1 - Linear truss analysis                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fea()

close all
clc

%--- Input file ----------------------------------------------------------%
assignment1_ex1_input_v2             % Input file

neqn = size(X,1)*size(X,2);         % Number of equations
ne = size(IX,1);                    % Number of elements

%--- Initialize arrays ---------------------------------------------------%
Kmatr=zeros(neqn,neqn);                 % Stiffness matrix
P=zeros(neqn,1);                        % Force vector
D=zeros(neqn,1);                        % Displacement vector
R=zeros(neqn,1);                        % Residual vector
strain=zeros(ne,1);                     % Element strain vector
stress=zeros(ne,1);                     % Element stress vector
force = zeros(ne,1);                    % Element force vector


%--- Calculate displacements ---------------------------------------------%

[P]=buildload(P,loads);                     % global load vector

[Kmatr]=buildstiff(X,IX,ne,mprop,Kmatr);    % global stiffness matrix

[Kmatr,P]=enforce(Kmatr,P,bound);           % Enforce boundary conditions

D = Kmatr\P;                                % Solve system of equations

[strain,stress,force]=recover(mprop,X,IX,D,ne,strain,stress,force); % Calculate element 
                                                        % stress, strain
                                                        % and force
disp('reaction forces')
disp(force(1))                              % F1x
disp(force(2))                              % F1y
disp(force(4))                              % F2y

disp('force and stress in element A - B')
disp(force(21))                           
disp(stress(21))

disp('horizontal and vertical displacements in node C')
disp(D(31))
disp(D(32))

disp('max stress')
disp(max(stress))

%--- Plot results --------------------------------------------------------%                                                        

PlotStructure(X,IX,ne,neqn,bound,loads,D,stress)        % Plot structure

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build global load vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P]=buildload(P,loads)

% This subroutine builds the global load vector

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

% This subroutine builds the global stiffness matrix from
% the local element stiffness matrices

function [K]=buildstiff(X,IX,ne,mprop,K)

for e=1:ne

    % compute undeformed lenght
    delx = X(IX(e,2),1) - X(IX(e,1),1);
    dely = X(IX(e,2),2) - X(IX(e,1),2);
    L0 = sqrt(delx^2 + dely^2);

    % compute linear displacement vector
    [Bo] = (1/L0^2).*[-delx -dely delx dely]';

    % compute E and A
    propno = IX(e, 3);
    E = mprop(propno,1);
    A = mprop(propno,2);

    % compute element stiffness matrix
    [ke] = E.*A.*L0.*Bo*Bo';
    
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

    % compute undeformed lenght
    delx = X(IX(e,2),1) - X(IX(e,1),1);
    dely = X(IX(e,2),2) - X(IX(e,1),2);
    L0 = sqrt(delx^2 + dely^2);

    % compute linear displacement vector
    [Bo] = (1/L0^2).*[-delx -dely delx dely]';

    % compute displacement vector
    ui = D(IX(e,1).*2 - 1);
    vi = D(IX(e,1).*2);
    uj = D(IX(e,2).*2 - 1);
    vj = D(IX(e,2).*2);
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

    % compute force matrix, minus sign to print the reaction forces
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
    elseif abs(elstress) < max(stress)*10^-5  
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

legend([h1 h2],{'Undeformed state','Deformed state'}, 'location','southeast','FontSize',15)

axis equal;

hold off

return
