%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Exercise 4 - Topology optimization                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fea()

close all
clc

tic

%--- Input files ----------------------------------------------------------%
nodes9by7_all

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
v = zeros(ne,1);                        % Element volume vector
rho = zeros(ne,1);                      % Element density vector
rho_new = zeros(ne,1);                  % New element density vector
fp=zeros(ne,1);                         % Gradient f
gp=zeros(ne,1);                         % Gradient g
C = [];                                 % Compliance
i = [];                                 % Iteration


V_con = 8;
rho_min = 10^-6;
max_iopt = 100;
p = 1;
eps = 10^-8;
eta = 0.99;

%--- Calculate displacements ---------------------------------------------%

[P]=buildload(P,loads);                     % global load vector

[v] = buildvolume(X, IX, mprop, v, V_con, ne);  % volume vector
rho = ones(ne,1)*V_con./sum(v);

for iopt = 1:max_iopt
    [Kmatr]=buildstiff(X,IX,ne,mprop,Kmatr,rho,p);  % global stiffness matrix
    [Kmatr,P]=enforce(Kmatr,P,bound);         % Enforce boundary conditions
    D = Kmatr\P;                              % Solve system of equations
    [strain,stress,force,fp,gp]=recover(mprop,...
        X,IX,D,ne,strain,stress,force,v,rho,p,fp,gp); % Find gradients
    [lambmid,rho_new] = Bisectfunction(eps,ne,rho,fp,gp,eta,v,V_con,rho_min,rho_new);  % Bisection method
    if norm(rho-rho_new) < eps*norm(rho_new)
        break
    end

    f = D'*Kmatr*D;                          % compute compliance 
    rho = rho_new;
    C = [C,f];
end

i = 1:length(C);

figure
plot(i,C)

%--- Plot results --------------------------------------------------------%                                                        
figure

PlotStructure(X,IX,ne,neqn,bound,loads,D,stress,rho,rho_min)        % Plot structure

toc

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
%%% Build volume  vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v] = buildvolume(X, IX, mprop, v, ~, ne)

% This subroutine builds the volume vector

for e=1:ne

    % compute undeformed length
    delx = X(IX(e,2),1) - X(IX(e,1),1);
    dely = X(IX(e,2),2) - X(IX(e,1),2);
    L0 = sqrt(delx^2 + dely^2);
    
    propno = IX(e, 3);
    A = mprop(propno,2);

    % calculate volume

    ve = L0*A;
    v(e) = ve;

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build global stiffness matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This subroutine builds the global stiffness matrix from
% the local element stiffness matrices

function [K]=buildstiff(X,IX,ne,mprop,K,rho,p)

K = K*0;

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
    [ke] = E.*A.*L0.*Bo*Bo'*(rho(e))^p;
    
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

function [strain,stress,force,fp,gp]=recover(mprop,X,IX,D,ne,strain,stress,force,v,rho,p,fp,gp)

% This subroutine recovers the element stress, element strain, 
% and nodal reaction forces
        
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
    [ke0] = E.*A.*L0.*Bo*Bo';

  
    % compute displacement vector
    ui = D(IX(e,1).*2 - 1);
    vi = D(IX(e,1).*2);
    uj = D(IX(e,2).*2 - 1);
    vj = D(IX(e,2).*2);
    d = [ui vi uj vj]';

    % find gradients

    fpe = -p*rho(e)^(p-1).*d'*ke0*d;

    fp(e) = fpe;
    gp(e) = v(e);
    
    % compute strain matrix
    elstrain = Bo'*d;
    strain(e) = elstrain;

    % compute stress matrix
    elstress = (rho(e))^p*E*elstrain;
    stress(e) = elstress;
    
    % compute force matrix, minus sign to print the reaction forces
    elforce = elstress.*A;
    force(e) = elforce;

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bisect function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lambmid,rho_new] = Bisectfunction(eps,ne,rho,fp,gp,eta,v,V_con,rho_min,rho_new)

% This subroutine that applies the bi-section method

lamb1 = 10^-10;
lamb2 = 10^10;

while (lamb2-lamb1)/(lamb1+lamb2) > eps
    lambmid = (lamb1+lamb2)/2;
    for e=1:ne
        rho_newe = rho(e)*(-fp(e)/(lambmid*gp(e)))^eta;
        if rho_newe <= rho_min
            rho_newe = rho_min;
        elseif rho_newe >= 1
            rho_newe = 1;
        else
            rho_newe = rho_newe;
        end
        rho_new(e) = rho_newe;
    end
    g = rho_new' * v - V_con;
    if g > 0
        lamb1 = lambmid;
    else 
        lamb2 = lambmid;
    end

end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotStructure(X,IX,ne,neqn,bound,loads,D,stress,rho, rho_min)

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
    %h1=plot(xx,yy,'k:','LineWidth',1.);
    edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)];
    xx = xx + D(edof(1:2:4));
    yy = yy + D(edof(2:2:4));

    % color
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

    % thickness
    tick = rho(e)*7;
    if tick > 10*rho_min
        h2=plot(xx,yy, 'Color', [0, 0, 1], 'LineWidth',tick); 
    end
    
end
plotsupports
plotloads

legend(h2,'Deformed state', 'location','southeast','FontSize',15)

axis equal;

hold off

return

