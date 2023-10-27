%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Material Nonlinearity                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fea()

close all
clc

%--- Input file ----------------------------------------------------------%
TrussExercise2_2023                     % Input file

neqn = size(X,1)*size(X,2);             % Number of equations
ne = size(IX,1);                        % Number of elements
disp(['Number of DOF ' sprintf('%d',neqn) ...
    ' Number of elements ' sprintf('%d',ne)]);

%--- All methods ---------------------------------------------------------%

[PE, DE] = Euler(X, IX, ne, neqn, mprop,loads, bound); % Runs Euler method
[PNR,DNR] = NewtonRhapson(X, IX, ne, neqn, i_max, eSTOP,...
    mprop,loads, bound);                      % Runs Newton Rhapson method
[PNRm,DNRm] = NewtonRhapsonModified(X, IX, ne, neqn, i_max, eSTOP,...
    mprop,loads, bound);                      % Runs modified NR method


plot(DE(48,:),PE(48,:))
hold on
plot(DNR(48,:),PNR(48,:))
plot(DNRm(48,:),PNRm(48,:),"--")
hold off

title('Force vs Displacement for Different FEA Methods','FontSize',14)
xlabel('Vertical Displacement','FontSize',12)
ylabel('Load','FontSize',12)
legend({'Euler','Newton Rhapson', 'Modified Newton Rhapson'},...
    'Location','northwest','FontSize',12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Euler Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PE, DE] = Euler(X, IX, ne, neqn, mprop,loads, bound)

Kmatr=zeros(neqn,neqn);                 % Stiffness matrix
P=zeros(neqn,1);                        % Force vector
delP=zeros(neqn,1);                     % Force increment vector
D=zeros(neqn,1);                        % Displacement vector
R=zeros(neqn,1);                        % Residual vector
Rint=zeros(neqn,1);
Darray = [];                            % Store D
Parray = [];                            % Store P 


%--- Main Euler code loop ------------------------------------------------%

[delP]=buildloadeuler(delP,loads);      % Build delta P vector

for n= 1:loads(4)+1                     
    Parray = [Parray,P];                % Save loads for plotting
    P = P + delP;                       % Update new load P
    [Kmatr]=buildstiff(X,IX,ne,mprop,Kmatr,D); % Build stiffness matrix
    delP1 = delP - R;                   % Subtract residual from delta P
    [Kmatr,delP1]=enforce(Kmatr,delP1,bound);  % Enforce boundary
    delD = Kmatr\(delP1);               % Solve equation for delta D
    Darray = [Darray,D];                % Save displacement for plotting
    D = D + delD;                       % Update dsplacement
    [R] = residual(X,IX,ne,mprop,P,D,R,Rint);  % Recalculate residual
end

[PE] = Parray;
[DE] = Darray;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Newton Rhapson Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PNR,DNR] = NewtonRhapson(X, IX, ne, neqn, i_max, eSTOP,...
    mprop,loads, bound)

Kmatr=zeros(neqn,neqn);                 % Stiffness matrix
P=zeros(neqn,1);                        % Force vector
Pfinal = zeros(neqn,1);                 % Final force vector
delP=zeros(neqn,1);                     % Force increment vector
D=zeros(neqn,1);                        % Displacement vector
R=zeros(neqn,1);                        % Residual vector
Rint=zeros(neqn,1);
Darray = [];                            % Store D
Parray = [];                            % Store P 

%--- Main NR code loop ---------------------------------------------------%

[delP,Pfinal]=buildloadNR(delP,Pfinal, loads);  % Build global load vectors

for n= 1:loads(4)+1
    Parray = [Parray,P];                   % Save loads for plotting
    P = P + delP;                          % Update new load P
    Darray = [Darray,D];                   % Save displacement for plotting
    for i = 0:i_max                         % Iterate on residual
        [R] = residual(X,IX,ne,mprop,P,D,R,Rint);  % Find R
        [R]=enforceR(R,bound);           % Enforce boundary conditions on R
        if norm(R) <= eSTOP*norm(Pfinal)       % Exit loop if R is small 
            break
        end
        [Kmatr]=buildstiff(X,IX,ne,mprop,Kmatr,D); % Build stiffness matrix
        [Kmatr,P]=enforce(Kmatr,P,bound);   % Enforce boundary on Kmatr
        [LM,UM] = lu(Kmatr);                % Factorize Kmatr
        delD = UM \ (LM\R);                 % Solve for delta D
        D = D - delD;                       % Update displacement
    end
end

[PNR] = Parray;
[DNR] = Darray;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Modified Newton Rhapson Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PNRm,DNRm] = NewtonRhapsonModified(X, IX, ne, neqn, i_max, ...
    eSTOP, mprop,loads, bound)  

Kmatr=zeros(neqn,neqn);                 % Stiffness matrix
P=zeros(neqn,1);                        % Force vector
Pfinal = zeros(neqn,1);                 % Final force vector
delP=zeros(neqn,1);                     % Force increment vector
D=zeros(neqn,1);                        % Displacement vector
R=zeros(neqn,1);                        % Residual vector
Rint=zeros(neqn,1);
Darray = [];                            % Store D
Parray = [];                            % Store P 

%--- Main Modified NR code loop ------------------------------------------%

[delP,Pfinal]=buildloadNR(delP,Pfinal, loads);  % Build global load vectors

for n= 1:loads(4)+1
    Parray = [Parray,P];                % Save loads for plotting
    P = P + delP;                       % Update new load P
    Darray = [Darray,D];                % Save displacement for plotting
    [Kmatr]=buildstiff(X,IX,ne,mprop,Kmatr,D);     % Build stiffness matrix
    [Kmatr,P]=enforce(Kmatr,P,bound);   % Enforce boundary on Kmatr
    [LM,UM] = lu(Kmatr);            % Factorize Kmatr
    for i = 0:i_max                      % Iterate on residual
        [R] = residual(X,IX,ne,mprop,P,D,R,Rint);      % Find R
        [R]=enforceR(R,bound);          % Enforce boundary conditions on R
        if norm(R) <= eSTOP*norm(Pfinal)   % Exit loop if R is small 
            break
        end
        delD = UM \ (LM\R);             % Solve for delta D
        D = D - delD;                   % Update displacement                                   
    end
end

[PNRm] = Parray;
[DNRm] = Darray;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build global load vector Euler %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delP]=buildloadeuler(delP,loads)

%This subroutine builds a vector containing load increments for the Euler
%method

for i=1:size(loads,1)
    % Extract information from load matrix
    node = loads(i,1);
    dir = loads(i, 2);
    loadstep = loads(i, 3)/loads(i, 4);
    
    % Compute corrisponding dof
    dof = node.*2 - 2 + dir;        

    % Insert value in matrix
    delP(dof) = loadstep;
end
return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build global stiffness matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K]=buildstiff(X,IX,ne,mprop,K,D)

% This subroutine builds the global stiffness matrix from the local 
% element stiffness matrices

%Zero K matrix
K = K*0;                                    

for e=1:ne

    % Find element nodes
    node1 = IX(e,1);
    node2 = IX(e,2);

    % Evaluate coordinates of element nodes
    node1x = X(node1,1);
    node1y = X(node1,2);
    node2x = X(node2,1);
    node2y = X(node2,2);

    % Compute undeformed length
    delx = node2x - node1x;
    dely = node2y - node1y;
    Lo = sqrt(delx^2 + dely^2);

    % Compute linear strain displacement vector
    [Bo] = (1/Lo^2).*[-delx -dely delx dely]';

    % Extract material properties
    propno = IX(e, 3);
    A = mprop(propno,5);
    c1 = mprop(IX(e,3),1);
    c2 = mprop(IX(e,3),2);
    c3 = mprop(IX(e,3),3);
    c4 = mprop(IX(e,3),4);

    % Build displacement vector
    ui = D(node1.*2 - 1);
    vi = D(node1.*2);
    uj = D(node2.*2 - 1);
    vj = D(node2.*2);
    d = [ui vi uj vj]';
    
    % Compute element strain, lambda, and E
    elstrain = Bo'*d;
    lambda = 1 + c4*elstrain;
    Et = c4*(c1*(1+2*lambda^-3)+3*c2*lambda^-4+3*c3* ...
        (-1+lambda^2-2*lambda^-3+2*lambda^-4));

    % Compute element stiffness matrix
    [ke] = Et.*A.*Lo.*Bo*Bo';
    
    % Convert local dof to global dof
    for i = 1:2
        edof(2.*i-1) = 2.*IX(e,i)-1;
        edof(2.*i) = 2.*IX(e,i);
    end

    % Assemble global stiffness matrix
    K(edof,edof) = K(edof,edof) + ke;

end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Enforce boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Kmatr,P]=enforce(Kmatr,P,bound)

% This subroutine enforces the support boundary conditions

for i=1:size(bound,1)

    % Compute corrisponding dof
    dof = bound(i,1).*2 -2 + bound(i,2);

    % Enforce row and column = 0 in stiffness matrix
    Kmatr(dof,:) = 0;
    Kmatr(:,dof) = 0;
    Kmatr(dof,dof) = 1;

    % Enforce row = 0 in force vector
    P(dof) = 0;

end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build residual vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R]=residual(X,IX,ne,mprop,P,D,R,Rint)

% This subroutine builds the residual matrix from the local element 
% residual matrices

%Zero Rint
Rint = Rint*0;                       

for e=1:ne

    % Find corresponding nodes
    node1 = IX(e,1);
    node2 = IX(e,2);

    % Evaluate coordinates of the nodes
    node1x = X(node1,1);
    node1y = X(node1,2);
    node2x = X(node2,1);
    node2y = X(node2,2);

    % Compute undeformed length
    delx = node2x - node1x;
    dely = node2y - node1y;
    Lo = sqrt(delx^2 + dely^2);

    % Compute linear strain displacement vector
    [Bo] = (1/Lo^2).*[-delx -dely delx dely]';

    % Extract material properties
    propno = IX(e, 3);
    A = mprop(propno,5);
    c1 = mprop(IX(e,3),1);
    c2 = mprop(IX(e,3),2);
    c3 = mprop(IX(e,3),3);
    c4 = mprop(IX(e,3),4);

    % Build displacement vector
    ui = D(node1.*2 - 1);
    vi = D(node1.*2);
    uj = D(node2.*2 - 1);
    vj = D(node2.*2);
    d = [ui vi uj vj]';

    % Compute element strain, lambda, and sigma
    elstrain = Bo'*d;
    lam = 1 + c4*elstrain;
    sig = 1.*(lam-lam.^-2)+50.*(1-lam.^-3)+ 0.1.*...
    (1-3.*lam+lam.^3-2.*lam.^-3+3.*lam.^-2);

    % Compute element rint
    rint = Bo*A*sig*Lo;

    for i = 1:2
       edof(2.*i-1) = 2.*IX(e,i)-1;
       edof(2.*i) = 2.*IX(e,i);
    end

    % Build global Rint
    Rint(edof,1) = Rint(edof,1) + rint;

end

%Subtract P from Rint
R = Rint - P;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build global load vector NR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delP,Pfinal]=buildloadNR(delP,Pfinal, loads)

%This subroutine builds a vector containing load increments and a vector 
% containing max load for the NR method

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
%%% Enforce boundary conditions on R %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R]=enforceR(R,bound)

% This subroutine enforces the support boundary conditions on R

for i=1:size(bound,1)

    % compute corrisponding dof
    dof = bound(i,1).*2 -2 + bound(i,2);

    % enforce row = 0 in residual vector
    R(dof) = 0;
end
