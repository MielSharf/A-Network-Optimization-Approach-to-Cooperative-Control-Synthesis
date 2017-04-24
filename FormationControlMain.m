function [] = FormationControlMain()

%% Graph
n = 4; %# of agents
Edges = {[1,2],[2,3],[3,4]}; %Connected edges
E = constructAdjacencyMatrix(n,Edges);

%% Agents - Physical Parameters
%%% Original Random Choice 
% omegas = 1+5*randn(n,1).^2; 
% bs = 0.2+0.2*randn(n,1).^2;
% Zerospots = 12*randi(n,1)-6;  %Centers of oscillators
% ws = omegas.^2.*Zerospots;  %Steady-States for zero control.

%%% Reconstruction of Paper Results
omegas = [15.54,5.13,7.89,4.29]';
bs = [1.66,1.22,4.62,1.23]';
Zerospots = [5,-2,1,0]'; %Centers of oscillators
ws = omegas.^2.*Zerospots;  %Steady-States for zero control.

%% Construction of Dynamical Systems
% Agents' dynamical systems are of the form:
% dx/dt = Ax+Bu+Lw
% y     = Cx
As = cell(n,1); 
Bs = cell(n,1);
Cs = cell(n,1);
Ls = cell(n,1);

% Linearly forced damped oscillator
for i=1:n
    As{i} = [0 , 1 ; -omegas(i)^2 , -bs(i)];
    Bs{i} = [0 ; 1];
    Cs{i} = [1 , 0];
    Ls{i} = [0 ; 1];
end

% Stacked Plant
A = blkdiagcell(As);
B = blkdiagcell(Bs);
C = blkdiagcell(Cs);
L = blkdiagcell(Ls);

%% Run the System for a First Time - Find zeta_0.
yout = [];
tout = [];
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
x0 = [0,0,0,0,0,0,0,0]';
eta0 = [0,0,0]';

y0 = [x0;eta0];

t0 = 0;
te = 1000;

zeta = [0;0;0]; %Unaugmented Plant
muzeta = [0;0;0]; %Unaugmented Plant
w = [zeta;muzeta];
[~, yout] = ode15s(@IntegratorFormationControlODE,[t0 te],y0,options, E,A,B,C,L,ws,w);
y0 = yout(end,[1,3,5,7])';
zeta_0 = E'*y0;


%% Specification of Wanted Formations and Computation of the Corresponding Support Vectors \mu_\zeta
NumberOfSwitches = 5; %Number of different formations considered
zetas = cell(NumberOfSwitches,1); 
% Different specified formation vectors.
zetas{1} = [0;0;0];
zetas{2} = [1;1;1];
zetas{3} = [2;2;2];
zetas{4} = [0;3;0];
zetas{5} = [3;0;-3];

%Compute different mu_zetas from plant and zeta.
muzeta_0 = FindSupportVectorsForLinearFormation(zeta_0,omegas,ws);
muzetas = cell(size(zetas));
for i=1:NumberOfSwitches
    muzetas{i} = FindSupportVectorsForLinearFormation(zetas{i},omegas,ws);
end

%% Run the System
yout = [];
tout = [];
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
%Initial Conditions - reconstruction
x0 = [6,-2,4,1,-2,6,0,1]';
eta0 = [0,0,0]';

y0 = [x0;eta0];

dt = 25; %Time between a different zeta is issued.
zeta0 = E'*Zerospots;
mu_zeta0 = zeros(n-1,1);

for i=1:NumberOfSwitches
te = i*dt;
t0 = (i-1)*dt;
zeta = zetas{i};
muzeta = muzetas{i};
w = [zeta-zeta_0;muzeta-muzeta_0];
[toutcurr, youtcurr] = ode15s(@IntegratorFormationControlODE,[t0 te],y0,options, E,A,B,C,L,ws,w);
tout = [tout;toutcurr];
yout = [yout;youtcurr];
y0 = yout(end,:);
end

figure;
plot(tout, yout(:,[1,3,5,7]), 'Linewidth',1);
grid;
legend('Agent1','Agent2','Agent3','Agent4');
% print -depsc UnifiedControlSchemePositionsFig

FormationVectors = zeros(length(tout),length(Edges));
for i=1:length(tout)
    FormationVector(i,:) = (E'*yout(i,[1,3,5,7])')';
end
figure;
plot(tout, FormationVector, 'Linewidth',1);
grid;
legend('\zeta_1','\zeta_2','\zeta_3');
% print -depsc UnifiedControlSchemeFormationFig

end

function E = constructAdjacencyMatrix(n,Edges)
m = length(Edges);
E = zeros(n,m);
for i=1:length(Edges)
    E(Edges{i}(1),i) = -1;
    E(Edges{i}(2),i) = 1;
end
end

function A = blkdiagcell(As)
l = length(As);
A = As{1};
for i=2:l
    A = blkdiag(A,As{i});
end
end

function mu = FindSupportVectorsForLinearFormation(zeta,omegas,ws)
% We need to solve the problem:
% min K^\star(y) s.t. E^Ty=zeta
% We know that K^\star(y) = 1/2 y^T W^2 y - y^T x where W=diag(omegas) and
% x = Zerospots. We note that ws = W^2*Zerospots, so we get:

% If the edges are i->i+1 in order, this y0 below satisfies E^Ty0 = zeta.
% This is the initial guess for y_zeta from the proof.
y0 = zeros(length(omegas),1); %Initial Guess vector.
y0(1) = 0;
for i=2:length(y0)
    y0(i) = y0(i-1)+zeta(i-1);
end

% The collection of all y's such that E^Ty=zeta is given by {y0+\beta*1}. 
% Find the optimal beta by deriving K\star(y). This can also be done using
% methods of convex optimization (like CVX).
OneVec = ones(size(y0));
W = diag(omegas);

%Analytic form of optimal beta.
beta = - (OneVec' * W^2 * y0) / (OneVec' * W^2 * OneVec) + (ws'*OneVec)/(OneVec' * W^2 * OneVec);
y_zeta = y0 + beta*OneVec;

% Now find the correpsonding mu. In this case, we are dealing with a tree,
% so the computation is easy (start from leaves and move inside). In
% general, one can use a spanning tree and have mu_k=0 at all edges not in
% the tree.
k_inv_y = W^2*y_zeta - ws; %k^{-1} of y_zeta.
mu = zeros(length(zeta),1);
mu(1) = k_inv_y(1);
mu(2) = mu(1) + k_inv_y(2);
mu(3) = mu(2) + k_inv_y(3);
end
