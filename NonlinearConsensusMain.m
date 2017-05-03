function [] = FormationControlMain()

%% Graph
n = 20; %# of agents
p = 0.1; %Probability of a random edge to appear. We use a very low probability to get a random sparse graph.
IsConnected = 0;
AttemptNumber = 0;
while(IsConnected == 0) % We need a connected graph, so we take a random graph until getting one.
    AttemptNumber = AttemptNumber+1;
    disp(sprintf('Attemp #%d',AttemptNumber));
    Edges = GetRandomEdges(n,p);
    [IsConnected,SecondLowestEigen] = PlotGraph(n,Edges);
end
E = constructAdjacencyMatrix(n,Edges);
SecondLowestEigen

%% Run the System
yout = [];
tout = [];
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
%Initial Conditions - randomized
x0 = 8*randn(n,1);
y0 = [x0];

te = 30;
t0 = 0;
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
[tout, yout] = ode23(@NonlinearConsensusODE,[t0 te],y0,options,E);

figure;
plot(tout, yout(:,:), 'Linewidth',1);
grid;
% print -depsc NonlinearControl
end

function [Edges] = GetRandomEdges(n,p)
Edges = {};
for i=1:n
    for j=(i+1):n
        if(rand<p)
            Edges{end+1} = [i,j];
        end
    end
end
end

function [IsConnected,SecondLowestEigen] = PlotGraph(n,Edges)
A = zeros(n,n);
for i=1:length(Edges)
    A(Edges{i}(1),Edges{i}(2))=1;
    A(Edges{i}(2),Edges{i}(1))=1;
end
G = graph(A);
Eig = eig(G.laplacian);
SecondLowestEigen = min(Eig(Eig > 0));
if max(G.conncomp) == 1
    IsConnected = 1;
else
    IsConnected = 0;
end
plot(G);
end

function E = constructAdjacencyMatrix(n,Edges)
m = length(Edges);
E = zeros(n,m);
for i=1:length(Edges)
    E(Edges{i}(1),i) = -1;
    E(Edges{i}(2),i) = 1;
end
end