function [dstate] = NonlinearConcensusODE(t,state,E)
% [na,ne] = size(E); %na = # of agents. ne = # of edges
% 
% x=state(1:na,1); %Each node has one state variable.

y = state; %Position measurement
zeta = E'*y; %Relative Positions.
mu = zeta.*exp(abs(zeta).^2); %Highly nonlinear controller.
u = -E*mu;
dstate = u;

end