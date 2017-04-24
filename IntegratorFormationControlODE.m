function [dstate] = IntegratorFormationControlODE(t,state, E, A,B,C,L,ws,omega)
[na,ne] = size(E); %na = # of agents. ne = # of edges

x=state(1:2*na,1); %Each node has two state variables - x1,x2 - for position and velocity.
eta=state(2*na+1:end,1); %Controller internal states

zeta_star_part = omega(1:ne); %\zeta^\star - \zeta_0
mu_zeta_star_part = omega(ne+1:end); % \mu_{\zeta^\star} - \mu_{\zeta_0}

y = C*x; %Position measurement
zeta = E'*y; %Relative positions
deta = -eta + (zeta-zeta_star_part); %Controller augmented dynamics
mu = tanh(eta) + mu_zeta_star_part; %Controller augmented measurement
u = -E*mu; %Input to nodes is (minus) the sum of the adjacent controllers' measurements.
dx = A*x + B*u + L*ws; %Linear Plant.

dstate=[dx; deta];

end