function [t, cg, we] = relative_transition_strengths(F1, mF1, F2, mF2)

% input: F, Mf, I, J of each level
% output: relative transition strengths... up to what factor?

% ---- enter these by hand for now,
% ---- first entry ground state, second entry excited state
I = 0;
F = [F1 F2];
Mf = [mF1 mF2];
J = [1/2 3/2];
L = [0 1];
S = 1/2;
k = 1;          % rank of tensor, = 1 for dipole transitions

% --- first find the Clebsch Gordan coefficient
j1 = F(2);
j  = F(1);
j2 = k;         % j2 = k = 1 for dipole transitions
m1 = Mf(2);
m  = Mf(1);
m2 = m - m1;    % such that m1 + m2 = m

cg = ClebschGordan(j1, j2, j, m1, m2, m);

% --- now the Wigner-Eckart bit

% part a
a1 = F(2) + J(1) + k + I;   % from F = I + J
a2 = J(2) + L(1) + k + S;   % from J = L + S
a = (-1)^(a1+a2);

% part b
b1 = (2*F(2)+1) * (2*J(1)+1);   % from F = I + J
b2 = (2*J(2)+1) * (2*L(1)+1);   % from J = L + S
b = sqrt(b1*b2);

% part c
c1 = Wigner6j(J(1), J(2), k, F(2), F(1), I);    % from F = I + J
c2 = Wigner6j(L(1), L(2), k, J(2), J(1), S);    % from J = L + S
c = c1*c2;

we = a*b*c;
we2 = a2*b2*c2;
we_sq = we2.^2

% --- and altogether...

t = we*cg;

% --- yay, success!

end