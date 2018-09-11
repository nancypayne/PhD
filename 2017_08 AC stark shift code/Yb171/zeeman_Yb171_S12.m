% function [V1,e] = zeeman_Yb171_S12(B)

% This function calculates the Zeeman splitting of the S level of  Yb+ (171) ion.
% It returns matrices V1 and e of the energies and eigenstates of the mf basis, respectively, for a
% magnetic field B. Energy values are the diagonal elements of the energy
% matrices. ordered low to high, i.e. energy(1,1) corresponds to mf = -1, etc.

B = 20;

A_hfs = 12.6428121185e9; % Hyperfine interaction constant In Hz
% as given by equation 6.10, pg 99 in Foot?
% or taken from experimental data, I assume makes more sense

% electron proton mass ratio
mass_ratio = 1836.15267389;     % 1/mass_ratio = 5.4462e-4;
me = 9.10938188e-31;            % Electron mass - Kg
u = 1.66053873e-27;             % Atomic mass unit - Kg
m = 173.9388621*u;              % mass of Yb 174

uB = 1.399624624e6;             % Bohr magneton in Hz/Gauss (uB/h)
% NIST says the above is 1.3996 245 042 e6
% compared to the above, 1.3996 246 24  e6

L = 0;      LL = L*(L+1);
J = 1/2;    JJ = J*(J+1);
I = 1/2;    II = I*(I+1);
S = 1/2;    SS = S*(S+1);

F  = [0  1 1 1];    % F of each level you want to calculate shift for
mF = [0 -1 0 1];    % mF of the corresponding levels above
FF = F.*(F+1);

%----------- gL, gS, gJ

% gL (electron orbital g-factor)
% gL = 1 - 1/M, where M is nucleus electron mass ratio
gL = 1 - me/m;   % close to 1

% gJ = 2.00226206;          % original code

% from Steck eqn 7.305:
gJ = gL + (gS-gL)*(JJ + SS - LL)/(2*JJ);
    
% ------- gI

% gI = 0.0002134779853*gJ;  % original code

% I'm pretty sure gI should be smaller by the proton/electron mass ratio,
% 1/mass_ratio = 0.0005, not sure what this 0.0002 number is. 
gI = gJ/mass_ratio; % approx 0.001    
    
% -------- gF - Steck eqn 7.309
idx = find(mF ~= 0);    % find the elements where mF = 0 (no shift)
gF = gJ + (gI-gJ)*(FF(idx) + II - JJ)./(2*FF(idx));


% --------------- splitting in weak field (Zeeman effect) --------------- %

E_hfs = zeros(size(mF));        % if mJ = 0, shift is zero
E_hfs(idx) = uB*gF.*mF(idx)*B;  % calculate shift for mJ ~= 0


% --------------- splitting in intermediate field ----------------------- %

% for details, see Steck section 7.4.1.3

mI = [-1/2 -1/2 1/2 1/2];
mJ = [-1/2 1/2 -1/2 1/2];

% the diagonal elements of the Hamiltonian have two contributions:
% < H_hfs >   = A_hfs*mI*mJ
% < H_B_hfs > = uB*(gJ*mJ + gI*mI)*B

for i = 1:4
    H_diag(i) = A_hfs*mI(i)*mJ(i) + uB*(gJ*mJ(i) + gI*mI(i))*B;
end

% the upper diagonal element of the Hamiltonian (Steck, eqn 7.321)
r = 2;  % row number
H_diag_u = A_hfs/2*sqrt((J-mJ(r)+1)*(J+mJ(r))*(I+mI(r)+1)*(I-mI(r)));

% the lower diagonal element of the Hamiltonian (Steck, eqn 7.321)
r = 3;
H_diag_l = A_hfs/2*sqrt((J+mJ(r)+1)*(J-mJ(r))*(I-mI(r)+1)*(I+mI(r)));

% put together into H;
H_hfs = diag(H_diag);
H_hfs(2,3) = H_diag_u;
H_hfs(3,2) = H_diag_l;

% find the eigenvalues (energies)
[v,e] = eig(H);