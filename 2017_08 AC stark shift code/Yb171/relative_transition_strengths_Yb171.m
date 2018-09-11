function [trd, tru] = relative_transition_strengths_Yb171

% for Yb 171
% calculate relative transition strengths from S1/2 F,mF = (0,0) and (1,0)
% to all P1/2 and P3/2 hyperfine levels, in the order:
% [P1/2 (0,0) (1,-1) (1,0) (1,+1),
%  P3/2 (1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2)]

% --- the output:
% the output is a matrix (tr), where each row gives the transition strengths from a
% particular state (in this case row 1 corresponds to transitions from the down state,
% and row 2 corresponds to transitions from the up state), to the upper hyperfine 
% levels in the following order:
% [P1/2 (0,0) (1,-1) (1,0) (1,+1), P3/2 (1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2)]
% In the same language as Hermann's code, his trd is my tr(1,:).

clear all

% inputs (these should be tidied, or more automatic... info is degenerate)
I = 1/2;
S = 1/2;
L_vec = [0, 1]; % L for all levels you want to calculate transitions between (here, S and P)
Jg = 1/2;       % J of ground states
Lg = L_vec(1);  % L of ground states  
Le = L_vec(2);  % L of excited states, fine for now as all excited states are L = 1
Fg = [0,1];     % F of ground states
mFg = 0;        % mF of ground states
k = 1;          % k = 1 for dipole transitions

% -------------------------------------------------------------------- %

% find vector of J (from low to high energy)
J_vec = [];

for i = 1 : length(L_vec)
    min_J = abs(L_vec(i) - S);
    max_J = L_vec(i) + S;
    J = min_J : 1 : max_J;
    J_vec = [J_vec J];
end

% find vector of Fs for levels (from low to high energy)
% (will Fs always be increasing? i.e., is F = 2 always higher in energy than F = 1?)
F_vec = [];
JJ_vec = [];

for i = 1 : length(J_vec)
    min_F = abs(I - J_vec(i));
    max_F = I + J_vec(i);
    F = min_F : 1 : max_F;
    F_vec = [F_vec F];
    
    % for every entry made in F_vec, record the corresponding J value
    Length = length(F);
    JJ_vec = [JJ_vec J_vec(i)*ones(1,Length)];
end

% now we can work out the relative transition strengths between all levels,
cg_vec = 0; % initialise some things for the loop
we_vec = 0;
count = 0;
tr = [];

for i = 1 : length(Fg) % i.e., for Fg = 0 (down state) and Fg = 1 (up state)
    count = 0;
    
    for n = 3 : length(F_vec) % i.e., only consider transitions to P1/2(0,0) and higher
        Fe = F_vec(n);  % note: g denotes ground state, e denotes excited state
        Je = JJ_vec(n);
        
        % find the Wigner-Eckart factor, breaking it up
        a1 = Je + Lg + k + S;       % from decomposition of F = I + J
        a2 = Fe + Jg + k + S;       % from decomposition of J = L + S
        a = (-1)^(a1+a2);
        b1 = (2*Fe+1) * (2*Jg+1);   % from F = I + J
        b2 = (2*Je+1) * (2*Lg+1);   % from J = L + S
        b = sqrt(b1*b2);
        c1 = Wigner6j(Jg, Je, k, Fe, Fg(i), I);    % from F = I + J
        c2 = Wigner6j(Lg, Le, k, Je, Jg, S);       % from J = L + S
        c = c1*c2;
        we = a*b*c; % combine everything together
        
        for mFe = -Fe : 1 : Fe % look at transitions to every mF level within each F
            count = count + 1;
            % find the Clebsch Gordan coefficient
            cg = ClebschGordan(Fe, 1, Fg(i), mFe, mFg-mFe, mFg);
            cg_vec(count) = cg;
            % create corresponding vector for the Wigner-Eckart factor
            we_vec(count) = we;
        end
    end
    rel_trans_strengths = cg_vec.^2 .* we_vec.^2;
    we_sq = we_vec.^2
    cg_sq = cg_vec.^2
    tr = [tr ; rel_trans_strengths];
end

% in the same wording as Hermann:
trd = tr(1,:);
tru = tr(2,:);

% exagerated energy level scheme
Evec = [ 0 0.8 1 1.2   3 3.8 4 4.2   6.8 7 7.2  8.3 8.5 8.7 8.9 9.1 ];
LL = length(trd);

figure(1), clf, hold on
for i = 1 : length(Evec)
    plot([-2 18], [Evec(i), Evec(i)], 'k-')
end
for j = 1 : LL
    plot([0.5*j 0.5*j], [0, Evec(j+4)], 'r-')
    text(0.5*j - 0.15, 1.1+0.4*j, num2str(trd(j)))
end
for k = 1 : LL
    plot([9+0.5*k 9+0.5*k], [Evec(3), Evec(k+4)], 'r-')
    text(8.85+0.5*k, 1.1+0.4*k, num2str(tru(k)))
end
hold off
xlim([-2 18])
ylim([-1 10])
shg

end
