function [trd, tru] = relative_transition_strengths_Yb174

% for Yb 174
% calculate relative transition strengths from S1/2 mJ = -1/2, +1,2
% to all P1/2 and P3/2 levels, in the order:
% [P1/2:  -1/2, +1,2, P3/2: -3/2, -1/2, +1,2, +3/2]

% --- the output:
% the output is a matrix (tr), where each row gives the transition strengths from a
% particular state (in this case row 1 corresponds to transitions from the down state,
% and row 2 corresponds to transitions from the up state), to the upper levels

% inputs:
I = 0;
S = 1/2;
L_vec = [0, 1];     % L for all levels you want to calculate transitions between (here, S and P)
Jg = 1/2;           % J of ground states
% mJg = [-1/2, 1/2];  % mJ of chosen ground states
Lg = L_vec(1);      % L of ground states  
Le = L_vec(2);      % L of excited states, fine for now as all excited states are L = 1
% Fg = [0,1];       % F of ground states
% mFg = 0;          % mF of ground states
k = 1;              % k = 1 for dipole transitions

% -------------------------------------------------------------------- %

% find vector of J (from low to high energy)
J_vec = [];

for i = 1 : length(L_vec)
    min_J = abs(L_vec(i) - S);
    max_J = L_vec(i) + S;
    J = min_J : 1 : max_J;
    J_vec = [J_vec J];
end

% now we can work out the relative transition strengths between all levels,
cg_vec = 0; % initialise some things for the loop
we_vec = 0;
count = 0;
tr = [];

for mJg = -Jg : 1 : Jg % i.e., for mJ = -1/2 (down state) and mJ = +1/2 (up state) 
    count = 0;
    
    for n = 2 : length(J_vec) % i.e., only consider transitions to P1/2 and higher
        Je = J_vec(n);
        
        % find the Wigner-Eckart factor, breaking it up
        % a1 = Fe + Jg + k + S;         % from decomposition of F = I + J
        a2 = Je + Lg + k + S;           % from decomposition of J = L + S
        a = (-1)^(a2);
        % b1 = (2*Fe+1) * (2*Jg+1);     % from F = I + J
        b2 = (2*Je+1) * (2*Lg+1);       % from J = L + S
        b = sqrt(b2);
        % c1 = Wigner6j(Jg, Je, k, Fe, Fg(i), I);   % from F = I + J
        c2 = Wigner6j(Lg, Le, k, Je, Jg, S);        % from J = L + S
        c = c2;
        we = a*b*c; % combine everything together
        
        for mJe = -Je : 1 : Je % look at transitions to every mJ level within each J
            count = count + 1;
            
            % find the Clebsch Gordan coefficient
            cg = ClebschGordan(Je, 1, Jg, mJe, mJg-mJe, mJg);
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
trd = tr(1,:)
tru = tr(2,:)

% exagerated energy level scheme
Evec = [ 0 0.3   3 3.3   6 6.3 6.6 6.9 ];
LL = length(trd);

figure(1), clf, hold on
for i = 1 : length(Evec)
    plot([-2 18], [Evec(i), Evec(i)], 'k-')
end
for j = 1 : LL
    plot([0.5*j 0.5*j], [0, Evec(j+2)], 'r-')
    text(0.5*j - 0.15, 0.35+0.4*j, num2str(trd(j)))
end
for k = 1 : LL
    plot([9+0.5*k 9+0.5*k], [Evec(2), Evec(k+2)], 'r-')
    text(8.85+0.5*k, 0.35+0.4*k, num2str(tru(k)))
end
hold off
xlim([-2 18])
ylim([-1 8])
shg

end
