% This code calculates the AC Stark shifts on the Zeeman qubit of Yb174+.
% This will then used to find the frequencies which null the differential 
% shift and determine appropriate operating conditions for entanglement gates.

clear all

% Physical constants
hbar = 1.0545718e-34;       % Planck's constant [J.s]
c = 2.99792458e8;           % speed of light
eps0 = 8.85418782e-12;      % permittivity of free space


% -------------- calculate dbar matrix elements ------------------------ %

% --- Ion data: S1/2 - P1/2 (cooling transition)
gamma12 = 1/8.07e-9/(2*pi);     % Natural line width of cooling transition [Hz] % ---------Ref?
lambda12 = 369.53e-9;           % Yb cooling wavelength between S1/2 -> P1/2
f12 = c/lambda12;               % Yb cooling frequency [Hz]
w12 = 2*pi*f12;                 % Yb cooling frequency in [rad/s]

Jg = 1/2;                       % Total angular momentum J of S1/2 ground states
Lg = 0;                         % Orbital angular momentum L of ground states
Le = 1;                         % Orbital angular momentum L of excited states
Sg = 1/2;                       % Spin angular momentum S of ground states
Je12 = 1/2;                     % Total angular momentum J of P1/2 excited state

% --- Calculate the double bar matrix element
% This formula comes from Steck (equation 7.242), where the <Jg||d||Je> element 
% has been decomposed into the <Lg||d||Le> matrix element (equation 7.274)
dbar12 = sqrt(3*pi*eps0*hbar*c^3*gamma12/w12^3/(2*Jg+1)/(2*Lg+1)/Wigner6j(Lg,Le,1,Je12,Jg,Sg));

% --- Ion data: S1/2 - P3/2 transition
gamma32 = 1/6.15e-9/(2*pi);     % Natural line width of S1/2 - P3/2 transition [Hz] % -----Ref?
lambda32 = 328e-9;              % Yb S1/2 - P3/2 wavelength
f32 = c/lambda32;               % Yb S1/2 - P3/2 frequency [Hz]
w32 = 2*pi*f32;                 % Yb S1/2 - P3/2 frequency in [rad/s]
Je32 = 3/2;                     % Total J angular momentum of P1/2 excited state

% Calculate the double bar matrix element as above
dbar32 = sqrt(3*pi*eps0*hbar*c^3*gamma32/w32^3/(2*Jg+1)/(2*Lg+1)/abs(Wigner6j(Lg,Le,1,Je32,Jg,Sg)));

% create muvec (length 12)
% (1:2 for S1/2--P1/2 transitions, 3:6 for S1/2--P3/2 transitions)
muvec = [dbar12*ones(1,2) dbar32*ones(1,4)];

% Note: the dbar matrix elements give us an *absolute* measure of the
% strength of the corresponding transitions. The relative transition
% strengths between different terms/levels will be found using the
% Clebsch-Gordan coefficients and the Wigner-Eckart theorem (for further
% details, see Steck (equation


% ---------- the counter-propagating beams ----------------------------- %

% Field strengths of counter-propagating beams
P = 20e-3;                      % Power in laser beam [W]
r = 50e-6;                      % Laser beam radius
A = pi*r^2;                     % Laser beam cross-sectional area
E0 = sqrt(2*P/eps0/c/A);        % Electric field strength

% --- Polarization vectors of counter-propagating fields
% Here I use the NIST jargon of red-raman and blue-raman for the lowest 
% energy and highest energy beams respectively

ss = 1;                       % for example
pp = sqrt(1-2*ss^2);            % for example

pr = [ss pp ss];                % [sig- pi sig+]
pr = pr/sqrt(pr*pr.');          % normalize to 1

pb = [ss pp ss];                % [sig- pi sig+]
pb = pb/sqrt(pb*pb.');          % normalize to 1


% ---------------- Zeeman shifts --------------------------------------- %

% B = 20;                                 % applied magnetic field in units of Gauss as required for code
% 
% % splitting of S1/2
% [vS12, eS12] = zeeman_Yb_I_12(B);       % calculate eigenvectors and eigenvalues
% deS12 = diag(eS12);                     % extract eigenvalues (energies)
% deS12 = deS12 + abs(deS12(1));          % Set lowest m_F = 0 level to zero energy
% 
% % splitting of P1/2
% [vP12, eP12] = zeeman_P12_Yb_I_12(B);
% deP12 = diag(eP12);
% deP12 = deP12 + abs(deP12(1));          % Set lowest m_F = 0 level to zero energy
% 
% % splitting of P3/2
% [vP32, eP32] = zeeman_P32_Yb_I_12(B);
% deP32 = diag(eP32);
% deP32 = deP32 + abs(deP32(2));          % Set lowest m_F = 0 level to zero energy,
%                                         % but be careful, that is eP32(2), not eP3291).
%                                         
% Evec = [deS12.' (f12+deP12).' (f32 + deP32).'];   % vector of absolute FREQUENCY [Hz] of shifted levels
% lvec = c./Evec(2:end);                            % vector of wavelengths

% without taking Zeeman splitting into account, for the levels:
% [ S1/2 mJ = -1/2, +1/2,  P1/2 mJ = -1/2, +1/2,  P3/2 mJ = -3/2, -1/2, +1/2, +3/2 ]

% make up splitting to test code
spl = 10e6;

Evec = [0 0 c/lambda12 c/lambda12+spl c/lambda32-2*spl c/lambda32-spl c/lambda32+spl c/lambda32+2*spl];
lvec = c./Evec(2:end);


% --------------- relative transition strenths ------------------------- %

% ------ Yb 174:

% transition strengths from up and down states S1/2 mJ = -1/2, +1/2 to upper levels, in the order:
% [P1/2:  -1/2, +1,2, P3/2: -3/2, -1/2, +1,2, +3/2]

% using the Matlab script:
% [trd, tru] = relative_transition_strengths_Yb174
trd = [8/9 16/9 12 8 4 0];
tru = [16/9 8/9 0 4 8 12];

% % --- sort fields by the transitions that they couple to:
polvec_dr = [pr(2) pr(3) pr 0];     % for transitions from down, red raman beam
polvec_db = [pb(2) pb(3) pb 0];     % for transitions from down, blue raman beam

polvec_ur = [pr(1) pr(2) 0 pr];     % for transitions from up, red
polvec_ub = [pb(1) pb(2) 0 pb];     % for transitions from up, blue


% ----------------- calculate AC Stark shift --------------------------- %

% --- Calculate the vector of Rabi frequencies for each transition
% (calculate transitions from up and down states separately)

% Rabi frequencies for down, with red and blue respectively 
Rabi_dr = E0*polvec_dr.*trd.*muvec/hbar;
Rabi_db = E0*polvec_db.*trd.*muvec/hbar;

% Rabi frequencies for up, with red and blue respectively 
Rabi_ur = E0*polvec_ur.*tru.*muvec/hbar;   
Rabi_ub = E0*polvec_ub.*tru.*muvec/hbar;

% --- Calculate AC Stark shfit for the up and down states, for range of detuning:

lambda_g = 369.535e-9;              % wavelength of gate beams [nm]
fg0 = c/lambda_g;                   % frequency of gate beams [Hz]
fg328 = c/lambda32;
tunevec = (-400:1e-1:400)*1e9;      % range of detuning [Hz]
fg = fg0 + tunevec;                 % ...displaced by gate frequency fg0 (i.e., around P1/2)

% some extra stuff:
% tunevec = (9:1e-5:13)*1e9;        % range of detuning [Hz], if you want
                                    % to see the splitting 
% fg = fg328 + tunevec;             % detining displaced by fg328 (i.e., around P3/2)

% to find the total AC Stark shift on each qubit level, sum over all contributions:
nd = length(tunevec);
Eac_dtot = zeros(1,nd);
Eac_utot = zeros(1,nd);

for jj=1:nd
    % qubit downstate detunings
    deltadown = 2*pi*(fg(jj) - Evec(3:end));            % rad/s

    % qubit upstate detunings
    deltaup = 2*pi*(fg(jj) + Evec(2) - Evec(3:end));    % rad/s

    % Stark shifts for down
    Eac_d = Rabi_dr.^2./(4*deltadown) + Rabi_db.^2./(4*deltadown);
    Eac_dtot(jj) = sum(Eac_d);

    % Stark shifts for up
    Eac_u = Rabi_ur.^2./(4*deltaup) + Rabi_ub.^2./(4*deltaup);
    Eac_utot(jj) = sum(Eac_u);
end    

figure(1), clf
plot(tunevec/1e9,Eac_dtot/(2*pi),'r-', 'Linewidth', 1.5), hold on
plot(tunevec/1e9,Eac_utot/(2*pi),'k-', 'Linewidth', 1.5)
hold off
xlim([tunevec(1)/1e9, tunevec(end)/1e9])
xlabel('Detuning [GHz]')
ylabel('Stark shift [Hz]')
legend('Stark shift on down','Stark shift on up')
set(gca, 'Fontsize', 15)

AC_diff_1 = abs(Eac_dtot(1) - Eac_utot(1))
AC_diff_end = abs(Eac_dtot(end) - Eac_utot(end))
abs_shift_down = Eac_dtot(end)
abs_shift_up = Eac_utot(end)