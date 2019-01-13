% This code calculates the AC Stark shifts on the |F=0, m_F=0> to
% |F=1, m_F=0> hyperfine qubit transition of Yb171+. This is then used to
% find the frequencies which null the differential shift and determine
% appropriate operating conditions for entanglement gates.

% Nancy update: also being edited to do Yb174 in places (in progress)

clear all

% Physical constants
hbar = 1.0545718e-34;       % Planck's constant [J.s]
c = 2.99792458e8;           % speed of light
eps0 = 8.85418782e-12;      % permittivity of free space


% -------------- calculate dbar matrix elements ------------------------ %

% --- Ion data: S1/2 - P1/2 (cooling transition)
gamma12 = 1/8.07e-9/(2*pi);     % Natural line width of cooling transition [Hz] Ref Li et al, J. Phys. B 32, 1731 (1999)
lambda12 = 369.53e-9;           % Yb cooling wavelength between S1/2, F=1 -> P1/2, F=0
f12 = c/lambda12;               % Yb cooling frequency [Hz]
w12 = 2*pi*f12;                 % Yb cooling frequency in [rad/s]

Jg = 1/2;                       % Total angular momentum J of S1/2 ground states
Lg = 0;                         % Orbital angular momentum L of ground states
Le = 1;                         % Orbital angular momentum L of excited states
Sg = 1/2;                       % Spin angular momentum S of ground states
Je12 = 1/2;                     % Total angular momentum J of P1/2 excited state

% --- Calculate the double bar matrix element
% This formula comes from Steck (equation 7.242), where the <Jg||d||Je> element 
% has been reduced into the <Lg||d||Le> matrix element (equation 7.252)
dbar12 = sqrt(3*pi*eps0*hbar*c^3*gamma12/w12^3/(2*Jg+1)/(2*Lg+1)/Wigner6j(Lg,Le,1,Je12,Jg,Sg));

% --- Ion data: S1/2 - P3/2 transition
gamma32 = 1/6.15e-9/(2*pi);     % Natural line width of S1/2 - P3/2 transition [Hz] Ref Li et al, J. Phys. B 32, 1731 (1999)
lambda32 = 328e-9;              % Yb S1/2 - P3/2 wavelength
f32 = c/lambda32;               % Yb S1/2 - P3/2 frequency [Hz]
w32 = 2*pi*f32;                 % Yb S1/2 - P3/2 frequency in [rad/s]
Je32 = 3/2;                     % Total J angular momentum of P1/2 excited state

% Calculate the double bar matrix element as above
dbar32 = sqrt(3*pi*eps0*hbar*c^3*gamma32/w32^3/(2*Jg+1)/(2*Lg+1)/abs(Wigner6j(Lg,Le,1,Je32,Jg,Sg)));

% create muvec (length 12)
% (1:4 for S1/2--P1/2 transitions, 5:12 for S1/2--P3/2 transitions)
muvec = [dbar12*ones(1,4) dbar32*ones(1,8)];  


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

B = 20;                                 % applied magnetic field in units of Gauss as required for code

% splitting of S1/2
[vS12, eS12] = zeeman_Yb_I_12(B);       % calculate eigenvectors and eigenvalues
deS12 = diag(eS12);                     % extract eigenvalues (energies)
deS12 = deS12 + abs(deS12(1));          % Set lowest m_F = 0 level to zero energy

% splitting of P1/2
[vP12, eP12] = zeeman_P12_Yb_I_12(B);
deP12 = diag(eP12);
deP12 = deP12 + abs(deP12(1));          % Set lowest m_F = 0 level to zero energy

% splitting of P3/2
[vP32, eP32] = zeeman_P32_Yb_I_12(B);
deP32 = diag(eP32);
deP32 = deP32 + abs(deP32(2));          % Set lowest m_F = 0 level to zero energy,
                                        % but be careful, that is eP32(2), not eP3291).

Evec = [deS12.' (f12+deP12).' (f32 + deP32).'];   % vector of absolute energies of shifted levels
lvec = c./Evec(2:end);                            % vector of wavelengths


% --------------- relative transition strenths ------------------------- %

% ------ Yb 171:
% transition strengths from up and down states S1/2 (1,0) and (0,0) to upper levels, in the order:
% [P1/2 (0,0) (1,-1) (1,0) (1,+1),
%  P3/2 (1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2)]

% using the Matlab script:
% [trd, tru] = relative_transition_strengths_Yb171;
trd = [0 1/9 1/9 1/9 2/9 2/9 2/9 0 0 0 0 0];
tru = [1/9 1/9 0 1/9 1/18 0 1/18 0 3/18 2/9 3/18 0];

% --- Sort fields by the transitions that they couple to:
% e.g., polvecr(1) = pr(2) = pi, i.e., transition from S1/2 (0,0) to P1/2 (0,0)
% e.g., polvecr(2:4) = pr = [sigma-, pi, sigma+], i.e., transitions from (0,0) to (1,-1),(1,0),(1,1)

polvecr = [pr(2) pr pr 0 pr 0];
polvecb = [pb(2) pb pb 0 pb 0];


% ----------------- calculate AC Stark shift --------------------------- %

% --- Calculate the vector of Rabi frequencies for each transition
% (calculate transitions from up and down states separately)

% Rabi frequencies for down 
Rabidr = E0*polvecr.*trd.*muvec/hbar;
Rabidb = E0*polvecb.*trd.*muvec/hbar;

% Rabi frequencies for up
Rabiur = E0*polvecr.*tru.*muvec/hbar;
Rabiub = E0*polvecb.*tru.*muvec/hbar;

% --- Calculate AC Stark shfit for the up and down states, for range of detuning:

lambda_g = 369.5300000001e-9;       % wavelength of gate beams [nm]
fg0 = c/lambda_g;                   % frequency of gate beams [Hz]
fg328 = c/lambda32;
tunevec = (-200:1e-1:200)*1e9;        % range of detuning... [Hz]
fg = fg0 + tunevec;                 % ...displaced by gate frequency fg0 (i.e., around P1/2)
% fg = fg328 + tunevec;             % ...or displaced by fg328 (i.e., around P3/2)

% to find the total AC Stark shift on each qubit level, sum over all contributions:
nd = length(tunevec);
Eacdtot = zeros(1,nd);
Eacutot = zeros(1,nd);

for jj = 1:nd
    % qubit downstate detunings
    deltadown = 2*pi*(fg(jj) - Evec(5:end));            % rad/s

    % qubit upstate detunings
    deltaup = 2*pi*(fg(jj) + Evec(3) - Evec(5:end));    % rad/s

    % Stark shifts for down
    Eac_d = Rabidr.^2./(4*deltadown) + Rabidb.^2./(4*deltadown);
    Eac_dtot(jj) = sum(Eac_d);

    % Stark shifts for up
    Eac_u = Rabiur.^2./(4*deltaup) + Rabiub.^2./(4*deltaup);
    Eac_utot(jj) = sum(Eac_u);
end    

figure(1)
plot(tunevec/1e9,Eac_dtot/(2*pi),'r', 'Linewidth', 1.5), hold on
plot(tunevec/1e9,Eac_utot/(2*pi),'k', 'Linewidth', 1.5)
hold off
xlim([tunevec(1)/1e9, tunevec(end)/1e9])
xlabel('Detuning [GHz]')
ylabel('Stark shift [Hz]')
legend('Stark shift on down','Stark shift on up', 'Location', 'Northeast')
set(gca, 'Fontsize', 15)

AC_diff_1 = Eac_dtot(1) - Eac_utot(1)
AC_diff_end = Eac_dtot(end) - Eac_utot(end)
abs_AC_shift_down = Eac_dtot(1) 
abs_AC_shift_up = Eac_utot(1)