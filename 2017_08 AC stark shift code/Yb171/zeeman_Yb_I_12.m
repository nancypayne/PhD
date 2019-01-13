function [V1,e] = zeeman_Yb_I_12(B)
% This function calculates the Zeeman splitting of the S level of  Yb+ (171) ion.
% It returns matrices V1 and e of the energies and eigenstates of the mf basis, respectively, for a
% magnetic field B. Energy values are the diagonal elements of the energy
% matrices. ordered low to high, i.e. energy(1,1) corresponds to mf = -1, etc.

Ahfs = 12.6428121185e9; % Hyperfine interaction constant In Hz
% as given by equation 6.10, pg 99 in Foot?
% or taken from experimental data, I assume makes more sense

L = 0;
J = 1/2;
I = 1/2;
S = 1/2;
F1 = 0;
F2 = 1;

meuB = 1.399624624e6; % In Hz/Gauss

gJ = 2.00226206;
gI = 0.0002134779853*gJ;

mi = [-1/2, 1/2];
mj = [-1/2, 1/2];

% writing down the Hamiltonian matrix:
H = zeros(length(mi)*length(mj));

for i1 = 1:length(mi)
    for j1 = 1:length(mj)
        for i2 = 1:length(mi)
            for j2 = 1:length(mj)  
                if mi(i1) == mi(i2)
                    if mj(j1) == mj(j2)
                        H((i1-1)*length(mj)+j1,(i2-1)*length(mj)+j2) = Ahfs*mi(i1)*mj(j1)+meuB*B*(gJ*mj(j1)+gI*mi(i1));
                    end
                end
                if i1 < length(mi)
                    if mi(i1+1) == mi(i2)
                        if j1 > 1
                            if mj(j1-1) == mj(j2)
                                H((i1-1)*length(mj)+j1,(i2-1)*length(mj)+j2) = Ahfs/2*sqrt((J+mj(j1))*(J-mj(j1)+1)*(I-mi(i1))*(I+mi(i1)+1));
                            end
                        end
                    end
                end
                if j1 < length(mj)
                    if mj(j1+1) == mj(j2)
                        if i1 > 1
                            if mi(i1-1) == mi(i2)
                                H((i1-1)*length(mj)+j1,(i2-1)*length(mj)+j2)=Ahfs/2*sqrt((J-mj(j1))*(J+mj(j1)+1)*(I+mi(i1))*(I-mi(i1)+1));
                            end
                        end
                    end
                end
            end
        end
    end
end
[V1,e] = eig(H);
              
