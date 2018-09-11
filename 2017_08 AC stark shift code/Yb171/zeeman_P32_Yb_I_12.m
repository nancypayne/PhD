function [V1,e] = zeeman_P32_Yb_I_12(B)

% Still to be adapted to ytterbium

Ahfs = 0.877e9; %In Hz [Ref - Berends and Maleki, J. Opt. Soc. Am. B]
L = 1;
J = 3/2;
I = 1/2;
S = 1/2;

meuB = 1.399624624e6;   % In Hz/Gauss
me = 9.10938188e-31;    % Electron mass - Kg
u = 1.66053873e-27;     % Atomic mass unit - Kg
m = 170.936323*u;       % Appx. mass of Be9+

gL = 1-me/m;
gS = 2.0023193043737;
gJ = gL*(J*(J+1)-S*(S+1)+L*(L+1))/(2*J*(J+1))+gS*(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1));
gI = 0.0002134779853*2.00226206; %Same value that we take for the 2S1/2.

mi = [-1/2, 1/2];
mj = [-3/2, -1/2, 1/2, 3/2];


% writing down the Hamiltonian matrix.
H = zeros(length(mi)*length(mj));

for i1 = 1:length(mi)
    for j1 = 1:length(mj)
      for i2 = 1:length(mi)
            for j2 = 1:length(mj)
                
                if mi(i1) == mi(i2) 
                    if mj(j1) == mj(j2)
                       H((i1-1)*length(mj)+j1,(i2-1)*length(mj)+j2)=Ahfs*mi(i1)*mj(j1)+meuB*B*(gJ*mj(j1)+gI*mi(i1));
                    end
                end
                if i1 < length(mi)
                    if mi(i1+1) == mi(i2)
                      if j1 > 1
                        if mj(j1-1) == mj(j2)
                        H((i1-1)*length(mj)+j1,(i2-1)*length(mj)+j2)=Ahfs/2*sqrt((J+mj(j1))*(J-mj(j1)+1)*(I-mi(i1))*(I+mi(i1)+1));
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


              
