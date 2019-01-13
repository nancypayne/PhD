% apply code to set of levels

% format of M = [ F1 Mf1 F2 Mf2 ];

M = [ 0 0 1  0 ;    % F = 0 to 1
    0  0 1  1 ;
    1 -1 1 -1 ;     % F = 1 to 1, pi
    1  1 1  1 ;
    1  -1 1 0 ;     % F = 1 to 1, sigma+
    1  0 1  1 ;
    1 -1 2 -1 ;     % F = 1 to 2, pi
    1  0 2  0 ;
    1  1 2  1 ;
    1 -1 2  0 ;     % F = 1 to 2, sigma+
    1  0 2  1 ;
    1  1 2  2 ];

t = zeros(size(M,1), 1);
cg = t;
we = t;

for i = 1 : size(M,1)
    [t(i), cg(i), we(i)] = relative_transition_strengths(M(i,1), M(i,2), M(i,3), M(i,4));
end

cg_sq = cg.^2;
cg_sq = cg_sq / min(cg_sq);
we_sq = we.^2;
we_sq = we_sq / min(we_sq);
rel_strengths = t.^2;
rel_strengths = rel_strengths / min(rel_strengths);

fprintf('\t\t F1 \t Mf1 \t  F2 \t  Mf2 \t\t cg^2 \t   we^2 \t Relative \n')
results = [M cg_sq we_sq rel_strengths]