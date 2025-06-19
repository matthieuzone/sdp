proba = zeros(8, 8);
for x = 1:8
    proba(x, 1) = 1;
end
ineq = readmatrix('interesting.txt');
ineq = ineq(1,:);
score(proba, ineq)

W = eye(8*8)/8;
M = {{{eye(4)/2, zeros(4,4)}, {eye(4)/2, zeros(4,4)}}, 
    {{eye(4)/2, zeros(4,4)}, {eye(4)/2, zeros(4,4)}}, 
    {{eye(4)/2, zeros(4,4)}, {eye(4)/2, zeros(4,4)}}};

probabis = calc_proba(W, M);
disp(probabis);
disp('Score:');
disp(score(probabis, ineq));