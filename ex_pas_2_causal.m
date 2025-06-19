dp = DepolarizingChannel(2,1);
id2 = eye(2);
idd = [1,0,0,1];
id_CJ = pure_to_mixed(idd);
psi = 1/sqrt(2)*[1,1];
ppsi = pure_to_mixed(psi);
p0 = pure_to_mixed([1,0]);
p1 = pure_to_mixed([0,1]);

d = [2,2,2,2,2,2,2,2,2];
p = {{1},{2,3},{4,5},{6,7},{[8,9]}};

Wabc = Tensor([0,1],psi,idd,idd,idd,[0,1]);
Wbca = Tensor([1,0],psi,idd,idd,idd,[1,0]);
Wbca = PermuteSystems(Wbca,[1,6,7,2,3,4,5,8,9],d);

w = Wabc + Wbca;
W = pure_to_mixed(w);

disp("W est valide : ");
disp((is_valid(W, d, p)) && is_PSD(W));

%A = 2, B = 3, C = 4
disp("W est compatible : B < C");
disp((is_compatible(W, p{3}, p{4}, d, p)))

%pas 2-Causal : rang 1 donc pas de décomposition convexe ?
%et compatible avec aucun des 3 premiers ou derniers

disp("W est compatible : A < {B,C}");
disp((is_compatible(W, p{2}, {p{3},p{4}}, d, p)))
disp("W est compatible : {B,C} < A");
disp((is_compatible(W, {p{3},p{4}}, p{2}, d, p)))
disp("W est compatible : B < {A,C}");
disp((is_compatible(W, p{3}, {p{2},p{4}}, d, p)))
disp("W est compatible : {A,C} < B");
disp((is_compatible(W, {p{2},p{4}}, p{3}, d, p)))
disp("W est compatible : C < {A,B}");
disp((is_compatible(W, p{4}, {p{2},p{3}}, d, p)))
disp("W est compatible : {A,B} < C");
disp((is_compatible(W, {p{2},p{3}}, p{4}, d, p)))

disp("W est 2-Causal")
disp(is_2causal(W, d, p))