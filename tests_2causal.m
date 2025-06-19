disp("switch 2-causal ?");
[W, d, p] = quantum_switch(2,2,false, Tensor(pure_to_mixed([1,1]), pure_to_mixed([1,0])));

ops = sdpsettings('solver','mosek','verbose',0);


r = sdpvar(1);
W2 = W + r*eye(length(W));
cst = is_2causal(W2, d, p);
optimize(cst, norm(r), ops);
disp("2causal : ");
r

disp("grenoble 2-causal ?");
[W, d, p] = Grenoble_process(true, Tensor(pure_to_mixed(c), pure_to_mixed([1,0])));
r = sdpvar(1);
W2 = (1-r)*W + r*eye(length(W));
optimize(is_2causal(W2, d, p), norm(r), ops)
disp("2causal : ");
r

dp = DepolarizingChannel(2,1);
id2 = eye(2);
idd = [1,0,0,1];
id_CJ = pure_to_mixed(idd);
psi = 1/sqrt(2)*[1,1];
ppsi = pure_to_mixed(psi);
p0 = pure_to_mixed([1,0]);
p1 = pure_to_mixed([0,1]);

disp("swich with C 2-causal ?");

d = [1,2,2,2,2,2,2,1,1];
p = {{1},{2,3},{4,5},{[6,7],8},{9}};

Wa = Tensor(psi,idd,idd,[1,0]);
Wb = Tensor(psi,idd,idd,[0,1]);
Wb = PermuteSystems(Wb,[1,4,5,2,3,6,7,8,9],d);
w = Wa + Wb;
W = pure_to_mixed(w);

disp("valide : ");
disp(is_valid(W, d, p) && is_PSD(W));
disp("{A,B} < C : ");
disp(is_compatible(W, {p{2},p{3}}, p{4}, d, p));
disp("2-causal : ");
r = sdpvar(1);
W2 = W + r*eye(length(W));
optimize(is_2causal(W2, d, p), norm(r), ops);
disp("2causal : ");
r

[W, d, p] = Lugano_process();
W2 = W + r*eye(length(W));
cst = is_2causal(W2, d, p);
optimize(cst, norm(r), ops);
disp("2causal : ");
r
