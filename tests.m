diary log
diary on

disp("grenoble");
c = 1/sqrt(3)*[1,1,1];
[W, d, p] = Grenoble_process(true, Tensor(pure_to_mixed(c), pure_to_mixed([1,0])));
disp("QCCC :");
disp(superop_random_robustness(W, d, p, 'QCCC'));
disp("2-causal :");
disp(superop_random_robustness(W, d, p, '2-causal'));
disp("some causality :");
disp(superop_random_robustness(W, d, p, 'some_causality'));


disp("swich with C");
id2 = eye(2);
idd = [1,0,0,1];
id_CJ = pure_to_mixed(idd);
psi = 1/sqrt(2)*[1,1];
ppsi = pure_to_mixed(psi);
p0 = pure_to_mixed([1,0]);
p1 = pure_to_mixed([0,1]);
d = [1,2,2,2,2,2,2,1,1];
p = {{1},{2,3},{4,5},{[6,7],8},{9}};
Wa = Tensor([1, 0],idd,idd,[1,0]);
Wb = Tensor([1, 0],idd,idd,[0,1]);
Wb = PermuteSystems(Wb,[1,4,5,2,3,6,7,8,9],d);
w = Wa + Wb;
W = pure_to_mixed(w);

disp("A < B : ")
disp(is_compatible(W, {p{2}}, {p{3}}, d, p))
disp("B < A : ")
disp(is_compatible(W, {p{3}}, {p{2}}, d, p))
disp("A < C : ")
disp(is_compatible(W, {p{2}}, {p{4}}, d, p))
disp("C < A : ")
disp(is_compatible(W, {p{4}}, {p{2}}, d, p))
disp("B < C : ")
disp(is_compatible(W, {p{3}}, {p{4}}, d, p))
disp("C < B : ")
disp(is_compatible(W, {p{4}}, {p{3}}, d, p))


disp("QCCC : ");
disp(superop_random_robustness(W, d, p, 'QCCC'));
disp("2-causal : ");
disp(superop_random_robustness(W, d, p, '2-causal'));
disp("some causality : ");
disp(superop_random_robustness(W, d, p, 'some_causality'));

disp("Lugano");
[W, d, p] = Lugano_process();
disp("QCCC : ");
disp(superop_random_robustness(W, d, p, 'QCCC'));
disp("2-causal : ");
disp(superop_random_robustness(W, d, p, '2-causal'));
disp("some causality : ");
disp(superop_random_robustness(W, d, p, 'some_causality'));

disp("example")
d = [1,2,2,2,2,2,2,2,2];
p = {{1},{2,3},{4,5},{6,7},{[8,9]}};
Wa = Tensor(psi, idd, idd, idd, [1,0])
Wb = Tensor(psi, idd, idd, idd, [0,1]);
Wb = PermuteSystems(Wb, [1,4,5,6,7,2,3,8,9], d);
w = Wa + Wb;
W = pure_to_mixed(w);

disp("A < B : ")
disp(is_compatible(W, {p{2}}, {p{3}}, d, p))
disp("B < A : ")
disp(is_compatible(W, {p{3}}, {p{2}}, d, p))
disp("A < C : ")
disp(is_compatible(W, {p{2}}, {p{4}}, d, p))
disp("C < A : ")
disp(is_compatible(W, {p{4}}, {p{2}}, d, p))
disp("B < C : ")
disp(is_compatible(W, {p{3}}, {p{4}}, d, p))
disp("C < B : ")
disp(is_compatible(W, {p{4}}, {p{3}}, d, p))

disp("QCCC : ");
disp(superop_random_robustness(W, d, p, 'QCCC'));
disp("2-causal : ");
disp(superop_random_robustness(W, d, p, '2-causal'));
disp("some causality : ");
disp(superop_random_robustness(W, d, p, 'some_causality'));

diary off