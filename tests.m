ineqs = readmatrix('data/interesting.txt');

% best found instruments in dim 2*2
idd = [1,0,0,1];
id_CJ = pure_to_mixed(idd);
p0 = pure_to_mixed([1,0]);
p1 = pure_to_mixed([0,1]);
t01 = Tensor(p0, p1);
t10 = Tensor(p1, p0);
M = {{{id_CJ, zeros(4,4)}, {t01, t10}}, 
    {{id_CJ, zeros(4,4)}, {t01, t10}}, 
    {{id_CJ, zeros(4,4)}, {t01, t10}}};

% canonical instruments from paper
% change in max_inequality : d = [2, 4, 2, 4, 2, 4];
% Ma = {{Tensor(eye(4)/2, p0), zeros(8,8)}, {Tensor(Tensor(p0,p0),p0), Tensor(Tensor(p1,p1),p0)}};
% M = {Ma, Ma, Ma};

% simple instruments
% M = {{{eye(4)/2, zeros(4,4)}, {eye(4)/2, zeros(4,4)}}, 
%     {{eye(4)/2, zeros(4,4)}, {eye(4)/2, zeros(4,4)}}, 
%     {{eye(4)/2, zeros(4,4)}, {eye(4)/2, zeros(4,4)}}};

cls = {'causal', '2-causal', 'some-causality', 'valid'};


for i = 1:size(ineqs, 1)
    ineq = ineqs(i, :);
    disp(['Inequality ', num2str(i), ':']);
    for j = 1:length(cls)
        class = cls{j};
        disp(['Class: ', class]);
        [W, M] = max_inequality(ineq, class, 1, false, M);
        s = score(calc_proba(W, M), ineq);
        disp(['Score: ', num2str(s)]);
    end
end