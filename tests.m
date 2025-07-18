%teste les violations sur les inégalités dans le dossier data, qui sont prises sous le format donné par le repo python causal-inequalities

ineqs = readmatrix('data/interesting.txt');
% ineqs = readmatrix('data/3lazy_some-causality_repr.txt');

%canonical instruments from paper
d = [2, 4, 2, 4, 2, 4];
p0 = pure_to_mixed([1,0]);
p1 = pure_to_mixed([0,1]);
Ma = {{Tensor(eye(4)/2, p0), zeros(8,8)}, {Tensor(Tensor(p0,p0),p0), Tensor(Tensor(p1,p1),p0)}};
M = {Ma, Ma, Ma};

% simple instruments
% M = {{{eye(4)/2, zeros(4,4)}, {eye(4)/2, zeros(4,4)}}, 
%     {{eye(4)/2, zeros(4,4)}, {eye(4)/2, zeros(4,4)}}, 
%     {{eye(4)/2, zeros(4,4)}, {eye(4)/2, zeros(4,4)}}};

cls = {'2-causal', 'some-causality', 'valid'};

for i = 1:size(ineqs, 1)
    ineq = ineqs(i, :);
    disp(['Inequality ', num2str(i), ':']);
    for j = 1:length(cls)
        class = cls{j};
        disp(['Class: ', class]);
        [W, M] = max_inequality(ineq, class, 1, false, M, d);
        s = score(calc_proba(W, M), ineq);
        disp(['Score: ', num2str(s)]);
    end
end