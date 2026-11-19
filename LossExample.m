%% LossExample_v1.m  (SCRIPT)
% Runs PF using V_G = 1 + gamma*v1 (magnitudes only), prints Pg/Qg for table

clear; clc;
define_constants;

%% MATPOWER path (EDIT IF NEEDED)
%% Check MATPOWER installation
if isempty(which('runpf')) || isempty(which('loadcase')) || isempty(which('makeYbus'))
    error('MATPOWER not found. Install MATPOWER and add it to the MATLAB path.');
end
disp('MATPOWER found.');

%% Build 2-area case: IEEE 14 duplicated -> 28 buses
mpc1 = loadcase('case14');
mpc2 = mpc1;

BUS_I = 1; F_BUS = 1; T_BUS = 2;
shift = size(mpc1.bus,1);

mpc2.bus(:, BUS_I)            = mpc2.bus(:, BUS_I) + shift;
mpc2.gen(:, GEN_BUS)          = mpc2.gen(:, GEN_BUS) + shift;
mpc2.branch(:, F_BUS)         = mpc2.branch(:, F_BUS) + shift;
mpc2.branch(:, T_BUS)         = mpc2.branch(:, T_BUS) + shift;

mpc = mpc1;
mpc.bus    = [mpc1.bus;    mpc2.bus];
mpc.gen    = [mpc1.gen;    mpc2.gen];
mpc.branch = [mpc1.branch; mpc2.branch];
if isfield(mpc1,'gencost') && isfield(mpc2,'gencost')
    mpc.gencost = [mpc1.gencost; mpc2.gencost];
end

% Fix optional metadata fields with wrong size
if isfield(mpc,'bus_name')
    mpc = rmfield(mpc,'bus_name');
end

fprintf('Combined case: %d buses, %d generators\n', size(mpc.bus,1), size(mpc.gen,1));

%% Compute YGGM
[Ybus, ~, ~] = makeYbus(mpc);

gen_buses  = unique(mpc.gen(:, GEN_BUS)).';
all_buses  = mpc.bus(:, BUS_I).';
load_buses = setdiff(all_buses, gen_buses);

n = numel(gen_buses);
m = numel(load_buses);
fprintf('n = %d generator buses, m = %d load buses\n', n, m);

[~, gen_rows]  = ismember(gen_buses,  mpc.bus(:, BUS_I));
[~, load_rows] = ismember(load_buses, mpc.bus(:, BUS_I));

YGG  = Ybus(gen_rows,  gen_rows);
YGL  = Ybus(gen_rows,  load_rows);
YLG  = Ybus(load_rows, gen_rows);
YLL  = Ybus(load_rows, load_rows);

YGGM = full(YGG - YGL*(YLL\YLG));

%% Nullspace basis and pick v1
tol_null = 1e-1;
Vbasis = null(YGGM, tol_null);


% Numerical rank and column permutation
s = svd(YGGM);
r = sum(s > tol_null);

[~, R, p] = qr(YGGM, 'vector');
P = eye(n);
P = P(:, p);

% Approximate basis in the form of Theorem 2.1
W = -R(1:r,1:r) \ R(1:r,r+1:n);
Z = P * [W; eye(n-r)];

disp('Generator-bus ordering:');
disp(gen_buses);

disp('Column permutation:');
disp(p);

disp('Matrix W:');
disp(W);

fprintf('Numerical rank: %d\n', r);

for j = 1:size(Z,2)
    residual = norm(YGGM*Z(:,j),2) / norm(Z(:,j),2);
    fprintf('Basis %d normalized residual: %.12e\n', ...
            j, residual);
end


if size(Vbasis,2) < 2
    error('Nullspace dimension < 2. Increase tol_null slightly.');
end

v1 = Vbasis(:,1);
v1 = v1 / max(abs(v1));   % normalize

%% Build a feasible generator-voltage profile from v1
gamma = 0.10;                           % small mixing factor
VGprof = ones(n,1) + gamma*v1;          % complex profile
VGset  = abs(VGprof);                   % MATPOWER uses magnitudes

%% Apply VG magnitudes to PV/slack setpoints
mpc_v1 = mpc;
for k = 1:n
    bnum = gen_buses(k);
    bus_row = find(mpc_v1.bus(:, BUS_I) == bnum, 1);
    mpc_v1.bus(bus_row, VM) = VGset(k);     % bus voltage magnitude

    idx = find(mpc_v1.gen(:, GEN_BUS) == bnum);
    mpc_v1.gen(idx, VG) = VGset(k);         % gen voltage setpoint
end

%% Run PF
mpopt = mpoption('verbose',0,'out.all',0);
res = runpf(mpc_v1, mpopt);
if ~res.success
    error('runpf failed for v1 profile.');
end

Pg = res.gen(:, PG);
Qg = res.gen(:, QG);

disp('--- Generator buses (order used) ---'); disp(gen_buses);
disp('--- VG magnitudes used (p.u.) ---');    disp(VGset.');

fprintf('\nPg for v1 (MW):\n'); fprintf('%.2f ', Pg); fprintf('\nTotal Pg: %.2f\n', sum(Pg));
fprintf('\nQg for v1 (MVAr):\n'); fprintf('%.2f ', Qg); fprintf('\n');

save('case14x2_v1_results.mat','mpc','YGGM','gen_buses','VGset','v1','gamma','Pg','Qg','res');
disp('Saved to case14x2_v1_results.mat');
