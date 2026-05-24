clear
close all
clc

%% Raw data import
% Load historical Day-Ahead Market (DAM) and Real-Time Market (RTM) spread data
% data source NYTI
Input  = readtable('MATLAB_DRO_Inputs.xlsx','Sheet','spread data','VariableNamingRule','preserve');
Test   = readtable('NYISO_Market_Data_2026Jan.xlsx','Sheet','DAM_RTM_spread','VariableNamingRule','preserve');
tg_day = unique(Input.Target_Date);
tic
%% Parameter definition
zone   = 11;      % N: Number of geographic zones (e.g., NYISO zones)
interv = 24;      % T: 24 hours scheduling horizon
K      = 2;       % Index set for CVaR piecewise linear approximation
eps    = 1:2.5:50;% Wasserstein ball radius (epsilon) defining the ambiguity set
D      = 7;       % Number of historical data samples in the empirical distribution D_j
Lamb   = 3500;    % Lambda: Maximum bound of the spread uncertainty set (Box constraint)
rho    = 0:0.1:1; % Risk-aversion factor balancing expected loss and CVaR
alpha  = 0.1;     % Confidence level for CVaR (alpha = 10%)
L      = 400;     % Maximum bidding volume limit (MWh) across all buses

% Coefficients for the CVaR reformulation
a_1    = -rho;
a_2    = -rho-(1-rho)/alpha;
b_1    = 1-rho;
b_2    = (1-rho)*(1-1/alpha);
a_k    = [a_1;a_2];
b_k    = [b_1;b_2];

%%
% decision variable
q_t    = sdpvar(interv,zone,'full');        % Bidding quantity (MWh) at each node and time
tau    = sdpvar(1,1);                       % CVaR threshold variable
lamb   = sdpvar(1,1);                       % Dual variable for Wasserstein distance (lambda)
x_d    = sdpvar(D,1,'full');                % Auxiliary variable for worst-case expectation (x^d)
v_td   = sdpvar(zone,interv,D,K,'full');    % Auxiliary variable for uncertainty set mapping (v_{t,k}^d)
mu_kd  = sdpvar(D,K,'full');                % Auxiliary variable for box uncertainty bound (mu_k^d)

% Pre-extract all spread data
num_days = length(tg_day);
s_td_inputs = cell(num_days, 1);
s_td_tests  = cell(num_days, 1);
for i = 1:num_days
    raw_input = Input(Input.Target_Date == tg_day(i),:);
    s_td_inputs{i} = reshape(raw_input.Spread, zone, interv, D);
    
    raw_test = Test(dateshift(Test.Time, 'start', 'day') == tg_day(i),:);
    s_td_tests{i} = reshape(raw_test.Spread, zone, interv);
end

% Gurobi Solver Configuration
ops = sdpsettings('solver','gurobi','verbose', 0);
ops.gurobi.Method    = 2;                 % Use Barrier method for QP
ops.gurobi.Presolve  = 2;
ops.gurobi.Crossover = 0;
ops.gurobi.Threads   = 2;
ops.gurobi.BarConvTol = 1e-4;
ops.gurobi.FeasibilityTol = 1e-5;
ops.gurobi.OptimalityTol = 1e-5;

q_t_sens  = cell(length(rho), length(eps));
prof_sens = cell(length(rho), length(eps));

%% Main Optimization Loop
for r = 1:length(rho)
    fprintf('Building & Solving for rho = %.2f (%d/%d)...\n', rho(r), r, length(rho));
    a_k_tem = a_k(:,r);
    b_k_tem = b_k(:,r);
    
    yalmip('clear'); 
    q_t    = sdpvar(interv, zone, 'full');
    tau    = sdpvar(1, 1);
    lamb   = sdpvar(1, 1);
    x_d    = sdpvar(D, 1, 'full');
    v_td   = sdpvar(zone, interv, D, K, 'full');
    mu_kd  = sdpvar(D, K, 'full');
    
    s_td_p = sdpvar(zone, interv, D, 'full'); 
    eps_p  = sdpvar(1, 1);
    
    % Construct constraints based on the DRO-CVaR model (Eq. 13)
    q_tcell = num2cell(q_t, 2);

    % Constraint (13e): L1-norm bidding limit constraint
    consts_13e = cellfun(@(x) norm(x, 1), q_tcell, 'UniformOutput', false);
    consts = [consts_13e{:}] <= L;
    
    consts_inner = cell(D, K);
    for k = 1:K
        % Constraint (13b): Epigraph formulation combining CVaR and Wasserstein ambiguity
        consts = [consts; 
                  b_k_tem(k)*tau + a_k_tem(k)*squeeze(sum(s_td_p.*repmat(q_t.',[1 1 D]),[1 2])) + ...
                  squeeze(sum(v_td(:,:,:,k).*s_td_p,[1 2])) + Lamb*mu_kd(:,k) <= x_d];
        for d = 1:D
            v_td_temp  = num2cell(v_td(:,:,d,k), 1).';
            % Constraint (13c): Dual norm constraint derived from Wasserstein metric
            consts_13c = cellfun(@(x,y) norm(x.'-a_k_tem(k)*y, 2), v_td_temp, q_tcell, 'UniformOutput', false);
            % Constraint (13d): Auxiliary constraint bounding the Box Uncertainty Set
            consts_13d = cellfun(@(x) norm(x, 1), v_td_temp, 'UniformOutput', false);
            consts_inner{d,k} = [[consts_13d{:}] <= mu_kd(d,k);
                                 [consts_13c{:}] <= lamb];  
        end
    end
    consts = [consts; consts_inner{:}]; 
    
    % Objective (13a): Minimize worst-case expected loss + Wasserstein penalty
    obj = lamb * eps_p + sum(x_d) / D;
    
    % Compile the Optimizer Object (Inputs: Spread data & Epsilon; Output: Optimal Bids)
    DRO_solver = optimizer(consts, obj, ops, {s_td_p, eps_p}, q_t);
    
    % parfor may be able to spped up finishing the task 
    for e = 1:length(eps)
        eps_tem = eps(e);
        q_t_zone  = zeros(interv, zone, num_days);
        prof_zone = zeros(zone, num_days);
        
        for i = 1:num_days
            % Solve for the optimal bidding strategy (q_t) using historical training dataset
            q_t_sol = DRO_solver({s_td_inputs{i}, eps_tem});
            
            % Clean up numerical noise
            q_t_sol(abs(q_t_sol) <= 1e-3) = 0;
            q_t_sol = round(q_t_sol, 3);
            
            % Evaluate real-world profitability on the test dataset
            q_t_zone(:,:,i) = q_t_sol;
            prof_zone(:,i)  = sum(q_t_sol .* s_td_tests{i}.', 1).';
        end
        q_t_sens{r,e}  = q_t_zone;
        prof_sens{r,e} = prof_zone;
    end
end

%%
Calmar_matrix = zeros(length(rho), length(eps));
Return_matrix = zeros(length(rho), length(eps));
for r = 1:length(rho)

     for e = 1:length(eps)
     % Equity Curve
     daily_profits = sum(prof_sens{r,e});
     equity_curve  = cumsum(daily_profits);
     
     % everyday's maxmimum reward since the first day 
     running_max = cummax(equity_curve);

     % drawdown
     drawdowns = running_max - equity_curve;

     % maximum dawdown 
     MDD_abs      = max(drawdowns);
     max_drawdown = MDD_abs/max(running_max);
     
     % average reward
     R_abs = mean(daily_profits) * 365;

     % Calmar ratio
     if MDD_abs == 0
        C = 0; 
        else
        C = R_abs / MDD_abs;
     end
     Calmar_matrix(r, e) = C;
     Return_matrix(r, e) = sum(daily_profits);
     end
end
toc
%for e = 1:length(eps)
%    tem      = cell2num(q_t_sens(:,e));
%    sparsity = sum(tem>0);
%end 

%% visualize
% Calmar Ratio sensitivity analysis
[X,Y] = meshgrid(rho,eps);
figure
surface(X,Y,Calmar_matrix.')
view(125, 35) 
colormap(flip(othercolor('YlGnBu9')))
shading flat
shading interp
%[~, linear_index] = max(Calmar_matrix(:)); 
%[optm_rho, optm_eps] = ind2sub(size(Calmar_matrix), linear_index);
grid on;
box on;
xlim([min(rho), max(rho)]);
ylim([min(eps), max(eps)]);

hx = xlabel('Risk Aversion Level $\rho$', 'Interpreter', 'latex', 'FontSize', 12);
hy = ylabel('Wasserstein Radius $\epsilon$', 'Interpreter', 'latex', 'FontSize', 12);
hz = zlabel('Calmar Ratio', 'Interpreter', 'latex', 'FontSize', 12);
title('Risk-Adjusted Performance Surface (Calmar Ratio)');
hx.Rotation = 33;  
hy.Rotation = -16;
hz.Rotation = 90;  
hx.Position(2) = max(eps)-0.08*(max(eps)- min(eps));
hx.Position(1) = min(rho)+0.3*(max(rho)- min(rho));
hy.Position(1) = max(rho)-0.1*(max(rho)- min(rho));
hy.Position(2) = min(eps)+0.3*(max(eps)- min(eps));

ch = colorbar;
caxis([0 500])
ylabel(ch, 'Calmar Ratio', 'Interpreter', 'latex', 'FontSize', 12);

% Bidding Volume Distribution sensitivity
figure('Position', [100, 50, 700, 500]);
rho_tem = 5;
node_total_volume = zeros(length(eps),zone);
for e = 1:length(eps)
    node_total_volume(e,:) = sum(abs(q_t_sens{rho_tem,e}), [1 3]);    
end
imagesc(eps, 1:zone, node_total_volume.');
set(gca, 'YDir', 'normal');
colormap(flip(othercolor('YlGnBu9')))
xlim([min(eps) max(eps)])
ylim([1, zone]);

xlabel('Wasserstein ball radius $\epsilon$', 'Interpreter', 'latex');
ylabel('Zone index $n$', 'Interpreter', 'latex');
title('Bidding Volume Distribution under Different Wasserstein Radii');

ch = colorbar;
caxis([0 7e4])
ylabel(ch, 'Total Bidding Volume (MW)', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');