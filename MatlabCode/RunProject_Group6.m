%%%%%%%%%%%%%%%%%%%%%     SOLVENCY II Project      %%%%%%%%%%%%%%%%%%%%%
% Main code for the Insurance Project 2023/2024 - SII Project
% Group 6
%   -> Matteo Ceve
%   -> Taha Rafai
%   -> Marco Ronchetti
%   -> Nicolò Toia

% clear the work sapce
clc;
clear;
close all;

% fix the seed
rng(42); % The answer to life, the universe, and everything

%% Add path to the data folder

addpath('Data');

% start the timer
tic;
%% Data

today = '03-31-2024';   % set the initial date 31-Mar-2024

% Contract maturity
T = 50;

% initial invested capital and fund capital
C0 = 1e5; 
F0 = C0;    

% Assets

    equity_perc = 0.8;                  % equity percentage
    property_perc = 1 - equity_perc;    % property percentage
    
    equity = equity_perc * F0;          % initial equity value
    property = property_perc * F0;      % initial property value
    
    % null dividend yield both equity and property
    dividend = 0;
    
    % volatilities equity & property
    sigma_equity = 0.20;
    sigma_property = 0.10; 

% Laiabilities
    
    % benefits
        % penalties in case of lapse (€)
        penalties = 20;   
        % death payoff: max(C0, Fund_value)
    
    % others
        % regular deduction
        RD = 0.022;
        % commisions  to the distribution channels
        comm = 0.014;

%Time vector
time=(1:T)';

% age of the insured
age = 60;

% lapse: flat annual rates
l_t = 0.15*ones(1,T)';

% expanses (€)
expenses_0 = 50; 

% inflation: flat annual rate of 2%
inflation = 0.02;

% display the data
disp('Data uploaded:');
disp(' ');
disp('Today: 31-Mar-2024');
disp('Contract maturity: 50 years');
disp('Initial invested capital: €100,000');
disp('Equity percentage: 80%');
disp('Property percentage: 20%');
disp('Dividend yield: 0%');
disp('Equity volatility: 20%');
disp('Property volatility: 10%');
disp('Penalties in case of lapse: €20');
disp('Regular deduction: 2.2%');
disp('Commissions to the distribution channels: 1.4%');
disp('Age of the insured: 60 years');
disp('Lapse rate: 15%');
disp('Expenses (per year): €50 (+ inflation)');
disp('Inflation rate: 2%');
disp(' ');

%% Upload risk free rates from EIOPA_Term_Structures.xlsx and mortality probability from LifeTableMale2022.xlsx

% EIOPA risk free rates EU 31-Mar-20204 & ISTAT Life Tables 2022 (male)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNCOMMENT THIS SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % rates_base = xlsread ('EIOPA_Term_Structures.xlsx',3 , 'C11:C60');
% % 
% % % up scenario (without VA)
% % rates_up = xlsread ('EIOPA_Term_Structures.xlsx',5 , 'C11:C60');
% % 
% % % down scenario (without VA)
% % rates_down = xlsread ('EIOPA_Term_Structures.xlsx',6 , 'C11:C60');
% % 
% % % Upload from excel the mortality rates
% % qx = xlsread ("LifeTableMale2022.xlsx",'D63:D112')/1000; % From 60 y-o to 109 y-o
% % 
% % % display the updated data (both rates and mortality rates) in a table 
% % table = table(rates_base, rates_up, rates_down, qx, 'VariableNames', {'Base', 'Up', 'Down', 'Mortality'});
% % disp('Risk free rates and mortality rates uploaded:');
% % disp(' ');
% % disp(table);
% % disp(' ');

% save data above in a .mat file
% save('data.mat', 'rates_base', 'rates_up', 'rates_down', 'qx');
load('data.mat');

%% Dates

% check if today is a business date, where today is the settlement date
if ~isbusday(today)
    today = busdate(today);
end

% from the settlement date, find the dates of the next 50 years
dates = datetime(today, 'ConvertFrom', 'datenum') + calyears(1:T)';
% check and move dates to business dates
dates(~isbusday(dates,0)) = busdate(dates(~isbusday(dates,0)), "modifiedfollow", 0);
% convvert to datenum
dates = [datenum(today); datenum(dates)];

% Set the year convention
ACT_365 = 3;

% compute the deltas between t0 and t_i
delta_t0 = yearfrac(dates(1), dates, ACT_365);

% find the time difference for each time step
deltas = delta_t0(2:end) - delta_t0(1:end-1);

% Uncomment the following for deterministic computation
% or as an alternative to have integer values and same step

dt = (0:T)';
delta_t0 = dt;
deltas = delta_t0(2:end) - delta_t0(1:end-1);

% Plot rates curves and mortality rates in sub plots
figure;
subplot(1,2,1);
plot(dates(2:end) , rates_base, 'b', 'LineWidth', 1.5);
hold on;
plot(dates(2:end) , rates_up, 'r', 'LineWidth', 1.5);
plot(dates(2:end) , rates_down, 'g', 'LineWidth', 1.5);
title('Risk free rates curves');
% set boudaries in x-axis
xlim([dates(1) dates(end)]);
% set the x-axis in datestr onli each 10 years
set(gca, 'XTick', dates(1:10:end));
set(gca, 'XTickLabel', datestr(dates(1:10:end), 'yyyy'));
xlabel('Time');
ylabel('Risk free rates');
legend('Base scenario', 'Up scenario', 'Down scenario');
grid on;
hold off;

subplot(1,2,2);
plot(60:109, qx*1000, 'b', 'LineWidth', 1.5);
title('Mortality rates');
xlabel('Age');
ylabel('qx (per mille)');
grid on;

%% Forward rates & discounts

% rates shifted by 100bp
% rates_base = rates_base + 0.01;
% rates_up = rates_up + 0.01;
% rates_down = rates_down + 0.01;

% compute the zero rates
rates_base = log(1+rates_base);
rates_up = log(1+rates_up);
rates_down = log(1+rates_down);

% compute the discounts base scenario
discounts = exp(-[1; rates_base].*delta_t0);
fwd_discount = discounts./[1; discounts(1:end-1)];
% compute the forward rates base scenario
fwd_rates_base = -log(fwd_discount);

% compute the discounts up scenario
discounts_up = exp(-[1; rates_up].*delta_t0);
fwd_discount_up = discounts_up./[1; discounts_up(1:end-1)];
% compute the forward rates up scenario
fwd_rates_up = -log(fwd_discount_up);

% compute the discounts down scenario
discounts_down = exp(-[1; rates_down].*delta_t0);
fwd_discount_down = discounts_down./[1; discounts_down(1:end-1)];
% compute the forward rates down scenario
fwd_rates_down = -log(fwd_discount_down);

%% Martingale Test
% perform the martingale test to select an adequate number of simulations for MC

% vector of N simulations to test (number of elements <= 6)
N_sim = [10, 1e2, 1e3, 1e4 ,1e5, 1e6];

% loop over the number of simulations to test Equity
disp('The number of simulations that pass the Martingale test for equity are the following:');
disp(' ');

for i = 1:length(N_sim)
    % perform the martingale test for equity
    martingaleTest = mtgTest(equity ,delta_t0, deltas ,N_sim(i), rates_base, fwd_rates_base, dividend ,sigma_equity, 0, i, 2);  
    
    if max(martingaleTest) < 0.01*equity
        % display the results
        disp(['N = ', num2str(N_sim(i))]);
    end
end

% loop over the number of simulations to test Property
disp(' ');
disp('The number of simulations that pass the Martingale test for property are the following:');
disp(' ');

for i = 1:length(N_sim)
    % perform the martingale test for property
    martingaleTest = mtgTest(property ,delta_t0, deltas ,N_sim(i), rates_base, fwd_rates_base, dividend ,sigma_property, 0, i, 3);  
    
    if max(martingaleTest) < 0.01*property
        % display the results
        disp(['N = ', num2str(N_sim(i))]);
    end
end

%% Simulate asset and property in the three rates scenario
% simulate the fund value in the three scenarios via Geometric Brownian Motion and Monte Carlo simulation

N = 1e6 ; % number of MC simulations

% simulate GBM equity and property base scenario
equity_MC = MC_simulation(equity ,deltas ,N ,fwd_rates_base, dividend ,sigma_equity, RD);
property_MC = MC_simulation(property ,deltas ,N ,fwd_rates_base, dividend ,sigma_property, RD);

% simulate GBM equity and property up scenario
equity_MC_up = MC_simulation(equity ,deltas ,N ,fwd_rates_up, dividend ,sigma_equity, RD);
property_MC_up = MC_simulation(property ,deltas ,N ,fwd_rates_up, dividend ,sigma_property, RD);

% simulate GBM equity and property down scenario
equity_MC_down = MC_simulation(equity ,deltas ,N ,fwd_rates_down, dividend ,sigma_equity, RD);
property_MC_down = MC_simulation(property ,deltas ,N ,fwd_rates_down, dividend ,sigma_property, RD);

% Value of the fund for the three scenarios: base, up and down
Fund_base = equity_MC + property_MC;
Fund_up = equity_MC_up + property_MC_up;
Fund_down = equity_MC_down + property_MC_down;

% create a structure to store the results
Value_Fund = struct("Fund_t0", F0, "Fund_base", Fund_base, "Fund_up", Fund_up, "Fund_down", Fund_down);

%% Plot the fund value in the three scenarios and the initial value
% set on the x-axis the dates (in datestr) and on the y-axis the fund value

figure;
plot(dates , F0*ones(1,51), '-.k', 'LineWidth', 1);
hold on;
plot(dates , mean(Fund_base,1), 'b', 'LineWidth', 1.5);
plot(dates , mean(Fund_up,1), 'r', 'LineWidth', 1.5);
plot(dates , mean(Fund_down,1), 'g', 'LineWidth', 1.5);
title('Fund value in the three scenarios');
% set boudaries in x-axis
xlim([dates(1) dates(end)]);
% set the x-axis in datestr onli each 10 years
set(gca, 'XTick', dates(1:10:end));
set(gca, 'XTickLabel', datestr(dates(1:10:end), 'yyyy'));
xlabel('Time');
ylabel('Fund value');
% legend on the top left side
legend('Initial value', 'Base scenario', 'Up scenario', 'Down scenario', 'Location', 'NorthWest');
grid on;

%% DETERMINISTIC CHECK 
% uncomment to run the deterministic code as a check of robustness

% Fund_base = ones(1e6,51) * 1e5;
% Fund_up = ones(1e6,51)*1e5;
% Fund_down = ones(1e6,51)*1e5;
% 
% for i = 2:51
%     Fund_base(:,i) = Fund_base(:,i-1) * (1 - RD) / (fwd_discount(i));
%     Fund_up(:,i) = Fund_up(:,i-1) * (1 - RD) / (fwd_discount_up(i));
%     Fund_down(:,i) = Fund_down(:,i-1) * (1 - RD) / (fwd_discount_down(i));
% end
% 
% % display the results
% disp(' ');
% disp('The deterministic value of the fund at the end of the simulation is:');
% disp(['Base scenario: ', num2str(mean(Fund_base(:,end)))]);
% disp(['Up scenario: ', num2str(mean(Fund_up(:,end)))]);
% disp(['Down scenario: ', num2str(mean(Fund_down(:,end)))]);
% 
% Value_Fund = struct("Fund_t0", F0, "Fund_base", Fund_base, "Fund_up", Fund_up, "Fund_down", Fund_down);

%% Expenses

% Compute the vector of future expanses assumed a fixed infaltion rate
expenses = ones(T,1);
expenses(1,1) = expenses_0;

for i=2:T
    expenses(i,1) = expenses_0 * (1 + inflation)^(i-1);
end

%% Interest Rate Risk

% Base scenario
% Compute the liabilities, Macaulay duration and BEL components for the base scenario
[Liabilities_base, M_duration_base, Lapse_BEL_base, Death_BEL_base, Expenses_BEL_base, Commissions_BEL_base] = ...
    Liabilities(C0, Fund_base, discounts, delta_t0, l_t, qx, penalties, expenses, RD, comm);

% Compute Basic Own Funds (BOF) and the Delta BOF for the base scenario
BOF_base = F0 - Liabilities_base;
Delta_BOF_base = max(0, BOF_base - BOF_base);

% Up scenario
% Compute the liabilities, Macaulay duration and BEL components for the up scenario
[Liabilities_up, M_duration_up, Lapse_BEL_up, Death_BEL_up, Expenses_BEL_up, Commissions_BEL_up] = ...
    Liabilities(C0, Fund_up, discounts_up, delta_t0, l_t, qx, penalties, expenses, RD, comm);

% Compute Basic Own Funds (BOF) and the Delta BOF for the up scenario
BOF_up = F0 - Liabilities_up;
Delta_BOF_up = max(0, BOF_base - BOF_up);

% Down scenario
% Compute the liabilities, Macaulay duration and BEL components for the down scenario
[Liabilities_down, M_duration_down, Lapse_BEL_down, Death_BEL_down, Expenses_BEL_down, Commissions_BEL_down] = ...
    Liabilities(C0, Fund_down, discounts_down, delta_t0, l_t, qx, penalties, expenses, RD, comm);

% Compute Basic Own Funds (BOF) and the Delta BOF for the down scenario
BOF_down = F0 - Liabilities_down;
Delta_BOF_down = max(0, BOF_base - BOF_down);

% Compute the Solvency Capital Requirement (SCR) for the interest rate risk
SCR_IR = max(Delta_BOF_up, Delta_BOF_down);

%% Equity Risk

% shock equity by 39% plus 5.25% symmetry adjastment by Eiopa (consult data folder for more details)
equity_shock = 0.39;
symmetry_adj = 0.0525;
% equity stressed
E_stressed = equity*(1 - equity_shock - symmetry_adj); 
F0_stressed_E = E_stressed + property;

% simulate the equity value in the stressed scenario
equity_MC_stressed = MC_simulation(E_stressed ,deltas ,N ,fwd_rates_base, dividend ,sigma_equity, RD);

% Fund value in the equity stressed scenario
F_E_stressed = equity_MC_stressed + property_MC;

%%%%%%%%%% Uncomment for deterministic calculation %%%%%%%%%%%%%
% F_E_stressed = ones(1e6,51)*F0_stressed_E;
% for i = 2:51
%     F_E_stressed(:,i) = F_E_stressed(:,i-1) * (1 - RD) / (fwd_discount(i));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the liabilities, Macaulay duration and BEL components for the equity stressed scenario
[Liabilities_EQ, M_duration_EQ, Lapse_BEL_EQ, Death_BEL_EQ, Expenses_BEL_EQ, Commissions_BEL_EQ] = ...
        Liabilities(C0, F_E_stressed , discounts, delta_t0, l_t, qx, penalties, expenses, RD, comm);

% Compute Basic Own Funds (BOF) and the Delta BOF for the equity stressed scenario
BOF_EQ = F0_stressed_E - Liabilities_EQ;
Delta_BOF_EQ = max(0, BOF_base - BOF_EQ);

% Compute the Solvency Capital Requirement (SCR) for the equity risk
SCR_EQ = Delta_BOF_EQ;

%% Property Risk

% shock property by 25%
property_shock = 0.25;
% property stressed
P_stressed = property*(1 - property_shock);
F0_stressed_P = equity + P_stressed;

% simulate the property value in the stressed scenario
property_MC_stressed = MC_simulation(P_stressed ,deltas ,N ,fwd_rates_base, dividend ,sigma_property, RD);

% Fund value in the property stressed scenario
F_P_stressed = equity_MC + property_MC_stressed;

%%%%%%%%%% Uncomment for deterministic calculation %%%%%%%%%%%%%
% F_P_stressed = ones(1e6,51)*F0_stressed_P;
% for i = 2:51
%     F_P_stressed(:,i) = F_P_stressed(:,i-1) * (1 - RD) / (fwd_discount(i));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the liabilities, Macaulay duration and BEL components for the property stressed scenario
[Liabilities_PR, M_duration_PR, Lapse_BEL_PR, Death_BEL_PR, Expenses_BEL_PR, Commissions_BEL_PR] = ...
        Liabilities(C0, F_P_stressed , discounts, delta_t0, l_t, qx, penalties, expenses, RD, comm);

% Compute Basic Own Funds (BOF) and the Delta BOF for the property stressed scenario
BOF_PR = F0_stressed_P - Liabilities_PR;
Delta_BOF_PR = max(0, BOF_base - BOF_PR);

% Compute the Solvency Capital Requirement (SCR) for the property risk
SCR_PR = Delta_BOF_PR;

%% Mortality Risk

% Shock the mortality rate by 15% (check that the probability of death is not greater than 1)
qx_stressed = min(1,qx*(1+0.15));

% Compute the liabilities, Macaulay duration and BEL components for the mortality stressed scenario
[Liabilities_MT, M_duration_MT, Lapse_BEL_MT, Death_BEL_MT, Expenses_BEL_MT, Commissions_BEL_MT] = ...
        Liabilities(C0, Fund_base , discounts, delta_t0, l_t, qx_stressed, penalties, expenses, RD, comm);

% Compute Basic Own Funds (BOF) and the Delta BOF for the mortality stressed scenario
BOF_MT = F0 - Liabilities_MT;
Delta_BOF_MT = max(0, BOF_base - BOF_MT);

% Compute the Solvency Capital Requirement (SCR) for the mortality risk
SCR_MT = Delta_BOF_MT;

%% Lapse Risk 

% UP Scenario shok by 50% the probability of lapse (check that the probability of lapse is not greater than 1)
l_t_up = min(1.5*l_t(1),1)*ones(T,1);

% Compute the liabilities, Macaulay duration and BEL components for the lapse up scenario
[Liabilities_LTU, M_duration_LTU, Lapse_BEL_LTU, Death_BEL_LTU, Expenses_BEL_LTU, Commissions_BEL_LTU] = ...
        Liabilities(C0, Fund_base, discounts, delta_t0, l_t_up, qx, penalties, expenses, RD, comm);

% Compute Basic Own Funds (BOF) and the Delta BOF for the lapse up scenario
BOF_LTU = F0 - Liabilities_LTU;
Delta_BOF_LTU = max(0, BOF_base - BOF_LTU);

% DOWN Scenario shok by 50% the probability of lapse (check that the probability 
%       of lapse is not smaller than basic scenario minus 20%)            
l_t_dw = max(0.5*l_t(1),l_t(1)-0.2)*ones(T,1);

% Compute the liabilities, Macaulay duration and BEL components for the lapse down scenario
[Liabilities_LTD, M_duration_LTD, Lapse_BEL_LTD, Death_BEL_LTD, Expenses_BEL_LTD, Commissions_BEL_LTD] = ...
        Liabilities(C0, Fund_base, discounts, delta_t0, l_t_dw, qx, penalties, expenses, RD, comm);


% Compute Basic Own Funds (BOF) and the Delta BOF for the lapse down scenario
BOF_LTD = F0 - Liabilities_LTD;
Delta_BOF_LTD = max(0, BOF_base - BOF_LTD);

% Lapse MASS Risk shock the first year by 40%
l_t_mass = l_t;
l_t_mass(1) = (0.4 + l_t(1));

% Compute the liabilities, Macaulay duration and BEL components for the lapse mass scenario
[Liabilities_LTM, M_duration_LTM, Lapse_BEL_LTM, Death_BEL_LTM, Expenses_BEL_LTM, Commissions_BEL_LTM] = ...
        Liabilities(C0, Fund_base, discounts, delta_t0, l_t_mass, qx, penalties, expenses, RD, comm);

% Compute Basic Own Funds (BOF) and the Delta BOF for the lapse mass scenario
BOF_LTM = F0 - Liabilities_LTM;
Delta_BOF_LTM = max(0, BOF_base - BOF_LTM);

% Compute the Solvency Capital Requirement (SCR) for the lapse risk
SCR_LAP = max(max(Delta_BOF_LTM, Delta_BOF_LTD), Delta_BOF_LTU);

%% Catastrophe Risk

% Shock the mortality rate of the first year by 15 bps
qx_CAT = qx;
qx_CAT(1) = qx_CAT(1) + 0.0015;

% Compute the liabilities, Macaulay duration and BEL components for the catastrophe scenario
[Liabilities_CAT, M_duration_CAT, Lapse_BEL_CAT, Death_BEL_CAT, Expenses_BEL_CAT, Commissions_BEL_CAT] = ...
        Liabilities(C0, Fund_base , discounts, delta_t0, l_t, qx_CAT, penalties, expenses, RD, comm);

% Compute Basic Own Funds (BOF) and the Delta BOF for the catastrophe scenario
BOF_CAT = F0 - Liabilities_CAT;
Delta_BOF_CAT = max(0, BOF_base - BOF_CAT);

% Compute the Solvency Capital Requirement (SCR) for the catastrophe risk
SCR_CAT = Delta_BOF_CAT;

%% Expense Risk

% Shock the expenses by 10% and the inflation rate by 1%
expenses_stressed = zeros(T,1);
expenses_stressed(1) = expenses_0 * (1 + 0.1);

for i = 2:T
    expenses_stressed(i) = expenses_stressed(i-1) * (1+ inflation + 0.01);
end

% Compute the liabilities, Macaulay duration and BEL components for the expense stressed scenario
[Liabilities_EX, M_duration_EX, Lapse_BEL_EX, Death_BEL_EX, Expenses_BEL_EX, Commissions_BEL_EX] = ...
        Liabilities(C0, Fund_base , discounts, delta_t0, l_t, qx, penalties, expenses_stressed, RD, comm);

% Compute Basic Own Funds (BOF) and the Delta BOF for the expense stressed scenario
BOF_EX = F0 - Liabilities_EX;
Delta_BOF_EX = max(0, BOF_base - BOF_EX);

% Compute the Solvency Capital Requirement (SCR) for the expense risk
SCR_EX = Delta_BOF_EX;

%% Solvency Capital Requirement (SCR) and Basic Solvency Capital Requirement (BSCR)

% Aggregate the risks in:
% - Market Risk: Interest Rate Risk, Equity Risk, Property Risk
% - Life Risk: Mortality Risk, Lapse Risk, Catastrophe Risk
% build the SCR vectors

% Market Risk
SCR_market_Risk = [SCR_IR, SCR_EQ, SCR_PR]';
% Life Risk
SCR_life_Risk = [SCR_MT, SCR_LAP, SCR_EX, SCR_CAT]';

% Find the adeguate correlation matrices for the risks

% Choose the right market risk correlatio matrix based on 
%       the risk factors we are exposed to (eg. IR_up or IR_down)
if Delta_BOF_down > Delta_BOF_up    % Down scenario
    mkt_corr_matrix = [1 0.5 0.5;
                       0.5 1 0.75;
                       0.5 0.75 1];
else     % Up scenario
    mkt_corr_matrix = [1 0  0;
                       0 1 0.75;
                       0 0.75 1];
end

% Life correlation matrix
Life_corr_matrix = [1 0 0.25 0.25;
                    0 1 0.5 0.25;
                    0.25 0.5 1 0.25;
                    0.25 0.25 0.25 1];

% BSCR correlation matrix
BSCR_corr_matrix = [1 0.25;
                    0.25 1];

% Calculate the  aggregate SCR

% Market SCR
SCR_Market = sqrt(SCR_market_Risk' * mkt_corr_matrix * SCR_market_Risk);

% Life SCR
SCR_Life = sqrt(SCR_life_Risk' * Life_corr_matrix * SCR_life_Risk);

% build SCR vector
SCR = [SCR_Market, SCR_Life]';

% Calculate the Basic Solvency Capital Requirement (BSCR)
BSCR = sqrt(SCR' * BSCR_corr_matrix * SCR);

%% Save all the results in a structure

% First layer of the struct is DataStructResults
% insiede DataStructResults we have the following fields:
% - Rates scenarios
% - Equity and property shock
% - Mortality shock
% - Lapse shock
% - Catastrophe shock
% - Expense shock
% insiede each field we have the following sub_fields:
% - Results: (which contains the following fields)
            % - Assets (Fund value in t0)
            % - Liabilities
            % - BOF
            % - Delta BOF
            % - Macaulay Duration         - 
% - BEL components (Lapse, Death, Expenses, Commissions)

DataStructResults = struct("Rates", struct("Base", struct(), "Up", struct(), "Down", struct()), ...
                 "Equity_Property", struct("Equity", struct(), "Property", struct()), ...
                 "Mortality", struct(), ...
                 "Lapse", struct("Up", struct(), "Down", struct(), "Mass", struct()), ...
                 "Catastrophe", struct(), ...
                 "Expenses", struct());

% Rates scenarios
DataStructResults.Rates.Base.Results = struct("Assets", F0, "Liabilities", Liabilities_base, "BOF", BOF_base, "Delta_BOF", ...
    Delta_BOF_base, "Macaulay_Duration", M_duration_base);

DataStructResults.Rates.Base.BEL = struct("Lapse", Lapse_BEL_base, "Death", Death_BEL_base, "Expenses", Expenses_BEL_base, ...
    "Commissions", Commissions_BEL_base);

DataStructResults.Rates.Up.Results = struct("Assets", F0, "Liabilities", Liabilities_up, "BOF", BOF_up, "Delta_BOF", ...
    Delta_BOF_up, "Macaulay_Duration", M_duration_up);

DataStructResults.Rates.Up.BEL = struct("Lapse", Lapse_BEL_up, "Death", Death_BEL_up, "Expenses", Expenses_BEL_up, ...
    "Commissions", Commissions_BEL_up);

DataStructResults.Rates.Down.Results = struct("Assets", F0, "Liabilities", Liabilities_down, "BOF", BOF_down, "Delta_BOF", ...
    Delta_BOF_down, "Macaulay_Duration", M_duration_down);

DataStructResults.Rates.Down.BEL = struct("Lapse", Lapse_BEL_down, "Death", Death_BEL_down, "Expenses", Expenses_BEL_down, ...
    "Commissions", Commissions_BEL_down);

% Equity and property shock
DataStructResults.Equity_Property.Equity.Results = struct("Assets", F0_stressed_E, "Liabilities", Liabilities_EQ, "BOF", BOF_EQ, "Delta_BOF", ...
    Delta_BOF_EQ, "Macaulay_Duration", M_duration_EQ);

DataStructResults.Equity_Property.Equity.BEL = struct("Lapse", Lapse_BEL_EQ, "Death", Death_BEL_EQ, "Expenses", Expenses_BEL_EQ, ...
    "Commissions", Commissions_BEL_EQ);

DataStructResults.Equity_Property.Property.Results = struct("Assets", F0_stressed_P, "Liabilities", Liabilities_PR, "BOF", BOF_PR, "Delta_BOF", ...
    Delta_BOF_PR, "Macaulay_Duration", M_duration_PR);

DataStructResults.Equity_Property.Property.BEL = struct("Lapse", Lapse_BEL_PR, "Death", Death_BEL_PR, "Expenses", Expenses_BEL_PR, ...
    "Commissions", Commissions_BEL_PR);

% Mortality shock
DataStructResults.Mortality.Results = struct("Assets", F0, "Liabilities", Liabilities_MT, "BOF", BOF_MT, "Delta_BOF", ...
    Delta_BOF_MT, "Macaulay_Duration", M_duration_MT);

DataStructResults.Mortality.BEL = struct("Lapse", Lapse_BEL_MT, "Death", Death_BEL_MT, "Expenses", Expenses_BEL_MT, ...
    "Commissions", Commissions_BEL_MT);

% Lapse shock
DataStructResults.Lapse.Up.Results = struct("Assets", F0, "Liabilities", Liabilities_LTU, "BOF", BOF_LTU, "Delta_BOF", ...
    Delta_BOF_LTU, "Macaulay_Duration", M_duration_LTU);

DataStructResults.Lapse.Up.BEL = struct("Lapse", Lapse_BEL_LTU, "Death", Death_BEL_LTU, "Expenses", Expenses_BEL_LTU, ...
    "Commissions", Commissions_BEL_LTU);

DataStructResults.Lapse.Down.Results = struct("Assets", F0, "Liabilities", Liabilities_LTD, "BOF", BOF_LTD, "Delta_BOF", ...
    Delta_BOF_LTD, "Macaulay_Duration", M_duration_LTD);

DataStructResults.Lapse.Down.BEL = struct("Lapse", Lapse_BEL_LTD, "Death", Death_BEL_LTD, "Expenses", Expenses_BEL_LTD, ...
    "Commissions", Commissions_BEL_LTD);

DataStructResults.Lapse.Mass.Results = struct("Assets", F0, "Liabilities", Liabilities_LTM, "BOF", BOF_LTM, "Delta_BOF", ...
    Delta_BOF_LTM, "Macaulay_Duration", M_duration_LTM);

DataStructResults.Lapse.Mass.BEL = struct("Lapse", Lapse_BEL_LTM, "Death", Death_BEL_LTM, "Expenses", Expenses_BEL_LTM, ...
    "Commissions", Commissions_BEL_LTM);

% Catastrophe shock
DataStructResults.Catastrophe.Results = struct("Assets", F0, "Liabilities", Liabilities_CAT, "BOF", BOF_CAT, "Delta_BOF", ...
    Delta_BOF_CAT, "Macaulay_Duration", M_duration_CAT);

DataStructResults.Catastrophe.BEL = struct("Lapse", Lapse_BEL_CAT, "Death", Death_BEL_CAT, "Expenses", Expenses_BEL_CAT, ...
    "Commissions", Commissions_BEL_CAT);

% Expense shock
DataStructResults.Expenses.Results = struct("Assets", F0, "Liabilities", Liabilities_EX, "BOF", BOF_EX, "Delta_BOF", ...
    Delta_BOF_EX, "Macaulay_Duration", M_duration_EX);

DataStructResults.Expenses.BEL = struct("Lapse", Lapse_BEL_EX, "Death", Death_BEL_EX, "Expenses", Expenses_BEL_EX, ...
    "Commissions", Commissions_BEL_EX);


% Build a structure to store the results of SCR and BSCR
% First layer of the struct is SCR_BSCR
% insiede SCR_BSCR we have the following fields:
% - SCR_MKT_Risk (% insiede each field we have the following sub_fields: interest rate, equity, property)
% - SCR_Life_Risk (% insiede each field we have the following sub_fields: mortality, lapse, catastrophe, expenses)
% - SCR (% insiede each field we have the following sub_fields: market and life)
% - BSCR

SCR_BSCR = struct("SCR_MKT_Risk", struct("Interest", struct(), "Equity", struct(), "Property", struct()), ...
                  "SCR_Life_Risk", struct("Mortality", struct(), "Lapse", struct(), "Catastrophe", struct(), "Expenses", struct()), ...
                  "SCR", struct("Market", struct(), "Life", struct()), ...
                  "BSCR", struct());

SCR_BSCR.SCR_MKT_Risk.Interest = SCR_IR;
SCR_BSCR.SCR_MKT_Risk.Equity = SCR_EQ;
SCR_BSCR.SCR_MKT_Risk.Property = SCR_PR;

SCR_BSCR.SCR_Life_Risk.Mortality = SCR_MT;
SCR_BSCR.SCR_Life_Risk.Lapse = SCR_LAP;
SCR_BSCR.SCR_Life_Risk.Catastrophe = SCR_CAT;
SCR_BSCR.SCR_Life_Risk.Expenses = SCR_EX;

SCR_BSCR.SCR.Market = SCR_Market;
SCR_BSCR.SCR.Life = SCR_Life;
SCR_BSCR.BSCR = BSCR;


%% Display all the results

disp(' ')
disp('  --------------------------------------------------- SII Project Results --------------------------------------------------------------------')
disp(' ')

% Extracting .Results for each field
ratesBaseResults = DataStructResults.Rates.Base.Results;
ratesUpResults = DataStructResults.Rates.Up.Results;
ratesDownResults = DataStructResults.Rates.Down.Results;
equityResults = DataStructResults.Equity_Property.Equity.Results;
propertyResults = DataStructResults.Equity_Property.Property.Results;
mortalityResults = DataStructResults.Mortality.Results;
lapseUpResults = DataStructResults.Lapse.Up.Results;
lapseDownResults = DataStructResults.Lapse.Down.Results;
lapseMassResults = DataStructResults.Lapse.Mass.Results;
catastropheResults = DataStructResults.Catastrophe.Results;
expensesResults = DataStructResults.Expenses.Results;

% Formatting the numbers
formatNumber = @(x) sprintf('%.2f', x);

% Creating a table
tableResults = table({formatNumber(ratesBaseResults.Assets); formatNumber(ratesUpResults.Assets); formatNumber(ratesDownResults.Assets); formatNumber(equityResults.Assets); ...
                      formatNumber(propertyResults.Assets); formatNumber(mortalityResults.Assets); formatNumber(lapseUpResults.Assets); formatNumber(lapseDownResults.Assets); ...
                      formatNumber(lapseMassResults.Assets); formatNumber(catastropheResults.Assets); formatNumber(expensesResults.Assets)}, ...
                     {formatNumber(ratesBaseResults.Liabilities); formatNumber(ratesUpResults.Liabilities); formatNumber(ratesDownResults.Liabilities); formatNumber(equityResults.Liabilities); ...
                      formatNumber(propertyResults.Liabilities); formatNumber(mortalityResults.Liabilities); formatNumber(lapseUpResults.Liabilities); formatNumber(lapseDownResults.Liabilities); ...
                      formatNumber(lapseMassResults.Liabilities); formatNumber(catastropheResults.Liabilities); formatNumber(expensesResults.Liabilities)}, ...
                     {formatNumber(ratesBaseResults.BOF); formatNumber(ratesUpResults.BOF); formatNumber(ratesDownResults.BOF); formatNumber(equityResults.BOF); ...
                      formatNumber(propertyResults.BOF); formatNumber(mortalityResults.BOF); formatNumber(lapseUpResults.BOF); formatNumber(lapseDownResults.BOF); ...
                      formatNumber(lapseMassResults.BOF); formatNumber(catastropheResults.BOF); formatNumber(expensesResults.BOF)}, ...
                     {formatNumber(ratesBaseResults.Delta_BOF); formatNumber(ratesUpResults.Delta_BOF); formatNumber(ratesDownResults.Delta_BOF); formatNumber(equityResults.Delta_BOF); ...
                      formatNumber(propertyResults.Delta_BOF); formatNumber(mortalityResults.Delta_BOF); formatNumber(lapseUpResults.Delta_BOF); formatNumber(lapseDownResults.Delta_BOF); ...
                      formatNumber(lapseMassResults.Delta_BOF); formatNumber(catastropheResults.Delta_BOF); formatNumber(expensesResults.Delta_BOF)}, ...
                     {formatNumber(ratesBaseResults.Macaulay_Duration); formatNumber(ratesUpResults.Macaulay_Duration); formatNumber(ratesDownResults.Macaulay_Duration); ...
                      formatNumber(equityResults.Macaulay_Duration); formatNumber(propertyResults.Macaulay_Duration); formatNumber(mortalityResults.Macaulay_Duration); ...
                      formatNumber(lapseUpResults.Macaulay_Duration); formatNumber(lapseDownResults.Macaulay_Duration); formatNumber(lapseMassResults.Macaulay_Duration); ...
                      formatNumber(catastropheResults.Macaulay_Duration); formatNumber(expensesResults.Macaulay_Duration)}, ...
                     'VariableNames', {'Assets', 'Liabilities', 'BOF', 'Delta_BOF', 'Macaulay_Duration'}, ...
                     'RowNames', {'Rates_Base', 'Rates_Up', 'Rates_Down', 'Equity', 'Property', 'Mortality', 'Lapse_Up', ...
                                  'Lapse_Down', 'Lapse_Mass', 'Catastrophe', 'Expenses'});

% Displaying the table
disp(tableResults);
disp(' ')

%% Disp BEL Components

% Extracting .BEL for each field
ratesBaseBEL = DataStructResults.Rates.Base.BEL;
ratesUpBEL = DataStructResults.Rates.Up.BEL;
ratesDownBEL = DataStructResults.Rates.Down.BEL;
equityBEL = DataStructResults.Equity_Property.Equity.BEL;
propertyBEL = DataStructResults.Equity_Property.Property.BEL;
mortalityBEL = DataStructResults.Mortality.BEL;
lapseUpBEL = DataStructResults.Lapse.Up.BEL;
lapseDownBEL = DataStructResults.Lapse.Down.BEL;
lapseMassBEL = DataStructResults.Lapse.Mass.BEL;
catastropheBEL = DataStructResults.Catastrophe.BEL;
expensesBEL = DataStructResults.Expenses.BEL;

% Creating a table for .BEL
tableBEL = table({formatNumber(ratesBaseBEL.Lapse); formatNumber(ratesUpBEL.Lapse); formatNumber(ratesDownBEL.Lapse); formatNumber(equityBEL.Lapse); ...
                  formatNumber(propertyBEL.Lapse); formatNumber(mortalityBEL.Lapse); formatNumber(lapseUpBEL.Lapse); formatNumber(lapseDownBEL.Lapse); ...
                  formatNumber(lapseMassBEL.Lapse); formatNumber(catastropheBEL.Lapse); formatNumber(expensesBEL.Lapse)}, ...
                 {formatNumber(ratesBaseBEL.Death); formatNumber(ratesUpBEL.Death); formatNumber(ratesDownBEL.Death); formatNumber(equityBEL.Death); ...
                  formatNumber(propertyBEL.Death); formatNumber(mortalityBEL.Death); formatNumber(lapseUpBEL.Death); formatNumber(lapseDownBEL.Death); ...
                  formatNumber(lapseMassBEL.Death); formatNumber(catastropheBEL.Death); formatNumber(expensesBEL.Death)}, ...
                 {formatNumber(ratesBaseBEL.Expenses); formatNumber(ratesUpBEL.Expenses); formatNumber(ratesDownBEL.Expenses); formatNumber(equityBEL.Expenses); ...
                  formatNumber(propertyBEL.Expenses); formatNumber(mortalityBEL.Expenses); formatNumber(lapseUpBEL.Expenses); formatNumber(lapseDownBEL.Expenses); ...
                  formatNumber(lapseMassBEL.Expenses); formatNumber(catastropheBEL.Expenses); formatNumber(expensesBEL.Expenses)}, ...
                 {formatNumber(ratesBaseBEL.Commissions); formatNumber(ratesUpBEL.Commissions); formatNumber(ratesDownBEL.Commissions); formatNumber(equityBEL.Commissions); ...
                  formatNumber(propertyBEL.Commissions); formatNumber(mortalityBEL.Commissions); formatNumber(lapseUpBEL.Commissions); formatNumber(lapseDownBEL.Commissions); ...
                  formatNumber(lapseMassBEL.Commissions); formatNumber(catastropheBEL.Commissions); formatNumber(expensesBEL.Commissions)}, ...
                 'VariableNames', {'Lapse', 'Death', 'Expenses', 'Commissions'}, ...
                 'RowNames', {'Rates_Base', 'Rates_Up', 'Rates_Down', 'Equity', 'Property', 'Mortality', 'Lapse_Up', ...
                              'Lapse_Down', 'Lapse_Mass', 'Catastrophe', 'Expenses'});

% Displaying the table for .BEL
disp(tableBEL);

%% Displaying the Solvency Capital Requirement (SCR) and Basic Solvency Capital Requirement (BSCR)

disp(' ')
disp('  --------------------------------------------------- SCR & BSCR --------------------------------------------------------------------')
disp(' ')

% Display table for SCR_MKT_Risk
table_SCR_MKT_Risk = table(SCR_BSCR.SCR_MKT_Risk.Interest, SCR_BSCR.SCR_MKT_Risk.Equity, SCR_BSCR.SCR_MKT_Risk.Property, ...
    'VariableNames', {'Interest', 'Equity', 'Property'});

disp("SCR Market Risk:");
disp(table_SCR_MKT_Risk);

% Display table for SCR_Life_Risk
table_SCR_Life_Risk = table(SCR_BSCR.SCR_Life_Risk.Mortality, SCR_BSCR.SCR_Life_Risk.Lapse, SCR_BSCR.SCR_Life_Risk.Catastrophe, SCR_BSCR.SCR_Life_Risk.Expenses, ...
    'VariableNames', {'Mortality', 'Lapse', 'Catastrophe', 'Expenses'});

disp("SCR Life Risk:");
disp(table_SCR_Life_Risk);

% Display table for SCR (Market and Life)
table_SCR = table(SCR_BSCR.SCR.Market, SCR_BSCR.SCR.Life, ...
    'VariableNames', {'Market', 'Life'});

disp("SCR (Market and Life):");
disp(table_SCR);

% Display BSCR
disp("BSCR: " + formatNumber(SCR_BSCR.BSCR));
disp(' ')

toc; % stop the timer