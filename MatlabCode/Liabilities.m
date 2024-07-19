function [Liabilities, M_duration, Lapse_BEL, Death_BEL, Expenses_BEL, Commissions_BEL] = ...
            Liabilities(C0, F, discounts, time, lt, qx, penalties, expenses, RD, COMM)
% This function calculates the liabilities of a portfolio of contracts, 
% moreover it computes the Macaulay Duration of liabilities
% and the Best Estimate Liabilities components (BEL) such as Lapse Benefits,
% Death Benefits, Expenses and Commissions
%
% INPUTS:
% C0                    : initial fund value
% F                     : matrix of fund values
% discounts             : vector of discounts
% time                  : vector of time
% lt                    : vector of lapse rates
% qx                    : vector of mortality rates
% penalties            : penalties for early lapse
% expenses              : vector of expenses
% RD                    : regular deduction
% COMM                  : commission
%
% OUTPUTS:
% Liabilities           : total liabilities
% M_duration            : Macaulay Duration of liabilities
% Lapse_BEL             : Lapse BEL component
% Death_BEL             : Death BEL component
% Expenses_BEL          : Expenses BEL component
% Commissions_BEL       : Commissions BEL component
 
% Maturity of the contracts
T = round(time(end));

% compute the probability of remaining in the contract at each year
% cumprod function compute the cumuluative product of the elements of the vector
% we initialize the vector with 1 because the first year the probability of
% remaining in the contract is 1
P_remain_contract = cumprod([1; (1-qx(1:end-1)).*(1-lt(1:end-1))]); 

% initialize the vectors
V = zeros(T,1);

% initialize the vectors of BEL components
Lapse_ben=zeros(T,1);
Death_ben=zeros(T,1);
Expenses_ben=zeros(T,1);
Commissions_ben=zeros(T,1);

% loop over the years
for i = 1 : T

    % cash flows at each year for the benefits

    % death benefits, use the full matrix F since max is a non linear operator
    death_cf = (max(C0, F(:,i+1)))*qx(i);

    % after taking the max we can compute the mean for semplicity
    death_cf = mean(death_cf);  % mean of the death benefits
    F_mean = mean(F(:,i+1));    % mean of the fund value


    % Lapse benefits at each year
    if i == T   % condition because in the last year it's mandatory to lapse 
        % the contract (contract is expired)
        lapse_cf = (F_mean)*(1-qx(i));
    else
        lapse_cf = (F_mean-penalties)*lt(i)*(1-qx(i));
    end

    % commission  cash flow at each year
    commision = (F_mean / (1-RD)) * COMM;

    % Total Liabilities at each year pre discounting (cash flows at each year)
    V(i) = P_remain_contract(i) * ( lapse_cf + death_cf + expenses(i) + commision );  

    % Benefits at each year
    Lapse_ben(i)        = (P_remain_contract(i) * lapse_cf);
    Death_ben(i)        = (P_remain_contract(i) * death_cf);
    Expenses_ben(i)     = (P_remain_contract(i) * expenses(i));
    Commissions_ben(i)  = (P_remain_contract(i) * commision);
end

% liabilities as the sum of discounted expected cash flows
Liabilities = V' * discounts(2:end);

% Macaulay Duration of liabilities
M_duration = ( (V .* discounts(2:end) )' * time(2:end) ) / Liabilities;

% compute Best Estimate Liabilities components (BEL)
Lapse_BEL       = (Lapse_ben)' * discounts(2:end);
Death_BEL       = (Death_ben)' * discounts(2:end);
Expenses_BEL    = (Expenses_ben)' * discounts(2:end);
Commissions_BEL = (Commissions_ben)' * discounts(2:end);

end