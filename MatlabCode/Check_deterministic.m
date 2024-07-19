%% Check excel computations fund value deterministic
% % % Fund_det = zeros(T+1,1);
% % % Fund_det(1) = F0;
% % % 
% % % Fund_det_up = zeros(T+1,1);
% % % Fund_det_up(1) = F0;
% % % 
% % % Fund_det_down = zeros(T+1,1);
% % % Fund_det_down(1) = F0;
% % % 
% % % for i = 2:T+1
% % %     Fund_det(i) = Fund_det(i-1) * (1 - RD) / (fwd_discount(i));
% % % end
% % % 
% % % for i = 2:T+1
% % %     Fund_det_up(i) = Fund_det_up(i-1) * (1 - RD) / (fwd_discount_up(i));
% % % end
% % % 
% % % for i = 2:T+1
% % %     Fund_det_down(i) = Fund_det_down(i-1) * (1 - RD) / (fwd_discount_down(i));
% % % end
% % % 
% % % disp(['The deterministic value of the fund at the end of the simulation is: ', num2str(Fund_det(end))]);
% % % disp(['The deterministic value of the fund at the end of the simulation is: ', num2str(Fund_det_up(end))]);
% % % disp(['The deterministic value of the fund at the end of the simulation is: ', num2str(Fund_det_down(end))]);