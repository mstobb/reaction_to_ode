function [time,y] = example_ode(t_final,t_start)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 1 
%
%
% This file was created by issuing command: 
% 	python make_ode_m_file.py example_stoich.csv example_ode.m
%


if nargin == 1
     t_start = 0;  % Default start time is 0 
elseif nargin == 0
     t_start = 0; % Default start time is 0
     t_final = 1; % Default final time is 1
end


% Kinetic Parameters 
k_1 = 1; 
k_3 = 1; 
k_1 = 1; 
k_3 = 1; 
k_4 = 1; 

p = [ k_1, k_3, k_1, k_3, k_4];


% Initial Conditions 
A_IC = 0; 
B_IC = 0; 
C_IC = 0; 
D_IC = 0; 
E_IC = 0; 
F_IC = 0; 

init_cond = [ A_IC, B_IC, C_IC, D_IC, E_IC, F_IC];


options = odeset('RelTol',1e-12,'AbsTol',1e-23);


%------------------------------ Main Solve ---------------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,p), [t_start t_final], init_cond, options);
%---------------------------------------------------------------------%


% Rename solution components
A = y(:,1); 
B = y(:,2); 
C = y(:,3); 
D = y(:,4); 
E = y(:,5); 
F = y(:,6); 



%  
% Place plots or other calculations here
%   
% Example: 
% plot(time, A, 'k-o', 'LineWidth', 4, 'MarkerSize', 4); legend('A');


end



%-----------------------------------------------------%
%-------------------- RHS Function -------------------%
%-----------------------------------------------------%

function dy = RHS(t,y,p)

dy = zeros(6,1);


% Rename Variables 

A   = y(1); 
B   = y(2); 
C   = y(3); 
D   = y(4); 
E   = y(5); 
F   = y(6); 


% Rename Kinetic Parameters 
k_1 = p(1);  
k_3 = p(2);  
k_1 = p(3);  
k_3 = p(4);  
k_4 = p(5);  


% Reactions to be used 

%{ 

('0', 'k_1', 'k_3', 'k_1', 'k_3', 'k_4')
('A', '-1', '1', '0', '-1', '1')
('B', '-2', '2', '-1', '0', '0')
('C', '1', '-1', '0', '0', '0')
('D', '0', '0', '-1', '0', '0')
('E', '0', '0', '1', '-3', '3')
('F', '0', '0', '0', '1', '-1')


%} 



% ODEs from reaction equations 

% A
 dy(1)  =  -  k_1 * A * B^2  +  k_3 * C  -  k_3 * A * E^3  +  k_4 * F;

% B
 dy(2)  =  -  k_1 * A * B^2  +  k_3 * 2 * C  -  k_1 * B * D;

% C
 dy(3)  =  +  k_1 * A * B^2  -  k_3 * C;

% D
 dy(4)  =  -  k_1 * B * D;

% E
 dy(5)  =  +  k_1 * B * D  -  k_3 * A * E^3  +  k_4 * 3 * F;

% F
 dy(6)  =  +  k_3 * A * E^3  -  k_4 * F;





end