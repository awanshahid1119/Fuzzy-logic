

clear;
clc;

%% Load the dataset generated from 14 initial states.
load tb1.dat
load tb2.dat
load tb3.dat
load tb4.dat
load tb5.dat
load tb6.dat
load tb7.dat
load tb8.dat
load tb9.dat
load tb10.dat
load tb11.dat
load tb12.dat
load tb13.dat
load tb14.dat
tb = [tb1; tb2; tb3; tb4; tb5; tb6; tb7;tb8; tb9; tb10; tb11; tb12; tb13; tb14];
q = 4; % length of truck

%% Step 1: Divide the input and output spaces into fuzzy regions.
% Divide each domain interval into 2N+1 regions
% m_x, m_phi and m_theta

x = 0:0.1:20; % Universe of discourse of input x.
M_x(:,1) = trapmf(x,[0 0 1.5 7]);
M_x(:,2) = trimf(x,[4 7 10]);
M_x(:,3) = trimf(x,[9 10 11]);
M_x(:,4) = trimf(x,[10 13 16]);
M_x(:,5) = trapmf(x,[13 18.5 20 20]);

theta = -40:0.4:40; % Universe of discourse of output theta.
M_theta(:,1) = trimf(theta,[-40 -40 -20]);
M_theta(:,2) = trimf(theta,[-33 -20 -7]);
M_theta(:,3) = trimf(theta,[-14 -7 0]);
M_theta(:,4) = trimf(theta,[-4 0 4]);
M_theta(:,5) = trimf(theta,[0 7 14]);
M_theta(:,6) = trimf(theta,[7 20 33]);
M_theta(:,7) = trimf(theta,[20 40 40]);

phi = -115:0.1:295; % Universe of discourse of input phi.
M_phi(:,1) = trimf(phi,[-115 -65 -15]);
M_phi(:,2) = trimf(phi,[-45 0 45]);
M_phi(:,3) = trimf(phi,[15 52.5 90]);
M_phi(:,4) = trimf(phi,[80 90 100]);
M_phi(:,5) = trimf(phi,[90 127.5 165]);
M_phi(:,6) = trimf(phi,[135 180 225]);
M_phi(:,7) = trimf(phi,[195 245 295]);



%% Step 2: Generate Fuzzy Rules from Given Data Pairs.
% Calculate the Degree for given datasets

% Input x
MU_x_1 = trapmf(tb(:,1),[0 0 1.5 7]);
MU_x_2 = trimf(tb(:,1),[4 7 10]);
MU_x_3 = trimf(tb(:,1),[9 10 11]);
MU_x_4 = trimf(tb(:,1),[10 13 16]);
MU_x_5 = trapmf(tb(:,1),[13 18.5 20 20]);
MU_x = [MU_x_1 MU_x_2 MU_x_3 MU_x_4 MU_x_5];

%Output theta
MU_theta_1 = trimf(tb(:,3),[-40 -40 -20]);
MU_theta_2 = trimf(tb(:,3),[-33 -20 -7]);
MU_theta_3 = trimf(tb(:,3),[-14 -7 0]);
MU_theta_4 = trimf(tb(:,3),[-4 0 4]);
MU_theta_5 = trimf(tb(:,3),[0 7 14]);
MU_theta_6 = trimf(tb(:,3),[7 20 33]);
MU_theta_7 = trimf(tb(:,3),[20 40 40]);
MU_theta = [MU_theta_1 MU_theta_2 MU_theta_3 MU_theta_4 MU_theta_5 MU_theta_6 MU_theta_7];

% Input phi
MU_phi_1 = trimf(tb(:,2),[-115 -65 -15]);
MU_phi_2 = trimf(tb(:,2),[-45 0 45]);
MU_phi_3 = trimf(tb(:,2),[15 52.5 90]);
MU_phi_4 = trimf(tb(:,2),[80 90 100]);
MU_phi_5 = trimf(tb(:,2),[90 127.5 165]);
MU_phi_6 = trimf(tb(:,2),[135 180 225]);
MU_phi_7 = trimf(tb(:,2),[195 245 295]);
MU_phi = [MU_phi_1 MU_phi_2 MU_phi_3 MU_phi_4 MU_phi_5 MU_phi_6 MU_phi_7];

[deg_x,ind_x] = max(MU_x,[],2);
[deg_phi,ind_phi] = max(MU_phi,[],2);
[deg_theta,ind_theta] = max(MU_theta,[],2);

rule = [ind_x ind_phi ind_theta];
rule_deg = deg_x.*deg_phi.*deg_theta;

%% Step 3: Assign a Degree to Each Rule.
% For conflicting rule in rule base, assign the maximum degree for the
% conflict group.
temp = sortrows([rule rule_deg]);
[p q]=size(temp);
value=temp(1,4);
j=1;
m=zeros(p,4);
c=0;
for i=1:p-1
    if(isequal([temp(i,1) temp(i,2)],[temp(i+1,1) temp(i+1,2)]))
        if(value<temp(i+1,4)) 
            value=temp(i+1,4);
            j=i+1;
        end    
    else
      value=temp(j+1,4);  
      m(c+1,:)=temp(j,:);
      c=c+1;
    end
end
M=m(1:c,:);
M_final=[M;  5 7 7 0.3625];
rules = sortrows(M_final);

%% Step 4: Create a Combined Fuzzy Rule Base.
% Not required in this case. It requires only in the case where or
% connectors are used.

%% Step 5: Determine a Mapping Based on the Combined Fuzzy Rule Base.
% Defuzzification also present in this step. The defuzzification strategy
% is Centroid of Area (COA).
figure;
axis([0 20 -10 100]); % x axis and y-axis limit.
iterat = 70;
y_t = [-40 -20 -7 0 7 20 40]; % the points where membership value is 1.

%% Input 1
in_x = 9; % Input x (by user)
in_phi = 220; % Input phi (by user)
% theta_output = 0; % Output theta
sam_y = 2; % initial y position

% plot(x_input,y_sample,'r*');
% title(sprintf('Case 1: x = %s and phi = %s - Truck Trajectory', num2str(x_input), num2str(in_phi)));
% hold on;

fin_y = zeros(1,1);
fin_x = zeros(1,1);
fin_x(1,1) = in_x;
fin_y(1,1) = sam_y;
for i = 1:iterat % no. of iterations for trajectory tracking.
    
    % value of x.
    mu_x1 = trapmf(in_x,[0 0 1.5 7]);
    mu_x2 = trimf(in_x,[4 7 10]);
    mu_x3 = trimf(in_x,[9 10 11]);
    mu_x4 = trimf(in_x,[10 13 16]);
    mu_x5 = trapmf(in_x,[13 18.5 20 20]);
    
    % value of phi.
    mu_phi1 = trimf(in_phi,[-115 -65 -15]);
    mu_phi2 = trimf(in_phi,[-45 0 45]);
    mu_phi3 = trimf(in_phi,[15 52.5 90]);
    mu_phi4 = trimf(in_phi,[80 90 100]);
    mu_phi5 = trimf(in_phi,[90 127.5 165]);
    mu_phi6 = trimf(in_phi,[135 180 225]);
    mu_phi7 = trimf(in_phi,[195 245 295]);

    mu_x = [mu_x1 mu_x2 mu_x3 mu_x4 mu_x5];
    mu_p = [mu_phi1 mu_phi2 mu_phi3 mu_phi4 mu_phi5 mu_phi6 mu_phi7];
    mo = mu_x(rules(:,1)).*mu_p(rules(:,2)); % product operation to determine the degree of output control.
    y_bar = y_t(rules(:,3));
    out_theta = sum(mo.*y_bar)/sum(mo); % Value of theta from the rule base.
    
    % Approximate kinematics of truck backer upper control for calculating the next states.
    in_x = in_x + cosd(in_phi + out_theta) + sind(out_theta)*sind(in_phi);
    fin_x(i+1,1) = in_x;
    sam_y = sam_y + sind(in_phi + out_theta) - sind(out_theta)*cosd(in_phi); % displacement on y axis.
    fin_y(i+1,1) = sam_y;
    in_phi = in_phi - asind(2*sind(out_theta)/q);
    
end

% Plot the trajectory.
plot(fin_x,fin_y,'.-','MarkerSize',12);
set(gca,'Fontsize',10,'FontName','Times New Roman');
hold on;
xlabel('x','FontSize',10,'FontName','Times New Roman'); 

%% Input 2
in_x = 0; % Input x (by user)
in_phi = -30; % Input phi (by user)
% theta_output = 0; % Output theta
sam_y = 2; % initial y position

% plot(x_input,y_sample,'r*');
% title(sprintf('Case 1: x = %s and phi = %s - Truck Trajectory', num2str(x_input), num2str(in_phi)));
% hold on;

fin_y = zeros(1,1);
fin_x = zeros(1,1);
fin_x(1,1) = in_x;
fin_y(1,1) = sam_y;
for i = 1:iterat % no. of iterations for trajectory tracking.
    
    % value of x.
    mu_x1 = trapmf(in_x,[0 0 1.5 7]);
    mu_x2 = trimf(in_x,[4 7 10]);
    mu_x3 = trimf(in_x,[9 10 11]);
    mu_x4 = trimf(in_x,[10 13 16]);
    mu_x5 = trapmf(in_x,[13 18.5 20 20]);
    
    % value of phi.
    mu_phi1 = trimf(in_phi,[-115 -65 -15]);
    mu_phi2 = trimf(in_phi,[-45 0 45]);
    mu_phi3 = trimf(in_phi,[15 52.5 90]);
    mu_phi4 = trimf(in_phi,[80 90 100]);
    mu_phi5 = trimf(in_phi,[90 127.5 165]);
    mu_phi6 = trimf(in_phi,[135 180 225]);
    mu_phi7 = trimf(in_phi,[195 245 295]);

    mu_x = [mu_x1 mu_x2 mu_x3 mu_x4 mu_x5];
    mu_p = [mu_phi1 mu_phi2 mu_phi3 mu_phi4 mu_phi5 mu_phi6 mu_phi7];
    mo = mu_x(rules(:,1)).*mu_p(rules(:,2)); % product operation to determine the degree of output control.
    y_bar = y_t(rules(:,3));
    out_theta = sum(mo.*y_bar)/sum(mo); % Value of theta from the rule base.
    
    % Approximate kinematics of truck backer upper control for calculating the next states.
    in_x = in_x + cosd(in_phi + out_theta) + sind(out_theta)*sind(in_phi);
    fin_x(i+1,1) = in_x;
    sam_y = sam_y + sind(in_phi + out_theta) - sind(out_theta)*cosd(in_phi); % displacement on y axis.
    fin_y(i+1,1) = sam_y;
    in_phi = in_phi - asind(2*sind(out_theta)/q);
    
end

% Plot the trajectory.
plot(fin_x,fin_y,'.-','MarkerSize',12);
set(gca,'Fontsize',10,'FontName','Times New Roman');
hold on;
xlabel('x','FontSize',10,'FontName','Times New Roman');

%% Input 3
in_x = 14; % Input x (by user)
in_phi = 30; % Input phi (by user)
% theta_output = 0; % Output theta
sam_y = 2; % initial y position

% plot(x_input,y_sample,'r*');
% title(sprintf('Case 1: x = %s and phi = %s - Truck Trajectory', num2str(x_input), num2str(in_phi)));
% hold on;

fin_y = zeros(1,1);
fin_x = zeros(1,1);
fin_x(1,1) = in_x;
fin_y(1,1) = sam_y;
for i = 1:iterat % no. of iterations for trajectory tracking.
    
    % value of x.
    mu_x1 = trapmf(in_x,[0 0 1.5 7]);
    mu_x2 = trimf(in_x,[4 7 10]);
    mu_x3 = trimf(in_x,[9 10 11]);
    mu_x4 = trimf(in_x,[10 13 16]);
    mu_x5 = trapmf(in_x,[13 18.5 20 20]);
    
    % value of phi.
    mu_phi1 = trimf(in_phi,[-115 -65 -15]);
    mu_phi2 = trimf(in_phi,[-45 0 45]);
    mu_phi3 = trimf(in_phi,[15 52.5 90]);
    mu_phi4 = trimf(in_phi,[80 90 100]);
    mu_phi5 = trimf(in_phi,[90 127.5 165]);
    mu_phi6 = trimf(in_phi,[135 180 225]);
    mu_phi7 = trimf(in_phi,[195 245 295]);

    mu_x = [mu_x1 mu_x2 mu_x3 mu_x4 mu_x5];
    mu_p = [mu_phi1 mu_phi2 mu_phi3 mu_phi4 mu_phi5 mu_phi6 mu_phi7];
    mo = mu_x(rules(:,1)).*mu_p(rules(:,2)); % product operation to determine the degree of output control.
    y_bar = y_t(rules(:,3));
    out_theta = sum(mo.*y_bar)/sum(mo); % Value of theta from the rule base.
    
    % Approximate kinematics of truck backer upper control for calculating the next states.
    in_x = in_x + cosd(in_phi + out_theta) + sind(out_theta)*sind(in_phi);
    fin_x(i+1,1) = in_x;
    sam_y = sam_y + sind(in_phi + out_theta) - sind(out_theta)*cosd(in_phi); % displacement on y axis.
    fin_y(i+1,1) = sam_y;
    in_phi = in_phi - asind(2*sind(out_theta)/q);
    
end

% Plot the trajectory.
plot(fin_x,fin_y,'.-','MarkerSize',12);
set(gca,'Fontsize',10,'FontName','Times New Roman');
hold on;
xlabel('x','FontSize',10,'FontName','Times New Roman');




