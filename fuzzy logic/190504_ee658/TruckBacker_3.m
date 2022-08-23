%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

clear;
close all;
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
tb = [tb1; tb2; tb3; tb4; tb5; tb6; tb7; tb8; tb9; tb10; tb11; tb12; tb13; tb14;];
b = 4; % length of truck

%% Step 1: Divide the input and output spaces into fuzzy regions.
% Divide each domain interval into 2N+1 regions
% m_x, m_phi and m_theta

x = 0:0.1:20; % Universe of discourse of input x.
m_x(:,1) = trapmf(x,[0 0 1.5 7]);
m_x(:,2) = trimf(x,[4 7 10]);
m_x(:,3) = trimf(x,[9 10 11]);
m_x(:,4) = trimf(x,[10 13 16]);
m_x(:,5) = trapmf(x,[13 18.5 20 20]);

phi = -115:0.1:295; % Universe of discourse of input phi.
m_phi(:,1) = trimf(phi,[-115 -65 -15]);
m_phi(:,2) = trimf(phi,[-45 0 45]);
m_phi(:,3) = trimf(phi,[15 52.5 90]);
m_phi(:,4) = trimf(phi,[80 90 100]);
m_phi(:,5) = trimf(phi,[90 127.5 165]);
m_phi(:,6) = trimf(phi,[135 180 225]);
m_phi(:,7) = trimf(phi,[195 245 295]);

theta = -40:0.4:40; % Universe of discourse of output theta.
m_theta(:,1) = trimf(theta,[-40 -40 -20]);
m_theta(:,2) = trimf(theta,[-33 -20 -7]);
m_theta(:,3) = trimf(theta,[-14 -7 0]);
m_theta(:,4) = trimf(theta,[-4 0 4]);
m_theta(:,5) = trimf(theta,[0 7 14]);
m_theta(:,6) = trimf(theta,[7 20 33]);
m_theta(:,7) = trimf(theta,[20 40 40]);

%% Step 2: Generate Fuzzy Rules from Given Data Pairs.
% Calculate the Degree for given datasets

% Input x
mu_x1 = trapmf(tb(:,1),[0 0 1.5 7]);
mu_x2 = trimf(tb(:,1),[4 7 10]);
mu_x3 = trimf(tb(:,1),[9 10 11]);
mu_x4 = trimf(tb(:,1),[10 13 16]);
mu_x5 = trapmf(tb(:,1),[13 18.5 20 20]);
mu_x_test = [mu_x1 mu_x2 mu_x3 mu_x4 mu_x5];

% Input phi
mu_phi1 = trimf(tb(:,2),[-115 -65 -15]);
mu_phi2 = trimf(tb(:,2),[-45 0 45]);
mu_phi3 = trimf(tb(:,2),[15 52.5 90]);
mu_phi4 = trimf(tb(:,2),[80 90 100]);
mu_phi5 = trimf(tb(:,2),[90 127.5 165]);
mu_phi6 = trimf(tb(:,2),[135 180 225]);
mu_phi7 = trimf(tb(:,2),[195 245 295]);
mu_phi = [mu_phi1 mu_phi2 mu_phi3 mu_phi4 mu_phi5 mu_phi6 mu_phi7];

%Output theta
mu_theta1 = trimf(tb(:,3),[-40 -40 -20]);
mu_theta2 = trimf(tb(:,3),[-33 -20 -7]);
mu_theta3 = trimf(tb(:,3),[-14 -7 0]);
mu_theta4 = trimf(tb(:,3),[-4 0 4]);
mu_theta5 = trimf(tb(:,3),[0 7 14]);
mu_theta6 = trimf(tb(:,3),[7 20 33]);
mu_theta7 = trimf(tb(:,3),[20 40 40]);
mu_theta = [mu_theta1 mu_theta2 mu_theta3 mu_theta4 mu_theta5 mu_theta6 mu_theta7];

[avg_degree_x,avg_index_x] = max(mu_x_test,[],2);
[avg_degree_phi,avg_index_phi] = max(mu_phi,[],2);
[avg_degree_theta,avg_index_theta] = max(mu_theta,[],2);

rules = [avg_index_x avg_index_phi avg_index_theta];
rules_degree = avg_degree_x.*avg_degree_phi.*avg_degree_theta;

%% Step 3: Assign a Degree to Each Rule.
% For conflicting rule in rule base, assign the maximum degree for the
% conflict group.
r_temp = sortrows([rules rules_degree]);
% r_temp = r_temp(end:-1:1,:);
% rules_degree = rules_degree(end:-1:1,:);
% rules = rules(end:-1:1,:);
temp1(:,1:2) = r_temp(:,1:2);
% temp1(:,3) = r_temp(:,4);
[new_rule, rule_index] = unique(temp1,'rows','sorted');
Rule_temp = zeros(size(new_rule,1),4);
for i = 1:size(new_rule,1)
    if i==size(new_rule,1)
        [temp3, temp4] = max( r_temp(rule_index(i):size(r_temp,1),4) );
        Rule_temp(i,:) = r_temp(temp4+rule_index(i)-1,:);
%     elseif i==1
%         [temp3, temp4] = max( r_temp(rule_index(i):rule_index(i+1)-1,4) );
%         Rule_temp(i,:) = r_temp(temp4,:);
    else
        [temp3, temp4] = max( r_temp(rule_index(i):rule_index(i+1)-1,4) );
        Rule_temp(i,:) = r_temp(temp4+rule_index(i)-1,:);
    end
end
% new_rule = [new_rule r_temp(rule_index,3:4)];
% final_rules = new_rule;
final_rules = sortrows(Rule_temp);

%% Step 4: Create a Combined Fuzzy Rule Base.
% Not required in this case. It requires only in the case where or
% connectors are used.

%% Step 5: Determine a Mapping Based on the Combined Fuzzy Rule Base.
% Defuzzification also present in this step. The defuzzification strategy
% is Centroid of Area (COA).
figure;
axis([0 20 -10 100]); % x axis and y-axis limit.
iter = 5;
y_t = [-40 -20 -7 0 7 20 40]; % the points where membership value is 1.

%% Input 1
x_input = 10; % Input x (by user)
phi_input = 220; % Input phi (by user)
% theta_output = 0; % Output theta
y_sample = 2; % initial y position
final_rules = sortrows(Rule_temp(6:22,:));
% plot(x_input,y_sample,'r*');
% title(sprintf('Case 1: x = %s and phi = %s - Truck Trajectory', num2str(x_input), num2str(phi_input)));
% hold on;

y_final = zeros(1,1);
x_final = zeros(1,1);
x_final(1,1) = x_input;
y_final(1,1) = y_sample;
for i = 1:iter % no. of iterations for trajectory tracking.
    
    % value of x.
    mu_x1_test = trapmf(x_input,[0 0 1.5 7]);
    mu_x2_test = trimf(x_input,[4 7 10]);
    mu_x3_test = trimf(x_input,[9 10 11]);
    mu_x4_test = trimf(x_input,[10 13 16]);
    mu_x5_test = trapmf(x_input,[13 18.5 20 20]);
    
    % value of phi.
    mu_phi1_test = trimf(phi_input,[-115 -65 -15]);
    mu_phi2_test = trimf(phi_input,[-45 0 45]);
    mu_phi3_test = trimf(phi_input,[15 52.5 90]);
    mu_phi4_test = trimf(phi_input,[80 90 100]);
    mu_phi5_test = trimf(phi_input,[90 127.5 165]);
    mu_phi6_test = trimf(phi_input,[135 180 225]);
    mu_phi7_test = trimf(phi_input,[195 245 295]);

    mu_x_test = [mu_x1_test mu_x2_test mu_x3_test mu_x4_test mu_x5_test];
    mu_p_test = [mu_phi1_test mu_phi2_test mu_phi3_test mu_phi4_test mu_phi5_test mu_phi6_test mu_phi7_test];
    mo = mu_x_test(final_rules(:,1)).*mu_p_test(final_rules(:,2)); % product operation to determine the degree of output control.
    y_bar = y_t(final_rules(:,3));
    theta_output = sum(mo.*y_bar)/sum(mo); % Value of theta from the rule base.
    
    % Approximate kinematics of truck backer upper control for calculating the next states.
    x_input = x_input + cosd(phi_input + theta_output) + sind(theta_output)*sind(phi_input);
    x_final(i+1,1) = x_input;
    y_sample = y_sample + sind(phi_input + theta_output) - sind(theta_output)*cosd(phi_input); % displacement on y axis.
    y_final(i+1,1) = y_sample;
    phi_input = phi_input - asind(2*sind(theta_output)/b);
    
end

% Plot the trajectory.
plot(x_final,y_final,'.-','MarkerSize',12);
set(gca,'Fontsize',10,'FontName','Times New Roman');
hold on;
xlabel('x','FontSize',10,'FontName','Times New Roman'); 

%% Input 2
x_input = 3; % Input x (by user)
phi_input = -30; % Input phi (by user)
% theta_output = 0; % Output theta
y_sample = 2; % initial y position
final_rules = sortrows(Rule_temp(1:2,:));

% plot(x_input,y_sample,'r*');
% title(sprintf('Case 1: x = %s and phi = %s - Truck Trajectory', num2str(x_input), num2str(phi_input)));
% hold on;

y_final = zeros(1,1);
x_final = zeros(1,1);
x_final(1,1) = x_input;
y_final(1,1) = y_sample;
for i = 1:iter % no. of iterations for trajectory tracking.
    
    % value of x.
    mu_x1_test = trapmf(x_input,[0 0 1.5 7]);
    mu_x2_test = trimf(x_input,[4 7 10]);
    mu_x3_test = trimf(x_input,[9 10 11]);
    mu_x4_test = trimf(x_input,[10 13 16]);
    mu_x5_test = trapmf(x_input,[13 18.5 20 20]);
    
    % value of phi.
    mu_phi1_test = trimf(phi_input,[-115 -65 -15]);
    mu_phi2_test = trimf(phi_input,[-45 0 45]);
    mu_phi3_test = trimf(phi_input,[15 52.5 90]);
    mu_phi4_test = trimf(phi_input,[80 90 100]);
    mu_phi5_test = trimf(phi_input,[90 127.5 165]);
    mu_phi6_test = trimf(phi_input,[135 180 225]);
    mu_phi7_test = trimf(phi_input,[195 245 295]);

    mu_x_test = [mu_x1_test mu_x2_test mu_x3_test mu_x4_test mu_x5_test];
    mu_p_test = [mu_phi1_test mu_phi2_test mu_phi3_test mu_phi4_test mu_phi5_test mu_phi6_test mu_phi7_test];
    mo = mu_x_test(final_rules(:,1)).*mu_p_test(final_rules(:,2)); % product operation to determine the degree of output control.
    y_bar = y_t(final_rules(:,3));
    theta_output = sum(mo.*y_bar)/sum(mo); % Value of theta from the rule base.
    
    % Approximate kinematics of truck backer upper control for calculating the next states.
    x_input = x_input + cosd(phi_input + theta_output) + sind(theta_output)*sind(phi_input);
    x_final(i+1,1) = x_input;
    y_sample = y_sample + sind(phi_input + theta_output) - sind(theta_output)*cosd(phi_input); % displacement on y axis.
    y_final(i+1,1) = y_sample;
    phi_input = phi_input - asind(2*sind(theta_output)/b);
    
end

% Plot the trajectory.
plot(x_final,y_final,'.-','MarkerSize',12);
set(gca,'Fontsize',10,'FontName','Times New Roman');
hold on;
xlabel('x','FontSize',10,'FontName','Times New Roman');

%% Input 3
x_input = 13; % Input x (by user)
phi_input = 30; % Input phi (by user)
% theta_output = 0; % Output theta
y_sample = 2; % initial y position
final_rules = sortrows(Rule_temp(17:27,:));

% plot(x_input,y_sample,'r*');
% title(sprintf('Case 1: x = %s and phi = %s - Truck Trajectory', num2str(x_input), num2str(phi_input)));
% hold on;

y_final = zeros(1,1);
x_final = zeros(1,1);
x_final(1,1) = x_input;
y_final(1,1) = y_sample;
for i = 1:iter % no. of iterations for trajectory tracking.
    
    % value of x.
    mu_x1_test = trapmf(x_input,[0 0 1.5 7]);
    mu_x2_test = trimf(x_input,[4 7 10]);
    mu_x3_test = trimf(x_input,[9 10 11]);
    mu_x4_test = trimf(x_input,[10 13 16]);
    mu_x5_test = trapmf(x_input,[13 18.5 20 20]);
    
    % value of phi.
    mu_phi1_test = trimf(phi_input,[-115 -65 -15]);
    mu_phi2_test = trimf(phi_input,[-45 0 45]);
    mu_phi3_test = trimf(phi_input,[15 52.5 90]);
    mu_phi4_test = trimf(phi_input,[80 90 100]);
    mu_phi5_test = trimf(phi_input,[90 127.5 165]);
    mu_phi6_test = trimf(phi_input,[135 180 225]);
    mu_phi7_test = trimf(phi_input,[195 245 295]);

    mu_x_test = [mu_x1_test mu_x2_test mu_x3_test mu_x4_test mu_x5_test];
    mu_p_test = [mu_phi1_test mu_phi2_test mu_phi3_test mu_phi4_test mu_phi5_test mu_phi6_test mu_phi7_test];
    mo = mu_x_test(final_rules(:,1)).*mu_p_test(final_rules(:,2)); % product operation to determine the degree of output control.
    y_bar = y_t(final_rules(:,3));
    theta_output = sum(mo.*y_bar)/sum(mo); % Value of theta from the rule base.
    
    % Approximate kinematics of truck backer upper control for calculating the next states.
    x_input = x_input + cosd(phi_input + theta_output) + sind(theta_output)*sind(phi_input);
    x_final(i+1,1) = x_input;
    y_sample = y_sample + sind(phi_input + theta_output) - sind(theta_output)*cosd(phi_input); % displacement on y axis.
    y_final(i+1,1) = y_sample;
    phi_input = phi_input - asind(2*sind(theta_output)/b);
    
end

% Plot the trajectory.
plot(x_final,y_final,'.-','MarkerSize',12);
set(gca,'Fontsize',10,'FontName','Times New Roman');
hold on;
xlabel('x','FontSize',10,'FontName','Times New Roman');



