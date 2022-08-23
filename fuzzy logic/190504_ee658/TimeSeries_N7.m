
clear;
close all;
clc;

%% Input: Mackey-Glass Chaotic Time Series and Output: Predicted Series
load mgdata.dat; % Mackey-Glass Chaotic Time Series Data (inbuilt in matlab).

%% Input Output pairs formation
m = 9; % points taken for predicting one value (gievn in paper)
input_mgdata = zeros(1,1);
output_mgdata = zeros(1,1);
for ki = 1:(size(mgdata,1)-m) % creating (M-m) pairs of input and output (refer the paper).
   input_mgdata(ki,1:m) = mgdata(ki:ki+m-1,2)'; % taking only second column because first column has only indices to represent on x-axis.
   output_mgdata(ki,1) = mgdata(ki+m,2)'; % taking (m+1)th column value as output which we will get in next state.
end

train_data = [input_mgdata(1:700,:) output_mgdata(1:700,:)]; % according to paper training data has 700 samples.
test_data = input_mgdata(701:1000,:); % according to paper test data has remaining 300 samples.

% Plot of Mackey-Glass chaotic times series (data is already provided)
figure;
plot(mgdata(:,1),mgdata(:,2),'LineWidth',1.0);
set(gca,'FontSize',13,'FontName','Times New Roman');
xlabel('t','FontSize',20,'FontName','Times New Roman');
ylabel('x(t)','FontSize',20,'FontName','Times New Roman');
title(sprintf('Mackey-Glass Chaotic Time Series for tau = 17'));

%%  Step 1: Divide the input and output spaces into fuzzy regions.
% Divide each domain interval into 2N+1 regions
   
N = 7;
x = 0.2:0.01:1.6;
[x_R, R] = divide_fuzzyregions(N, x);
[y_t] = center_value(x_R, x);  % Find center value with membership value = 1
    
% Plot Membership Function
figure;
for i = 1:size(x_R,2)
    plot(x,x_R{1,i},'Linewidth',1.0);
    hold on;
end
set(gca,'FontSize',13,'FontName','Times New Roman');
title(sprintf('Membership Function for Chaotic Time Series Prediction (N = %d and R = %d)', N, R));
xlabel('x(t)','FontSize',20);
ylabel('\mu(x)','FontSize',20)
    
%% Step 2: Generate Fuzzy Rules from Given Data Pairs.
% Calculate the Degree for given input output pairs

Degree_Value = zeros(1,1);
Rule_Value = zeros(1,1);
for i = 1:size(train_data,1)
    for j = 1:size(train_data,2)
        [Degree_Value(i,j), Rule_Value(i,j)] = find_degree(train_data(i,j), x, x_R, R);
    end
end
clear i;

%% Step 3: Assign a Degree to Each Rule.
% For conflicting rule in rule base, assign the maximum degree for the
% conflict group.

Degree_Rule1 = prod(Degree_Value,2);
[tmp, index] = unique(Rule_Value,'rows','stable'); % unique rules. (max among the similar will be taken)
NewRule_Degree = Degree_Rule1(index); % Degree of unique rules obtained.
[a,b,c] = unique(tmp(:,1:m),'rows','stable'); % value of m is 9 as given in paper (see before step 1).

k = 1;
% final_matrix = zeros(1,1);
for i = 1:size(b,1)
        dup_rows_index = find(c==i);                        % Identify rows having same index
        if length(dup_rows_index) > 1                       % Check no. of matching
             [u,v] = max(NewRule_Degree(dup_rows_index));   % Find row having maximum rule degree
             final_matrix(k,:) = tmp(dup_rows_index(v),:);  % Keep that rule and put in final matrix
        else
             final_matrix(k,:) = tmp(dup_rows_index,:) ;  
        end
        k = k + 1;
end

%% Step 4: Create a Combined Fuzzy Rule Base.
% Not required in this case. It requires only in the case where or
% connectors are used.

fuzzy_rule_base = final_matrix;
    
%% Step 5: Determine a Mapping Based on the Combined Fuzzy Rule Base.
% Defuzzification also present in this step. The defuzzification strategy
% is Centroid of Area (COA).

Degree_test = zeros(1,1);
in_mf_prod = zeros(1,1);
y_bar = zeros(1,1);
test_output = zeros(1,1);
for m = 1:size(test_data,1)
    sample = test_data(m,:); 
    for i = 1:size(fuzzy_rule_base,1)
        for j = 1:size(test_data,2)
            val = find_degree_test(sample(j),x,x_R,R); 
            Degree_test(1,j) = val(fuzzy_rule_base(i,j));
        end
        in_mf_prod(i,:) = prod(Degree_test);
        y_bar(i,:) = y_t(fuzzy_rule_base(i,10)) ;  
    end
    test_output(m,1) = sum(in_mf_prod.*y_bar)/sum(in_mf_prod) ;
end

%% Plot the result
figure;
plot(1:300,output_mgdata(701:1000),'LineWidth',1.0);
set(gca,'FontSize',13,'FontName','Times New Roman');
hold on;
plot(1:300,test_output,'-r','LineWidth',1.0);
title(sprintf('Chaotic Time Series Prediction using %d regions',R));
xlabel('x(t)','FontSize',20);
ylabel('\mu(x)','FontSize',20);
legend('Actual Value','Predicted Value');
legend boxoff;






