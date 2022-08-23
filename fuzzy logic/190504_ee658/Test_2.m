clc;
clear;
close all;

%% Initializing the input variable xvar1 and xvar2 as random variable
xvar1 = 20.*rand(1,20);
xvar2 = 10.*rand(1,20);

%% Generation of Rule-I
x11(:,1) = trimf(xvar1,[0 0 20]);
x21(:,1) = trimf(xvar2,[5.5 10 10]);
y1(:,1) = -0.14*xvar1 + 0.65*xvar2 + 27.6;

%% Generation of Rule-II
x12(:,2) = trimf(xvar1,[0 0 20]);
x22(:,2) = trimf(xvar2,[0 0 10]);
y2(:,2) = 0.91*xvar1 + 2.06*xvar2 + 2.3;

%% Generation of Membership values corresponding to inputs
A = [x11(:,1) x12(:,2)];
B = [x21(:,1) x22(:,2)];

%% Application of min operation for membership values corresponding to each rule for finding weights
C = min(A,B);

%% Summation of weight parameters corresponding to Rule
D = sum(C(1:20,:));

%% Set of output parameters corresponding to each rule
y = [y1(:,1) y2(:,2)];

%% Normalization of weight vector
for i = 1:20
    for j = 1:2
        E(i,j) = C(i,j)/D(1,j);
        P(i,:) = [1,xvar1(:,i),xvar2(:,i)];
    end
end

%% Generation of normalized output parameters Y by multiplication of weights with output parameters corresponding to each rule
N = (E.*y);
Y = N(:,1)+N(:,2);

%% Generation of X matrices by multiplying the normalized weights with input parameters
m = 1;
for k = 1:3
    for l = 1:2
        X(:,m) = [E(:,l).*P(:,k)];
        m = m+1;
    end
end

%% Calculation of parameters by least square method with
P = inv(X'*X)*X'*Y;
M = X\Y;

