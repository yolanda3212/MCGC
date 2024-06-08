close all; clear all£» clc
load('NH.mat');
maxIters = 30;
alpha = 0.01; 
d = 45;
gamma = 0.1;
lambda = 0.000001;
result = MCGC(X, alpha, lambda, d, gamma, maxIters, gt)
