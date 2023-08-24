clc
clear
addpath(genpath('funs'))
NITER=20;

load('MSRC');
for i = 1 :length(X)
    X{i} = full((X{i} - mean(X{i}, 2)) ./ repmat(std(X{i}, [], 2), 1, size(X{i}, 2)));
end
c = length(unique(Y));
Avd = constructA_vd(X, 5, 6);
y = main(Avd,c,NITER,28, 1000);
result = ClusteringMeasure_new(Y, y');

