clear;close all;clc;
addpath(genpath('./'));
dbstop if error;
%% Nausea
options.dataType = 1;
dataFull = generateRealTest(options);
load('RealNausea.mat');
fprintf('\nNausea\n');
for i = 1:2
    res = resFull{i};
    res_rev = resFull_rev{i};
    data = dataFull{i};
    % printStat(data,res);
    % printStat(data,res_rev);
    printStat(data,res,res_rev);
end

%% Bladder
options.dataType = 2;
dataFull = generateRealTest(options);
load('RealBladder.mat');
fprintf('\nBladder\n');
for i = 1:2
    res = resFull{i};
    res_rev = resFull_rev{i};
    data = dataFull{i};
    printStat(data,res,res_rev);
end

%% Skin cancer
options.dataType = 3;
dataFull = generateRealTest(options);
load('RealSkin.mat');
fprintf('\nSkin\n');
for i = 1:4
    res = resFull{i};
    res_rev = resFull_rev{i};
    data = dataFull{i};
    printStat(data,res,res_rev);
end
