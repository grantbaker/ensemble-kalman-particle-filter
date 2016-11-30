clear all
close all
clc

deltaT = 10^-5;
numPaths = 1;
tFin = 5;
length = round(tFin/deltaT);
InitialCond = [1, 0, 0];

% Solutions has vectors of solutions
tic
[Time, Paths] = SDESolver(deltaT, numPaths, tFin, InitialCond);
toc
