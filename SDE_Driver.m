clear all
close all
clc

deltaT = 10^-5;
numPaths = 1;
tFin = 5;
length = round(tFin/deltaT);

% Solutions has vectors of solutions
tic
[Time, Paths] = SDESolver(deltaT, numPaths, tFin);
toc
