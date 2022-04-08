function [ dudt ] = ddt_heat(t, u)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    global A B D

    dudt = D * (A * u + B);
    
end

