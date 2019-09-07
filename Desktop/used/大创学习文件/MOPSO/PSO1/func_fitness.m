function [ result ] = func_fitness( pop )
%OBJFUNCTION 求适应度，最小值
%   待优化目标函数
% x:　种群或者个体
% result : 种群适应度
objValue =  func_objValue(pop);
result  = 4001 - objValue ;
end

