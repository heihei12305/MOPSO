% 标准粒子群(实数编码),求最小值
% author zhaoyuqiang
clear all ;
close all ;
clc ;
N = 100 ; %　种群规模
D = 10 ; % 粒子维度
T = 100 ; %　迭代次数
Xmax = 20 ;
Xmin = -20 ;
C1 = 1.5 ; %　学习因子1
C2 = 1.5 ; %　学习因子2
W = 0.8 ; %　惯性权重
Vmax = 10 ; %　最大飞行速度
Vmin = -10 ; %　最小飞行速度
popx = rand(N,D)*(Xmax-Xmin)+Xmin ; % 初始化粒子群的位置(粒子位置是一个D维向量)
popv = rand(N,D)*(Vmax-Vmin)+Vmin ; % 初始化粒子群的速度(粒子速度是一个D维度向量) 
% 初始化每个历史最优粒子
pBest = popx ; 
pBestValue = func_fitness(pBest) ; 
%初始化全局历史最优粒子
[gBestValue,index] = max(func_fitness(popx)) ;
gBest = popx(index,:) ;
for t=1:T
    for i=1:N
        % 更新个体的位置和速度
        popv(i,:) = W*popv(i,:)+C1*rand*(pBest(i,:)-popx(i,:))+C2*rand*(gBest-popx(i,:)) ;
        popx(i,:) = popx(i,:)+popv(i,:) ;
        % 边界处理，超过定义域范围就取该范围极值
        index = find(popv(i,:)>Vmax | popv(i,:)<Vmin);
        popv(i,index) = rand*(Vmax-Vmin)+Vmin ; %#ok<*FNDSB>
        index = find(popx(i,:)>Xmax | popx(i,:)<Xmin);
        popx(i,index) = rand*(Xmax-Xmin)+Xmin ;
        % 更新粒子历史最优
        if func_fitness(popx(i,:))>pBestValue(i)    
           pBest(i,:) = popx(i,:) ;
           pBestValue(i) = func_fitness(popx(i,:));
        elseif pBestValue(i) > gBestValue
            gBest = pBest(i,:) ;
            gBestValue = pBestValue(i) ;
        end
    end
    % 每代最优解对应的目标函数值
    tBest(t) = func_objValue(gBest); %#ok<*SAGROW>
end
figure
plot(tBest);
xlabel('迭代次数') ;
ylabel('适应度值') ;
title('适应度进化曲线') ;