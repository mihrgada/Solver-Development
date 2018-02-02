clear 
close all
clc 


ug = xlsread('u_g.csv');
vg = xlsread('v_g.csv');
Pg = xlsread('P_g.csv');
% X = xlsread('X.csv');
% Y = xlsread('Y.csv');
% uav = xlsread('u_avg.csv');
% vav = xlsread('v_avg.csv');


for j = 1:1:61
    for i = 12:1:351  
    
    Pg(i,j) = NaN;
    vg(i,j) = NaN;
    ug(i,j) = NaN;
    
    end
end

for j = 82:1:142
    for i = 12:1:351  
    
    Pg(i,j) = NaN;
    vg(i,j) = NaN;
    ug(i,j) = NaN;
    
    end
end

for j = 1:141:142
    for i = 1:1:362 
    
    Pg(i,j) = NaN;
    vg(i,j) = NaN;
    ug(i,j) = NaN;
    
    end
end

for j = 1:1:142
    for i = 1:361:362 
    
    Pg(i,j) = NaN;
    vg(i,j) = NaN;
    ug(i,j) = NaN;
    
    end
end
figure()
%contourf(X',Y',P')
% contourf(X',Y',uav')
contourf(ug',20)
% ylim([0,3.6])
title('ug')
contourcbar

figure()
%contourf(X',Y',P')
% contourf(X',Y',uav')
contourf(vg',20)
% ylim([0,3.6])
title('vg')
contourcbar

figure()
%contourf(X',Y',P')
% contourf(X',Y',uav')
contourf(Pg',20)
% ylim([0,3.6])
title('Pg')
contourcbar