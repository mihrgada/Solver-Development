clc 
clear all
itr = 1300

% num =0:1250:2500;
k=1;
for i=0:100:itr
P{k} = xlsread(['P_n_',num2str(i), '.csv']);k=k+1;
end

l=1;
for j=0:100:itr
u_avg{l} = xlsread(['u_avg_',num2str(j), '.csv']);
l=l+1;
end

m=1;
for z=0:100:itr
v_avg{m} = xlsread(['v_avg_',num2str(z), '.csv']);
m=m+1;
end

u=1;
for z=0:100:itr
vel_res{u} = sqrt((u_avg{u}.*u_avg{u})+(v_avg{u}.*v_avg{u}));
u=u+1;
end


% Y=1:1:360;
% Z=1:1:140;
% for Y=1:1:360
%     for Z=1:1:140
% %    iden(Y,Z)=     
%     end
% end
% iden = ones(360,140)
% A=P{1};
% contourf(A)
%plot(P{1}(:,1),P{1}(:,2));


gif('Presuure.gif','DelayTime',1)
for k = 1:1:(itr/100)+1
   %  axis tight
    hold all
    contourf(vel_res{k}')
    contourcbar
%     plot(iden,P{k});
    gif
    
end
% contourf(P{13}')
%  contourcbar