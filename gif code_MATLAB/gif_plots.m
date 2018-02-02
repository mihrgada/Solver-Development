clc 
close
clear all

% itr = 100000;
% step = 200;

% num =0:1250:2500;

% a=1;
% for i=0:step:itr
% P{a} = xlsread(['pressure\P_n_',num2str(i),'.csv']);
% a = a+1
% end
% 
% b=1;
% for j=0:step:itr
% u_avg{b} = xlsread(['u\u_avg_',num2str(j),'.csv']);
% b=b+1
% end
% 
% c=1;
% for k=0:step:itr
% v_avg{c} = xlsread(['v\v_avg_',num2str(k),'.csv']);
% c=c+1
% end
% 
% d=1;
% for l=0:step:itr
% c_avg{d} = xlsread(['conc\c_n_',num2str(l),'.csv']);
% d=d+1
% end
% 
% u=1;
% for z=0:step:itr
% vel_res{u} = sqrt((u_avg{u}.*u_avg{u})+(v_avg{u}.*v_avg{u}));
% u=u+1
% end


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

% %% Residual Velocity
% gif('Res_vel.gif')
% for k = 1:1:(itr/step)+1
%    %  axis tight
%     hold all
%     contourf(vel_res{k}')
% %     contourcbar
% %     plot(iden,P{k});
%     gif    
% end
% 
% contourf(P{13}')
% contourcbar

%%
ss = 100000;
W = 1e-1; % mm
dx = W/5;
dy = dx;
L = 36*W;
B = 14*W;
M = L/dx;
N = B/dy;

x = linspace(dx/2,L-dx/2,M);
y = linspace(dy/2,B-dy/2,N);
[X,Y] = meshgrid(x,y);

% P = 1e3*xlsread(['pressure\P_n_',num2str(ss),'.csv']);
% for j=1:1:70
%   for i=1:1:180
%       if P(i,j) == 0
%           P(i,j) = NaN;
%       end
%   end
%  end
% % contourf(P')
% figure()
% contourf(X,Y,P')
% contourcbar
% xlabel('x (mm)')
% ylabel('y (mm)')
% title('Pressure Contours (Pa)')
% xlim([0,L])
% ylim([0,B])
% 
% u = xlsread(['u\u_avg_',num2str(ss),'.csv']);
% v = xlsread(['v\v_avg_',num2str(ss),'.csv']);
% v_res = 1e3*sqrt(u.*u + v.*v);
% figure()
% % contourf(X,Y,c',50,'LineColor','none')
% contourf(X,Y,v_res',10,'LineColor','none')
% contourcbar
% xlabel('x (mm)')
% ylabel('y (mm)')
% title('Resultant Velocity (mm/s)')
% xlim([0,L])
% ylim([0,B])

ss = 194400;

c = xlsread(['conc\c_n_',num2str(ss),'.csv']);
for j=1:1:70
  for i=1:1:180
      if c(i,j) == 0
          c(i,j) = NaN;
      end
  end
 end
figure()
contourf(X,Y,c',50,'LineColor','none')
contourcbar
xlabel('x (mm)')
ylabel('y (mm)')
title('Concentration (mol/m^3)')
xlim([0,L])
ylim([0,B])