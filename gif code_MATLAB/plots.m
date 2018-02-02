clc 
clear all
itr = 190000;
step = 1000;

% num =0:1250:2500;
k=1;
for i=0:step:itr
C{k} = xlsread(['conc\c_n_',num2str(i), '.csv']);
 A = C{k};
 for j=1:1:70
  for i=1:1:180
      if A(i,j) == 0
          A(i,j) = NaN;
      end
  end
 end
C{k}=A;
k=k+1;
end

gif('Conc.gif','Delay Time',1/500)
for k = 1:1:(itr/step)+1
   %  axis tight
%     hold all
    contourf(C{k}',50,'LineColor','none')
%     contourcbar
%     plot(iden,P{k});
    gif    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itr = 70000;
step = 400;

% figure();
% k=1;
% for i=0:step:itr
% P{k} = xlsread(['pressure\P_n_',num2str(i), '.csv']);
%  B = P{k};
%  for j=1:1:70
%   for i=1:1:180
%       if B(i,j) == 0
%           B(i,j) = NaN;
%       end
%   end
%  end
% P{k}=B;
% k=k+1;
% end
% 
% gif('Press.gif','Delay Time')
% for k = 1:1:(itr/step)+1
%    %  axis tight
% %     hold all
%     contourf(P{k}')
% %     contourcbar
% %     plot(iden,P{k});
%     gif
%     
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure();
% l=1;
% for j=0:step:itr
% u_avg{l} = xlsread(['u\u_avg_',num2str(j), '.csv']);
% l=l+1;
% end
% 
% m=1;
% for z=0:step:itr
% v_avg{m} = xlsread(['v\v_avg_',num2str(z), '.csv']);
% m=m+1;
% end
% 
% u=1;
% for z=0:step:itr
% vel_res{u} = sqrt((u_avg{u}.*u_avg{u})+(v_avg{u}.*v_avg{u}));
% u=u+1;
% end
% 
% gif('Vel.gif','Delay Time')
% for k = 1:1:(itr/step)+1
%    %  axis tight
% %     hold all
%     contourf(vel_res{k}')
% %     contourcbar
% %     plot(iden,P{k});
%     gif
%     
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% k=1;
% for i=0:step:itr
% C{k} = xlsread(['conc\c_n_',num2str(i), '.csv']);
%  A = C{k};
%  for j=1:1:70
%   for i=1:1:180
%       if A(i,j) == 0
%           A(i,j) = NaN;
%       end
%   end
%  end
% C{k}=A;
% k=k+1;
% end

% l=1;
% for j=0:step:itr
% u_avg{l} = xlsread(['u_avg_',num2str(j), '.csv']);
% l=l+1;
% end
% 
% m=1;
% for z=0:step:itr
% v_avg{m} = xlsread(['v_avg_',num2str(z), '.csv']);
% m=m+1;
% end
% 
% u=1;
% for z=0:step:itr
% vel_res{u} = sqrt((u_avg{u}.*u_avg{u})+(v_avg{u}.*v_avg{u}));
% u=u+1;
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




% gif('Press.gif')
% for k = 1:1:(itr/step)+1
%    %  axis tight
%     hold all
%     contourf(C{k}')
%     contourcbar
% %     plot(iden,P{k});
%     gif
%     
% end

% contourf(P{13}')
%  contourcbar