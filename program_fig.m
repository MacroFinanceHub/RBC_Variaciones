%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare RBC04A.mod 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'y','c','innv','lab','w','r','kap','g'};
names  ={'PBI','Consumo', 'Inversión','Empleo','Salario real','Tasa de interés','Capital','Gasto público'};
[nper,junk1] = size(y_e_z);
nvar = length(varble);
dates = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'e_g';
size_shock = 1;
resp_mat1 = [];
for ii=1:nvar
    eval(['y1=',char(varble(ii)),'_',char(shock),';']);
    y1= y1*size_shock;
    y1=[0;y1];
    resp_mat1 = [resp_mat1 y1];
end
save('Model01','resp_mat1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare RBC04B.mod 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'y','c','innv','lab','w','r','kap','g'};
names  ={'PBI','Consumo', 'Inversión','Empleo','Salario real','Tasa de interés','Capital','Gasto público'};
[nper,junk1] = size(y_e_z);
nvar = length(varble);
dates = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'e_g';
size_shock = 1;
resp_mat2 = [];
for ii=1:nvar
    eval(['y2=',char(varble(ii)),'_',char(shock),';']);
    y2= y2*size_shock;
    y2=[0;y2];
    resp_mat2 = [resp_mat2 y2];
end
save('Model02','resp_mat2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare RBC04C.mod 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varble = {'y','c','innv','lab','w','r','kap','g'};
names  ={'PBI','Consumo', 'Inversión','Empleo','Salario real','Tasa de interés','Capital','Gasto público'};
[nper,junk1] = size(y_e_z);
nvar = length(varble);
dates = (0:1:nper)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shock  = 'e_g';
size_shock = 1;
resp_mat3 = [];
for ii=1:nvar
    eval(['y3=',char(varble(ii)),'_',char(shock),';']);
    y3= y3*size_shock;
    y3=[0;y3];
    resp_mat3 = [resp_mat3 y3];
end
save('Model03','resp_mat3');

load Model01
load Model02
load Model03

figure(1)
for ii=1:nvar
    subplot(2,4,ii);
    plot(dates,resp_mat1(:,ii)*100,'-r',dates,resp_mat2(:,ii)*100,'--b',dates,resp_mat3(:,ii)*100,'--g','LineWidth',1.5); 
    grid on; xlim([1 40]);
    hold on; 
    plot([0 40],[0 0],'-k','LineWidth',1.5)
    hold off;
    if ii>8
        xlabel('Trimestres','Fontsize',8)
    end
    ylabel('Desv. % EE','Fontsize',8)
    title(names(ii),'Interpreter','none','Fontsize',10);
   if ii==nvar
        legend('\sigma=1','\sigma=0.5', '\sigma=2');
   end
end
