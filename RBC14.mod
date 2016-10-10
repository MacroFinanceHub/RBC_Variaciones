%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Teoría Macrodinámica %%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%    FIEECS - UNI      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo RBC basico loglinealizado (expansion de
// Taylor de 1er orden). Las CPO se presentan de forma no lineal.
// El Estado Estacionario es obtenido manualmente.
// La funcion de utilidad corresponde a la forma logaritmico-lineal.
// Se incluye trabajo indivisble.
// Se incluye dinero en la función de utilidad.
// Se replica el modelo del capítulo 9 del McCandless.
// (c) Carlos Rojas Quiroz 

var r w c kap lab y M z g p m;
predetermined_variables kap;
varexo e_z e_g;
parameters betta delta theta B D rho_z rho_g g_ss z_ss r_ss w_ss c_ss kap_ss lab_ss y_ss m_ss p_ss;
betta=0.99;
delta=0.025;
theta=0.36;
B=-2.5805;
rho_z=0.95;
rho_g=0.50;
g_ss=1;
z_ss=1;
D=0.01;
r_ss=1/betta-1+delta;
w_ss=(1-theta)*(theta/r_ss)^(theta/(1-theta));
c_ss=-w_ss/B;
kap_ss=c_ss/((r_ss*(1-theta)/(w_ss*theta))^(1-theta)-delta);
lab_ss=r_ss*(1-theta)*kap_ss/(w_ss*theta);
y_ss=c_ss+delta*kap_ss;
m_ss=D*g_ss*c_ss/(g_ss-betta);
p_ss=1;

model;
1/exp(c)=betta*exp(p)/(exp(c(+1))*exp(p(+1)))+D*exp(p)/exp(M);
1/exp(c)=betta*1/exp(c(+1))*(exp(r(+1))+1-delta);
1/exp(c)=-B/exp(w);
exp(y)=exp(z)*exp(kap)^theta*exp(lab)^(1-theta);
exp(kap(+1))+exp(M)/exp(p)+exp(c)=(1-delta)*exp(kap)+exp(w)*exp(lab)+exp(r)*exp(kap)+exp(M(-1))/exp(p)+(exp(g)-1)*exp(M)/exp(p); 
exp(r)=theta*exp(z)*exp(kap)^(theta-1)*exp(lab)^(1-theta);
exp(w)=(1-theta)*exp(z)*exp(kap)^(theta)*exp(lab)^(-theta);
z=(1-rho_z)*log(z_ss)+rho_z*z(-1)+e_z;
g=(1-rho_g)*log(g_ss)+rho_g*g(-1)+e_g;
exp(M)=exp(g)*exp(M(-1));
exp(m)=exp(M)/exp(p);
end;

steady_state_model;
r   =log(r_ss); 
w   =log(w_ss); 
c   =log(c_ss); 
kap =log(kap_ss); 
lab =log(lab_ss);
p   =log(p_ss); 
y   =log(y_ss); 
M   =log(m_ss); 
z   =log(z_ss); 
g   =log(g_ss);
m   =log(m_ss/p_ss);
end;

shocks;
var e_z; stderr 0.01;
var e_g; stderr 0.01;
end;
 
resid;
steady;
check;

stoch_simul(order = 1, nograph);
