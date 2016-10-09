%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Teoría Macrodinámica %%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%    FIEECS - UNI      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo RBC basico loglinealizado (expansion de
// Taylor de 1er orden). Las CPO se presentan de forma no lineal.
// El Estado Estacionario es obtenido manualmente.
// La funcion de utilidad corresponde a la forma logaritmico-lineal.
// Se incluye trabajo indivisible.
// Se incluye una restricción Cash-in-Advance.
// Se replica el modelo del capítulo 8 del McCandless (sin señoreaje)
// (c) Carlos Rojas Quiroz 

var r w c p lab kap m y g z pic;
predetermined_variables kap;
varexo e_z e_g;
parameters betta delta theta B rho_z rho_g r_ss w_ss c_ss p_ss lab_ss kap_ss y_ss g_ss z_ss pic_ss m_ss;
betta  =0.990;
delta  =0.025; 
theta  =0.360;
B      =-2.5805;
rho_z  =0.750;
rho_g  =0.500;
g_ss   =1.200;
r_ss   =1/betta-(1-delta);
w_ss   =(1-theta)*(r_ss/theta)^(theta/(theta-1));
c_ss   =-betta*w_ss/(g_ss*B);
p_ss   =1/c_ss;
kap_ss =c_ss/(r_ss/theta-delta);
lab_ss =(r_ss/theta)^(1/(1-theta))*kap_ss;
y_ss   =c_ss+delta*kap_ss;
z_ss   =1;
pic_ss =1;
m_ss   =1;

model;
1/betta=exp(w)/exp(w(+1))*(1-delta+exp(r(+1)));
B/(exp(w)*exp(p))=-betta/(exp(p(+1))*exp(c(+1))*exp(g(+1)));
exp(p)*exp(c)=(exp(m(-1))+exp(g)-1)/exp(g);
exp(kap(+1))+exp(m)/exp(p)=(1-delta)*exp(kap)+exp(w)*exp(lab)+exp(r)*exp(kap); 
exp(w)=(1-theta)*exp(z)*(exp(kap)/exp(lab))^theta;
exp(r)=theta*exp(z)*(exp(kap)/exp(lab))^(theta-1);
exp(y)=exp(c)+delta*exp(kap);
exp(m)=1;
g=(1-rho_g)*log(g_ss)+rho_g*g(-1)+e_g;
z=(1-rho_z)*log(z_ss)+rho_z*z(-1)+e_z;
exp(pic)=exp(p)/exp(p(-1));
end;

steady_state_model;
r   =log(r_ss); 
w   =log(w_ss);
c   =log(c_ss);
p   =log(p_ss);
lab =log(lab_ss); 
kap =log(kap_ss); 
y   =log(y_ss);
g   =log(g_ss);
z   =log(z_ss);
pic =log(pic_ss);
m   =log(m_ss);
end;

shocks;
var e_z; stderr 0.01;
var e_g; stderr 0.01;
end;
 
resid;
steady;
check;

stoch_simul(order = 1, nograph);
