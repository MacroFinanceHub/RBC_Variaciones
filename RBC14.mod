%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Teoría Macrodinámica %%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%    FIEECS - UNI      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo RBC basico loglinealizado (expansion de
// Taylor de 1er orden). Las CPO se presentan de forma no lineal.
// El Estado Estacionario es obtenido manualmente.
// La funcion de utilidad corresponde a la forma logaritmico-lineal.
// Se incluye dinero con una restricción Cash-in-Advance
// (c) Carlos Rojas Quiroz 

var lab c w r y kap innv z g lambda mu m inom pic;
predetermined_variables kap;
varexo e_z e_g e_m;
parameters alpha delta betta theta rho_z rho_g rho_m
z_ss lab_ss r_ss kap_ss w_ss y_ss c_ss inv_ss g_ss pi_star pic_ss lambda_ss mu_ss m_ss inom_ss C_Y I_Y G_Y;

alpha       = 1-0.33;
delta       = 0.023;
betta       = 0.99;
theta       = 1/2.75;
rho_z       = 0.95;
rho_g       = 0.75;
rho_m       = 0.75;
z_ss        = 1;
G_Y         = 0.155;
pi_star     = 0.00;
pic_ss      = pi_star;
lab_ss      = 1/((1-theta)*(1+pic_ss)/(alpha*theta*betta)*((1-betta+alpha*betta*delta)/(1-betta+betta*delta)-G_Y)+1);
y_ss        = z_ss*(((1-alpha)*betta/(1-betta+betta*delta))^((1-alpha)/alpha))*lab_ss;
w_ss        = alpha*y_ss/lab_ss;
kap_ss      = (1-alpha)*betta/(1-betta+betta*delta)*y_ss;
inv_ss      = delta*kap_ss;
r_ss        = (1-alpha)*y_ss/kap_ss-delta;
c_ss        = ((1-betta+alpha*betta*delta)/(1-betta+betta*delta)-G_Y)*y_ss;
g_ss        = G_Y*y_ss;
C_Y         = c_ss/y_ss;
I_Y         = inv_ss/y_ss;
mu_ss       = theta/c_ss*(1+pic_ss-betta)/(1+pic_ss);
lambda_ss   = theta/c_ss-mu_ss;
m_ss        = c_ss;
inom_ss     = (1/betta)*(1+pic_ss)-1;

model;
exp(m)          =exp(c);
theta/exp(c)    =exp(lambda)+exp(mu); 
exp(lambda)*exp(w)=(1-theta)/(1-exp(lab));
exp(lambda)     =betta*exp(lambda(+1))*(1+exp(r(+1)));
exp(lambda)     =betta*exp(lambda(+1))*(1+exp(inom))/(1+pic(+1));
exp(lambda)     =betta/(1+pic(+1))*(exp(mu(+1))+exp(lambda(+1)));
exp(y)          =exp(c)+exp(innv)+exp(g);
exp(kap(+1))    =(1-delta)*exp(kap)+exp(innv);
exp(r)          =(1-alpha)*exp(y)/exp(kap)-delta;
exp(w)          =alpha*exp(y)/exp(lab);
exp(y)          =exp(z)*exp(kap)^(1-alpha)*exp(lab)^alpha;
z               =(1-rho_z)*log(z_ss) + rho_z*z(-1) + e_z;
g               =(1-rho_g)*log(g_ss) + rho_g*g(-1) + e_g;
m-m(-1)         =(1-rho_m)*pi_star-pic+rho_m*pic(-1)+rho_m*(m(-1)-m(-2)) + e_m;
end;

steady_state_model;
lambda  =log(lambda_ss);
mu      =log(mu_ss);
m       =log(m_ss);
lab     =log(lab_ss);
c       =log(c_ss); 
w       =log(w_ss); 
r       =log(r_ss); 
y       =log(y_ss); 
kap     =log(kap_ss); 
innv    =log(inv_ss);
z       =log(z_ss);
g       =log(g_ss);
pic     =pic_ss;
inom    =log(inom_ss);
end;

shocks;

var e_z; stderr 0.01;
var e_g; stderr 0.01;
var e_m; stderr 0.01;
end;
 
resid;
steady;
check;

stoch_simul(order = 1, nograph);
