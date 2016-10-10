%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Teoría Macrodinámica %%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%    FIEECS - UNI      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo RBC basico loglinealizado (expansion de
// Taylor de 1er orden). Las CPO se presentan de forma no lineal.
// El Estado Estacionario es obtenido manualmente.
// La funcion de utilidad corresponde a la forma logaritmico-lineal.
// Incluye impuestos distorsivos al consumo, al trabajo y al capital.
// (c) Carlos Rojas Quiroz 

var lab c w r y kap innv z g tauL tauK;
predetermined_variables kap;
varexo e_z e_g e_tauL e_tauK;
parameters alpha delta betta theta rho_z rho_g rho_tauL rho_tauK
z_ss lab_ss r_ss kap_ss w_ss y_ss c_ss inv_ss g_ss tauL_ss tauK_ss C_Y I_Y G_Y;

alpha   = 1-0.33;
delta   = 0.023;
betta   = 0.99;
theta   = 1/2.75;
rho_z   = 0.95;
rho_g   = 0.75;
rho_tauL= 0.90;
rho_tauK= 0.90;
z_ss    = 1;
G_Y     = 0.155;
lab_ss  = 1/((1-theta)/(alpha*theta*z_ss)*((1-betta+alpha*betta*delta)/(1-betta+betta*delta)-G_Y)+1);
y_ss    = z_ss*(((1-alpha)*betta/(1-betta+betta*delta))^((1-alpha)/alpha))*lab_ss;
w_ss    = alpha*y_ss/lab_ss;
kap_ss  = (1-alpha)*betta/(1-betta+betta*delta)*y_ss;
inv_ss  = delta*kap_ss;
r_ss    = (1-alpha)*y_ss/kap_ss-delta;
c_ss    = ((1-betta+alpha*betta*delta)/(1-betta+betta*delta)-G_Y)*y_ss;
g_ss    = G_Y*y_ss;
C_Y     = c_ss/y_ss;
I_Y     = inv_ss/y_ss;
tauL_ss = 0.2;
tauK_ss = 0.1;

model;
theta/exp(c) =(1-theta)/((1-exp(lab))*(1-exp(tauL))*exp(w));
1/exp(c)     =betta*1/exp(c(+1))*(1+exp(r(+1)))*(1-exp(tauK(+1)));
exp(w)       =alpha*exp(y)/exp(lab);
exp(r)+delta =(1-alpha)*exp(y)/exp(kap);
exp(y)       =exp(c)+exp(innv)+exp(g);
exp(kap(+1)) =(1-delta)*exp(kap)+exp(innv);
exp(y)       =exp(z)*exp(kap)^(1-alpha)*exp(lab)^alpha;
z            =(1-rho_z)*log(z_ss) + rho_z*z(-1) + e_z;
g            =(1-rho_g)*log(g_ss) + rho_g*g(-1) + e_g;
tauL         =(1-rho_tauL)*log(tauL_ss) + rho_tauL*tauL(-1) + e_tauL;
tauK         =(1-rho_tauK)*log(tauK_ss) + rho_tauK*tauK(-1) + e_tauK;
end;

steady_state_model;
lab =log(lab_ss);
c   =log(c_ss); 
w   =log(w_ss); 
r   =log(r_ss); 
y   =log(y_ss); 
kap =log(kap_ss); 
innv=log(inv_ss); 
z   =log(z_ss);
g   =log(g_ss);
tauL=log(tauL_ss);
tauK=log(tauK_ss);
end;

shocks;

var e_z; stderr 0.01;
var e_g; stderr 0.01;
var e_tauL; stderr 0.01;
var e_tauK; stderr 0.01;
end;
 
resid;
steady;
check;

stoch_simul(order = 1, nograph);
