%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Teoría Macrodinámica %%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%    FIEECS - UNI      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo RBC basico loglinealizado (expansion de
// Taylor de 1er orden). Las CPO se presentan de forma no lineal.
// El Estado Estacionario es obtenido manualmente.
// La funcion de utilidad corresponde a la forma logaritmico-lineal.
// Incluye un choque de mark-up derivado de la presencia de un sector en
// competencia imperfecta. 
// (c) Carlos Rojas Quiroz 

var lab c w r y kap innv z g psi;
predetermined_variables kap;
varexo e_z e_g e_psi;
parameters alpha delta betta theta nu rho_z rho_g rho_psi
z_ss lab_ss r_ss kap_ss w_ss y_ss c_ss inv_ss g_ss psi_ss C_Y I_Y G_Y;

alpha  = 1-0.33;
delta  = 0.023;
betta  = 0.99;
theta  = 1/2.75;
rho_z  = 0.95;
rho_g  = 0.75;
rho_psi= 0.50;
z_ss   = 1;
nu     = 6;
psi_ss = nu/(nu-1);
G_Y    = 0.155;
lab_ss = 1/(1+(1-delta*betta*(1-alpha)*(nu-1)/((1-betta+betta*delta)*nu)-G_Y)*((1-theta)/(alpha*theta))*(nu/(nu-1)));
y_ss   = z_ss*(((1-alpha)*betta*(nu-1)/((1-betta+betta*delta)*nu))^((1-alpha)/alpha))*lab_ss;
w_ss   = (nu-1)/nu*alpha*y_ss/lab_ss;
kap_ss = (1-alpha)*betta/(1-betta+betta*delta)*y_ss*(nu-1)/nu;
inv_ss = delta*kap_ss;
r_ss   = (nu-1)/nu*(1-alpha)*y_ss/kap_ss-delta;
c_ss   = (1-delta*betta*(1-alpha)*(nu-1)/((1-betta+betta*delta)*nu)-G_Y)*y_ss;
g_ss   = G_Y*y_ss;
C_Y    = c_ss/y_ss;
I_Y    = inv_ss/y_ss;

model;
theta/exp(c) =(1-theta)/((1-exp(lab))*exp(w));
1/exp(c)     =betta*1/exp(c(+1))*(1+exp(r(+1)));
exp(w)       =1/exp(psi)*alpha*exp(y)/exp(lab);
exp(r)+delta =1/exp(psi)*(1-alpha)*exp(y)/exp(kap);
exp(y)       =exp(c)+exp(innv)+exp(g);
exp(kap(+1)) =(1-delta)*exp(kap)+exp(innv);
exp(y)       =exp(z)*exp(kap)^(1-alpha)*exp(lab)^alpha;
z            =(1-rho_z)*log(z_ss) + rho_z*z(-1) + e_z;
g            =(1-rho_g)*log(g_ss) + rho_g*g(-1) + e_g;
psi          =(1-rho_psi)*log(psi_ss) + rho_psi*psi(-1) + e_psi;
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
psi =log(psi_ss);
end;

shocks; 

var e_z; stderr 0.01;
var e_g; stderr 0.01;
var e_psi; stderr 0.01;
end;
  
resid;
steady;
check;

stoch_simul(order = 1, nograph);
