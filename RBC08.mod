%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Teoría Macrodinámica %%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%    FIEECS - UNI      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo RBC basico loglinealizado (expansion de
// Taylor de 1er orden). Las CPO se presentan de forma no lineal.
// El Estado Estacionario es obtenido manualmente.
// La funcion de utilidad corresponde a la forma logaritmico-lineal.
// Se introduce un ratio de uso de capital en forma de costo. 
// (c) Carlos Rojas Quiroz 

var lab c w r y kap innv u z g;
predetermined_variables kap;
varexo e_z e_g;
parameters alpha delta betta theta phi1 phi2 rho_z rho_g 
z_ss lab_ss r_ss  kap_ss w_ss y_ss c_ss inv_ss u_ss g_ss C_Y I_Y G_Y;

alpha  = 1-0.33;
delta  = 0.023;
betta  = 0.99;
phi1   = 1/betta-1;
phi2   = 1.00;
I_Y    = 0.229298138541349;
rho_z  = 0.95;
rho_g  = 0.75;
z_ss   = 1;
G_Y    = 0.155;
lab_ss = 0.385844633250606;
u_ss   = (delta*(1-alpha)*betta/I_Y-(1-betta*(1+phi1+phi2/2)))/(betta*(delta+phi1+phi2/2));
theta  = 1/(1+(1-lab_ss)/lab_ss*alpha*z_ss/(1-I_Y-G_Y));
y_ss   = z_ss^(1/alpha)*(((1-alpha)*betta/(1-betta+betta*(delta*u_ss+phi1*(u_ss-1)+phi2/2*(u_ss-1)^2)))^((1-alpha)/alpha))*lab_ss;
w_ss   = alpha*y_ss/lab_ss;
kap_ss = (1-alpha)*betta/(1-betta+betta*(delta*u_ss+phi1*(u_ss-1)+phi2/2*(u_ss-1)^2))*y_ss;
inv_ss = I_Y*y_ss;
r_ss   = (1-alpha)*y_ss/(u_ss*kap_ss)-delta;
c_ss   = (1-I_Y-G_Y)*y_ss;
g_ss   = G_Y*y_ss;
C_Y    = c_ss/y_ss;

model;
theta/exp(c) =(1-theta)/((1-exp(lab))*exp(w));
1/exp(c)     =betta*1/exp(c(+1))*(1+exp(r(+1))*exp(u(+1))-(phi1*(exp(u(+1))-1) + phi2/2*(exp(u(+1))-1)^2));
exp(w)       =alpha*exp(y)/exp(lab);
exp(r)+delta =(1-alpha)*exp(y)/(exp(kap)*exp(u));
exp(y)       =exp(c)+exp(innv)+exp(g);
exp(kap(+1)) =(1-delta)*exp(kap)+exp(innv);
exp(y)       =exp(z)*(exp(kap)*exp(u))^(1-alpha)*exp(lab)^alpha;
exp(r)       =phi1 + phi2*(exp(u)-1);
z            =(1-rho_z)*log(z_ss) + rho_z*z(-1) + e_z;
g            =(1-rho_g)*log(g_ss) + rho_g*g(-1) + e_g;
end;

steady_state_model;
lab =log(lab_ss);
c   =log(c_ss); 
w   =log(w_ss); 
r   =log(r_ss); 
y   =log(y_ss); 
kap =log(kap_ss); 
innv=log(inv_ss);
u   =log(u_ss);
z   =log(z_ss);
g   =log(g_ss);
end;

shocks;

var e_z; stderr 0.01;
var e_g; stderr 0.01;
end;
 
resid;
steady;
check;

stoch_simul(order = 1, nograph);