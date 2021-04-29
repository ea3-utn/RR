##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------
## _______     ___       __   __   __  |
##|   ____|   /   \     |  | |  | |  | | CATEDRA ESTRUCTURAS AERONAUTICAS III
##|  |__     /  ^  \    |  | |  | |  | |
##|   __|   /  /_\  \   |  | |  | |  | | METODO DE RAYLEIGH-RITZ (VIGA DE NAVIER) 
##|  |____ /  _____  \  |  | |  | |  | |
##|_______/__/     \__\ |__| |__| |__| | RR: script principal
##                                     |
##---------CICLO LECTIVO 2020----------------------------------------------------

## CONFIGURACION

clear;

pkg load symbolic; # Carga de paquete que me permite hacer operaciones algebraicas

warning ('off','OctSymPy:sym:rationalapprox');

## DECLARACIONES

syms x;

markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

color=["k","r","g","b","m","c","k","r","g","b","m","c"];

######### FUNCIONES ################

function d=D(z)
  
  syms x;

  d=diff(z,x,1);

endfunction

function d=D2(z)
  
  syms x;

  d=diff(z,x,2);

endfunction

######## DATOS DE LA VIGA #####################

RelA=.75;

E=7.311924e9; #Kgf/m2

Ixx=10*6.19768592718e-4;

xca=16.23;

CAVIG=[.5*Ixx,0,E,0,5.07;Ixx,0,E,5.07,33.02-4.5;.5*Ixx,0,E,33.02-4.5,33.02]; ## Ix,Ax,E,Xi,Xf

######## RESORTES #################################

KL=[0,0]; # K,X

KT=[0,0]; # K,X

######## CARGAS APLICADAS #########################

[w0,Nl,Pm,Tal,Mt,Mct,Macn,Nt,xca,Lf]=equilibrioEstatico();

q=(w0)*sin(pi*x/Lf);

Qdis=[-q,0,Lf]; ## Qx,Xi,Xf

P=[-Pm,xca;Nl,xca;Nt,Lf]; ## CARGA,X

M=[Mt,xca;-Mct,xca;Macn,xca];

######## POLINOMIO DE INTERPOLACION ###############

function d=forN(grados)
 
  d=sym(zeros(1,grados));
  
  for i=1:grados

    syms x;
    
    d(i)=(1000*x^(i+1)+1000+1*x); # Aca se carga el polinomio de interpolacion
       
  endfor
  
endfunction


##//////// CONFIGURACION //////////////////////////////

gl=1 ;# Grados de libertad iniciales 

gli=gl;

maxIteraciones=12;

CRITERIO=5; # Criterio de convergencia en %

##//////// DECLARACIONES //////////////////////////////

CONVG=0;

N=sym(ones(1,1));

xInicial=double(CAVIG(1,4));

xFinal=double(CAVIG(end,5));

abscisas=linspace(xInicial,xFinal,50);

#############################################################
######### VIGA DE NAVIER ####################################
#############################################################

figure(1)

clf;

hold on;

grid on;

title ('APROXIMACION METODO DE RAYLEIGHT-RITZ')

while (CONVG<=0 && gl<=maxIteraciones)
  
  
  #### DEFINICIONES PREVIAS #####
  
  N=forN(gl);

  NT=symTr(N);

  DN=D(N);

  DNT=symTr(DN);

  D2N=D2(N);

  D2NT=symTr(D2N);

  Kv=zeros(gl,gl);

  Kl=Kv;Kt=Kv;
  
  fq=zeros(gl,1);

  fm=zeros(gl,1);

  fp=zeros(gl,1);
  
######### VECTORES y MATRICES ############################

### ENERGIAS DE DEFORMACION ###
  
  for u=1:size(CAVIG,1)

    Iactual=CAVIG(u,1);
    
    Eactual=CAVIG(u,3);

    XiActual=double(CAVIG(u,4));

    XfActual=double(CAVIG(u,5));

    KNX=D2NT*Eactual*Iactual*D2N;

    Kv=Kv+[integrador(KNX,XiActual,XfActual)];  
    
  endfor

   
  for u=1:size(KL,1)

    K=KL(u,1);
    
    X=KL(u,2);

    Kl=function_handle(NT*K*N);

    Kv=Kv+Kl(X);  
    
  endfor

  for u=1:size(KT,1)

    K=KT(u,1);
    
    X=KT(u,2);

    Kt=function_handle(DNT*K*DN);

    Kv=Kv+Kt(X);  
    
  endfor
  
### TRABAJO DE LAS CARGAS ###
  
  for u=1:size(Qdis,1)  ## CARGA DISTRIBUIDA

    XiActual=double(Qdis(u,2));

    XfActual=double(Qdis(u,3));

    QN=NT*Qdis(u,1);

    fq=fq+[integrador(QN,XiActual,XfActual)];  

  endfor

  for u=1:size(P,1)  ## CARGA 

    X=P(u,2);

    Pu=function_handle(NT);

    fp=fp+Pu(X)*P(u,1);  

  endfor

  for u=1:size(M,1)  ## MOMENTO 

    X=M(u,2);

    Mu=function_handle(DNT);

    fm=fm+Mu(X)*-M(u,1);  

  endfor
    
  
######### RESOLUCION ######################################

  C=inv(Kv)*(fq+fp+fm);

######### PLOTEO #########################################
  
  RRf=function_handle(N*C);

  Norma=function_handle(abs(N*C));
    
  AreaActual=trapz(abscisas,Norma(abscisas));

  try
  
    AreaPrevia=trapz(abscisas,NormaPrevia);

  catch

    AreaPrevia=2*AreaActual;

  end_try_catch

  NormaPrevia=Norma(abscisas);
  
  Error=abs((AreaPrevia-AreaActual)/AreaPrevia)*100;
  
  if  (Error<=CRITERIO || gl>=maxIteraciones)

    CONVG=1;

  endif

  RESULTADO=N*C;
  
  rr=RRf(abscisas);
  
  plot(abscisas,rr,["--" markStyle(gl) color(gl) ";V(X) GL= " num2str(gl) ";"]);

  clear Kv fp fq fm N NT DN DNT D2N D2NT;

  N=sym(ones(1,1));

  gl++;
  
endwhile

hold off;

printf("\n\n\n RR %d GL -->  ERROR RELATIVO= %d%%\n",(gl-1),Error);

rrSymbolico=RESULTADO;

############# DIAGRAMA DE CARACTERISTICAS ###################


for u=1:size(CAVIG,1)

  Iactual=CAVIG(u,1);
  
  Eactual=CAVIG(u,3);

  absAct=abscisas;

  
  Mx=function_handle(-Iactual*Eactual*diff(RESULTADO,x,2));

  Qx=function_handle(Iactual*Eactual*diff(RESULTADO,x,3));

  Sx=function_handle(-Eactual*diff(RESULTADO,x,2));
  
  
  XiActual=double(CAVIG(u,4));

  XfActual=double(CAVIG(u,5));
  
  condiciones=absAct(:,:)<XiActual | absAct(:,:)>=XfActual;
  
  absAct(:,condiciones)=[];
  
  try

    mx=[mx,Mx(absAct)];

    qx=[qx,Qx(absAct)];

    sx=[sx,Sx(absAct)];

  catch

    
    mx=[Mx(absAct)];

    qx=[Qx(absAct)];

    sx=[Sx(absAct)];

  end_try_catch  

endfor


mx=[mx,Mx(XfActual)];

qx=[qx,Qx(XfActual)];

sx=[sx,Sx(XfActual)];

figure (2);clf;hold on;grid on;

title ('DIAGRAMA DE CARACTERISTICAS')

plot(abscisas,mx,["--" markStyle(3) color(3) ";Mx;"]);

plot(abscisas,qx,["--" markStyle(4) color(4) ";Qx;"]);

hold off




