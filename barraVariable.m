##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------
## _______     ___       __   __   __  |
##|   ____|   /   \     |  | |  | |  | | CATEDRA ESTRUCTURAS AERONAUTICAS III
##|  |__     /  ^  \    |  | |  | |  | |
##|   __|   /  /_\  \   |  | |  | |  | | METODO DE RAYLEIGH-RITZ (BARRA VARIABLE) 
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

######## DATOS DE LA BARRA #####################

E=1;
Ao=1;
p=100;
L=10;



CABAR=[Ao*(1-x/(2*L)),0,E,0,L]; ## Ax,NoUtil,E,Xi,Xf

######## RESORTES #################################

KL=[0,0]; # K,X

######## CARGAS APLICADAS #########################


Qdis=[p,0,L]; ## Qx,Xi,Xf

P=[0,0]; ## CARGA,X


######## POLINOMIO DE INTERPOLACION ###############

function d=forN(grados)
 
  d=sym(zeros(1,grados));
  
  for i=1:grados

    syms x;
    
    d(i)=x^i; # Aca se carga el polinomio de interpolacion
       
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

xInicial=double(CABAR(1,4));

xFinal=double(CABAR(end,5));

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

  Kb=zeros(gl,gl);

  Kl=Kb;Kt=Kb;
  
  fq=zeros(gl,1);

  fm=zeros(gl,1);

  fp=zeros(gl,1);
  
######### VECTORES y MATRICES ############################

### ENERGIAS DE DEFORMACION ###
  
  for u=1:size(CABAR,1)

    Aactual=CABAR(u,1);
    
    Eactual=CABAR(u,3);

    XiActual=double(CABAR(u,4));

    XfActual=double(CABAR(u,5));

    KNX=DNT*Eactual*Aactual*DN;

    Kb=Kb+[integrador(KNX,XiActual,XfActual)];  
    
  endfor

   
  for u=1:size(KL,1)

    K=KL(u,1);
    
    X=KL(u,2);

    Kl=function_handle(NT*K*N);

    Kb=Kb+Kl(X);  
    
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
    
  
######### RESOLUCION ######################################

  C=inv(Kb)*(fq+fp+fm);

######### PLOTEO #########################################
  
  RRf=function_handle(N*C);
  
  AreaActual=trapz(abscisas,RRf(abscisas));

  try
  
    AreaPrevia=trapz(abscisas,rr);

  catch

    AreaPrevia=2*AreaActual;

  end_try_catch
  
  Error=abs((AreaPrevia-AreaActual)/AreaPrevia)*100;
  
  if  (Error<=CRITERIO || gl>=maxIteraciones)

    CONVG=1;

  endif

  RESULTADO=N*C;
  
  rr=RRf(abscisas);
  
  plot(abscisas,rr,["--" markStyle(gl) color(gl) ";V(X) GL= " num2str(gl) ";"]);

  clear Kb fp fq fm N NT DN DNT D2N D2NT;

  N=sym(ones(1,1));

  gl++;
  
endwhile

hold off;

printf("\n\n\n RR %d GL -->  ERROR RELATIVO= %d%%\n",(gl-1),Error);

rrSymbolico=RESULTADO;

## ############# DIAGRAMA DE CARACTERISTICAS ###################


## for u=1:size(CABAR,1)

##   Aactual=CABAR(u,1);
  
##   Eactual=CABAR(u,3);

##   absAct=abscisas;

  
##   Mx=function_handle(-Aactual*Eactual*diff(RESULTADO,x,2));

##   Qx=function_handle(Aactual*Eactual*diff(RESULTADO,x,3));

##   Sx=function_handle(-Eactual*diff(RESULTADO,x,2));
  
  
##   XiActual=double(CABAR(u,4));

##   XfActual=double(CABAR(u,5));
  
##   condiciones=absAct(:,:)<XiActual | absAct(:,:)>=XfActual;
  
##   absAct(:,condiciones)=[];
  
##   try

##     mx=[mx,Mx(absAct)];

##     qx=[qx,Qx(absAct)];

##     sx=[sx,Sx(absAct)];

##   catch

    
##     mx=[Mx(absAct)];

##     qx=[Qx(absAct)];

##     sx=[Sx(absAct)];

##   end_try_catch  

## endfor


## mx=[mx,Mx(XfActual)];

## qx=[qx,Qx(XfActual)];

## sx=[sx,Sx(XfActual)];

## figure (2);clf;hold on;grid on;

## title ('DIAGRAMA DE CARACTERISTICAS')

## plot(abscisas,mx,["--" markStyle(3) color(3) ";Mx;"]);

## plot(abscisas,qx,["--" markStyle(4) color(4) ";Qx;"]);

## hold off




