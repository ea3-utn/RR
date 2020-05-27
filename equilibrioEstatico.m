##-----UTN FACULTAD REGIONAL HAEDO----------------* - Octave - *-----------------
## _______     ___       __   __   __  |
##|   ____|   /   \     |  | |  | |  | | CATEDRA ESTRUCTURAS AERONAUTICAS III
##|  |__     /  ^  \    |  | |  | |  | |
##|   __|   /  /_\  \   |  | |  | |  | |    METODO DE RAYLEIGH-RITZ
##|  |____ /  _____  \  |  | |  | |  | |
##|_______/__/     \__\ |__| |__| |__| |  DLR ASTRA A320 - CASO DE ESTUDIO TC1000 
##                                     |
##----------Unidades en MKS-----------------------------------------------------

function [w0,Nl,Pm,Tal,Mt,Mct,Mca,Nt,xca,Lf]=equilibrioEstatico()

  ## DECLARACIONES

  clear;

  syms x Mac deltaW;

  markStyle=["+","o","*",".","x","s","d","^","v",">","<","p","h"];

  color=["k","r","g","b","m","c","k","r","g","b","m","c"];


######### INGRESO DE DATOS

  g=9.81; # Aceleracion de la gravedad

  ## Geometria de la Aeronave

  Lf=33.02; # Longitud del fuselaje

  hca=1.46; # Distancia en Y --> CA,EJE NEUTRO

  hm=1.36; # Distancia en Y --> eje traccion,CA

  xca=16.23; # Centro aerodinamico desde origen

  Sref=122.4; # Superficie alar de referencia

  ## Parametros aerodinamicos

  alfa=10; # Angulo de ataque en grados

  V=67.166; # Velocidad

  rho=0.9555; # Densidad

  Cl=2; # # "Aerodinamic design of transport aircraft"

  CdClRatio=.2/.8; # Landing Config. "Nummerical max lift predictions ..."

  ## Cargas

  Pmotor=2404; # Peso individual motor

  CantidadMotores=2;

  AGW=50200; # Aircraft gross weight

########### Calculos auxiliares

  Cd=Cl*CdClRatio;

  alpha=alfa*pi/180;

  L=.5*rho*(V^2)*Sref*Cl/g;

  D=.5*rho*(V^2)*Sref*Cd/g;

  Nl=L*cos(alpha)+D*sin(alpha);

  Ct=D*cos(alpha)-L*sin(alpha);

  T=D/(cos(alpha));

  Pm=Pmotor*CantidadMotores*cos(alpha);

######### CARGA DISTRIBUIDA

  AGWestruc=AGW*cos(alpha)-Pm;

				#q=AGWestruc/Lf;

  w0=AGWestruc*pi/(2*Lf);

  q=(w0)*sin(pi*x/Lf);
  
######### EQUILIBRIO DE FUERZAS

  Tal=T*sin(alpha);
  
  Lt=AGW-L-Tal;  # Por sumatoria de fuerzas en Y

  Nt=Lt*cos(alpha);

######### EQUILIBRIO DE MOMENTOS en Y

  Mt=T*(hm+hca); # Momento propulsor

  Mct=Ct*hca;

  Mnt=Nt*(Lf-xca);

  Mq=int(q*(xca-x),x,0,Lf);

  Mca=Mt-Mct-Mnt-Mq+Mac;

  [Macn]=solve(Mca==0,Mac);

  Mca=double(Macn);
  
  ## MOMENTO FUNCION DE X

  M0A=function_handle(-int(q*x,x));

  MAF=function_handle(-int(q*x,x)+Macn+Mt-Mct+(Nl-Pm)*(x-xca));
  
				# ARMADO DE LA FUNCION

  abs1=linspace(0,xca,50);abs1(:,end)=[];

  abs2=linspace(xca,Lf,50);

  abs=[abs1,abs2];

  M=[M0A(abs1),MAF(abs2)];

  figure(3);clf;

  hold on;

  grid on;

  title ('DIAGRAMA DE CARACTERISTICAS');

  plot(abs,M,["--" markStyle(1) color(1) ";M(x)  Kgf.m " num2str(3) ";"]);

  ## CORTE FUNCION DE X

  Q0A=function_handle(-int(q,x));

  QAF=function_handle(-int(q,x)+Nl-Pm);
  
  

				# ARMADO DE LA FUNCION

  Q=[Q0A(abs1),QAF(abs2)];

  plot(abs,Q,["--" markStyle(2) color(2) ";Q(x) Kgf " num2str(3) ";"]);

  ## IMPRESION DE RESULTADOS

  printf ("L=%d Kgf\n",double(L));
  printf ("D=%d Kgf\n",double(D));
  printf ("T=%d Kgf\n",double(T));
  printf ("Lt=%d Kgf\n",double(Lt));
  printf ("Mac=%d Kgf.m\n\n Estructurales:\n",double(Macn));


  printf ("Q=%d Kgf\n",double(int(q,x,0,Lf)));
  printf ("Nl=%d Kgf\n",double(Nl));
  printf ("Pm=%d Kgf\n",double(Pm));
  printf ("T=%d Kgf\n\n",double(T*sin(alpha)));

  printf ("Mt=%d Kgf m\n",double(Mt));
  printf ("Mct=%d Kgf m\n",double(-Mct));
  printf ("Mca=%d Kgf m\n",double(Macn));

  printf ("Nl=%d Kgf\n",double(Nl));
  printf ("Pm=%d Kgf\n",double(Pm));
  printf ("Nt=%d Kgf\n",double(Nt));
  printf ("Wx=%d Kgf\n",double(int(q,x,0,Lf)));

endfunction
