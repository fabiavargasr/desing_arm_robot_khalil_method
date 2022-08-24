close all;
clear all;
clc;

syms t1 t2 t3  R4 d3
% hace una grafica del engowrist cuando varian los angulos 0  a 90 

%% Tabla de parametros


% gm al dj th rj 

matrix_p = [0   0     0  t1  0;
            0   pi/2  0  t2  0;
            0   -pi/2  d3 t3  0
            0   0     0  0   R4];

%% Se extraen los datos de la tabla.

for j = 1:4
    
al = matrix_p(j,2);
di = matrix_p(j,3);
th = matrix_p(j,4); 
rj = matrix_p(j,5);

if(j == 1)
M0T1 = [ cos(th)           -sin(th)         0          di;
         cos(al)*sin(th) cos(al)*cos(th) -sin(al) -rj*sin(al);
         sin(al)*sin(th) sin(al)*cos(th) cos(al)   rj*cos(al);
         0                      0           0          1]; 
elseif(j == 2)
M1T2 = [ cos(th)           -sin(th)         0          di;
         cos(al)*sin(th) cos(al)*cos(th) -sin(al) -rj*sin(al);
         sin(al)*sin(th) sin(al)*cos(th) cos(al)   rj*cos(al);
         0                      0           0          1]; 
elseif(j == 3)
M2T3 = [ cos(th)           -sin(th)         0          di;
         cos(al)*sin(th) cos(al)*cos(th) -sin(al) -rj*sin(al);
         sin(al)*sin(th) sin(al)*cos(th) cos(al)   rj*cos(al);
         0                      0           0          1];  

elseif(j == 4)
M3T4 = [ cos(th)           -sin(th)         0          di;
         cos(al)*sin(th) cos(al)*cos(th) -sin(al) -rj*sin(al);
         sin(al)*sin(th) sin(al)*cos(th) cos(al)   rj*cos(al);
         0                      0           0          1];  
end
end

%%

M0T1
M1T2
M2T3
M3T4

%% MGD

% Asiganación de valores
% [th1, th2, th4 ]            = deal(0, -90, 0);
% [t1, t2, t4, R1, R2, r3, R4] = deal(deg2rad(th1), deg2rad(th2), deg2rad(th4), 0.35, 0.25, 0.35, 0.1);

% Posiciones cartesianas para artilación 1
PX1  = M0T1(1,4)
PY1  = M0T1(2,4)
PZ1  = M0T1(3,4)

% Matrices de transformación respecto a la base
M0T2 = M0T1*M1T2;
M0T3 = M0T2*M2T3;
M0T4 = M0T3*M3T4;

M1T3 = M1T2*M2T3;


% Posiciones cartesianas para artilación 2
PX2  = M0T2(1,4)
PY2  = M0T2(2,4)
PZ2  = M0T2(3,4)
% Posiciones cartesianas para artilación 3
PX3  = M0T3(1,4)
PY3  = M0T3(2,4)
PZ3  = M0T3(3,4)
% Posiciones cartesianas para artilación 4
PX4  = M0T4(1,4)
PY4  = M0T4(2,4)
PZ4  = M0T4(3,4)

%% CARGAR DATOS 
t1=0;  t2=0; t3=0; d3= 0.009; R4=0.004


 for t1=0:pi/30:pi/2 
 t2=t1;
 t3=t1;
 M0T1_e=subs(M0T1);
 M1T2_e=subs(M1T2);
 M2T3_e=subs(M2T3);
 M3T4_e=subs(M3T4);

 %  mgd eval

 % MGD

% Asiganación de valores
% [th1, th2, th4 ]            = deal(0, -90, 0);
% [t1, t2, t4, R1, R2, r3, R4] = deal(deg2rad(th1), deg2rad(th2), deg2rad(th4), 0.35, 0.25, 0.35, 0.1);

% Posiciones cartesianas para artilación 1
PX1  = M0T1_e(1,4);
PY1  = M0T1_e(2,4);
PZ1  = M0T1_e(3,4);

%PX1=0.01


% Matrices de transformación respecto a la base
M0T2_e = M0T1_e*M1T2_e;
M0T3_e = M0T2_e*M2T3_e;
M0T4_e = M0T3_e*M3T4_e;

M1T3_e = M1T2_e*M2T3_e;


% Posiciones cartesianas para artilación 2
PX2  = M0T2_e(1,4);
PY2  = M0T2_e(2,4);
PZ2  = M0T2_e(3,4);
% Posiciones cartesianas para artilación 3
PX3  = M0T3_e(1,4);
PY3  = M0T3_e(2,4);
PZ3  = M0T3_e(3,4);
% Posiciones cartesianas para artilación 4
PX4  = M0T4_e(1,4);
PY4  = M0T4_e(2,4);
PZ4  = M0T4_e(3,4);



i=1;

% graficar posiciones endo
  % Gráficos
 
   PX0(i) = 0; PY0(i) = 0; PZ0(i) = 0;
   Ax = [PX0(i) PX1(i) PX2(i) PX3(i) PX4(i)];
   By = [PY0(i) PY1(i) PY2(i) PY3(i) PY4(i)];
   Cz = [PZ0(i) PZ1(i) PZ2(i) PZ3(i) PZ4(i)];
   plot3 (Ax, By, Cz, '-+', 'DisplayName', strcat('Punto ', i+48), 'LineWidth', 3);
   hold on;


   pos4=[PX4 PY4 PZ4]
  end

%% Calculo del MGI

syms sx sy sz nx ny nz nz ax ay az px py pz 
 
for j=1:3
al=matrix_p(j,2);
di=matrix_p(j,3);
th=matrix_p(j,4); 
rj=matrix_p(j,5);
if(j == 1)

M1T0 = [cos(th) cos(al)*sin(th) sin(al)*sin(th) -di*cos(th);
       -sin(th) cos(al)*cos(th) sin(al)*cos(th)  di*sin(th);
          0        -sin(al)          cos(al)        -rj;
          0           0                   0          1];
      
elseif(j == 2)
    
M2T1 = [cos(th) cos(al)*sin(th) sin(al)*sin(th) -di*cos(th);
       -sin(th) cos(al)*cos(th) sin(al)*cos(th)  di*sin(th);
          0        -sin(al)          cos(al)        -rj;
          0           0                   0          1];

elseif(j == 3)
    
M3T2 = [cos(th) cos(al)*sin(th) sin(al)*sin(th) -di*cos(th);
       -sin(th) cos(al)*cos(th) sin(al)*cos(th)  di*sin(th);
          0        -sin(al)          cos(al)        -rj;
          0           0                   0          1];
end
end

M1T4 = M1T3*M3T4;
% posicion y orientacion deseada

Uo = [sx	nx	ax	px;
      sy	ny	ay	py;
      sz	nz	az	pz;
       0	0	0	1];
%%  con esto se compara las matrices y se hace interaciones de paul, hay que buscar la ecuacion mas sencilla para encontrar 
% los angulos la solicion hay que encontrarla en la ecuacion mas facil de
% resolver

% si las ecuaciones estan muy largas se busca en el libro de khalil
% apendice 1 los tipos de ecuaciones y las soliciones posibles
U1 = M1T0*Uo;
a = U1 == M1T4  % esta linea saca la comparacion de las matrices con las ecuaciones para despejar posibles soluciones

% de la primera iteracion la mejor es py*cos(t1) - px*sin(t1) == 0

% PARA T1
%syms x y
%solve('py*cos(t1) == px*sin(t1)',t1)
% para t2 la ecuacion mejor es  px*cos(t1) + py*sin(t1) == d3*cos(t2) - R4*sin(t2)


t1 = atan2(py,px)


%PARA T2
%la solicion se la busca asi, dado que es una ecuacion del tipo 
% Y COS(TH)  + X SEN(TH)=Z
X = -R4;  Y = d3;  Z = px*cos(t1) + py*sin(t1);
e = -1;
SQ1 = (X*Z + e*Y*sqrt(X^2 + Y^2 - Z^2))/(X^2 + Y^2);
CQ1 = (Y*Z - e*X*sqrt(X^2 + Y^2 - Z^2))/(X^2 + Y^2);

t2 = atan2(SQ1,CQ1);

% PARA T3  LA MEJOR ES nx*cos(t1) + ny*sin(t1) == -cos(t2)*sin(t3)


sin(t3)=(nx*cos(t1) + ny*sin(t1))/ -cos(t2);
cos(t3)=(ny*cos(t1) - nx*sin(t1)); 

t3= atan2(sin(t3),cos(t3));


%%


%%
%----Movimiento circular-----------
%-------------------------------
sx = 1;	nx = 0;  ax = 0;
sy = 0;	ny = 1;  ay = 0;
sz = 0;	nz = 0;  az = 1;
%

% Ubicación de los ejes
l = 0.015;
Ax = [0 l];  By = [0 0];  Cz = [0 0];
plot3(Ax, By, Cz, 'r', 'Linewidth', 2);  hold on;
Ax = [0 0];  By = [0 l];  Cz = [0 0];
plot3(Ax, By, Cz, 'g', 'Linewidth', 2);  hold on;
Ax = [0 0];  By = [0 0];  Cz = [0 l];
plot3(Ax, By, Cz, 'b', 'Linewidth', 2);  hold on;
legend('Eje x', 'Eje y', 'Eje z');  grid on;
%axis([-0.05 0.55 -0.3 0.2 -0.1 0.4]);  pause(1);
%% Animación

%----Movimiento circular-----------
%-------------------------------

% Tiempo de muestreo y duración final de la trayectoria:

Tfinal=1.0;
Tem=0.01;

% Cálculo del número de muestras:

nbech=(Tfinal/Tem)+1;
if ((round(nbech)-nbech) == 0)
instant=[0:Tem:Tfinal]';
else
nbech=nbech+1;
instant=[0:Tem:Tfinal+Tem]';
end
t=0;

radio=0.008

for h=1:1:nbech
t=t+Tem;
x1(h)=radio*sin(2*pi*1/Tfinal*t);
y1(h)=radio*cos(2*pi*1/Tfinal*t);
end
x1=x1';
y1=y1';
%--------------------
% se suma la posicion inicial
cons1= 0 + x1;
cons2= 0 + y1;
cons3= 0.006*ones(101,1);


%plot3 (cons1, cons2, cons3, '-k', 'DisplayName', 'Trayectoria', 'LineWidth', 1);
%hold on;
%%

for i = 1 : 50
   px = cons1(i);
   py = cons2(i);
   pz = cons3(i);
   
   % MGI
   
    
        t1_c = atan2(py,px);
        t1=t1_c(1);
        
        X = -R4;  
        Y = d3;  
        Z = px*cos(t1) + py*sin(t1);
        e = 1;
        SQ1 = (X*Z + e*Y*sqrt(X^2 + Y^2 - Z^2))/(X^2 + Y^2);
        CQ1 = (Y*Z - e*X*sqrt(X^2 + Y^2 - Z^2))/(X^2 + Y^2);
        
        t2 = atan2(SQ1,CQ1);
        
        sin_t3=(nx*cos(t1) + ny*sin(t1))/ -cos(t2);
        cos_t3=(ny*cos(t1) - nx*sin(t1)); 
        
        t3= atan2(sin_t3,cos_t3);


   % MGD
         px0(i)=  0    ; 
		 py0(i)=    0  ;   
		 pz0(i)=    0;
         px1(i)=  0   ;  
		 py1(i)=    0  ;   
		 pz1(i)=    0   ;  
		 px2(i)=    0    ; 
		 py2(i)=    0    ;
		 pz2(i)=    0     ; 
		 px3(i)=    d3*cos(t1)*cos(t2) ;   
		 py3(i)=    d3*cos(t2)*sin(t1) ;  
		 pz3(i)=    d3*sin(t2) ;    
		 px4(i)=    d3*cos(t1)*cos(t2) - R4*cos(t1)*sin(t2);     
		 py4(i)=    d3*cos(t2)*sin(t1) - R4*sin(t1)*sin(t2);   
		 pz4(i)=    R4*cos(t2) + d3*sin(t2);
   
   % Gráficos
   Ax = [px0(i) px1(i) px2(i) px3(i) px4(i)];
   By = [py0(i) py1(i) py2(i) py3(i) py4(i)];
   Cz = [pz0(i) pz1(i) pz2(i) pz3(i) pz4(i)];
   plot3 (Ax, By, Cz, '-+', 'DisplayName', strcat('Punto ', i+48), 'LineWidth', 3);
   hold on;
   plot3 (cons1, cons2, cons3, '-k', 'DisplayName', 'Trayectoria', 'LineWidth', 1);
   %axis([-0.05 0.55 -0.3 0.2 -0.1 0.4]);  
   
 % pause(0.0001); 
 % hold off;
end
hold on
% Trayectoria Realizada
%Ax = [px(1) px(5)]; By = [py(1) py(5)]; Cz = [pz(1) pz(5)];
% plot3 (cons1, cons2, cons3, '-k', 'DisplayName', 'Trayectoria', 'LineWidth', 5);
%axis([-0.05 0.55 -0.3 0.2 -0.1 0.4]);  hold off;

%% Modelo cinematico directo

% Columna 1 de matriz jacobiana, K=1
j04_1 = [-M1T4(2,4)*M0T1(1:3,1) + M1T4(1,4)*M0T1(1:3,2); M0T1(1:3,3)]; 
% Columna 2 de matriz jacobiana, K=2
j04_2 = [-M2T4(2,4)*M0T2(1:3,1) + M2T4(1,4)*M0T2(1:3,2); M0T2(1:3,3)];
% Columna 3 de matriz jacobiana, K=3
j04_3 = [M0T3(1:3,3); 0*M0T3(1:3,3)];
% Columna 3 de matriz jacobiana, K=3
j04_4 = [0;0;0; M0T4(1:3,3)];
% Matriz Jacobiana
J = [j04_1 j04_2 j04_3 j04_4];

