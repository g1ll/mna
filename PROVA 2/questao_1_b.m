%DISCIPLINA DE MÉTODOS NUMÉRICOS APLICADOS
%PROFESSOR: WILIAM C. MARQUES - 28/04/2015
close all; clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PROGRAMA PARA RESOLVER A EQUAÇÃO DIFERENCIAL EM 2 DIMENSÕES
%DU/DT = ni( D²U/DX² + D²U/DY²) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOMINIO NUMÉRICO
alfa = 1;%LIMITE DO DOMINIO EM X
beta = 1;%LIMITE DO DOMINIO EM Y
deltax = 0.1;%VARICAO ESPACIAL EM X
deltay = 0.1;%VARIACAO ESPACIAL EM Y
four = 0.25;%PARÂMETRO DE ESTABILIDADE - NUMERO DE FOURIER
deltat = four*(deltax^2);%PASSO DE TEMPO
ni = 1; %COEFICIENTE DE DIFUSÃO
timesim=50;
y=[0:deltax:beta]; %DOMINIO DO TEMPO
x=[0:deltay:alfa]; %DOMINIO ESPACIAL
u=zeros(length(y),length(x)); %MATRIZ DA PROPRIEDADE FISICA
U=zeros(length(y),length(x)); %MATRIZ DA PROPRIEDADE FISICA
equ=(length(y)-2)*(length(x)-2);
ti(1:timesim) = (1:timesim)*270*0.1+270;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONDIÇÃO INICIAL E DE CONTORNO
U(:,:)=270; % 1226.85º C
U(1,:) = 270; %% 1226.85º C
U(length(y),:) = 270; %% 1226.85º C
U(:,1)=270; %% 0º C
U(:,length(x)) = 270; %% 0º C
USAI(:,:,1)=U;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTRUTURA PARA ORDENAR AS EQUAÇÕES DE UM DOMINIO (N x M)
t=1;
for i = 2:length(x) - 1
    for j = 2: length(y) - 1
        if(equ <= 6)
        index(1,1:7,t) = [j j j j j j-1 j+1];
        index(2,1:7,t) = [i i i i-1 i+1 i i];
        else
        index(1,1:7,t) = [j j j j j j-1 j+1];
        index(2,1:7,t) = [i i i i-1 i+1 i i];
        index(1,8:equ,t)=0;
        index(2,8:equ,t)=0;
        end
        t=t+1;
    end
end
t=1;
for i = 2:length(x) - 1
    for j = 2: length(y) - 1
        order(1,t) = j;
        order(2,t) = i;
        order(3,t) = t;
        t=t+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTRUTURA PARA MONTAR O SISTEMA DE UM DOMINIO (N x M) E RESOLVER AO LONGO
% DO TEMPO
[a,b,c]=size(index);
A=zeros(equ,equ);
%
for tml = 1:timesim %LOOP DO TEMPO
f(t) = -1;%    
CC=0;
      for n = 1:(equ)
          for m = 1:(equ)
             if(index(1,m,n) == order(1,n) & index(2,m,n) == order(2,n))
                 if (n == order(3,n))
                   A(n,n) = 1;
                 else
                   A(n,n) = 0;
                 end
             end
          end
      end
%
      for n=1:(equ)
          t=1;
          for m=1:(equ)
             if(index(1,m,n) ~= order(1,n) | index(2,m,n) ~= order(2,n) & n == order(3,n))
                 if(index(1,m,n) ~= 0 & index(2,m,n) ~= 0)
                   indI(t) = index(2,m,n); indJ(t) = index(1,m,n);
                   t=t+1;  
                 end
             end
            
          end
          CC = (1 - 4*ni*four)*U(index(1,1,n),index(2,1,n)) + (ni*four)*( U(indJ(1),indI(1)) + U(indJ(2),indI(2)) + ...
          U(indJ(3),indI(3)) + U(indJ(4),indI(4)) );
          B(n,1) = CC;CC=0;
      end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SOLUÇÃO DO SISTEMA PELO MÉTODO DIRETO DE ELIMINAÇÃO DE GAUSS
t1=2;
[a1 b1]=size(A);
%
MATAMP(1:a1,1:b1+1)=[A(:,:),B(1:b1,1)]; %MONTAGEM DA MATRIZ AMPLIADA
[aa,bb]=size(MATAMP);
%
%disp('MATRIZ AMPLIADA')
%MATAMP
%pause
%
N=1;
for n=1:bb
%       
  for m=t1:aa
%      
  if (MATAMP(m-1,n) ~= 0)    
  COEF(m-1) = MATAMP(m,n)/MATAMP(N,n);% CÁLCULO DOS COEFICIENTES DE ESCALONAMENTO
  else
  COEF(m-1) = 0;    
  end 
  end
%  
  for m=1:bb
  for mm=t1:aa
%      
   if (COEF(mm-1) ~= 0)    
   L(mm,m)= MATAMP(mm,m) - COEF(mm-1)*MATAMP(N,m);% ATUALIZAÇÃO DAS LINHAS DA MATRIZ AMPLIADA
   else
   L(mm,m) = MATAMP(mm,m);    
   end
%   
  end
  end
%  
  for mm=t1:aa
  MATAMP(mm,:)=L(mm,:); % ATUALIZAÇÃO DA MATRIZ AMPLIADA
  end
  t1=t1+1;N=N+1;   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SOLUÇÃO DO SISTEMA Ax = B
a_tr = MATAMP(:,1:bb-1); % MATRIZ DE COEFICIENTES (A)
b_tr = MATAMP(:,bb); % MATRIZ DE CONSTANTES (B)
[n,m]=size(a_tr);% TAMANHO DO SISTEMA  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SOLUÇÃO UTILIZANDO TRIANGULARIZAÇÃO SUPERIOR 
x_tr(n)=b_tr(n)/a_tr(n,n);
for k=n-1:-1:1
    sum=0;
    for j=k+1:n
      sum = sum + (a_tr(k,j)*x_tr(j)); % SOMATÓRIA DA RELAÇÃO GERAL 
    end
    x_tr(k) = (b_tr(k) - sum)/a_tr(k,k); % RELAÇÃO GERAL
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('A MATRIZ DE COEFICIENTES')
%
%a_tr
%
%disp('A MATRIZ DE CONSTANTES')
%
%b_tr
%
%disp('A MATRIZ DE RESULTADOS')
%x_tr'
%pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATUALIZAÇÃO DOS RESULTADOS
u(1,:) = 0;F(1,:)=0;
u(length(y),:) = 0;F(length(y),:) = 0;
u(:,1)=0;F(:,1)=0;
u(:,length(x)) = 0;F(:,length(x)) = 0;
Ut=x_tr';
t=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 2:length(x) - 1    
    for n = 2:length(y) -1 
        u(n,m) = Ut(t);
        t=t+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ATUALIZAÇÃO NO TEMPO
U = u;

U(1,:) = 270;
U(length(y),:) = 270;
U(:,length(x)) = ti(tml);
U(:,1)= ti(tml);
USAI(:,:,tml+1)=U;
U
disp('TEMPO SIMULADO (S), ITERAÇÃO'), [tml*deltat tml U(1,1) U(length(y),1)]
end% FIM DO LOOP DO TEMPO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANÁLISE DOS RESULTADOS
y=[beta:-deltay:0];
[X,Y]=meshgrid(x,y);
%
% for j=1:timesim
% pcolor(X,Y,squeeze(USAI(:,:,j)))
% shading interp
% hold on
% plot(X,Y,'xk')
% xlabel('COMPRIMENTO (m)')
% ylabel('COMPRIMENTO (m)')
% colorbar
% M(j) = getframe;
% end 

A=min(min(min(USAI)));
B=max(max(max(USAI)));
for j=1:timesim
set(gcf,'Visible', 'off'); 
pcolor(X,Y,USAI(:,:,j))
shading interp
hold on
plot(X,Y,'xk')
hold on
xlabel('COMPRIMENTO (m)')
ylabel('COMPRIMENTO (m)')
title(sprintf('TEMP. LATERAIS: X_0 %.2f | X_1 %.2f | TEMP. CENTRO: %.2f',USAI(floor(length(y)/2),1,j),USAI(floor(length(y)/2),length(x),j),USAI(floor(length(x)/2),floor(length(y)/2),j)),'fontsize',10)
caxis([A B])
colorbar
print(strcat('fig_q1_b_',num2str(100+j)),'-dpng','-r300')
close(gcf)
%M(j) = getframe;
end
!convert -delay 50 -loop 0 fig_q1_b_*.png questao_1_b.gif
!rm fig_q1_b_*.png