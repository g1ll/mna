%DISCIPLINA DE MÉTODOS NUMÉRICOS APLICADOS
%PROFESSOR: WILIAM C. MARQUES - 28/04/2015
close all; clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PROGRAMA PARA RESOLVER A EQUAÇÃO DIFERENCIAL EM 2 DIMENSÕES
%DC/DT + UDC/DX + VDC/DY = mi(D²C/DX² + D²C/DY²)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOMINIO NUMÉRICO
alfa = 1;%LIMITE DO DOMINIO EM X
beta = 1;%LIMITE DO DOMINIO EM Y
deltax = 0.1;%VARICAO ESPACIAL EM X
deltay = 0.1;%VARIACAO ESPACIAL EM Y
deltat = 0.001;%PASSO DE TEMPO

Ka = deltat/deltax;%PARÂMETRO DE ESTABILIDADE
Kaa = deltat/(deltax^2);%PARÂMETRO DE ESTABILIDADE
mi = 1;%COEFICIENTE DE DIFUSÃO
timesim=50;%TEMPO DE SIMULAÇÃO

%Propriedade escalar
cp = 10;

%uadv(1:timesim) = linspace(2,4,timesim);%VELOCIDADE DE ADVECÇÃO EM X
%vadv(1:timesim) = linspace(4,2,timesim);%VELOCIDADE DE ADVECÇÃO EM Y

uadv(1:timesim) = linspace(4,1,timesim);%VELOCIDADE DE ADVECÇÃO EM X
vadv(1:timesim) = linspace(4,1,timesim);%VELOCIDADE DE ADVECÇÃO EM Y

y=[0:deltax:beta]; %DOMINIO DO TEMPO
x=[0:deltay:alfa]; %DOMINIO ESPACIAL
u=zeros(length(y),length(x)); %MATRIZ DA PROPRIEDADE FISICA
U=zeros(length(y),length(x)); %MATRIZ DA PROPRIEDADE FISICA
equ=(length(y)-2)*(length(x)-2);%MATRIZ DA PROPRIEDADE FISICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONDIÇÃO INICIAL E DE CONTORNO
U(:,:)= 0;
U(1,:) = 0;
%U(length(y),:) = 10;
U(length(y),5:7) = cp;
U(:,1)=0;
U(:,length(x)) = 0;
USAI(:,:,1)=U;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTRUTURA PARA ORDENAR AS EQUAÇÕES DE UM DOMINIO (N x M)
t=1;
for i = 2:length(x) - 1
    for j = 2: length(y) - 1
        if(equ <= 6)
        index(1,1:8,t) = [j j j j+1 j j j-1 j+1];
        index(2,1:8,t) = [i i i+1 i i-1 i+1 i i];
        else
        index(1,1:8,t) = [j j j j+1 j j j-1 j+1];
        index(2,1:8,t) = [i i i+1 i i-1 i+1 i i];
        index(1,9:equ,t)=0;
        index(2,9:equ,t)=0;
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
for tml = 1:timesim %LOOP DO TEMPO
    
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
          CC = (1 + uadv(tml)*Ka + vadv(tml)*Ka - 4*mi*Kaa)*U(index(1,1,n),index(2,1,n)) - Ka*( U(indJ(1),indI(1)) + ...
          U(indJ(4),indI(4)) ) + mi*Kaa*(U(indJ(1),indI(1)) + U(indJ(2),indI(2)) + U(indJ(3),indI(3)) + ...
          U(indJ(4),indI(4)) );
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
    for n = 2:length(y) -1 %length(y) - 1:-1:2
        u(n,m) = Ut(t);
        t=t+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ATUALIZAÇÃO NO TEMPO
U = u;
U(1,:) =0;
U(length(y),5:7) = cp;
U(:,1)=0;
U(:,length(x)) = 0;
USAI(:,:,tml+1)=U;
USAI
end% FIM DO LOOP DO TEMPO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANÁLISE DOS RESULTADOS
maxcp = 0;
for j=1:timesim+1
    mt = max(max(USAI(:,:,j)));
    if maxcp<mt
        maxcp=mt;
    end
end
maxcp
if maxcp <= cp
%
% for j=1:timesim
% pcolor(X,Y,squeeze(USAI(:,:,j)))
% shading interp
% hold on
% plot(X,Y,'xk')
% hold on
% quiver(X,Y,vu,vv,'k')
% xlabel('COMPRIMENTO (m)')
% ylabel('COMPRIMENTO (m)')
% caxis([0 10])
% colorbar
% M(j) = getframe;
% end
A=min(min(min(USAI)));
B=max(max(max(USAI)));
for j=1:timesim
    
y=[beta:-deltay:0];
[X,Y]=meshgrid(x,y);
[a,b]=size(X);
vu=uadv(j)*ones(a,b);vv=vadv(j)*ones(a,b);    
    
pcolor(X,Y,USAI(:,:,j))
set(gcf,'Visible', 'off'); 
shading interp
hold on
plot(X,Y,'xk')
hold on
quiver(X,Y,vu,vv,'k')
xlabel('COMPRIMENTO (m)')
ylabel('COMPRIMENTO (m)')
title(sprintf('Simulação no tempo: %d Vel X: %.2f Vel Y: %.2f',j,uadv(j),vadv(j)),'fontsize',12)
caxis([A B])
colorbar
print(strcat('fig_b_',num2str(100+j)),'-dpng','-r300');
if (j==1)||(j==ceil(timesim/2))||(j==timesim)
    print(strcat('plot_b_',num2str(j),'_cd_',num2str(mi),'_cp_',num2str(cp),'_vx_',num2str(uadv(j)),'_vy_',num2str(vadv(j))),'-djpeg','-r600')
end
close(gcf)
%M(j) = getframe;
end
!convert -delay 60 -loop 0 fig_b_*.png dispersao_2_b.gif
!rm fig_b_*.png
 
f1 = figure(1);
set(gcf,'Visible', 'off'); 
subplot(1,1,1)
plot(1:1:timesim,uadv,'r',1:1:timesim,vadv,'k--','LineWidth',2)
axis([0 timesim 0 5]);
legend('Velocidade em X','Velociade em Y')
xlabel('TEMPO (Segundos)','fontsize',12,'fontweight','b')
ylabel('VELOCIDADE','fontsize',12,'fontweight','b')
title('VELOCIDADES VARIÁVEIS','fontsize',12,'fontweight','b')
set(gca,'fontweight','b')       
print('plot_q_2_b','-dpng','-r600');
fprintf('\nQuestão 2 B: Simulação finalizada\n');
ls
else
fprintf('\nErro de estabilidade\n');
end
