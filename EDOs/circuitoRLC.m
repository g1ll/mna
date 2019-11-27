function [Qa  psi psirk2] = circuitoRLC(V,R,L,C,h,Tf)
    clc
    %APRESENTAÇÃO
    fprintf('\n        MODELO COMPUTACIONAL PARA SOLUÇÃO DO SISTEMA DINÂMICO');
    fprintf('\n\tCIRCUITO ELÉTRICO RLC EM SÉRIE (CIRCUITO RESSONANTE)\n');
    fprintf('\n\t\t\tAMORTECIMENTO FRACO\n');
    fprintf('\n\t\t\tR² - 4*L*1/C < 0\n');
    fprintf('\n\t\t           R (Resistor)');
    fprintf('\n\t\t   +-----/\\/\\/\\------+');
    fprintf('\n\t\t   |                 |');
    fprintf('\n\t\t   |                 |');
    fprintf('\n\t\t   +                === ');
    fprintf('\n\t\t  (~) V (Fonte)     === L (Indutor)');
    fprintf('\n\t\t   _                ===');
    fprintf('\n\t\t   |                 |');
    fprintf('\n\t\t   |                 |');
    fprintf('\n\t\t   +-------||--------+');
    fprintf('\n\t\t           C (Capacitor)\n');
    fprintf('\n\t\t\t MODELO MATEMÁTICO\n');
    fprintf('\n\t\t L*d²Q/dt² + R*dQ/dt + Q/C = V\n');
    fprintf('\n\tQ(T) = CARGA NO CAPACITOR EM FUNÇÃO DO TEMPO\n');
                
    %PARÂMETROS
    fprintf('\nPARÂMETROS:\n');
    fprintf('\n\tTENSÃO (Volts): \t%.2E',V);
    fprintf('\n\tRESISTÊNCIA (Ohms): \t%.2E',R);
    fprintf('\n\tINDUTÂNCIA (Henrys): \t%.2E',L);
    fprintf('\n\tCAPACITÂNCIA (Farads): \t%.2E',C);
      
    fprintf('\n\nMALHA DO TEMPO:\n');
    fprintf('\n\tPASSO DE TEMPO: \t%.2E',h);
    fprintf('\n\tLIMITE DE TEMPO (s): \t%.2E',Tf);
    %SOLUÇÃO DO SISTEMA DINÂMICO 
    %CIRCUITO ELÉTRICO RLC EM SÉRIE (CIRCUITO RESSONANTE)
    % AMORTECIMENTO FRACO
    % R² - 4*L*1/C

    %            R (Resistor)
    %   +-----/\/\/\------+
    %   |                 |
    %   |                 |
    %   +                === 
    %  (~) V (Fonte)     === L (Indutor)
    %   _                ===
    %   |                 |
    %   |                 |
    %   +-------||--------+
    %           C (Capacitor)

    % MODELO MATEMÁTICO

    % LQ"(t) + RQ'(t) + (1/C)Q(t) = V(t)

    % Q(T) = CARGA NO CAPACITOR EM FUNÇÃO DO TEMPO

    %PARÂMETROS

    % L = Indutância (Henrys)
    % R = Resistência (Ohms)
    % C = Capacitancia (Faradays)
    % V = Tensão (Volts)

    % CONDIÇÕES INICIAIS (t = tempo inicial | i = corrente inicial)
    % Q(t) = q
    % Q'(t) = I(t) = i

    % PROBLEMA:
    % PARÂMETROS
    %   L = 1, R = 10³ C = 10⁻⁶, V = 12

    % DETERMINAR CARGA Q(T) Coulombs em T = 1 Segundo | Q(T) =  ?

    % CONDIÇÕES INICIAIS 
    % Q(0) = 0
    % Q'(0) = I(0) = 0

    % PARÂMETROS DO PROBLEMA

    % V = 12;     % Tensão (Volts)
    % R = 10^3;    % Resistência (Ohms)
    % L = 1;      % Indutância (Henrys)
    % C = 10^-6;  % Capacitancia (Faradays)

    %SOLUÇÃO ANALÍTICA para L = 1, R = 10³ C = 10⁻⁶
    % RAIZES DE x² +10³x + 1/10⁻⁶ = 0
    % x1,2  = -500 +- 500*sqrt(3)i
    % Q(T) =  12/10^6*( 1- (e^(-500*T))*(  cos(500*sqrt(3)*T) + 1/(sqrt(3))*sen(500*sqrt(3)*T) ) )

    if R^2 - 4*L*1/C < 0

        % L*Q"(t) + R*Q'(t) + (1/C)*Q(t) = V.

        % SOLUÇÃO CARACTERÍSTICA
        %   Lx² +Rx + 1/C = 0

        % Raízes:

        x = roots([L R 1/C]); %retorna as raizes do polinômio

        % MALHA DO TEMPO
        % h = PASSO DE CÁLCULO
        T=0.00001:h:Tf; %DOMÍNIO NUMÉRICO (TEMPO) Tf = limite domínio
        
        % CALCULANDO SOLUÇÕES

        % SOLUÇÃO ANALÍTICA PARA RAÍZES COMPLEXAS
        r = real(x(1,1)); %Parte real
        i = imag(x(2,1)); %Parte imaginária
        Qa = (V/(1/C)).*(1-exp(r.*T).*(cos(i.*T)-(r/i).*sin(i.*T)));
        
        % SOLUÇÃO NUMÉRICA DO PROBLEMA - MÉTODO DE EULER
        
        %   L*d²Q/dt² + R*dQ/dt + Q/C = V
        %
        %
        %   CONDIÇÕES INICIAIS DO SISTEMA ORIGINAL
        %
        Qe(1)=0; 
        dQe(1)=0;
        %
        %TROCA DE VARIÁVEIS PARA MONTAR O SISTEM TRANSFORMADO
        %
        psi(1)=Qe(1);
        phi(1)=dQe(1);
        %
        %
        % CÁLCULO DO SISTEMA TRANSFORMADO
        %
        %       phi = dQ 
        %       phid = ddQ
        %
        %       psi = Q
        %       dpsi = dQ = phi
        %
        % L*phid + R*phi + psi/C = V --> phid = (V - ( R*phi - psi/C))/L
        %
        %RESOLVENDO UTILIZANDO MÉTODO DE EULER
        %
        %psi(n) = psi(n-1) + h*phi(n-1)
        %
        %phi(n) = phi(n-1) +h*(V - ( R*phi - psi/C))/L        
        %
        for n=2:length(T)    
            %SISTEMA TRANSFORMADO
            psi(n)=psi(n-1) + h*phi(n-1); %CALCULO DA CARGA
            phi(n)=phi(n-1) + h*(V - R*phi(n-1) - psi(n-1)/C)/L; %CALCULO DA CORRENTE
        end
       
        % SOLUÇÃO NUMÉRICA DO PROBLEMA - MÉTODO DE RUGE-KUTTA-2
        
        %   L*d²Q/dt² + R*dQ/dt + Q/C = V
        
        psirk2(1)=Qe(1);
        phirk2(1)=dQe(1);
        %
        for n=2:length(T)
            %SISTEMA TRANSFORMADO
            %
            k1 = phirk2(n-1);
            phirk2t=phirk2(n-1) + h*(V - R*phirk2(n-1) - psirk2(n-1)/C)/L;
            k2 = phirk2t;
                   
            psirk2(n) = psirk2(n-1) + (h/2)*(k1 + k2); %CALCULO DA CARGA

            k1 = (V - R*phirk2(n-1) - psirk2(n-1)/C)/L;
            Cj = phirk2(n-1) + h*k1;
            psirk2t=psirk2(n-1) + h*phirk2(n-1);
            k2 = (V-R*Cj -(1/C)*psirk2t)*1/L;
            phirk2(n)=phirk2(n-1) + (h/2)*(k1 + k2); %CALCULO DA CORRENTE
        end
        
        %--------------------IMPRESSÃO DE RESULTADOS-----------------------
       

        % RESULTADOS DA SOLUÇÃO ANALÍTICA
   
        fprintf('\n\nSOLUÇÃO ANALÍTICA\n');
        fprintf('\nQ(%.3d) = %.2E C coulombs \nQ(%d) = %.2E C coulombs',T(1),Qa(1)*10^6,T(length(T)),Qa(length(T)));
        fprintf('\nPico de Carga Qmax =  %.2E C coulombs',max(Qa)); 
        fprintf('\nCarga Mínima Qmin = %.2E C coulombs\n',min(Qa)); 
        
        % RESULTADOS DA SOLUÇÃO POR EULER
       
        fprintf('\n\nSOLUÇÃO NUMÉRICA EULER\n');
        fprintf('\nQ(%.3d) = %.2E C coulombs \nQ(%d) = %.2E C coulombs',T(1),psi(1),T(length(T)),psi(length(T)));
        fprintf('\nPico de Carga Qmax =  %.2E C coulombs',max(psi)); 
        fprintf('\nCarga Mínima Qmin = %.2E C coulombs\n',min(psi)); 
        fprintf('\nErro em Pico de Carga Qmax =  %.2E C coulombs',max(Qa)-max(psi)); 
        fprintf('\nErro em Carga Mínima Qmin = %.2E C coulombs\n',min(Qa)-min(psi));
        fprintf('\nErro  =  %.2E\n',mean(abs(Qa)-abs(psi))); 
        
        % RESULTADOS DA SOLUÇÃO POR RUGE-KUTTA
       
        fprintf('\n\nSOLUÇÃO NUMÉRICA RUGE-KUTTA\n');
        fprintf('\nQ(%.3d) = %.2E C coulombs \nQ(%d) = %.2E C coulombs',T(1),psirk2(1),T(length(T)),psirk2(length(T)));
        fprintf('\nPico de Carga Qmax =  %.2E C coulombs',max(psirk2)); 
        fprintf('\nCarga Mínima Qmin = %.2E C coulombs\n',min(psirk2)); 
        fprintf('\nErro em Pico de Carga Qmax =  %.2E C coulombs',max(Qa)-max(psirk2)); 
        fprintf('\nErro em Carga Mínima Qmin = %.2E C coulombs\n',min(Qa)-min(psirk2)); 
        fprintf('\nErro  =  %.2E\n',mean(abs(Qa)-abs(psi))); 
        
        % PLOTAGEM DE RESULTADOS E COMPARAÇÃO DE MODELOS
        
       
        subplot(2,2,1)
        grid on
        plot(T,Qa,'r','LineWidth',2)
        legend('CARGA ANALÍTICA')
        xlabel('TEMPO (Segundos)','fontsize',10,'fontweight','b')
        ylabel('CARGA (Coulombs)','fontsize',10,'fontweight','b')
        title('CARGA CIRCUITO RLC EM SÉRIE','fontsize',12,'fontweight','b')
        set(gca,'fontweight','b')
        
        subplot(2,2,2)
        grid on
        plot(T,psi,'k','LineWidth',2)
        legend('CARGA EULER')
        xlabel('TEMPO (Segundos)','fontsize',10,'fontweight','b')
        ylabel('CARGA (Coulombs)','fontsize',10,'fontweight','b')
        title('CARGA CIRCUITO RLC EM SÉRIE','fontsize',12,'fontweight','b')
        set(gca,'fontweight','b')
        
        subplot(2,2,3)
        grid on
        plot(T,psirk2,'b','LineWidth',2)
        legend('CARGA RUGE-KUTTA')
        xlabel('TEMPO (Segundos)','fontsize',10,'fontweight','b')
        ylabel('CARGA (Coulombs)','fontsize',10,'fontweight','b')
        title('CARGA CIRCUITO RLC EM SÉRIE','fontsize',12,'fontweight','b')
        set(gca,'fontweight','b')
        
        subplot(2,2,4)
        grid on
        plot(T,Qa,'r',T,psi,'k--',T,psirk2,'b--','LineWidth',2)
        legend('CARGA ANALÍTICA','CARGA EULER','CARGA RK2')
        xlabel('TEMPO (Segundos)','fontsize',10,'fontweight','b')
        ylabel('CARGA (Coulombs)','fontsize',10,'fontweight','b')
        title('CARGA CIRCUITO RLC EM SÉRIE','fontsize',12,'fontweight','b')
        set(gca,'fontweight','b')
        
        f1 = figure(2);
        subplot(1,1,1)
        plot(T,Qa,'r',T,psi,'k--',T,psirk2,'b--','LineWidth',2)
        legend('CARGA ANALÍTICA','CARGA EULER','CARGA RK2')
        xlabel('TEMPO (Segundos)','fontsize',12,'fontweight','b')
        ylabel('CARGA (Coulombs)','fontsize',12,'fontweight','b')
        title('CARGA CIRCUITO RLC EM SÉRIE','fontsize',12,'fontweight','b')
        set(gca,'fontweight','b')       
        
        f2 = figure(3);
        
        subplot(1,2,1)
        plot(T,abs(psi - Qa),'k')
        grid on
        legend('DIFERENÇA EULER - ANALÍTICA')
        xlabel('TEMPO (segundos)','fontsize',14,'fontweight','b')
        ylabel('DIFERENÇA (C)','fontsize',14,'fontweight','b')
        set(gca,'fontweight','b')
        
        subplot(1,2,2)
        plot(T,abs(psirk2 - Qa),'k')
        grid on
        legend('DIFERENÇA RK2 - ANALÍTICA')
        xlabel('TEMPO (segundos)','fontsize',14,'fontweight','b')
        ylabel('DIFERENÇA (C)','fontsize',14,'fontweight','b')
        set(gca,'fontweight','b')
        
        saveas(f1, 'plot1.png');
        saveas(f2, 'plot2.png');
    else
         fprintf('\n\n PARÂMETROS NÃO SUPORTADOS:\n'); 
         fprintf('\n Este Modelo Computacional suporta apenas circuitos RLC do tipo AMORTECIMENTO FRACO !\n\tOnde:'); 
         fprintf('\n\t\t R² - 4* 1/C < 0 e L > 0\n\t Porém: \n\t\t  %f ² - 4 x %f = %f >= 0 e %f <= 0\n ',R,C,(R^2 - 4*1/C),L); 
         fprintf('\n Informe valores de Resistência (R) e Capacitância (C) que respeitem a condição R² - 4* 1/C < 0 assim como L > 0\n'); 
         fprintf('\n\t Exemplo:\n\t\t R = 10³ , C = 10⁶\n'); 
    end
end
