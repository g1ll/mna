function [Qa  psi psirk2] = circuito_RLC(V,R,L,C,h,Tf)
    clc
    if R^2 - 4*L*1/C < 0
        x = roots([L R 1/C]) %retorna as raizes do polinômio
        T=0.00001:h:Tf; %DOMÍNIO NUMÉRICO (TEMPO) Tf = limite domínio
        % SOLUÇÃO ANALÍTICA PARA RAÍZES COMPLEXAS
        r = real(x(1,1)) %Parte real
        i = imag(x(1,1)) %Parte imaginária
        Qa = (V/(1/C)).*(1-exp(r.*T).*(cos(i.*T)-(r/i).*sin(i.*T)));
        
        % SOLUÇÃO NUMÉRICA DO PROBLEMA - MÉTODO DE EULER
        Qe(1)=0; 
        dQe(1)=0;
        psi(1)=Qe(1);
        phi(1)=dQe(1);
      
        for n=2:length(T)    
            %SISTEMA TRANSFORMADO
            psi(n)=psi(n-1) + h*phi(n-1); %CALCULO DA CARGA
            phi(n)=phi(n-1) + h*(V - R*phi(n-1) - psi(n-1)/C)/L; %CORRENTE
        end
       
        % SOLUÇÃO NUMÉRICA DO PROBLEMA - MÉTODO DE RUGE-KUTTA-2
        psirk2(1)=Qe(1);
        phirk2(1)=dQe(1);
        
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
    end
end
