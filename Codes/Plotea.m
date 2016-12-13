%Grafica los resultados del mapa log�stico

%% Reseta Matlab

clear
close all
clc

%% Par�metros Generales

Mapa = 'LogisticoB2';

%% Carga los datos

load(Mapa);

%% Define variables para dar formato a las figuras
PosicPlot=[50 50 800 600];
FuenteSize=26;
SubplotFuenteSize=20;
MarcaSize=10;
PuntoSize=30;
LineaAncho1=4;
LineaAncho2=2;

%% Plotea 'Hval', 'Qval', 'Cval', 'Hbp', 'Qbp', 'Cbp', 'Hbpw', 'Qbpw', 'Cbpw' y 'MP' vs P

Cuantis = {'Hval', 'Qval', 'Cval', 'Hbp', 'Qbp', 'Cbp', 'Hbpw', 'Qbpw', 'Cbpw', 'MP'};
NombresCuantis = {'H_{val}', 'Q_{val}', 'C_{val}', 'H_{BP}', 'Q_{BP}', 'C_{BP}', 'H_{BPW}', 'Q_{BPW}', 'C_{BPW}', 'MP'};
MaxCuanti = [1, 1, .5, 1, 1, .5, 1, 1, .5, factorial(D)];

for i_cuanti = 1:length(Cuantis);
    
    Cuanti = char(Cuantis(i_cuanti)); %Carga el nombre del cuantificador y lo pasa a string
    eval(['Datos = ' Cuanti ';']); %Carga los datos a plotear
    meanDatos = mean(Datos);
    
    NombreCuanti = char(NombresCuantis(i_cuanti)); %Carga el nombre del cuantificador y lo pasa a string
    
    figure('position',PosicPlot,'visible','off'); %Abro una figura con tama�o y posici�n
    hold on; grid on; box on
    %title('rand-ellip; roja=pasa-altos; azul=pasa-bajos; negro=rand-matlab')
    xlabel('P'); ylabel(NombreCuanti)
    ejes=newplot; %newplot returns handle of current axes y abre una figura
    set(ejes,'FontName','Arial','FontWeight','bold','FontSize',FuenteSize); %cambiamos las letras de los ejes de la figura
    set(gcf,'DefaultLineLineWidth',LineaAncho2); %fijamos el tamano de linea por default grosor 2
    axis([0 Precisiones(NPrecisiones) 0 MaxCuanti(i_cuanti)]);
    
    plot(Precisiones,meanDatos,'.:k','LineWidth', LineaAncho1, 'MarkerSize', PuntoSize)
    plot(Precisiones,Datos,'.r', 'MarkerSize', MarcaSize)    

    set(gcf, 'PaperType', 'e', 'PaperOrientation', 'Landscape', 'PaperUnits', 'Normalized', 'PaperPosition', [0 0 1 1]); %Papertype 'e' es el m�s cuadradito que encontr�, normalized normaliza las paperposition
    saveas(gcf,[Cuanti '_' Mapa],'pdf')
    close
    
end

%% Plotea per�odo

    Datos = Period;
    meanDatos = mean(Datos);
    
    NombreCuanti = 'Period'; %Carga el nombre del cuantificador y lo pasa a string
    
    figure('position',PosicPlot,'visible','off'); %Abro una figura con tama�o y posici�n
    hold on; grid on; box on
    %title('rand-ellip; roja=pasa-altos; azul=pasa-bajos; negro=rand-matlab')
    xlabel('P'); ylabel(NombreCuanti)
    ejes=newplot; %newplot returns handle of current axes y abre una figura
    set(ejes,'FontName','Arial','FontWeight','bold','FontSize',FuenteSize, 'YScale', 'log'); %cambiamos las letras de los ejes de la figura
    set(gcf,'DefaultLineLineWidth',LineaAncho2); %fijamos el tamano de linea por default grosor 2
    axis([0 Precisiones(NPrecisiones) .5 NMap]);
    
    semilogy(Precisiones,meanDatos, '.:k', 'LineWidth', LineaAncho1, 'MarkerSize', PuntoSize)
    semilogy(Precisiones,Datos,'.r', 'MarkerSize', MarcaSize)  
    
    set(gcf, 'PaperType', 'e', 'PaperOrientation', 'Landscape', 'PaperUnits', 'Normalized', 'PaperPosition', [0 0 1 1]); %Papertype 'e' es el m�s cuadradito que encontr�, normalized normaliza las paperposition
    saveas(gcf,['Period_' Mapa],'pdf')
    close
    
%% Plotea plano Hval-Hbp

  	Datos1 = Hval; %Carga el nombre del cuantificador y lo pasa a string
    Datos2 = Hbp;
    
    meanDatos1 = mean(Datos1);
    meanDatos2 = mean(Datos2);
    
    NombreCuanti1 = 'H_{Val}'; %Carga el nombre del cuantificador y lo pasa a string
    NombreCuanti2 = 'H_{BP}';
    
    figure('position',PosicPlot,'visible','off'); %Abro una figura con tama�o y posici�n
    hold on; grid on; box on
    %title('rand-ellip; roja=pasa-altos; azul=pasa-bajos; negro=rand-matlab')
    xlabel(NombreCuanti1); ylabel(NombreCuanti2)
    ejes=newplot; %newplot returns handle of current axes y abre una figura
    set(ejes,'FontName','Arial','FontWeight','bold','FontSize',FuenteSize); %cambiamos las letras de los ejes de la figura
    set(gcf,'DefaultLineLineWidth',LineaAncho2); %fijamos el tamano de linea por default grosor 2
    axis([0 1 0 1]);

    plot(meanDatos1,meanDatos2, '.:k', 'LineWidth', LineaAncho1, 'MarkerSize', PuntoSize)
    plot(Datos1,Datos2,'.r', 'MarkerSize', MarcaSize)
    
    set(gcf, 'PaperType', 'e', 'PaperOrientation', 'Landscape', 'PaperUnits', 'Normalized', 'PaperPosition', [0 0 1 1]); %Papertype 'e' es el m�s cuadradito que encontr�, normalized normaliza las paperposition
    saveas(gcf,['HbpHval_' Mapa],'pdf')
    close

    %% Plotea plano Hval-Hbpw

  	Datos1 = Hval; %Carga el nombre del cuantificador y lo pasa a string
    Datos2 = Hbpw;
    
    meanDatos1 = mean(Datos1);
    meanDatos2 = mean(Datos2);
    
    NombreCuanti1 = 'H_{Val}'; %Carga el nombre del cuantificador y lo pasa a string
    NombreCuanti2 = 'H_{BPW}';
    
    figure('position',PosicPlot,'visible','off'); %Abro una figura con tama�o y posici�n
    hold on; grid on; box on
    %title('rand-ellip; roja=pasa-altos; azul=pasa-bajos; negro=rand-matlab')
    xlabel(NombreCuanti1); ylabel(NombreCuanti2)
    ejes=newplot; %newplot returns handle of current axes y abre una figura
    set(ejes,'FontName','Arial','FontWeight','bold','FontSize',FuenteSize); %cambiamos las letras de los ejes de la figura
    set(gcf,'DefaultLineLineWidth',LineaAncho2); %fijamos el tamano de linea por default grosor 2
    axis([0 1 0 1]);

    plot(meanDatos1,meanDatos2, '.:k', 'LineWidth', LineaAncho1, 'MarkerSize', PuntoSize)
    plot(Datos1,Datos2,'.r', 'MarkerSize', MarcaSize)
    
    set(gcf, 'PaperType', 'e', 'PaperOrientation', 'Landscape', 'PaperUnits', 'Normalized', 'PaperPosition', [0 0 1 1]); %Papertype 'e' es el m�s cuadradito que encontr�, normalized normaliza las paperposition
    saveas(gcf,['HbpwHval_' Mapa],'pdf')
    close
    
%% Plotea plano HBP-CBP

  	Datos1 = Hbp; %Carga el nombre del cuantificador y lo pasa a string
    Datos2 = Cbp;
    
    meanDatos1 = mean(Datos1);
    meanDatos2 = mean(Datos2);
    
    NombreCuanti1 = 'H_{BP}'; %Carga el nombre del cuantificador y lo pasa a string
    NombreCuanti2 = 'C_{BP}';
    
    figure('position',PosicPlot,'visible','off'); %Abro una figura con tama�o y posici�n
    hold on; grid on; box on
    %title('rand-ellip; roja=pasa-altos; azul=pasa-bajos; negro=rand-matlab')
    xlabel(NombreCuanti1); ylabel(NombreCuanti2)
    ejes=newplot; %newplot returns handle of current axes y abre una figura
    set(ejes,'FontName','Arial','FontWeight','bold','FontSize',FuenteSize); %cambiamos las letras de los ejes de la figura
    set(gcf,'DefaultLineLineWidth',LineaAncho2); %fijamos el tamano de linea por default grosor 2
    axis([0 1 0 .5]);

    plot(meanDatos1,meanDatos2, '.:k', 'LineWidth', LineaAncho1, 'MarkerSize', PuntoSize)
    plot(Datos1,Datos2,'.r', 'MarkerSize', MarcaSize)
    plot(CotasInferiores(:,1), CotasInferiores(:,2), '-.', 'LineWidth', LineaAncho2, 'Color', [.5 .5 1])
   	plot(CotasSuperiores(:,1), CotasSuperiores(:,2), '-.', 'LineWidth', LineaAncho2, 'Color', [.5 .5 1])
    
    set(gcf, 'PaperType', 'e', 'PaperOrientation', 'Landscape', 'PaperUnits', 'Normalized', 'PaperPosition', [0 0 1 1]); %Papertype 'e' es el m�s cuadradito que encontr�, normalized normaliza las paperposition
    saveas(gcf,['CbpHbp_' Mapa],'pdf')
    close
    
%% Plotea plano HBPW-CBPW

  	Datos1 = Hbpw; %Carga el nombre del cuantificador y lo pasa a string
    Datos2 = Cbpw;
    
    meanDatos1 = mean(Datos1);
    meanDatos2 = mean(Datos2);
    
    NombreCuanti1 = 'H_{BPW}'; %Carga el nombre del cuantificador y lo pasa a string
    NombreCuanti2 = 'C_{BPW}';
    
    figure('position',PosicPlot,'visible','off'); %Abro una figura con tama�o y posici�n
    hold on; grid on; box on
    %title('rand-ellip; roja=pasa-altos; azul=pasa-bajos; negro=rand-matlab')
    xlabel(NombreCuanti1); ylabel(NombreCuanti2)
    ejes=newplot; %newplot returns handle of current axes y abre una figura
    set(ejes,'FontName','Arial','FontWeight','bold','FontSize',FuenteSize); %cambiamos las letras de los ejes de la figura
    set(gcf,'DefaultLineLineWidth',LineaAncho2); %fijamos el tamano de linea por default grosor 2
    axis([0 1 0 .5]);

    plot(meanDatos1,meanDatos2, '.:k', 'LineWidth', LineaAncho1, 'MarkerSize', PuntoSize)
    plot(Datos1,Datos2,'.r', 'MarkerSize', MarcaSize)
    plot(CotasInferiores(:,1), CotasInferiores(:,2), '-.', 'LineWidth', LineaAncho2, 'Color', [.5 .5 1])
   	plot(CotasSuperiores(:,1), CotasSuperiores(:,2), '-.', 'LineWidth', LineaAncho2, 'Color', [.5 .5 1])
    
    set(gcf, 'PaperType', 'e', 'PaperOrientation', 'Landscape', 'PaperUnits', 'Normalized', 'PaperPosition', [0 0 1 1]); %Papertype 'e' es el m�s cuadradito que encontr�, normalized normaliza las paperposition
    saveas(gcf,['CbpwHbpw_' Mapa],'pdf')
    close