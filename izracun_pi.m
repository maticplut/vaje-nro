function [] = izracun_pi(stTock)
  
    % Kliči funkcijo za oceno π s trenutnim številom točk
    [zKrog, zKvad] = mcc_pi(stTock);
    [ocenjenoPi, napaka] = area_pi(stTock);
    [devijacija, stTock2, rezultati] = calc_pi();

    % Izpis rezultata
    disp(['Ocenjeno π: ', num2str(ocenjenoPi,5)]);
    disp(['Napaka: ', num2str(napaka,5)]);

    % Nariši podgraf s točkami
    subplot(1, 2, 1)
    graf();
    hold on
    plot(zKvad(:,1),zKvad(:,2),'r.');
    hold on
    plot(zKrog(:,1),zKrog(:,2),'b.');
    hold on
    
    % Nariši drug podgraf
    subplot(1, 2, 2)
    plot(stTock2, rezultati(:, 2), 'b.');
    hold on;
    plot(stTock2, pi * ones(size(stTock2)), 'r--');
    hold on 
    plot(stTock2, devijacija,"g-",'LineWidth', 1.5);
    hold on

end

function [devijacija, stTock2, rezultati] = calc_pi()
    % Naraščajoče število naključnih točk
    stTock2 = 500:500:50000;
    
    rezultati = zeros(length(stTock2), 3);
    rez = zeros(length(stTock2),20);
    devijacija = zeros(1, length(stTock2));
    
    for j = 1:20
        for i = 1:length(stTock2)
            % Kliči funkcijo za oceno π s trenutnim številom točk
            [ocenjenoPi, napaka] = area_pi(stTock2(i));
    
            % Shrani rezultate
            rezultati(i, 1) = stTock2(i);
            rezultati(i, 2) = ocenjenoPi;
            rezultati(i, 3) = napaka;
            rez(i,j)= ocenjenoPi;
        end
    end
  
    for i = 1:length(stTock2)
        devijacija(i) =sqrt(sum((rez(i,:)-pi()).^2)/20) + pi(); 
    end

end




function [ocenjenoPi, napaka] = area_pi(stTock)
    
    [zKrog,~] = mcc_pi(stTock);

    % Izračun ocenjenega π in napake
    ocenjenoPi = 4 * (nnz(zKrog) / 2) / stTock;
    napaka = abs(ocenjenoPi - pi);
end

function [zKrog, zKvad] = mcc_pi(stTock)


    zKrog = zeros(stTock,2);
    zKvad = zeros(stTock,2);
    
    % Generiranje naključnih točk
    for i = 1:stTock
        x = 2 * rand() - 1;
        y = 2 * rand() - 1;

        if x^2 + y^2 <= 1
            zKrog(i,:) = [x,y];
        else
            zKvad(i,:) = [x,y];
        end
    end  
end

function [] = graf()
    % Izris korznice
    theta = linspace(0, 2 * pi, 1000);

    x2 = cos(theta);
    y2 = sin(theta);

    subplot(1, 2, 1)
    plot(x2, y2, 'g', 'LineWidth', 1);
    axis equal;
    hold on

    %Izris kvadrata
    T0=[1,1];
    T1=[1,-1];
    T2=[-1,-1];
    T3=[-1,1];
   
    plot(T1,T0,"k",'LineWidth', 1)
    hold on
    plot(T1,T2,"k",'LineWidth', 1)
    hold on
    plot(T0,T3,"k",'LineWidth', 1)
    hold on
    plot(T2,T3,"k",'LineWidth', 1)
    hold on
end