clc; clear; close all;
Folder = 'Input_data/';
Spacing = 1;
Config = {'Chern number $C=0$','Chern number $C=+3$'};
DataFile = {'C3_ChernNumber_0_01';
            'C3_ChernNumber_p3_02';
            'SBG_tau1.00_0.00_tpp_0.70_t2_0.10_IP_-2.60_-0.20_u0.00_v0.70_w0.00';
            'SBG_tau1.00_0.00_tpp_0.70_t2_0.10_IP_-2.60_-0.20_u0.00_v0.70_w0.60'};
%%
figure('outerposition',1E3*[-0.0062,0.0338,1.5504,0.8384]);
SubfigPos = [0.200,0.675,0.3,0.290;
             0.500,0.675,0.3,0.290;
             0.200,0.001,0.3,0.669;
             0.500,0.001,0.3,0.669];
% Shifting the 2nd column horizontally
dX_2ndColumn = -0.100;
SubfigPos([2,4],1) = SubfigPos([2,4],1)+dX_2ndColumn;
% Adjust height between the band & A(k) plots:
dYrelative = 0.0223;
SubfigPos(1:2,2) = SubfigPos(1:2,2)-dYrelative;
SubfigPos(1:2,4) = SubfigPos(1:2,4)+dYrelative;
SubfigPos(3:4,4) = SubfigPos(3:4,4)-dYrelative;
%% BZ meshing
% Determining the ka range:
aCC = 0.142; % [nm]
dGG = 0.335; % [nm]
a   = aCC*sqrt(3);
b01 = 4*pi/(sqrt(3)*a)*[cos(pi/3); +sin(pi/3)];
b02 = 4*pi/(sqrt(3)*a)*[cos(pi/3); -sin(pi/3)];
K01 = (2*b01+b02)/3;
% krange = [kxrange; kyrange]:
krange = [0,2*K01(1); 0,b01(2)]; % Rhombus BZ
Nk = [151; 171]*8; % Nk = [Nkx; Nky]
    dk = (krange(:,2)-krange(:,1))./(Nk-1); % dk = [dkx; dky]
    kx1 = (krange(1,1):dk(1):krange(1,2));
    ky1 = (krange(2,1):dk(2):krange(2,2));
    [kx,ky] = meshgrid(kx1,ky1);
    kx = kx+ky*sqrt(3)/3; % Rhombus BZ only
clear kx1 ky1;
%% Plot the band structures
% The energy range plotted:
Elim = [-1,1];
dE_threshold = 5E-3; % for finding the touching points
%
Eplot = zeros(Nk(2),Nk(1),4);
for i_config = 1:2
    Ek = load([Folder,DataFile{i_config},'.txt']);
    ik  = 0;
    for ikx = 1:Nk(1)
        for iky = 1:Nk(2)
            ik = ik+1;
            Eplot(iky,ikx,:) = Ek(ik,:)';
        end
    end
    clear Ek;
    axes('Position',SubfigPos(i_config,:));
    hold on; axis off;
    for i = 1:4
	    Ei = Eplot(:,:,i);
        Ei(Ei>Elim(2)) = NaN; Ei(Ei<Elim(1)) = NaN;
        Eplot(:,:,i) = Ei;
	    h = mesh(kx,ky,Ei);
        if(i == 2)||(i == 3)
            h.FaceAlpha = 0.1;       % transparency between 0 (fully transparent) and 1 (opaque)
            h.EdgeAlpha = 0.1;       % transparency for mesh edges
        end
        if(i >= 2)
            for j = 1:i-1
                dE = Eplot(:,:,i)-Eplot(:,:,j);
                [row,col] = find(dE < dE_threshold);
                for k = 1:size(row,1)
                    kxj = kx(row(k),col(k));
                    kyj = ky(row(k),col(k));
                    plot3(kxj,kyj,Eplot(row(k),col(k),j),'.r','MarkerSize',10);
                end
            end
        end
    end
    % Plot the boundary of the rhombus BZ:
        xBZ = [krange(1,:),krange(1,2)+krange(2,2)*sqrt(3)/3,krange(1,1)+krange(2,2)*sqrt(3)/3,0];
        yBZ = krange(2,2)*[0,0,1,1,0];
        plot3(xBZ,yBZ,Elim(1)*ones(1,5),'--k','LineWidth',1);
        plot3(xBZ,yBZ,Elim(2)*ones(1,5),'--k','LineWidth',1);
    % Mark the K point:
    K02 = 2*K01;
        plot3(K01(1),K01(2),Elim(1),'.m','MarkerSize',8);
        plot3(K02(1),K02(2),Elim(1),'.m','MarkerSize',8);
        plot3(K01(1),K01(2),Elim(2),'.m','MarkerSize',8);
        plot3(K02(1),K02(2),Elim(2),'.m','MarkerSize',8);
        plot3(K01(1)*[1,1],K01(2)*[1,1],Elim,'--m','LineWidth',1);
        plot3(K02(1)*[1,1],K02(2)*[1,1],Elim,'--m','LineWidth',1);
        text(K01(1),K01(2),Elim(2),'$\mathbf{K}$','VerticalAlignment','bottom',...
             'HorizontalAlignment','center','Color','b','FontSize',17,'Interpreter','latex');
        text(K02(1),K02(2),Elim(2),'$\mathbf{K^\prime}$','VerticalAlignment','bottom',...
             'HorizontalAlignment','center','Color','b','FontSize',17,'Interpreter','latex');
    view(0,15); daspect([1 1 0.11]);
    xlim([0,3*K01(1)]); ylim(krange(2,:)); zlim([Elim(1),Elim(2)+0.08]);
    title(Config{i_config},'FontSize',18,'interpreter','latex');
end
clear kx ky Eplot Ei dE;
%return;
%% Plot the A(k) field
MaxNorm = 2.8;
X = [0,b01(1),3*b01(1),2*b01(1)];
Y = [0,b01(2),b01(2),0];
NVertices = [0,0,2,5];
Vertices = zeros(5,3,4);
% Vertices(i,:,k) is [X_vertex,Y_vertex,sign(B(k))] of i-th vertex in k-th A(k) subfig. where:
% sign(B(k)) = 1 -> anti-clockwise vertex (blue), sign(B(k)) = 2 -> clockwise vertex (cyan)
Bsign = 'bc'; % color code for sign(B(k))
Vertices(1:2,:,3) = [norm(b01)/2,b01(2)*1/3,2;  % The VB2 for C =  0
                     norm(b01)  ,b01(2)*2/3,1];
Vertices(:,:,4) = [norm(b01)/2,b01(2)*1/3,2;    % The VB2 for C = +3
                   norm(b01)  ,b01(2)*2/3,1;
                   27.7485    ,17.9810,2;
                   31.2485    ,17.9810,2;
                   norm(b01)  ,15.0397,2];
for i_config = 3:4
    A = load([Folder,DataFile{i_config},'.txt']);
    k = A(:,1:2);
    A = A(:,3:6);
    Nk = size(k,1);
    Ak = zeros(2,2,Nk);
    NormAk = zeros(Nk,2);
    for ik = 1:Nk
        Aik = reshape(A(ik,:),2,2); % Ak(ik) = [Ak_band1, Ak_band2] where Ak_bandi = [Ak_x; Ak_y]
        Ak(:,:,ik) = Aik;
        NormAk(ik,:) = [norm(Aik(:,1)),norm(Aik(:,2))];
    end
    Ak = Ak*MaxNorm/max(NormAk,[],"all");
    clear A NormAk;
    axes('Position',SubfigPos(i_config,:));
    hold on; axis equal; axis off;
    xlim([0,3*b01(1)]); ylim([0,b01(2)+krange(2,2)+0.5]);
    plot(b01(1)*[1,1],[0,b01(2)+krange(2,2)+0.5],'--','Color','m','LineWidth',1);
    plot(2*b01(1)*[1,1],[0,b01(2)+krange(2,2)+0.5],'--','Color','m','LineWidth',1);
    for iband = [1,2]
        Y = Y+(krange(2,2)+0.5)*(iband-1);
        fill(X,Y,[0.5,0.5,0.5],'FaceAlpha',0.1,'LineStyle','none');
        An = Ak(:,iband,:);
        k(:,2) = k(:,2)+(krange(2,2)+0.5)*(iband-1);
        for ik = 1:Nk
            quiver(k(ik,1),k(ik,2),An(1,1,ik),An(2,1,ik),0,...
                    'MaxHeadSize',0.5,'Color','r','LineWidth',0.2);
        end
        text(b01(1)*3/2,b01(2)+(krange(2,2)+0.5)*(iband-1)-1,['$\mathcal{M}_{V_',num2str(iband),'}$'],...
        'VerticalAlignment','top','HorizontalAlignment','center','Color','k','FontSize',16,'Interpreter','latex');
    end
    if(iband == 2)
        for i = 1:NVertices(i_config)
            Vertex = Vertices(i,:,i_config);
            Vertex(2) = Vertex(2)+(krange(2,2)+0.5)*(iband-1);
            plot(Vertex(1),Vertex(2),'.','Color',Bsign(Vertex(3)),'MarkerSize',18);
        end
    end
end
%%
print(gcf,'Output_figures/Fig_5.png','-dpng','-r700');
close(gcf);