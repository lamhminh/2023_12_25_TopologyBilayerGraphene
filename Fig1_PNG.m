clc; clear; close all;
aCC = 0.142; % [nm]
dGG = 0.335; % [nm]
a   = aCC*sqrt(3);
%% Coordinates of basis vectors in the Cartesian basis:
a1    = a*[sqrt(3)/2; +1/2; 0];
a2    = a*[sqrt(3)/2; -1/2; 0];
d1    = (a1 + a2)/3;
d     = dGG*[0; 0; 1];
tau_x = d1;
tau_y = (a1 - a2)/2;
% Coordinates of layer-1 atoms in the basis {a1,a2,d1,d}:
c_atom = [-1,-1, 0.5,0;
             0,-1,-0.5,0;
             0,-1,+0.5,0;
             0, 0,-0.5,0;
            -1, 0,+0.5,0;
            -1, 0,-0.5,0;
             1,-1,-0.5,0;
             1,-1,+0.5,0;
             1, 0,-0.5,0;
             0, 0,+0.5,0;
             0, 1,-0.5,0;
            -1, 1,+0.5,0;
            -1, 1,-0.5,0;
             1, 0,+0.5,0;
             1, 1,-0.5,0;
             0, 1,+0.5,0]';
% Coordinates of layer-1 atoms in the Cartesian basis:
R_MLG = [a1,a2,d1,d]*c_atom;
% Atom bonds:
Bonds = [1,2; 2,3; 3,4; 4,5; 5,6; 6,1;
         3,7; 7,8; 8,9; 9,10; 10,11; 11,12; 12,13; 13,5;
         10,4; 9,14; 14,15; 15,16; 16,11];
% Coordinates of sliding vector in the basis {tau_x,tau_y}:
c_tau = [0,0; 0.5,0; 1,0; 0,0.5; 0,1]';
% Coordinates of sliding vector in the Cartesian basis:
tau = [tau_x,tau_y]*c_tau;
%% Plot the structure of 5 SL configurations
Na_MLG    = size(R_MLG,2);
Rl        = zeros(Na_MLG,3);
Color     = {'b','r'};
tau_head  = [1, 1, 1, 4, 4];
HeadSize  = [0, 1.5, 0.8, 2, 0.9];
tau_label = [0,0; -0.34,0.04; -0.35,0.04; -0.045,0.01; -0.045,0.01];
SubfigPos = [2, 0.17, 0.32,0.15,0.4;         % C1 lattice
             3, 0.35,0.51,0.15,0.4;   % C2 lattice
             5, 0.62,0.51,0.15,0.4;   % C3 lattice
             9, 0.35, 0.13,0.15,0.4;   % C4 lattice
             11,0.62, 0.13,0.15,0.4];  % C5 lattice
SubfigPos(:,4) = 0.16;
SubfigPos(:,5) = 0.35;
Sublabel  = {'$\mathbf{C1}$','$\mathbf{C2}$','$\mathbf{C3}$','$\mathbf{C4}$','$\mathbf{C5}$'};
figure('outerposition',1E3*[-0.0062,0.0338,1.5504,0.8384]);
radius = [0.017,0.006]; % [ball,stick]
for i_config = 1:5
    subplot('Position',SubfigPos(i_config,2:5));
    hold on; axis equal; axis off;
    % Label the subfigure
    text(0,0.42,Sublabel{i_config},'VerticalAlignment','middle','HorizontalAlignment','center',...
                                       'Color','k','FontSize',18,'Interpreter','latex');
    % Unit cell
    X_Unitcell = [-1,0,1,0,-1]*aCC*3/2;
    Y_Unitcell = [0,1,0,-1,0]*aCC*sqrt(3)/2;
    Z_Unitcell = [1,1,1,1,1]*dGG;
    line(X_Unitcell,Y_Unitcell,Z_Unitcell,'Color',[0,0.6,0],'LineWidth',1.5);
    %
    for l = 1:2 % layer index
        % Atom lattice
        for j = 1:Na_MLG
            Rl(j,:) = (R_MLG(:,j) + (-1)^l*0.5*(d + tau(:,i_config)))';
        end
        ball_and_stick_model(Rl,Bonds,Color{l},radius);
        if(i_config == 1)
            % The two Cartesian axes
            quiver(-0.41,0,0.86,0,0,'MaxHeadSize',0.15,'Color','k','LineWidth',1.1);
            quiver(0,-0.32,0,0.64,0,'MaxHeadSize',0.2,'Color','k','LineWidth',1.1);
            % Label the two Cartesian axes
            text(0.44,0.04,'$x$','VerticalAlignment','middle','HorizontalAlignment','center','Color','k','FontSize',16,'Interpreter','latex');
            text(0.04,0.32,'$y$','VerticalAlignment','middle','HorizontalAlignment','center','Color','k','FontSize',16,'Interpreter','latex');
            % The two unit vectors
            quiver3(-aCC*3/2,0,dGG*1.1,a1(1),a1(2),0,0,'MaxHeadSize',0.5,'Color','b','LineWidth',2);
            quiver3(-aCC*3/2,0,dGG*1.1,a2(1),a2(2),0,0,'MaxHeadSize',0.5,'Color','b','LineWidth',2);
            % Label the two unit vectors
            text(-0.03,+aCC*0.7,'$\mathbf{a}_1$','VerticalAlignment','bottom','HorizontalAlignment','right','Color','b','FontSize',15,'Interpreter','latex');
            text(-0.03,-aCC*0.7,'$\mathbf{a}_2$','VerticalAlignment','top','HorizontalAlignment','right','Color','b','FontSize',15,'Interpreter','latex');
            % The vector connecting 2 sublattices
            quiver3(Rl(4,1),Rl(4,2),dGG,d1(1),d1(2),0,0,'MaxHeadSize',0.8,'Color','k','LineWidth',2.5);
            % Label the vector connecting 2 sublattices
            text(0.04,0.05,'$\mathbf{r}_1$','VerticalAlignment','middle','HorizontalAlignment','center','Color','k','FontSize',15,'Interpreter','latex');
        elseif(l == 1)
            % The sliding vector
            H = Rl(tau_head(i_config),:);
            quiver3(H(1),H(2),dGG,tau(1,i_config),tau(2,i_config),0,0,'MaxHeadSize',HeadSize(i_config),'Color','m','LineWidth',2);
            % Label the sliding vector
            text(tau_label(i_config,1),tau_label(i_config,2),'𝛕','VerticalAlignment','middle','HorizontalAlignment','center','Color','m','FontSize',16);
        end
        if(i_config == 3)
            text(Rl(2,1)+0.002,1.1*Rl(2,2),['$\mathbf{B}_',num2str(l),'$'],'VerticalAlignment','bottom','HorizontalAlignment','right','Color',Color{l},'FontSize',14,'Interpreter','latex');
            text(Rl(1,1)+0.018,Rl(1,2)-0.015,['$\mathbf{A}_',num2str(l),'$'],'VerticalAlignment','top','HorizontalAlignment','right','Color',Color{l},'FontSize',14,'Interpreter','latex');
        end
    end
    xlim([-0.45,0.45]);
    ylim([-0.35,0.35]);
end
%% Plot the band structures
% Determining the ka range:
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
% The energy range plotted:
Elim = [-0.4,0.3];
dE_threshold = 5E-3;
%
SubfigPos = [1,0.06, 0.31,0.15,0.4;  % C1 bands
                    4, 0.51,0.51,0.15,0.4;   % C2 bands
                    6, 0.78,0.51,0.15,0.4;   % C3 bands
                    10, 0.51, 0.11,0.15,0.4;  % C4 bands
                    12,0.78, 0.11,0.15,0.4];  % C5 bands
SubfigPos(:,4) = 0.11;
SubfigPos(:,5) = 0.40;
Folder = 'Input_data/';
EnergyDataFile = {'RealNNN_01_03';
                  'RealNNN_05_09';
                  'RealNNN_12_09';
                  'RealNNN_19_03';
                  'RealNNN_22_03'};
Eplot = zeros(Nk(2),Nk(1),4);
for i_config = 1:5
    Ek = load([Folder,EnergyDataFile{i_config},'.txt']);
    i  = 0;
    for ikx = 1:Nk(1)
        for iky = 1:Nk(2)
            i = i+1;
            Eplot(iky,ikx,:) = Ek(i,:)';
        end
    end
    clear Ek;
    subplot('Position',SubfigPos(i_config,2:5)); % 2 row, 6 columns, index i_config
    hold on; axis equal; axis off;
    for i = 1:4
	    Ei = Eplot(:,:,i);
        Ei(Ei>Elim(2)) = NaN; Ei(Ei<Elim(1)) = NaN;
        Eplot(:,:,i) = Ei;
	    h = mesh(kx,ky,Ei);
        if(i == 2)
            h.FaceAlpha = 0.2;       % transparency between 0 (fully transparent) and 1 (opaque)
            h.EdgeAlpha = 0.2;       % transparency for mesh edges
        elseif(i == 3)
            h.FaceAlpha = 0.2;
            h.EdgeAlpha = 0.2;
        elseif(i == 4)
            h.FaceAlpha = 0.6;
            h.EdgeAlpha = 0.6;
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
    % Mark the K point:
    quiver3(K01(1),K01(2),Elim(1)-0.05,0,0,Elim(2)-Elim(1)+0.15,0,'MaxHeadSize',0.4,'Color','r','LineWidth',1.5);
    text(K01(1)-0.1,K01(2)-0.1,Elim(2)+0.1,'$E$','VerticalAlignment','middle','HorizontalAlignment','right','Color','r','FontSize',16,'Interpreter','latex');
    daspect([1 1 0.08]);
    if(i_config <= 3)
        quiver3(K01(1)-2.5,K01(2),0,5,0,0,0,'MaxHeadSize',0.5,'Color','r','LineWidth',1.5);
        text(K01(1)+2.5,K01(2),0,'$k_x$','VerticalAlignment','bottom','HorizontalAlignment','center','Color','r','FontSize',16,'Interpreter','latex');
        view(0,20);
    else
        quiver3(K01(1),K01(2)-2.5,0,0,5,0,0,'MaxHeadSize',0.5,'Color','r','LineWidth',1.5);
        text(K01(1),K01(2)+2.5,0,'$k_y$','VerticalAlignment','bottom','HorizontalAlignment','center','Color','r','FontSize',16,'Interpreter','latex');
        view(90,20);
    end
    xlim(K01(1)+3*[-1,1]); ylim(K01(2)+3.7*[-1,1]); zlim([Elim(1)-0.05,Elim(2)+0.1]);
end
%% Plot 2 arrows connecting the C1, C2 & C4 configs
annotation('arrow',[0.18,0.252]+0.1,[0.6,0.69],'LineWidth',2.5);
annotation('arrow',[0.18,0.252]+0.1,[0.4,0.31],'LineWidth',2.5);
annotation('textbox',[0.165+0.1,0.67,0.1,0.05],'String','𝛕 = 𝜉𝛕_1','EdgeColor','none','HorizontalAlignment','center', ...
           'FontWeight','bold','Color','m','FontSize',18);
annotation('textbox',[0.165+0.1,0.29,0.1,0.05],'String','𝛕 = 𝜁𝛕_2','EdgeColor','none','HorizontalAlignment','center', ...
           'FontWeight','bold','Color','m','FontSize',18);
%% Print the figure to PNG
print(gcf,'Output_figures/Fig_1.png','-dpng','-r500');
close(gcf);
%%
function ball_and_stick_model(points,connections,Color,radius)
    % === Parameters ===
    ball_radius = radius(1);
    stick_radius = radius(2);
    ball_color = Color;
    stick_color = Color;
    % === Plot Setup ===
    camlight;
    lighting gouraud;
    % === Plot Balls ===
    [XS, YS, ZS] = sphere(20); % resolution = 20
    for i = 1:size(points, 1)
        c = points(i, :);
        surf(ball_radius*XS + c(1), ...
             ball_radius*YS + c(2), ...
             ball_radius*ZS + c(3), ...
             'FaceColor', ball_color, ...
             'EdgeColor', 'none');
    end
    % === Plot Sticks ===
    for i = 1:size(connections, 1)
        idx1 = connections(i,1);
        idx2 = connections(i,2);
        p1 = points(idx1, :);
        p2 = points(idx2, :);
        draw_cylinder_between(p1, p2, stick_radius, stick_color);
    end
end
%%
function draw_cylinder_between(p1, p2, radius, color)
    % Create unit cylinder along z-axis
    [X, Y, Z] = cylinder(radius, 20); % resolution = 20
    height = norm(p2 - p1);
    Z = Z * height;
    % Compute rotation
    dir = (p2 - p1) / height;
    z_axis = [0 0 1];
    if norm(cross(z_axis, dir)) < 1e-10
        R = eye(3); % already aligned
    else
        rot_axis = cross(z_axis, dir);
        rot_axis = rot_axis / norm(rot_axis);
        angle = acos(dot(z_axis, dir));
        R = axang2rotm([rot_axis angle]);
    end
    % Apply rotation and translation
    n_pts = numel(X);
    for i = 1:n_pts
        pt = R * [X(i); Y(i); Z(i)];
        X(i) = pt(1) + p1(1);
        Y(i) = pt(2) + p1(2);
        Z(i) = pt(3) + p1(3);
    end
    surf(X, Y, Z, 'FaceColor', color, 'EdgeColor', 'none');
end
%%
function R = axang2rotm(axang)
    % Convert axis-angle to rotation matrix
    u = axang(1:3);
    theta = axang(4);
    ux = u(1); uy = u(2); uz = u(3);
    ct = cos(theta); st = sin(theta); vt = 1 - ct;
    R = [ct+ux^2*vt,    ux*uy*vt-uz*st, ux*uz*vt+uy*st;
         uy*ux*vt+uz*st, ct+uy^2*vt,    uy*uz*vt-ux*st;
         uz*ux*vt-uy*st, uz*uy*vt+ux*st, ct+uz^2*vt];
end