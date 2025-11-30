clc; clear; close all;
Folder = 'Input_data/';
Spacing = 1;
cl = {'#8D0E25','#FFA500','k','g','b','w','r','m','c','y','#CAAEC2'};
% C =    -5        -4     -3  -2  -1   0  +1  +2  +3  +4     +5
FigureFilename = {'PPDiagram_Indiv_SBG_tau1.00_0.00_tpp_0.70_t2_0.10_M0_0.00_u0.08_v0.70_w0.60';
                  'PPDiagram_Indiv_SBG_tau0.50_0.00_tpp_0.70_t2_0.56_M0_0.00_u0.90_v0.50_w0.50';
                  'PPDiagram_Indiv_SBG_tau0.00_0.50_tpp_0.70_t2_0.30_M0_0.00_u0.90_v0.70_w0.70';
                  'PPDiagram_Total_SBG_tau0.00_1.00_tpp_0.70_t2_0.30_M0_0.00_u0.50_v0.20_w0.20'};
Title = {'C3','C2','C4','C5'};
%%
SubfigPos = [0.056,0.169,0.208,0.700;
             0.296,0.169,0.208,0.700;
             0.536,0.169,0.208,0.700;
             0.776,0.169,0.208,0.700];
dX = 0.005;
SubfigPos(1:2,1) = SubfigPos(1:2,1)+dX;
SubfigPos(3:4,1) = SubfigPos(3:4,1)-dX;
dX = 0.010;
SubfigPos(1,1) = SubfigPos(1,1)+dX;
SubfigPos(4,1) = SubfigPos(4,1)-dX;
dX = -0.002;
SubfigPos(:,1) = SubfigPos(:,1)+dX;
dY = 0.002;
SubfigPos(:,2) = SubfigPos(:,2)+dY;
%%
figure('outerposition',1E3*[-0.0062,0.0338,1.5504,0.8384]);
for i = 1:4
    Data = load([Folder,FigureFilename{i},'.txt']);
    Ntp = size(Data,1);
    Np = sqrt(Ntp);
    X = Data(:,1); Y = Data(:,2);
    if(i <= 3) % C3, C2 & C4 configs.
        C = Data(:,5);
        for k = 1:Np
            C((Np-1)*k+1) = 0; % Remove the line phi_1 = -phi_2
            if(k <= Np/2)      % Correct the line phi_1 = phi_2
                C(Np*(k-1)+k) =  2;
            else
                C(Np*(k-1)+k) = -2;
            end
        end
    else
        C = Data(:,3);
    end
    clear Data;
    subplot('Position',SubfigPos(i,:)); % 1 row, 4 columns
    hold on; grid off; axis equal; box on; set(gca,'FontSize',17);
    xlim([-pi,pi]); ylim([-pi,pi]);
    %
    for iP = 1:Ntp
        CiP = C(iP);
        if(abs(CiP) > 0)&&(abs(CiP) <= 5)
            plot(X(iP),Y(iP),'.','Color',cl{CiP+6},'MarkerSize',2.8);
        end
    end
    %
    if(i == 4) % Correct the line phi_1 = phi_2 of C5 config.
        line(3.1*[-1,1],3.1*[-1,1],[1,1],'Color','w','LineWidth',1.7);
        line([-2.865,-0.305],[-2.865,-0.305],[2,2],'Color','m','LineWidth',1.7);
        line([0.305,2.865],[0.305,2.865],[2,2],'Color','g','LineWidth',1.7);
    end
    set(gca,'TickLength',[0.02,0.02]);
    set(gca,'LineWidth',1);
    xlabel('$\phi_1$','FontSize',20,'interpreter','latex');
    xticks([-pi,-2,0,2,pi]); xticklabels({'-\pi','-2','0','2','\pi'});
    yticks([-pi,-2,0,2,pi]);
    if(i == 1)
        yticklabels({'-\pi','-2','0','2','\pi'});
        ylabel('$\phi_2$','FontSize',20,'interpreter','latex');
    else
        yticklabels({});
    end
    title(['$\mathbf{',Title{i},'}$'],'FontSize',21,'interpreter','latex');
    rectangle('Position',[-pi,-pi,2*pi,2*pi],'LineWidth',1);
end
set(gcf,'PaperSize',[15.07,4.38]);                % Set page size
print(gcf,'Output_figures/Fig_3.pdf','-dpdf','-vector');
close(gcf);