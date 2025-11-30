clc; clear; close all;
Folder = 'Input_data/';
cl = {'#8D0E25','#FFA500','k','g','b','w','r','m','c','y','#CAAEC2'};
% C =    -5        -4     -3  -2  -1   0  +1  +2  +3  +4     +5
FigureFilename = {'PPDiagram_Indiv_SBG_tau1.00_0.00_tpp_0.18_t2_0.1006_M0_0.00_u0.08_v0.18_w0.08';
                  'PPDiagram_Indiv_SBG_tau0.50_0.00_tpp_0.18_t2_0.1006_M0_0.00_u0.14_v0.14_w0.10';
                  'PPDiagram_Indiv_SBG_tau0.00_0.00_tpp_0.18_t2_0.1006_M0_0.00_u0.18_v0.08_w0.08';
                  'PPDiagram_Indiv_SBG_tau0.00_0.50_tpp_0.18_t2_0.1006_M0_0.00_u0.15_v0.12_w0.12';
                  'PPDiagram_Indiv_SBG_tau0.00_1.00_tpp_0.18_t2_0.1006_M0_0.00_u0.10_v0.14_w0.14';
                  'PPDiagram_Indiv_SBG_tau1.00_0.00_tpp_0.50_t2_0.1006_M0_0.00_u0.22_v0.50_w0.22';
                  'PPDiagram_Indiv_SBG_tau0.50_0.00_tpp_0.50_t2_0.1006_M0_0.00_u0.40_v0.40_w0.27';
                  'PPDiagram_Indiv_SBG_tau0.00_0.00_tpp_0.50_t2_0.1006_M0_0.00_u0.50_v0.22_w0.22';
                  'PPDiagram_Indiv_SBG_tau0.00_0.50_tpp_0.50_t2_0.1006_M0_0.00_u0.43_v0.35_w0.35';
                  'PPDiagram_Indiv_SBG_tau0.00_1.00_tpp_0.50_t2_0.1006_M0_0.00_u0.27_v0.40_w0.40'};
Title = {'C3','C2','C1','C4','C5'};
Params = {'t_{12}=1.0','t_{12}=2.8'};
%%
figure('outerposition',1E3*[-0.0062,0.0338,1.5504,0.8384]);
SubfigPos = [0.13,0.51,0.16,0.4;
             0.30,0.51,0.16,0.4;
             0.48,0.51,0.16,0.4;
             0.65,0.51,0.16,0.4;
             0.83,0.51,0.16,0.4;
             0.13,0.13,0.16,0.4;
             0.30,0.13,0.16,0.4;
             0.48,0.13,0.16,0.4;
             0.65,0.13,0.16,0.4;
             0.83,0.13,0.16,0.4];
SubfigPos([1,6],1) = 0.15;
SubfigPos([2,7],1) = 0.31;
SubfigPos([3,8],1) = 0.47;
SubfigPos([4,9],1) = 0.63;
SubfigPos([5,10],1) = 0.79;
SubfigPos(1:5,2) = 0.477;
SubfigPos(6:10,2) = 0.160;
SubfigPos(:,3) = 0.142;
SubfigPos(:,4) = 0.4;
%
dX = -0.026; dY = -0.004;
SubfigPos(:,1) = SubfigPos(:,1)+dX;
SubfigPos(:,2) = SubfigPos(:,2)+dY;
%
for i = 1:10
    Data = load([Folder,FigureFilename{i},'.txt']);
    Ntp = size(Data,1);
    X = Data(:,1); Y = Data(:,2); C = Data(:,5); clear Data;
    Np = sqrt(Ntp);
    for k = 1:Np
        C((Np-1)*k+1) = 0; % Remove the line phi_1 = -phi_2
    end
    if(i==3)||(i==8)
        for k = 1:Np
            if(abs(C(Np*(k-1)+k)) ~= 2)
                C(Np*(k-1)+k) = 0; % Remove part of the line phi_1 = phi_2 for C1 config.
            end
        end
    else
        for k = 1:Np
            if(k <= Np/2)          % Correct the line phi_1 = phi_2 for other configs.
                C(Np*(k-1)+k) =  2;
            else
                C(Np*(k-1)+k) = -2;
            end
        end
    end
    axes('Position',SubfigPos(i,:));
    hold on; grid off; axis equal; box on; set(gca,'FontSize',13);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    xlim([-pi,pi]); ylim([-pi,pi]);
    yticks([-pi,-2,0,2,pi]);
    xticks([-pi,-2,0,2,pi]);
    if(i == 1)||(i == 6)
        yticklabels({'-\pi','-2','0','2','\pi'});
        ylabel('$\phi_2$','FontSize',16,'interpreter','latex');
    end
    if(i >= 6)
        xticklabels({'-\pi','-2','0','2','\pi'});
        xlabel('$\phi_1$','FontSize',16,'interpreter','latex');
    end
    set(gca,'TickLength',[0.02,0.02]);
    set(gca,'LineWidth',1);
    if(i <= 5)
        title(['$\mathbf{',Title{i},'}$'],'FontSize',18,'interpreter','latex');
    end
    %
    for iP = 1:Ntp
        CiP = C(iP);
        if(abs(CiP) > 0)&&(abs(CiP) <= 5)
            plot(X(iP),Y(iP),'.','Color',cl{CiP+6},'MarkerSize',3);
        end
    end
    
    if(i == 3)
        line([-pi,-2.8],[-pi,-2.8],'Color','w','LineWidth',2);
        line([2.8,pi],[2.8,pi],'Color','w','LineWidth',2);
        line([-0.35,0.35],[-0.35,0.35],'Color','w','LineWidth',2);
    end
    if(i == 8)
        line([-pi,-1.88],[-pi,-1.88],'Color','w','LineWidth',2);
        line([1.88,pi],[1.88,pi],'Color','w','LineWidth',2);
        line([-1.26,1.26],[-1.26,1.26],'Color','w','LineWidth',2);
    end
    rectangle('Position',[-pi,-pi,2*pi,2*pi],'LineWidth',1);
end
set(gcf,'PaperSize',[13.1,5.60]);                % Set page size
print(gcf,'Output_figures/Fig_2.pdf','-dpdf','-vector');
close(gcf);