clc; clear; close all;
Folder = 'Input_data/';
Spacing = 1;
cl = {'#8D0E25','#FFA500','k','g','b','w','r','m','c','y','#CAAEC2'};
% C =    -5        -4     -3  -2  -1   0  +1  +2  +3  +4     +5
FigureFilename = {'PPDiagram_Indiv_SBG_tau0.50_0.00_tpp_0.70_t2_0.56_M0_0.00_u0.90_v0.50_w0.50';
                  'PPDiagram_Indiv_SBG_tau0.50_0.00_tpp_0.70_t2_0.56_M0_0.00_u0.90_v0.50_w0.50';
                  'PPDiagram_VB1_SBG_tau0.00_1.00_tpp_0.70_t2_0.30_M0_0.00_u0.50_v0.20_w0.20';
                  'PPDiagram_VB2_SBG_tau0.00_1.00_tpp_0.70_t2_0.30_M0_0.00_u0.50_v0.20_w0.20';
                  'PPDiagram_Indiv_SBG_tau1.00_0.00_tpp_0.70_t2_0.10_M0_0.00_u0.00_v0.70_w0.00';
                  'PPDiagram_Indiv_SBG_tau1.00_0.00_tpp_0.70_t2_0.10_M0_0.00_u0.00_v0.70_w0.60'};
Title = {'\mathcal{M}_{V_1}','\mathcal{M}_{V_2}'};
Config = {'C2','','C5','','C3'};
%%
figure('outerposition',1E3*[-0.0062,0.0338,1.5504,0.8384]);
SubfigPos = [0.305,0.695,0.2,0.3;
             0.519,0.695,0.2,0.3;
             0.305,0.405,0.2,0.3;
             0.519,0.405,0.2,0.3;
             0.305,0.066,0.2,0.3;
             0.519,0.066,0.2,0.3];
SubfigPos(:,3) = 0.139;
SubfigPos(:,4) = 0.268;
SubfigPos([2,4,6],1) = SubfigPos([2,4,6],1)-0.068; % x-coordinates of the column-2 subfigs 
SubfigPos(:,1) = SubfigPos(:,1)+0.0680; % x-coordinates of all subfigs
SubfigPos(1:4,2) = SubfigPos(1:4,2)-0.013; % y-coordinates of the first 4 subfigs
SubfigPos(:,2) = SubfigPos(:,2)+0.003; % y-coordinates of all subfigs
for i = 1:6
    if(ismember(i,[1,3,4,5,6]))
        Data = load([Folder,FigureFilename{i},'.txt']);
        X = Data(:,1); Y = Data(:,2);
    end
    if(i == 1)||(i == 2)
        C = Data(:,i+2);
    elseif(i == 3)||(i == 4)
        C = Data(:,3);
        clear Data;
    else
        C = Data(:,5);
        clear Data;
    end
    Ntp = size(C,1);
    Np = sqrt(Ntp);
    if(ismember(i,[1,2,5,6]))
        for k = 1:Np
            C((Np-1)*k+1) = 0; % Remove the line phi_1 = -phi_2
        end
    end
    subplot('Position',SubfigPos(i,:)); % 2 rows, 2 columns
    hold on; grid off; axis equal; box on; set(gca,'FontSize',12);
    xlim([-pi,pi]); ylim([-pi,pi]);
    if(i == 5)
        title('$\mathcal{M}_V$','FontSize',14,'interpreter','latex');
    end
    if(i == 6)
        title('$\mathcal{M}_V$','FontSize',14,'interpreter','latex');
    end
    %
    for iP = 1:Ntp
        CiP = C(iP);
        if(abs(CiP) > 0)&&(abs(CiP) <= 5)
            plot(X(iP),Y(iP),'.','Color',cl{CiP+6},'MarkerSize',1.3);
        end
    end
    %
    if(i == 1)||(i == 3)||(i == 5)
        yticks([-pi,-2,0,2,pi]); yticklabels({'-\pi','-2','0','2','\pi'});
        ylabel('$\phi_2$','FontSize',14,'interpreter','latex');
        text(-4.4,2.8,['$\mathbf{',Config{i},'}$'],'VerticalAlignment','middle','HorizontalAlignment','center',...
                                                 'Color','k','FontSize',14,'Interpreter','latex');
    end
    if(i == 1)
        text(-4.43,3.8,'$\mathbf{(a)}$','VerticalAlignment','middle','HorizontalAlignment','center',...
                                                 'Color','k','FontSize',14,'Interpreter','latex');
    elseif(i == 5)
        text(-4.43,3.8,'$\mathbf{(b)}$','VerticalAlignment','middle','HorizontalAlignment','center',...
                                                 'Color','k','FontSize',14,'Interpreter','latex');
    end
    if(i == 5)||(i == 6)
        xticks([-pi,-2,0,2,pi]); xticklabels({'-\pi','-2','0','2','\pi'});
        xlabel('$\phi_1$','FontSize',14,'interpreter','latex');
    end
    set(gca,'TickLength',[0.02,0.02]);
    set(gca,'LineWidth',0.9);
    if(i == 1)||(i == 2)
        title(['$',Title{i},'$'],'FontSize',14,'interpreter','latex');
    end
    if(ismember(i,[1,2,3,4]))
        xticklabels({});
    end
    if(i == 2)||(i == 4)||(i == 6)
        yticklabels({});
    end
    rectangle('Position',[-pi,-pi,2*pi,2*pi],'LineWidth',0.9);
end
set(gcf,'PaperSize',[5.03,7.71]);                % Set page size
print(gcf,'Output_figures/Fig_4.pdf','-dpdf','-vector');
close(gcf);