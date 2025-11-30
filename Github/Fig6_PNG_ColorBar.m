clc; clear; close all;
figure('outerposition',1E3*[-0.0062,0.0338,1.5504,0.8384]);
hold on; axis off; axis equal;
xlim([-2,3]); ylim([-1,1]);
Range = [0,1];
Np = 201;
dp = (Range(2)-Range(1))/(Np-1);
for i = 1:Np
    p = Range(1)+dp*(i-1);
    RGB = zeros(1,3);
    RGB(3) = 1.0;
    if(p <= 0.5)
        RGB(1:2) = [0,2*(0.5-p)*1.0];
    else
        RGB(1:2) = [2*(p-0.5)*1.0,0];
    end
    if(ismember(i,[1,101,Np]))
        line([p,p],[-0.015,0.03],'Color',RGB,'LineWidth',1.5);
    else
        line([p,p],[0,0.03],'Color',RGB,'LineWidth',1);
    end
end
text(-0.04,-0.030,'(Left edge)','VerticalAlignment','bottom','HorizontalAlignment','right','Color','k',...
         'FontSize',14,'Interpreter','latex');
text(1.04,-0.030,'(Right edge)','VerticalAlignment','bottom','HorizontalAlignment','left','Color','k',...
         'FontSize',14,'Interpreter','latex');
text(0.5,0.04,'$X_\psi$','VerticalAlignment','bottom','HorizontalAlignment','center','Color','k',...
         'FontSize',15,'Interpreter','latex');
%
text(0,-0.02,'$0$','VerticalAlignment','top','HorizontalAlignment','center','Color','k',...
         'FontSize',14,'Interpreter','latex');
text(1,-0.02,'$1$','VerticalAlignment','top','HorizontalAlignment','center','Color','k',...
         'FontSize',14,'Interpreter','latex');
text(0.5,-0.02,'$0.5$','VerticalAlignment','top','HorizontalAlignment','center','Color','k',...
         'FontSize',14,'Interpreter','latex');
print(gcf,'Output_figures/Fig_6_ColorBar.png','-dpng','-r1000');
close(gcf);