clc; clear; close all;
Folder = 'Input_data/';
% Main data:
Band_data = {'C2_04_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.72_phi2_-0.43_ZNR_Energy';
             'C2_04_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.72_phi2_-0.43_ANR_Energy';
             'C2_07_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.34_phi2_-0.85_ZNR_Energy';
             'C2_07_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.34_phi2_-0.85_ANR_Energy';
             'C5_01_BLH_PP_t2_0.30_t12_3.94_M0_0.00_phi1_-1.80_phi2_-0.64_ZNR_Energy';
             'C5_01_BLH_PP_t2_0.30_t12_3.94_M0_0.00_phi1_-1.80_phi2_-0.64_ANR_Energy'};
Localization_data = {'C2_04_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.72_phi2_-0.43_ZNR_StateLocalization';
                     'C2_04_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.72_phi2_-0.43_ANR_StateLocalization';
                     'C2_07_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.34_phi2_-0.85_ZNR_StateLocalization';
                     'C2_07_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.34_phi2_-0.85_ANR_StateLocalization';
                     'C5_01_BLH_PP_t2_0.30_t12_3.94_M0_0.00_phi1_-1.80_phi2_-0.64_ZNR_StateLocalization';
                     'C5_01_BLH_PP_t2_0.30_t12_3.94_M0_0.00_phi1_-1.80_phi2_-0.64_ANR_StateLocalization'};
% Inset data:
Inset_Band_data = {'C2_03_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.72_phi2_-0.43_ZNR_Energy';
                   'None';
                   'C2_06_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.34_phi2_-0.85_ZNR_Energy';
                   'C2_06_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.34_phi2_-0.85_ANR_Energy';
                   'None';
                   'None'};
Inset_Localization_data = {'C2_03_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.72_phi2_-0.43_ZNR_StateLocalization';
                           'None';
                           'C2_06_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.34_phi2_-0.85_ZNR_StateLocalization';
                           'C2_06_BLH_PP_t2_0.56_t12_3.94_M0_0.00_phi1_-2.34_phi2_-0.85_ANR_StateLocalization';
                           'None';
                           'None'};
%% The data files are structured as follows:
%  1) Each row of Band_data files:
%        ka_i  E_i_1  E_i_2  ...  E_i_NEi
%     where the number of levels (NEi) is different between rows such that Emin < E_i_n < Emax.
%  2) Each row of Localization_data files:
%        ka_i  (X_1 X_2 X_3)_i_1  (X_1 X_2 X_3)_i_2  ...  (X_1 X_2 X_3)_i_NEi
%     where NEi is identical to the value in the corresponding Band_data file, and
%           (X_1 X_2 X_3) are the values of layers 1, 2 and their total, respectively.
%     Thus, the degree of localization of a state on n-th band at ka for the l-th layer is the
%     (3*in+ilayer-2)-th value in the ika-th row.
%%
ilayer = 1; % [1,2,3] = [Total,1,2]
Title = {'ZE','AE'};
Config = {'C2','','C2','','C5'};
Ylim = [-0.49,0.55; -0.11,0.21;  % 1st row
        -0.54,1.38; -0.15,0.69;  % 2nd row
         -0.9,0.44; -0.83,0.54]; % 3rd row
%% List of y-tick labels
NYticks = [7,7,6,5,8,7];
% 1st row
Yticks(1,1:NYticks(1)) = [-0.6,-0.4,-0.2,0,0.2,0.4,0.6]+0.045;
Yticklabels(1,1:NYticks(1)) = {'-0.6','-0.4','-0.2','0','0.2','0.4','0.6'};
Yticks(2,1:NYticks(2)) = [-0.15,-0.1,-0.05,0,0.05,0.1,0.15]+0.045;
Yticklabels(2,1:NYticks(2)) = {'-0.15','-0.1','-0.05','0','0.05','0.1','0.15'};
% 2nd row
Yticks(3,1:NYticks(3)) = [-0.8,-0.4,0,0.4,0.8,1.2]+0.026;
Yticklabels(3,1:NYticks(3)) = {'-0.8','-0.4','0','0.4','0.8','1.2'};
Yticks(4,1:NYticks(4)) = [-0.2,0,0.2,0.4,0.6]+0.026;
Yticklabels(4,1:NYticks(4)) = {'-0.2','0','0.2','0.4','0.6'};
% 3rd row
Yticks(5,1:NYticks(5)) = [-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6]+0.02;
Yticklabels(5,1:NYticks(5)) = {'-0.8','-0.6','-0.4','-0.2','0','0.2','0.4','0.6'};
Yticks(6,1:NYticks(6)) = [-0.8,-0.6,-0.4,-0.2,0,0.2,0.4]+0.02;
Yticklabels(6,1:NYticks(6)) = {'-0.8','-0.6','-0.4','-0.2','0','0.2','0.4'};
%%
figure('outerposition',1E3*[-0.0062,0.0338,1.5504,0.8384]);
SubfigPos = [0.305,0.544,0.2,0.4;
             0.547,0.544,0.2,0.4;
             0.305,0.544,0.2,0.4;
             0.547,0.544,0.2,0.4;
             0.305,0.088,0.2,0.4;
             0.547,0.088,0.2,0.4];
InsetAxes = [ 0.368,0.577,0.049, 0.10;
              0.612,0.577,0.045, 0.09;
             0.4205,0.631,0.042,0.146; % C2_ZNR: C=+4
              0.606,0.665,0.058, 0.10; % C2_ANR: C=+4
              0.370,0.121,0.045, 0.09;
              0.606,0.121,0.057, 0.09];
InsetAxes(:,1) = InsetAxes(:,1)+0.0125;
InsetAxes(:,2) = InsetAxes(:,2)-0.015;
InsetAxes([1,2,5,6],2) = InsetAxes([1,2,5,6],2)-0.008;
InsetRanges = [-0.65,0.65, 0.035,0.054;
                -0.2, 0.2, 0.035,0.055;
                2.85, 3.2, 0.006,0.051;    % C2_ZNR: C=+4
                2.45,3.85, 0.013,0.041; % C2_ANR: C=+4
                -0.3, 0.3,-0.033,0.048;
               -0.55,0.55,-0.033,0.048];
% Width & height of each subfigure:
SubfigPos(:,3) = 0.140; % Width
SubfigPos(:,4) = 0.270; % Height
% Shifing each row of subfigs vertically
dY_rows = [0.144,-0.173,-0.032]; % Row 1, 2, 3, respectively
SubfigPos(1:2,2) = SubfigPos(1:2,2)+dY_rows(1);
InsetAxes(1:2,2) = InsetAxes(1:2,2)+dY_rows(1);
SubfigPos(3:4,2) = SubfigPos(3:4,2)+dY_rows(2);
InsetAxes(3:4,2) = InsetAxes(3:4,2)+dY_rows(2);
SubfigPos(5:6,2) = SubfigPos(5:6,2)+dY_rows(3);
InsetAxes(5:6,2) = InsetAxes(5:6,2)+dY_rows(3);
% Shifing all subfigs vertically
dY_allsubfigs = +0.007;
SubfigPos(:,2) = SubfigPos(:,2)+dY_allsubfigs;
InsetAxes(:,2) = InsetAxes(:,2)+dY_allsubfigs;
% Shifting the 2nd column horizontally
dX_2ndColumn = -0.072;
SubfigPos([2,4,6],1) = SubfigPos([2,4,6],1)+dX_2ndColumn;
InsetAxes([2,4,6],1) = InsetAxes([2,4,6],1)+dX_2ndColumn;
% Shifting all the insets horizontally
dX_insets = -0.03;
InsetAxes(:,1) = InsetAxes(:,1)+dX_insets;
% Shifting the whole figure horizontally
dX_all = 0.057;
SubfigPos(:,1) = SubfigPos(:,1)+dX_all;
InsetAxes(:,1) = InsetAxes(:,1)+dX_all;
% Shifting the whole figure vertically
dY_all = -0.0013;
SubfigPos(:,2) = SubfigPos(:,2)+dY_all;
InsetAxes(:,2) = InsetAxes(:,2)+dY_all;
% Font sizes:
FontSizes = [19,19,15,12]-4; % [Title_&_ConfigNames,axes_labels,axes_ticklabels,inset_ticklabels]
%
kaRange = [-pi,pi; -pi,pi; 0,2*pi; 0,2*pi; -pi,pi; -pi,pi];
CNumbers = {'$C=+3$';'$C=+3$';'$C=+4$';'$C=+4$';'$C=+3$';'$C=+3$'};
%%
for i = 1:6
    for j = 1:2
        if(j == 1) % Main subplots
            fid = fopen([Folder,Band_data{i},'.txt'],'r');
            Band = textscan(fid,'%s','Delimiter','\n');
            fclose(fid);
            fid = fopen([Folder,Localization_data{i},'.txt'],'r');
            Edge = textscan(fid,'%s','Delimiter','\n');
            fclose(fid);
            Nk = length(Band{1});
            %
            Subfig = subplot('Position',SubfigPos(i,:)); % 2 row, 2 columns
            %annotation('rectangle',SubfigPos(i,:),'Color','b');
            hold on; grid on; box on; orient(gcf,'landscape'); set(gca,'FontSize',FontSizes(3));
            xlim(kaRange(i,:)); ylim(Ylim(i,:));
            if(ismember(i,[1,2,5,6]))
                xticks([-pi,-2,0,2,pi]);
                xticklabels({'-\pi','-2','0','2','\pi'});
            else
                xticks([0,pi,2*pi]);
                xticklabels({'0','\pi','2\pi'});
            end
            if(i == 1)||(i == 2)
                title(['$\mathbf{',Title{i},'}$'],'FontSize',FontSizes(1),'interpreter','latex');
            end
            if(i == 1)||(i == 3)||(i == 5)
                ylabel('$E\,[t_1]$','FontSize',FontSizes(2),'interpreter','latex');
                text(kaRange(i,1)-1.56,Ylim(i,2)-0.05,['$\mathbf{',Config{i},'}$'],...
                'VerticalAlignment','middle','HorizontalAlignment','center','Color','k',...
                'FontSize',FontSizes(1),'Interpreter','latex');
            end
            if(i == 5)||(i == 6)
                xlabel('$ka$','FontSize',FontSizes(2),'interpreter','latex');
            end
            set(gca,'TickLength',[0.02,0.02]);
            set(gca,'LineWidth',1);
            yticks(Yticks(i,1:NYticks(i)));
            yticklabels(Yticklabels(i,1:NYticks(i)));
            text(sum(kaRange(i,:))/2-0.95,Ylim(i,2)-0.09*(Ylim(i,2)-Ylim(i,1)),1,CNumbers{i},...
            'BackgroundColor','white','Margin',1.3,'EdgeColor','none','Interpreter','latex','FontSize',14);
            Thickness = [0.2,3.7]; % Thickness of the bands for i = [4,5,6]
            if(i == 2)||(i == 1)
                Thickness = [0.1,3.6]; % Thickness of the bands for i = [1,2]
            elseif(i == 3)
                Thickness = [0.05,3.7]; % Thickness of the bands for i = 3
            end
        else % Inset plot
            axes('Position',InsetAxes(i,:));
            hold on; grid on; box on; orient(gcf,'landscape');
            set(gca,'FontSize',FontSizes(4));
            set(gca, 'XColor', 'r', 'YColor', 'r', 'ZColor', 'r');
            xlim(InsetRanges(i,1:2)); ylim(InsetRanges(i,3:4));
            set(gca,'TickLength',[0.015,0.015]);
            set(gca,'LineWidth',0.9);
            if(ismember(i,[1,2,5,6]))
                xticks([InsetRanges(i,1),0,InsetRanges(i,2)]);
            end
            if(i == 1)||(i == 2)
                yticks((-0.01:0.005:0.01)+0.045);
                yticklabels({'-0.01','-0.005','0','0.005','0.01'});
                Thickness = [1.4,2.5]; % Thickness of the bands for i = 2
                if(i == 1)
                    Thickness = [0.05,2.3]; % Thickness of the bands for i = 1
                    fid = fopen([Folder,Inset_Band_data{i},'.txt'],'r');
                    Band = textscan(fid,'%s','Delimiter','\n');
                    fclose(fid);
                    fid = fopen([Folder,Inset_Localization_data{i},'.txt'],'r');
                    Edge = textscan(fid,'%s','Delimiter','\n');
                    fclose(fid);
                    Nk = length(Band{1});
                end
            elseif(i == 3)||(i == 4)
                yticks((-0.02:0.01:0.02)+0.026);
                yticklabels({'-0.02','-0.01','0','0.01','0.02'});
                Thickness = [0.02,2.1]; % Thickness of the bands for i = 3
                if(i == 4)
                    Thickness = [0.3,2.3]; % Thickness of the bands for i = 4
                end
                fid = fopen([Folder,Inset_Band_data{i},'.txt'],'r');
                Band = textscan(fid,'%s','Delimiter','\n');
                fclose(fid);
                fid = fopen([Folder,Inset_Localization_data{i},'.txt'],'r');
                Edge = textscan(fid,'%s','Delimiter','\n');
                fclose(fid);
                Nk = length(Band{1});
            else
                yticks((-0.04:0.02:0.04)+0.02);
                yticklabels({'-0.04','-0.02','0','0.02','0.04'});
                Thickness = [0.3,2.8]; % Thickness of the bands
            end
        end
        for ika = 1:Nk
            Band_k  = sscanf(Band{1}{ika},'%f')';
            ka = Band_k(1);
            Ek = Band_k(2:end);
            NE = size(Ek,2);
            Edge_k  = sscanf(Edge{1}{ika},'%f')';
            Edgeness = zeros(1,NE);
            RGB = zeros(NE,3);
            RGB(:,3) = 1.0;
            for ie = 1:NE
                Edgeness(ie) = Edge_k(3*ie+ilayer-2);
                if(Edgeness(ie) <= 0.5)
                    RGB(ie,1:2) = [0,2*(0.5-Edgeness(ie))*1.0];
                else
            	    RGB(ie,1:2) = [2*(Edgeness(ie)-0.5)*1.0,0];
                end
            end
            scatter(ka*ones(1,NE),Ek,Thickness(1)+Thickness(2)*abs(Edgeness-0.5)/0.5,RGB,'filled');
        end
        if(j == 1)
            rectangle('Position',[kaRange(i,1),Ylim(i,1),kaRange(i,2)-kaRange(1),Ylim(i,2)-Ylim(i,1)],...
                      'LineWidth',1);
        else
            rectangle('Position',[InsetRanges(i,1),InsetRanges(i,3),InsetRanges(i,2)-InsetRanges(i,1),...
                                  InsetRanges(i,4)-InsetRanges(i,3)],'EdgeColor','r','LineWidth',0.9);
        end
    end
    clear Band Edge Band_k Edge_k;
end
print(gcf,'Output_figures/Fig_6.png','-dpng','-r500');
close(gcf);