function C_simVC_theoryFigAutogen(X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5,state)
%CREATEFIGURE(X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5)
%  X1:  vector of x data
%  Y1:  vector of y data
%  X2:  vector of x data
%  Y2:  vector of y data
%  X3:  vector of x data
%  Y3:  vector of y data
%  X4:  vector of x data
%  Y4:  vector of y data
%  X5:  vector of x data
%  Y5:  vector of y data

%  Auto-generated by MATLAB on 26-Mar-2020 10:11:04

% Create figure
figure('OuterPosition',[839 306 576 567]);

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create plot
%plot(X1,Y1,'Marker','.','LineStyle','none','Color',[0 0 0]);


if state == 'L'

    plot(X4,Y4,'Marker','.','LineStyle','none',...
        'Color',[0.270588235294118 0.254901960784314 0.254901960784314]);

    % Create plot
    plot(X2,Y2,'Marker','.','LineStyle','none','Color',[0.8 0.8 0.8]);
else
        % Create plot
    plot(X2,Y2,'Marker','.','LineStyle','none','Color',[0.8 0.8 0.8]);
    plot(X4,Y4,'Marker','.','LineStyle','none',...
        'Color',[0.270588235294118 0.254901960784314 0.254901960784314]);
end
    

% Create plot
plot(X3,Y3,'Marker','.','LineStyle','none',...
    'Color',[0.850980392156863 0.325490196078431 0.098039215686274]);




% Create plot
plot(X5,Y5,'LineWidth',2,...
    'Color',[0 0 0]);

% Create ylabel
ylabel('C_{theory}');

% Create xlabel
xlabel('C_{sim}');

% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[-0.05 0.08]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[-0.05 0.08]);
axis(axes1,'square');
% Set the remaining axes properties
set(axes1,'FontSize',20,'LineWidth',2);