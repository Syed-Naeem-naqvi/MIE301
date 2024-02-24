%% MIE301 Lab 1 Case 1

close all; % closes all figures
clear all; % clears all variables from memory
clc;       % clears all calculations from the Matlab workspace

% Plot Parameters: these will be used to set the axis limits on the figure
% using axis()
xmin= -15;  % leftmost window edge
xmax=15;   % rightmost window edge
ymin= -15; % bottom window edge
ymax= 15;  % top window edge
%% Define the link and motion Parameters

stepsize = 30 *pi/180;                          % number of configuration steps to calculate along mechanism rotation
max_rotation_theta2 = 320 *pi/180;              % rotation limit of theta2, radians
theta2 = 0 : stepsize : max_rotation_theta2;    % link 2 rotation into 'increments' number of angles
length2 = 8;                                    % link #2 length
length3 = 4;                                    % link #3 length
%% set up figure

figure(1);                         %create new figure
set(1,'WindowStyle','Docked')      %dock the figure
%% calculate mechanism motion and plot it

for i=1:length(theta2)                            % step through motion of the mechanism
    hold off;                                     % erase the previously plotted data each cycle
    Ax(i) = 0;                                    % pivot point of link 2 position
    Ay(i) = 0;                                    % pivot point of link 2 position
    Bx(i) = length2*cos( theta2(i) );             % point B position
    By(i) = length2*sin( theta2(i) );
    Cx(i) = Bx(i);                                % Point C position   
    Cy(i) = By(i) - length3;
    
    %hold on; % For Testing purposes, COMMENT BEFORE SUBMISSION 

    % point B position
    plot( [Ax(i) Bx(i)], [Ay(i), By(i)],'Color','r','LineWidth',3 ) %draw the link from A to B for the current configuration
    
    hold on;
    % Point C Position
    plot( [Bx(i) Cx(i)], [By(i), Cy(i)],'Color','r','LineWidth',3 ) %draw the link from B to C for the current configuration

    
    hold on;
    
    basePivotSize = 1;                                 %size of drawn base pivot
    plot([0,basePivotSize],[0,-basePivotSize],'r');     %draw base pivot
    plot([0,-basePivotSize],[0,-basePivotSize],'r');    %draw base pivot
    plot(0,0,'ro','MarkerFaceColor','w');               %draw a small circle at the base pivot point
    plot(Bx(i), By(i), 'bo','MarkerFaceColor','w');     %draw a small circle at B

    text(Ax(i)+0.6, Ay(i), 'A','color','r');            %label point A
    text(Bx(i)+0.4, By(i), 'B','color','r');            %label point B
    text(Cx(i)+1.4, Cy(i), 'C','color','r');            %label point C
    text(Bx(i)/2+0.5, By(i)/2+0.5, '2','color','r');    %Label link 2
    text(Cx(i)+0.5, Cy(i)/1.2, '3','color','r');        %Label link 3



    dx = 1;                                                    % Helper value
    rectangle('position', [Cx(i)-dx, Cy(i)-2*dx, 2*dx, 2*dx]); % Make rectangle
    text(Cx(i)-dx/2, Cy(i)-dx, 'm','color','r');               % label mass


    hold on;                            % draw on top of the current figure
    xlabel('x (cm)', 'fontsize', 15);   % axis label
    ylabel('y (cm)', 'fontsize', 15);   % axil label
    grid on;                            % add a grid to the figure
    axis equal;                         % make sure the figure is not stretched
    title('Lab1 - Problem 1');          % add a title to the figure
    axis( [xmin xmax ymin ymax] );      % set figure axis limits
    
%%
    pause(0.1);                                                                    % wait to proceed to next configuration, in seconds
    
end

plot(Bx(:),By(:),'--g','linewidth',1);                                             % draw a trace of the path of point B
plot(Cx(:),Cy(:),'--r','linewidth',1);                                             % draw a trace of the path of point C

%% 
%