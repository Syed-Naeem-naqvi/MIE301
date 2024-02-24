close all; % closes all figures
clear all; % clears all variables from memory
clc;       % clears all calculations from the Matlab workspace

% Plot Parameters: these will be used to set the axis limits on the figure
%xmin= -10; % leftmost window edge
%xmax= 140;  % rightmost window edge
%ymin= -50; % bottom window edge
%ymax= 20;  % top window edge

% Link Parameters
increments = 50; % number of theta2 configuration steps to calculate along mechanism rotation %%%% YOU MAY WANT TO CHANGE THIS

max_rotation_theta2 = 360 *pi/180; % rotation limit of theta2, radians
theta2 = linspace(0,max_rotation_theta2,increments); % link 2 rotation into 'increments' number of angles

%minimum link length must be able to accomodate foot pedals
inc = 5;     %value to increase link length between each checked value
r1_increments = 115;   % link 1 length, cm
%r2 = [10 12 14 16 18 20 22 24 26];           % link 2 length, cm
r2_increments = 1:inc:50;  
%r3 = [106 108 110  112 114 116 118 120 122]
r3_increments = 107;            % link 3 length, cm
r4_increments = 5:inc:100;  
%r4 = 34;           % link 4 length, cm
offset_increments = 30:1:30;                    %height of offset (cm)
%footDist = 0.1:0.05:0.9

%variables to save stride length measurements
min_stride_length = 110; %should be set to largest value permitted for stride length
max_stride_length = 0;
%min_ratio = 1000;
%max_ratio = 0;

for v=1:length(r1_increments)
    r1=r1_increments(v);
    disp(r1)
    %disp(length(r1_increments))
    for o=1:length(offset_increments)
        offset = offset_increments(o);
        angOffset = asin(offset/r1);    %angle of offset
        R = [cos(angOffset) -sin(angOffset); sin(angOffset) cos(angOffset)];  %rotation matrix
        for w=1:length(r2_increments)
            for x=1:length(r3_increments)
                for y=1:length(r4_increments)
                    
                    r2=r2_increments(w);
                    r3=r3_increments(x);
                    r4=r4_increments(y); 
    
                    %parameters to make sure the link is a crank-rocker:
                    %to be constructed: r1<r2+r3+r4
                    if(r2+r3+r3>r1)
                        %disp('test')
                        list_of_lengths = [r1,r2,r3,r4];
                        [s,i]=min(list_of_lengths);
                        list_of_lengths(i) = [];
                        [l,j]=max(list_of_lengths);
                        list_of_lengths(j) = [];
                        p=list_of_lengths(1);
                        q=list_of_lengths(2);
    
                        if (((p+q)>(s+l))&&(s==r2))
                    %p+q>s+l and s must be input
                            footDist = 0.7*r3;     % distance from B to D where the foot is placed on link 3
                            h1 = r1/r2;
                            h2 = r1/r3;
                            h3 = r1/r4;
                            h4 =(-r1^2-r2^2-r3^2+r4^2)/(2*r2*r3);
                            h5 =(r1^2+r2^2-r3^2+r4^2)/(2*r2*r4);
                            %disp("Ran!")
                       
                            for i=1:increments 
                            %  geometric calculations (book eq. 4.3-56 to 4.3-62):
                                d = -h1 +(1-h3)*cos(theta2(i)) +h5;
                                b = -2*sin(theta2(i));
                                e = h1 -(1+h3)*cos(theta2(i)) +h5;
                                a_a = -h1 +(1+h2)*cos(theta2(i)) +h4;
                                c = h1 -(1-h2)*cos(theta2(i)) +h4;
                        
                                theta3_1(i) = 2*atan(((-b+(b^2-4*a_a*c)^0.5)/(2*a_a)));  %calculate angle of link 3 (eq. 4.3-64)
                                theta4_1(i) = 2*atan(((-b+(b^2-4*d*e)^0.5)/(2*d)));      %angle of theta 4  
                        
                                % Link Coordinates calculations:   
                                Ax(i) = 0;                                      % pivot point of link 2 position
                                Ay(i) = 0;                                      % pivot point of link 2 position
                                Bx(i) = r2*cos( theta2(i) );                    % point B position
                                By(i) = r2*sin( theta2(i) );                    % point B position
                                Cx(i) = Bx(i) + (r3 )*cos( theta3_1(i) );     % point C position. We're going to store these for each extension value (indexed by j)
                                Cy(i) = By(i) + (r3 )*sin( theta3_1(i) );     % point C position
                                Dx(i) = r1 + r4*cos( theta4_1(i) );             % point D position
                                Dy(i) = r4*sin( theta4_1(i) );                  % point D position   
                        
                                Fx(i) = (Bx(i)+Dx(i))*(footDist/r3);          % point F position
                                Fy(i) = (By(i)+Dy(i))*(footDist/r3);          % point F posistion  
      
                            end  
                            [Ax,Ay]=rotate(Ax,Ay,offset,r1);
                            [Bx,By]=rotate(Bx,By,offset,r1);
                            [Cx,Cy]=rotate(Cx,Cy,offset,r1);
                            [Dx,Dy]=rotate(Dx,Dy,offset,r1);
                            [Fx,Fy]=rotate(Fx,Fy,offset,r1);
         
                            %calculate slope, make sure selected links have
                            %apropriate slope
                            if (isreal(Ax)&&isreal(Ay)&&isreal(Bx)&&isreal(By)&&isreal(Cx)&&isreal(Cy)&&isreal(Dx)&&isreal(Dy))
                                stride_length = abs(max(Fx) - min(Fx));
                                stride_height = abs(max(Fy) - min(Fy));
                            
                                if ((stride_length>max_stride_length))
                                    max_stride_length = stride_length;
                                    max_stride_index = [v,w,x,y,o];
                                end
                                if ((abs(stride_length) <abs(min_stride_length))&&(stride_height < 18))
                                    min_stride_length = stride_length;
                                    min_stride_index = [v,w,x,y,o];
                                
                                %select link lengths for min stride length 
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%plot longest and shortest stride length
disp(max_stride_length)
disp(min_stride_length)
% set up figure
clear vars Ax Ay Bx By Cx Cy Dx Dy Fx Fy
figure(1);                         %create new figure
set(1,'WindowStyle','Docked')      %dock the figure
for i=1:increments 
    hold on;
                    r1=r1_increments(max_stride_index(1));
                    r2=r2_increments(max_stride_index(2));
                    r3=r3_increments(max_stride_index(3));
                    r4=r4_increments(max_stride_index(4)); 
                    offset = offset_increments(max_stride_index(5));
                    
                    h1 = r1/r2;
                    h2 = r1/r3;
                    h3 = r1/r4;
                    h4 =(-r1^2-r2^2-r3^2+r4^2)/(2*r2*r3);
                    h5 =(r1^2+r2^2-r3^2+r4^2)/(2*r2*r4);
                    %  geometric calculations (book eq. 4.3-56 to 4.3-62):
                        d = -h1 +(1-h3)*cos(theta2(i)) +h5;
                        b = -2*sin(theta2(i));
                        e = h1 -(1+h3)*cos(theta2(i)) +h5;
                        a_a = -h1 +(1+h2)*cos(theta2(i)) +h4;
                        c = h1 -(1-h2)*cos(theta2(i)) +h4;
                
                        theta3_1(i) = 2*atan(((-b+(b^2-4*a_a*c)^0.5)/(2*a_a)));  %calculate angle of link 3 (eq. 4.3-64)
                        theta4_1(i) = 2*atan(((-b+(b^2-4*d*e)^0.5)/(2*d)));      %angle of theta 4  
                
                        % Link Coordinates calculations:   
                        Ax(i) = 0;                                      % pivot point of link 2 position
                        Ay(i) = 0;                                      % pivot point of link 2 position
                        Bx(i) = r2*cos( theta2(i) );                    % point B position
                        By(i) = r2*sin( theta2(i) );                    % point B position
                        Cx(i) = Bx(i) + (r3 )*cos( theta3_1(i) );     % point C position. We're going to store these for each extension value (indexed by j)
                        Cy(i) = By(i) + (r3 )*sin( theta3_1(i) );     % point C position
                        Dx(i) = r1 + r4*cos( theta4_1(i) );             % point D position
                        Dy(i) = r4*sin( theta4_1(i) );                  % point D position   
                
                        Fx(i) = (Bx(i)+Dx(i))*(footDist/r3);          % point F position
                        Fy(i) = (By(i)+Dy(i))*(footDist/r3);          % point F posistion
                           
                        [Ax(i),Ay(i)]=rotate(Ax(i),Ay(i),offset,r1);
                        [Bx(i),By(i)]=rotate(Bx(i),By(i),offset,r1);
                        [Cx(i),Cy(i)]=rotate(Cx(i),Cy(i),offset,r1);
                        [Dx(i),Dy(i)]=rotate(Dx(i),Dy(i),offset,r1);
                        [Fx(i),Fy(i)]=rotate(Fx(i),Fy(i),offset,r1);

                        [Px,Py]=rotate(r1,0,offset,r1);
    if (true) % put true here to draw the mechanism each time, false to skip drawing that for increased simulation speed
                            plot( [Ax(i) Bx(i)], [Ay(i), By(i)],'Color','r','LineWidth',3 ); % draw link2
                            %hold on;            
                            plot( [Bx(i), Dx(i)], [ By(i), Dy(i)],'Color','b','LineWidth',3 ); % draw link3
                            plot( [Dx(i), Cx(i)], [ Dy(i), Cy(i)],'Color','b','LineWidth',3 ); % draw link3 extension to point C
                            plot( [Px, Dx(i)], [ Py, Dy(i)],'Color','m','LineWidth',3 ); % draw link4
                            
                            % Draw Base Pivots:
                            recsz = 2.5;                                        % size of drawn base pivot
                            plot([0,recsz],[0,-recsz],'r');                     % draw base pivot for link2
                            plot([0,-recsz],[0,-recsz],'r');                    % draw base pivot for link2
                            plot(0,0,'ro','MarkerFaceColor','w');               % draw a small circle at the base pivot point
                            plot(Bx(i), By(i), 'bo','MarkerFaceColor','w');     % draw a small circle at B
                            text(Bx(i)+0.9, By(i), 'B','color','b');            % label point B
                            %plot(Fx(i), Fy(i), 'bo','MarkerFaceColor','w');    % draw a small circle at F
                            text(Fx(i)+0.9, Fy(i)+0.9, 'F','color','b');        % label point F
                
                            plot([Px,Px+recsz],[Py,-recsz],'r');                 % draw base pivot for link4
                            plot([Px,Px-recsz],[Py,-recsz],'r');                 % draw base pivot for link4
                            plot(Px,Py,'ro','MarkerFaceColor','w');              % draw a small circle at the base pivot point
                            plot(Dx(i), Dy(i), 'bo','MarkerFaceColor','w');     % draw a small circle at D
                            text(Dx(i)+0.9, Dy(i), 'D','color','b');            % label point D
                            text(Cx(i)+0.9, Cy(i), 'C','color','b');        % label point C
                
                            plot(Fx(i), Fy(i), 'bo','MarkerFaceColor','w');     % draw a small circle at D                 % draw base pivot for link4

                            legend({['Length 1: ' , num2str(r1)],['Length 2: ' , num2str(r2)], ['Length 3: ', num2str(r3)], ['Length 4: ', num2str(r4)],['Offset: ', num2str(offset)]})
                            
                            xlabel('x (cm)', 'fontsize', 15);   % axis label
                            ylabel('y (cm)', 'fontsize', 15);   % axil label 
                            title('Elliptical Footpath');       % add a title to the figure
                            axis equal;                         % make sure the figure is not stretched
                            grid on;
                            %axis( [xmin xmax ymin ymax] );      % figure axis limits            

                         %pause(0.01);                            % wait to proceed to next configuration, seconds
    end
    %plot(Fx(:),Fy(:));
    hold off;
end

% set up figure
figure(2);                         %create new figure
%set(2,'WindowStyle','Docked')      %dock the figure
hold on;
plot(Fx, Fy);
title("Path Shape for Maximum Stride Length")
hold off;

clear vars Ax Ay Bx By Cx Cy Dx Dy Fx Fy
figure(3);                         %create new figure
%set(3,'WindowStyle','Docked')      %dock the figure
hold on;
title("Path Shape for Maximum Stride Length")
hold off;

for i=1:increments 
    hold on;
                    r1=r1_increments(min_stride_index(1));
                    r2=r2_increments(min_stride_index(2));
                    r3=r3_increments(min_stride_index(3));
                    r4=r4_increments(min_stride_index(4)); 
                    offset = offset_increments(min_stride_index(5));
                    
                    h1 = r1/r2;
                    h2 = r1/r3;
                    h3 = r1/r4;
                    h4 =(-r1^2-r2^2-r3^2+r4^2)/(2*r2*r3);
                    h5 =(r1^2+r2^2-r3^2+r4^2)/(2*r2*r4);
                    %  geometric calculations (book eq. 4.3-56 to 4.3-62):
                        d = -h1 +(1-h3)*cos(theta2(i)) +h5;
                        b = -2*sin(theta2(i));
                        e = h1 -(1+h3)*cos(theta2(i)) +h5;
                        a_a = -h1 +(1+h2)*cos(theta2(i)) +h4;
                        c = h1 -(1-h2)*cos(theta2(i)) +h4;
                
                        theta3_1(i) = 2*atan(((-b+(b^2-4*a_a*c)^0.5)/(2*a_a)));  %calculate angle of link 3 (eq. 4.3-64)
                        theta4_1(i) = 2*atan(((-b+(b^2-4*d*e)^0.5)/(2*d)));      %angle of theta 4  
                
                        % Link Coordinates calculations:   
                        Ax(i) = 0;                                      % pivot point of link 2 position
                        Ay(i) = 0;                                      % pivot point of link 2 position
                        Bx(i) = r2*cos( theta2(i) );                    % point B position
                        By(i) = r2*sin( theta2(i) );                    % point B position
                        Cx(i) = Bx(i) + (r3 )*cos( theta3_1(i) );     % point C position. We're going to store these for each extension value (indexed by j)
                        Cy(i) = By(i) + (r3 )*sin( theta3_1(i) );     % point C position
                        Dx(i) = r1 + r4*cos( theta4_1(i) );             % point D position
                        Dy(i) = r4*sin( theta4_1(i) );                  % point D position   
                
                        Fx(i) = (Bx(i)+Dx(i))*(footDist/r3);          % point F position
                        Fy(i) = (By(i)+Dy(i))*(footDist/r3);          % point F posistion
                           
                        [Ax(i),Ay(i)]=rotate(Ax(i),Ay(i),offset,r1);
                        [Bx(i),By(i)]=rotate(Bx(i),By(i),offset,r1);
                        [Cx(i),Cy(i)]=rotate(Cx(i),Cy(i),offset,r1);
                        [Dx(i),Dy(i)]=rotate(Dx(i),Dy(i),offset,r1);
                        [Fx(i),Fy(i)]=rotate(Fx(i),Fy(i),offset,r1);

                        [Px,Py]=rotate(r1,0,offset,r1);
    if (true) % put true here to draw the mechanism each time, false to skip drawing that for increased simulation speed
                            plot( [Ax(i) Bx(i)], [Ay(i), By(i)],'Color','r','LineWidth',3 ); % draw link2
                            %hold on;            
                            plot( [Bx(i), Dx(i)], [ By(i), Dy(i)],'Color','b','LineWidth',3 ); % draw link3
                            plot( [Dx(i), Cx(i)], [ Dy(i), Cy(i)],'Color','b','LineWidth',3 ); % draw link3 extension to point C
                            plot( [Px, Dx(i)], [ Py, Dy(i)],'Color','m','LineWidth',3 ); % draw link4
                            
                            % Draw Base Pivots:
                            recsz = 2.5;                                        % size of drawn base pivot
                            plot([0,recsz],[0,-recsz],'r');                     % draw base pivot for link2
                            plot([0,-recsz],[0,-recsz],'r');                    % draw base pivot for link2
                            plot(0,0,'ro','MarkerFaceColor','w');               % draw a small circle at the base pivot point
                            plot(Bx(i), By(i), 'bo','MarkerFaceColor','w');     % draw a small circle at B
                            text(Bx(i)+0.9, By(i), 'B','color','b');            % label point B
                            %plot(Fx(i), Fy(i), 'bo','MarkerFaceColor','w');    % draw a small circle at F
                            text(Fx(i)+0.9, Fy(i)+0.9, 'F','color','b');        % label point F
                
                            plot([Px,Px+recsz],[Py,-recsz],'r');                 % draw base pivot for link4
                            plot([Px,Px-recsz],[Py,-recsz],'r');                 % draw base pivot for link4
                            plot(Px,Py,'ro','MarkerFaceColor','w');              % draw a small circle at the base pivot point
                            plot(Dx(i), Dy(i), 'bo','MarkerFaceColor','w');     % draw a small circle at D
                            text(Dx(i)+0.9, Dy(i), 'D','color','b');            % label point D
                            text(Cx(i)+0.9, Cy(i), 'C','color','b');        % label point C
                
                            plot(Fx(i), Fy(i), 'bo','MarkerFaceColor','w');     % draw a small circle at D                 % draw base pivot for link4

                            legend({['Length 1: ' , num2str(r1)],['Length 2: ' , num2str(r2)], ['Length 3: ', num2str(r3)], ['Length 4: ', num2str(r4)],['Offset: ', num2str(offset)]})
                            
                            xlabel('x (cm)', 'fontsize', 15);   % axis label
                            ylabel('y (cm)', 'fontsize', 15);   % axil label 
                            title('Elliptical Footpath');       % add a title to the figure
                            axis equal;                         % make sure the figure is not stretched
                            grid on;
                            %axis( [xmin xmax ymin ymax] );      % figure axis limits            

                         %pause(0.01);                            % wait to proceed to next configuration, seconds
    end
    %plot(Fx(:),Fy(:));
    hold off;
end

%create a function to rotate the calculated coordinates
% set up figure
figure(4);                         %create new figure
%set(4,'WindowStyle','Docked')      %dock the figure
hold on;
plot(Fx, Fy);
title("Path Shape for Maximum Stride Length")
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity and Acceleration Analysis: Naeem
% To use this code, make Fx and Fy the x and y coordinates of the path
% shape
% NOTE: the set() function causes an error on my computer (no idea why), so I commented
% it out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define time stuff
theta2_dot = 30;                        % rotation rate, rpm
t_rev = 60/theta2_dot;                  % rotation time limit, seconds
time = linspace(0,t_rev,increments);    % link 2 rotation into 'increments' number of times. Increments is defined above
dt = time(2)                            % time values are linearly spaced, so pick the second one to obtain that spacing (first one's 0s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ok now find velocity using the optimized path (Fx, Fy)
velocity_vector = [];                                  % Store Velocities here
for i=2:length(Fx)-1
    velocity_value = abs(Fx(i+1) - Fx(i-1))/(2*dt);    % Aprroximate velocity
    velocity_vector = [velocity_vector; velocity_value];
end
%disp(velocity_vector)

% Now, acceleration

acceleration_vector = [];                                  % Store accelerations here
for i=2:length(velocity_vector)-1
    acceleration_value = abs(velocity_vector(i+1) - velocity_vector(i-1))/(2*dt);    % Aprroximate velocity
    acceleration_vector = [acceleration_vector; acceleration_value];
end
%disp(acceleration_vector)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graph the results as functions of time

velocity_times = time(2:end-1);               % For an Fx vector of size n, the velocity vector is of size (n-2)
acceleration_times = velocity_times(2:end-1); % For an Fx vector of size n, the acceleration vector is of size (n-4)
figure(5)
hold on;
plot(velocity_times, velocity_vector, 'r')
plot(acceleration_times, acceleration_vector, 'b')
xlim([0 t_rev])
xlabel('Time (sec)')
ylabel('Velocity (cm/s), Acceleration (cm/s^2)')
title('Velocity vs time and Acceleratio vs time of footpath')
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Find the min and max velocities and accelerations

maximum_velocity = max(velocity_vector)
minimum_velocity = min(velocity_vector)
maximum_acceleration = max(acceleration_vector)
minimum_acceleration = min(acceleration_vector)




%function to rotate link 4
%%
function [X,Y] = rotate(x,y,offset,r1)
    for i=1:length(x)          %iterate through length of the points that make up the footpath
        angOffset = asin(offset/r1);    %angle of offset
        R = [cos(angOffset) -sin(angOffset); sin(angOffset) cos(angOffset)];  %rotation matrix
        v = [x(i);y(i)];
        n = R*v;                %new  Fx and Fy are equal to the original Fx and Fy multimplied by the rotation matrix
        X(i)= n(1,:);
        Y(i) = n(2,:);
    end
end