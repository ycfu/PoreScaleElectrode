clc
clear all
close all

addpath(genpath('.\Functions'))
rng(1)

xlen = 0.2e-3;  %battery depth [m];
ylen = 0.1e-3;    %battery width [m];
zlen = 0.2e-3;   % battery height [m];
r = 3.55e-6;         % electrode radius = diameter [m]
theta1_range = 10;    % start angle range;
theta2_range = 30;    % arc angle range

ap = 200000;   % fiber specific area
% density = 0.5E12; % hole density  [number per m^2]

Boundary_Restriction = 0;  % control if the rod can penetrate the battery domain 0-> allow penetration;  1-> no penetration
Balls_Flag = 0;  % Defined Fiber Type: 0: traditional fiber;    1: Balls on fiber

%%
count = 0;

for k = 1:10
    disp(['Generating ', num2str(k),'th electrode'])
    w = ylen*0.95; % m
    theta1 = theta1_range *rand(); % degree
    theta2 = theta2_range *rand();
    
    xt = 0 ;   % xt
    zt = 0;    % zt
    yt = 0;
    rot = 0;
    
    R = w/(sind(theta1 + theta2) - sind(theta1));
    
    [Org, Trans] = tubegeneration(R, r ,theta1, theta2, rot, xt, yt, zt);
  
    
    if Boundary_Restriction == 1
        while checkboundary(Trans, xlen,zlen, r) ==0  & theta2>0
            theta2 = theta2-0.1;
            [Org, Trans] = tubegeneration(R, r ,theta1, theta2, rot, xt , yt, zt);
        end
    end
      count = count +1;
    
%     if theta2<0
%         continue;
%     else
%         count = count +1;
%         if count >10
%             break;
%         end
%     end
%     
    [Trans.Xm,Trans.Ym,Trans.Zm] = tubeplot([Trans.xvector; Trans.yvector; Trans.zvector],r,8);
    Electrode(count).Trans = Trans;
    Electrode(count).Org = Org;
    Electrode(count).theta1 = theta1;
    Electrode(count).theta2 = theta2;
    Electrode(count).R = R;
    Electrode(count).r = r;
    Electrode(count).rot = rot;
    Electrode(count).w = w;
    Electrode(count).xt = xt;
    Electrode(count).zt = zt;
    Electrode(count).area = 4*pi^2*  Electrode(count).R*  Electrode(count).r * (  Electrode(count).theta2/360);
    
%     r_bc = 1.75e-6;
%     r_bo = r;
%     ball_num = round(Electrode(count).area*0.49*0.85/(pi*r_bc^2));
%     
%     templength = 2*pi*Electrode(1).R*Electrode(1).theta2/360*0.85;
%     tempwidth = 2*pi*Electrode(1).r;
%     
%     n_width = round(sqrt(ball_num/(templength/tempwidth)));
%     n_length = round(n_width * (templength/tempwidth));
%     
%     circles =[];
%     for k = 1:ball_num
%         %         k
%         [circles flag] = add_circles(circles, r_bc, tempwidth, templength,2.2);
%         [k,flag]
%     end
%     
%     
%     theta_b1 = theta1+circles(:,2)/templength*theta2*0.85 +theta2*0.1;
%     theta_b2 = circles(:,1)/tempwidth*360;
%     %     theta_b2 = ones(ball_num,1)*90;
%     x_bo = Org.xc  + R*cosd(theta_b1);
%     y_bo = Org.yc  + R*sind(theta_b1);
%     x_bc = x_bo + r_bo*cosd(theta_b2).*cosd(theta_b1);
%     y_bc = y_bo + r_bo*cosd(theta_b2).*sind(theta_b1);
%     z_bc =        r_bo*sind(theta_b2);
%     
%     Balls(count).r_bc = r_bc;   %um
%     Balls(count).r_bo = r_bo;
%     Balls(count).theta_b1 =theta_b1;
%     Balls(count).theta_b2 = theta_b2;
%     Balls(count).x_bo = x_bo;
%     Balls(count).y_bo = y_bo;
%     
%     Balls(count).x_bc = x_bc;
%     Balls(count).y_bc = y_bc;
%     Balls(count).z_bc = z_bc;
%     
%     Balls(count).x_bct = xt + x_bc*cosd(rot) + z_bc*sind(rot);
%     Balls(count).y_bct = yt + y_bc;
%     Balls(count).z_bct = zt - x_bc*sind(rot) + z_bc*cosd(rot);
end

%% Shift Rod for Distribution
area = 0;
areatarget = xlen*ylen*zlen*ap;
rod_num = round(areatarget/Electrode(1).area);
rod_loc =[];

for k = 1:rod_num
    %     if k == 64
    %         k
    %     end
    
    [rod_loc flag] = add_circles(rod_loc, Electrode(1).r, xlen-2*r, zlen-2*r, 0.5);
    
    xt = rod_loc(k,1)+Electrode(1).r;
    zt = rod_loc(k,2)+Electrode(1).r;
    yt = 0;
    rot = 360*rand();
    [Org, Trans] = tubegeneration(Electrode(1).R, Electrode(1).r ,Electrode(1).theta1, Electrode(1).theta2, rot, xt, yt, zt);
    
    flag_overlapping = 1;
    %     while checkboundary(Trans, xlen,zlen, r) ==1  |  flag_overlapping == 1   % either outside boundary (1) or  tube overlapping(1) exists
    while checkboundary(Trans, xlen,zlen, r) ==1
        rod_loc(end,:) = [];
        [rod_loc flag] = add_circles(rod_loc, Electrode(1).r, xlen-2*r, zlen-2*r, 3);
        xt = rod_loc(k,1)+Electrode(1).r;
        zt = rod_loc(k,2)+Electrode(1).r;
        yt = 0;
        rot = 360*rand();
        [Org, Trans] = tubegeneration(R, r ,theta1, theta2, rot, xt , yt, zt);
        %         if k == 1
        %             flag_overlapping = 0;
        %         else
        %             flag_overlapping = checkoverlapping(Trans, ElectrodeB);
        %         end
    end
    [Trans.Xm,Trans.Ym,Trans.Zm] = tubeplot([Trans.xvector; Trans.yvector; Trans.zvector],r,8);
    ElectrodeB(k).Trans = Trans;
    ElectrodeB(k).Org = Org;
    ElectrodeB(k).theta1 = theta1;
    ElectrodeB(k).theta2 = theta2;
    ElectrodeB(k).R = R;
    ElectrodeB(k).r = r;
    ElectrodeB(k).rot = rot;
    ElectrodeB(k).w = w;
    ElectrodeB(k).xt = xt;
    ElectrodeB(k).zt = zt;
    ElectrodeB(k).area = 4*pi^2*  ElectrodeB(k).R*  ElectrodeB(k).r * (  ElectrodeB(k).theta2/360);
    
end


%%
figure,hold on,
for k = 1:length(Electrode)
    tubeplot([ Electrode(k).Trans.xvector;  Electrode(k).Trans.yvector;  Electrode(k).Trans.zvector], Electrode(k).r,8);
    plot3(Electrode(k).Trans.xvector,  Electrode(k).Trans.yvector,  Electrode(k).Trans.zvector,'r.')
end

% for kk = 1
%     for k = 1:length(Balls(k).x_bct)
%         k
%         %     plot3(Balls(kk).x_bct(k),Balls(kk).y_bct(k),Balls(kk).z_bct(k),'ko')
%         %     plot3(Balls(kk).x_bc(k),Balls(kk).y_bc(k),Balls(kk).z_bc(k),'bs')
%         %     plot3(Balls(kk).x_bo(k),Balls(kk).y_bo(k),Balls(kk).y_bo(k)*0,'gs')
%         plot_sphere( [Balls(kk).x_bct(k) Balls(kk).y_bct(k) Balls(kk).z_bct(k)],r_bc, [220 220 220]/255,1)
%     end
% end
axis equal
daspect([1,1,1]);
camlight;
view([90 90])
xlabel('x'),ylabel('y'),zlabel('z')


%%
figure,hold on,
for k = 1:length(ElectrodeB)
    tubeplot([ ElectrodeB(k).Trans.xvector;  ElectrodeB(k).Trans.yvector;  ElectrodeB(k).Trans.zvector], ElectrodeB(k).r,8);
    plot3(ElectrodeB(k).Trans.xvector,  ElectrodeB(k).Trans.yvector,  ElectrodeB(k).Trans.zvector,'r.')
    text(ElectrodeB(k).Trans.xvector(1),  ElectrodeB(k).Trans.yvector(1),  ElectrodeB(k).Trans.zvector(1),num2str(k))
end

axis equal
daspect([1,1,1]);
camlight;
view([90 90])
xlabel('x'),ylabel('y'),zlabel('z')
checkoverlapping( ElectrodeB(64).Trans, ElectrodeB)

%% Output Each single rod stl fi
% X =[]; Y =[]; Z = [];
% for k = 1:length(Electrode)
%     X = [Electrode(k).Trans.Xm];
%     Y = [Electrode(k).Trans.Ym];
%     Z = [Electrode(k).Trans.Zm];
%     surf2stl(['.\Single_Rods\STL_Rod_',num2str(k),'.stl'],X,Y,Z, 'ascii' );
% end
%% Save file for rod location
output = []; output1 = [];
for k = 1:length(Electrode)
    output1 = [   Electrode(k).w  Electrode(k).r Electrode(k).R Electrode(k).theta1 Electrode(k).theta2 ...
        0 0 0 ...   % translation xt zt and rot Electrode(k).xt Electrode(k).zt Electrode(k).rot
        Electrode(k).Org.x1 Electrode(k).Org.x2 Electrode(k).Org.xc ...
        Electrode(k).Org.y1 Electrode(k).Org.y2 Electrode(k).Org.yc];
    output(k,:) = output1;
end

output_scaled = output;
output_scaled(:,[1 2 3 6 7 9:14]) =  output_scaled(:,[1 2 3 6 7 9:14]) *1E4;

outputtable = array2table(output_scaled, 'VariableNames',{'w','r','R','theta1','theta2','xt','zt','rot', 'x1', 'x2', 'xc', 'y1', 'y2', 'yc'});
delete RodDistribution.csv
writetable(outputtable,'RodDistribution.csv')

%%
% output = []; output1 = [];
% for k = 1:length(ElectrodeB)
%     output1 = [ElectrodeB(k).xt ElectrodeB(k).zt ElectrodeB(k).rot];
%     output(k,:) = output1;
% end

VarName ={'sx','sy','sz',...
    'rdeg1','tx','ty','tz'...
    'rdeg2','rx','ry','rz',...
    'shift_a','shift_b','shift_c',...
    };
for k = 1:length(ElectrodeB)
    output = [ElectrodeB(k).xt ElectrodeB(k).zt ElectrodeB(k).rot];
    
    rn = rand();
    
    if rn <0.5
        rodpara(k,:) =[1    1      1   output(3) 0 1 0      0   0 1 0  output(1)*1E4     0   output(2)*1E4 ];
    elseif rn<0.75
        rodpara(k,:) =[1 xlen/ylen 1   output(3) 0 1 0     -90 0 0 1  xlen*0.02*1E4 output(1)*ylen/xlen*1E4 output(2)*1E4];
    elseif rn <1
        rodpara(k,:) =[1 zlen/ylen 1   output(3) 0 1 0     90 1 0 0  output(1)*1E4 output(2)*ylen/zlen*1E4  zlen*0.02*1E4];
    end
end
outputtable = array2table(rodpara, 'VariableNames',VarName);
delete Rod_Translation.csv
writetable(outputtable,'Rod_Translation.csv')


%%
% for k = 1:100
%     test(k) = int8(randsrc(1,1,[1,2,3;0.5,0.25,0.25]));
% end
% figure,histogram(test)
%%

% for ball l
output = []; output1 = [];
for k = 1
    output1 = [ Balls(k).x_bct Balls(k).y_bct Balls(k).z_bct  Balls(k).r_bc*ones(length(Balls(k).x_bct),1)];
    output = [output; output1];
end

output_scaled = output*1E4;
outputtable = array2table(output_scaled, 'VariableNames',{'xbct','ybct','zbct','rbc'});
delete Balls.csv
writetable(outputtable,'Balls.csv')

%%
function [Org, Trans] = tubegeneration(R, r ,theta1, theta2, rot, xt ,yt, zt)

Org.x1 = 0;
Org.y1 = 0;
Org.z1 = 0;
Org.x2 = -R*cosd(theta1) + R*cosd(theta1+theta2);
Org.y2 = -R*sind(theta1) + R*sind(theta1+theta2);
Org.z2 = 0;
Org.xc = -R*cosd(theta1);
Org.yc = -R*sind(theta1);

temptheta = linspace(theta1,theta1+theta2,100);
Org.xvector =  -R*cosd(theta1) + R*cosd(temptheta);
Org.yvector =  -R*sind(theta1) + R*sind(temptheta);
Org.zvector =  0 * ones(size(temptheta));

Trans.x1t = xt;    % translation vector
Trans.y1t = yt;
Trans.z1t = zt;

Trans.x2t = xt + Org.x2*cosd(rot) + Org.z2*sind(rot);  % translated tube end location
Trans.y2t = yt + Org.y2;
Trans.z2t = zt - Org.x2 *sind(rot)+ Org.z2*cosd(rot);

Trans.xvector = xt + Org.xvector*cosd(rot);
Trans.yvector = yt + Org.yvector;
Trans.zvector = zt - Org.xvector*sind(rot);
end


%% Checking Functions

function flag = checkboundary(Trans, xlen,zlen, r)
if (Trans.xvector+3*r)<xlen &  (Trans.xvector-3*r)> 0 &  (Trans.zvector+3*r)<zlen & (Trans.zvector-3*r) >0
    flag = 0;
else
    flag = 1;
end
end

function flag_overlapping = checkoverlapping(Trans, Electrode)

flag_overlapping =0; % start with no overlapping
D1 = [Trans.xvector' Trans.yvector' Trans.zvector'];

for k = 1:length(Electrode)
    D2 = [Electrode(k).Trans.xvector' Electrode(k).Trans.yvector' Electrode(k).Trans.zvector'];
    D = pdist2(D1,D2);
    D = diag(D);
    if min(D)<Electrode(k).r*3
        flag_overlapping =1;
    end
    
end


end


%% Plot Sphere
function plot_sphere(org,radius, facecolor, alpha)

[X,Y,Z] = sphere;
X2 = X * radius;
Y2 = Y * radius;
Z2 = Z * radius;

h = surf(X2+org(1),Y2+org(2),Z2+org(3));
set(h,'FaceColor',facecolor,'Facealpha',alpha,'EdgeColor','none','AmbientStrength',.5);
hold on
hold on
end

%%
function [circles flag] = add_circles(circles, r, Width, Length, seperationcof)


newCircleFound = false;
c =0;

while c<100000 & newCircleFound ==false
    
    x = r + (Width)*rand(1);
    y = r + (Length)*rand(1);
    
    
    if isempty(circles)
        newCircleFound = true;
        circles(1,:) = [x y];
    else
        %calculates distances from previous drawn circles
        prevCirclesY = [circles(:,2); circles(:,2); circles(:,2)];
        prevCirclesX = [circles(:,1); circles(:,1)+Width; circles(:,1)-Width];
        
        
        distFromPrevCircles = ((prevCirclesX-x).^2+(prevCirclesY-y).^2).^0.5;
        
        %if the distance is not to small - adds the new circle to the list
        if sum(distFromPrevCircles<=seperationcof*r)==0
            newCircleFound = true;
            circles(end+1,:) = [x y];
        end
    end
    c = c+1;
end
flag = newCircleFound;

end

