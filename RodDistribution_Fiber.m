clc
clear all
close all

addpath(genpath('.\Functions'))
rng(1)

xlen = 0.135e-3;  %battery depth;
ylen = 0.3e-3;    %battery width;
zlen = 0.54e-3;   % battery height;
r = 4e-6;         % m, electrode diameter
theta1_range = 30;    % start angle range;
theta2_range = 30;    % arc angle range

Boundary_Restriction = 1;  % control if the rod can penetrate the battery domain

%% Generate One Typical Fiber
count = 0;
for k = 1:100
    disp(['Generating ', num2str(k),'th electrode'])
    w = ylen*0.95; % m
    theta1 = theta1_range *rand(); % degree
    theta2 = theta2_range *rand();
    
    xt = r + (xlen-2*r)*rand() ;   % xt
    zt = r + (zlen-2*r)*rand();    % zt
    yt = 0;
    
    rot = 360*rand();
    R = w/(sind(theta1 + theta2) - sind(theta1));
    [Org, Trans] = tubegeneration(R, r ,theta1, theta2, rot, xt, yt, zt);
    
    while checkboundary(Trans, xlen,zlen, r) ==0  & theta2>0
        theta2 = theta2-0.1;
        [Org, Trans] = tubegeneration(R, r ,theta1, theta2, rot, xt , yt, zt);
    end
    
    if theta2<0
        continue;
    else
        count = count +1;
%         if count >1
%             break;
%         end
    end
    
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
end

figure,hold on,
for k = 1:length(Electrode)
    tubeplot([ Electrode(k).Trans.xvector;  Electrode(k).Trans.yvector;  Electrode(k).Trans.zvector], Electrode(k).r,8);
end
axis equal
daspect([1,1,1]);
camlight;
view([124 39])



%% Output Each single rod stl fi
X =[]; Y =[]; Z = [];
for k = 1:length(Electrode)
    X = [Electrode(k).Trans.Xm];
    Y = [Electrode(k).Trans.Ym];
    Z = [Electrode(k).Trans.Zm];
    surf2stl(['.\Single_Rods\STL_Rod_',num2str(k),'.stl'],X,Y,Z, 'ascii' );
end




%% Save file

for k = 1:length(Electrode)
    output1 = [   Electrode(k).w  Electrode(k).r Electrode(k).R Electrode(k).theta1 Electrode(k).theta2 ...
                  Electrode(k).xt Electrode(k).zt Electrode(k).rot ...
                  Electrode(k).Org.x1 Electrode(k).Org.x2 Electrode(k).Org.xc ...
                  Electrode(k).Org.y1 Electrode(k).Org.y2 Electrode(k).Org.yc];
    output(k,:) = output1;
end

outputtable = array2table(output, 'VariableNames',{'w','r','R','theta1','theta2','xt','zt','rot', 'x1', 'x2', 'xc', 'y1', 'y2', 'yc'})
delete RodDistribution.csv
writetable(outputtable,'RodDistribution.csv')

%%
function [Org, Trans] = tubegeneration(R, r ,theta1, theta2, rot, xt ,yt, zt)

Org.x1 = 0;
Org.y1 = 0;
Org.x2 = -R*cosd(theta1) + R*cosd(theta1+theta2);
Org.y2 = -R*sind(theta1) + R*sind(theta1+theta2);
Org.xc = -R*cosd(theta1);
Org.yc = -R*sind(theta1);

temptheta = linspace(theta1,theta1+theta2,100);
Org.xvector =  -R*cosd(theta1) + R*cosd(temptheta);
Org.yvector =  -R*sind(theta1) + R*sind(temptheta);
Org.zvector =  0 * ones(size(temptheta));

Trans.x1t = xt;
Trans.y1t = yt;
Trans.z1t = zt;

Trans.x2t = xt + Org.x2*cosd(rot);
Trans.y2t = yt + Org.y2;
Trans.z2t = zt - Org.x2 *sind(rot);

Trans.xvector = xt + Org.xvector*cosd(rot);
Trans.yvector = yt + Org.yvector;
Trans.zvector = zt - Org.xvector*sind(rot);
end

function flag = checkboundary(Trans, xlen,zlen, r)

if (Trans.x2t+1.5*r)<xlen &  (Trans.x2t-1.5*r)> 0 &  (Trans.z2t+1.5*r)<zlen & (Trans.z2t-1.5*r) >0
    flag = 1;
else
    flag =0;
end


end
