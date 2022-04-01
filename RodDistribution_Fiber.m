clc
clear all
close all

addpath(genpath('.\Functions'))
rng(1)

Para.xlen = 1e-3;  %battery depth [m];
Para.ylen = 4e-3;    %battery width [m];
Para.zlen = 1e-3;   % battery height [m];
Para.r = 6e-6;         % m, electrode radius
Para.theta1_range = 10;    % start angle range;
Para.theta2_range = 20;    % arc angle range
Para.ap = 34800; % specific area [m2/m3]
Para.porosity = 0.67; % porosity [-]
Para.Boundary_Restriction = 1;  % control if the rod can penetrate the battery domain
Para.volume = Para.xlen * Para.ylen * Para.zlen;
%% Generate One Typical Fiber
count = 0;
Electrode = [];
tempap = 0;
tempporosity =0;
for k = 1:1395
    disp(['Generating ', num2str(k),'th electrode'])
    %     if k < 50
    [Electrode] = ElectrodeGeneration(Electrode, Para, 'y');
    %     elseif k < 100
    %     [Electrode] = ElectrodeGeneration(Electrode, Para, 'z')
    %     else
    %     [Electrode] = ElectrodeGeneration(Electrode, Para, 'x')
    %     end
    tempap = tempap + Electrode(k).area/Para.volume;
    tempporosity = tempporosity + Electrode(k).volume/Para.volume;
end
tempap
1-tempporosity


%%
figure,hold on,
for k = 1:length(Electrode)
    tubeplot([ Electrode(k).Trans.xvector;  Electrode(k).Trans.yvector;  Electrode(k).Trans.zvector], Electrode(k).r,8);
end
axis equal
daspect([1,1,1]);
camlight;
view([124 39])



%% Output Each single rod stl fi
% X =[]; Y =[]; Z = [];
% for k = 1:length(Electrode)
%     
%     X = [Electrode(k).Trans.Xm];
%     Y = [Electrode(k).Trans.Ym];
%     Z = [Electrode(k).Trans.Zm];
%     surf2stl(['.\Single_Rods\STL_Rod_',num2str(k),'.stl'],X,Y,Z, 'ascii' );
% end




%% Save file

for k = 1:length(Electrode)
    output1 = [   Electrode(k).w  Electrode(k).r Electrode(k).R Electrode(k).theta1 Electrode(k).theta2 ...
        Electrode(k).xt Electrode(k).zt Electrode(k).rot ...
        Electrode(k).Org.x1 Electrode(k).Org.x2 Electrode(k).Org.xc ...
        Electrode(k).Org.y1 Electrode(k).Org.y2 Electrode(k).Org.yc];
    output(k,:) = output1;
end

output_scaled = output;
output_scaled(:,[1 2 3 6 7 9:14]) =  output_scaled(:,[1 2 3 6 7 9:14]) *1E2;

outputtable = array2table(output_scaled, 'VariableNames',{'w','r','R','theta1','theta2','xt','zt','rot', 'x1', 'x2', 'xc', 'y1', 'y2', 'yc'})
delete RodDistribution.csv
writetable(outputtable,'RodDistribution.csv')

%%
function [Electrode] = ElectrodeGeneration(Electrode, Para, GrowDirection)
if strcmp(GrowDirection,'y')
    thickness = Para.ylen;
    lenA = Para.xlen;
    lenB = Para.zlen;
    r = Para.r;
    theta1_range = Para.theta1_range;
    theta2_range = Para.theta2_range;
elseif strcmp(GrowDirection,'z')
    thickness = Para.zlen;
    lenA = Para.ylen;
    lenB = Para.xlen;
    r = Para.r;
    theta1_range = Para.theta1_range;
    theta2_range = Para.theta2_range;
elseif strcmp(GrowDirection,'x')
    thickness = Para.xlen;
    lenA = Para.zlen;
    lenB = Para.ylen;
    r = Para.r;
    theta1_range = Para.theta1_range;
    theta2_range = Para.theta2_range;
end

Org = [];
Trans =[];
while isempty(Trans)
    w = thickness*0.95; % m
    theta1 = theta1_range *rand(); % degree
    theta2 = theta2_range *rand();
    
    % tranlation
    xt = r + (lenA-2*r)*rand() ;   % xt
    zt = r + (lenB-2*r)*rand();    % zt
    yt = 0;
    
    rot = 360*rand();
    R = w/(sind(theta1 + theta2) - sind(theta1));
    [Org, Trans] = tubegeneration(R, r ,theta1, theta2, rot, xt, yt, zt);
    
    while checkboundary(Trans, lenA, lenB, r) ==0 & theta2>0
        theta2 = theta2-0.1;
        if theta2 > 0
            [Org, Trans] = tubegeneration(R, r ,theta1, theta2, rot, xt , yt, zt);
        else
            Trans =[];
        end
    end
end


if strcmp(GrowDirection,'z')
    temp = Trans;
    temp.xvector = Trans.zvector;
    temp.yvector = Trans.xvector;
    temp.zvector = Trans.yvector;
    Trans = temp;
elseif strcmp(GrowDirection,'x')
    temp = Trans;
    temp.xvector = Trans.yvector;
    temp.yvector = Trans.zvector;
    temp.zvector = Trans.xvector;
    Trans = temp;
end

[Trans.Xm,Trans.Ym,Trans.Zm] = tubeplot([Trans.xvector; Trans.yvector; Trans.zvector],r,8);
Electrode(end+1).Trans = Trans;
Electrode(end).Org = Org;
Electrode(end).theta1 = theta1;
Electrode(end).theta2 = theta2;
Electrode(end).R = R;
Electrode(end).r = r;
Electrode(end).rot = rot;
Electrode(end).w = w;
Electrode(end).xt = xt;
Electrode(end).zt = zt;
Electrode(end).area = 4*pi^2*  Electrode(end).R*  Electrode(end).r * (  Electrode(end).theta2/360);
Electrode(end).volume = pi^2*  Electrode(end).r^2 * (2*pi*Electrode(end).R)*(  Electrode(end).theta2/360);
end

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
if isempty(Trans)
    flag = 0;
elseif (Trans.x2t+1.5*r)<xlen &  (Trans.x2t-1.5*r)> 0 &  (Trans.z2t+1.5*r)<zlen & (Trans.z2t-1.5*r) >0
    flag = 1;
else
    flag = 0;
end


end
