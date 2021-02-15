close all;clear; clc;rng(0102);
snum  = 200;
gdist = 0;
nsig  = 0.1;

g1    = [-pi/2;-gdist]+nsig*randn(2,snum);
g4    = [pi/2;gdist]+nsig*randn(2,snum);

data  = [g1,g4];
lbls  = [ones(1,snum),2*ones(1,snum)];

%% Sensors locations
Sensor_num = 12;
height_num = 1;

LL = zeros(Sensor_num*height_num,3);
R = 1.1;

LL(1,1) = 0;
LL(1,2) = 0;
LL(1,3) = R;

LL(2,1) = 0;
LL(2,2) = 1/2*R;
LL(2,3) = sqrt(3)/2*R;

LL(3,1) = 0;
LL(3,2) = sqrt(3)/2*R;
LL(3,3) = 1/2*R;


LL(4,1) = 0;
LL(4,2) = R;
LL(4,3) = 0;

LL(5,1) = 0;
LL(5,2) =  sqrt(3)/2*R;
LL(5,3) = -1/2*R;

LL(6,1) = 0;
LL(6,2) = 1/2*R;
LL(6,3) = - sqrt(3)/2*R;


LL(7,1) = 0;
LL(7,2) = 0;
LL(7,3) = -R;

LL(8,1) = 0;
LL(8,2) = -1/2*R;
LL(8,3) = - sqrt(3)/2*R;

LL(9,1) = 0;
LL(9,2) =  -sqrt(3)/2*R;
LL(9,3) = -1/2*R;

LL(10,1) = 0;
LL(10,2) = -R;
LL(10,3) = 0;

LL(11,1) = 0;
LL(11,2) = -sqrt(3)/2*R;
LL(11,3) = 1/2*R;

LL(12,1) = 0;
LL(12,2) =- 1/2*R;
LL(12,3) = sqrt(3)/2*R;


%% Project 2D trajectories on a 3D semi-sphere
r        = 1;
loc(1,:) = r * cos(data(2,:)) .* sin(data(1,:));
loc(2,:) = r * sin(data(2,:)) .* sin(data(1,:));
loc(3,:) = r * cos(data(1,:));
%% center at two class
ord = zeros(2,1);
for oo = 1:2
    [~,ord(oo)] = min(sum((loc(:,lbls==oo) - mean(loc(:,lbls==oo),2)).^2));
end

figure
[Xv,Yv,Zv] = sphere;
surf(Xv,Yv,Zv,ones(length(Xv))); set(gca,'FontSize',20)
colormap([0.7,0.7,0.7])
alpha 0.7
hold on
DispRange = 1:size(loc,2);
hPlotHigh1 = scatter3(loc(1,lbls==1),loc(2,lbls==1),loc(3,lbls==1),'b','fill');
hPlotHigh2 = scatter3(loc(1,lbls==2),loc(2,lbls==2),loc(3,lbls==2),'r','fill');
hPlotHigh3 = scatter3(loc(1,lbls==3),loc(2,lbls==3),loc(3,lbls==3),'k','fill');
hPlotHigh4 = scatter3(loc(1,lbls==4),loc(2,lbls==4),loc(3,lbls==4),'y','fill');
loc_ = LL;
for hh = 1:length(LL)
    
        scatter3(LL(hh,1),LL(hh,2),LL(hh,3),200,'g','h','fill'); 
        text(loc_(hh,1),loc_(hh,2),loc_(hh,3),num2str(hh),'Fontsize',20);
end
set(gca,'View',[112,8]);axis equal
axis off
%% parameters
eps_num = 15;eps = logspace(-0.5,2,eps_num);
fac1    = 20;
nrel    = 100;
%%
repeat_time = 200;Q_prob = zeros(eps_num,Sensor_num);
mu = zeros(Sensor_num,1);
sigma1 = triu(abs(normrnd(0,0.5,Sensor_num,Sensor_num)),1);
sigma1 = sigma1 + sigma1';
sigma2 = triu(abs(normrnd(0,0.5,Sensor_num,Sensor_num)),1);
sigma2 = sigma2 + sigma2';


%%Realization

rate = zeros(snum*2,Sensor_num);
for s_rate = 1:Sensor_num*height_num
    rate(:,s_rate) = fac1*exp(-sqrt((loc(1,:) - loc_(s_rate,1)).^2+(loc(2,:) - loc_(s_rate,2)).^2+(loc(3,:) - loc_(s_rate,3)).^2));
end
y = zeros(nrel*Sensor_num,snum*2);
measure_y = cell(1,1);
for pos = 1:snum*2
    if(pos<snum+1)
        sigma1(logical(eye(Sensor_num))) = rate(pos,:);
        measure_y{pos} = mvnrnd(mu,sigma1,100);
        y(:,pos) = reshape(measure_y{pos},[],1);
    else
        sigma2(logical(eye(Sensor_num))) = rate(pos,:);
        measure_y{pos} = mvnrnd(mu,sigma1,100);
        y(:,pos) = reshape(measure_y{pos},[],1);
    end
end


plt_sensorM = 1;
if plt_sensorM == 1
    % Collect trajectories
    figure
    yyaxis left
    imagesc(y);colormap(jet);title('Sensor measurement, $$f_{j}(x_{i})$$','interpreter','latex');colorbar;
    set(gca,'Xtick',[100,300],'Xticklabel',{'Blue','Red'});
    set(gca,'Ytick',50:100:1150,'Yticklabel',...
        {'Sensor 1','Sensor 2','Sensor 3','Sensor 4','Sensor 5',...
        'Sensor 6','Sensor 7','Sensor 8','Sensor 9','Sensor 10',...
        'Sensor 11','Sensor 12','Sensor 13','Sensor 14','Sensor 15',...
        'Sensor 16','Sensor 17','Sensor 18','Sensor 19','Sensor 20'},...
        'Ycolor','k','fontsize',20);
    xlabel('Region'); 
    set(gca,'Fontsize',20)
    yyaxis right
    h = imagesc(y);
    ax = ancestor(h, 'axes');
    line([200,200],[0,length(y)+4],'linewidth',2,'color','k')
    for ss = 1:Sensor_num*height_num-1
        line([0,501],[100*ss+1,100*ss+1],'linewidth',2,'color','k')
    end
    ax.YTick = 50:100:1150;
    ax.YTickLabel =  mat2cell(repmat(100,fac1,1),ones(fac1,1));
    ax.YColor = 'k';
    ylabel('Dimension');
    c = colorbar;
    c.Label.String = 'Measurement value';
    set(gca,'Fontsize',20)
end
Psi = cell(1,1);
dist = squareform(pdist(y.'));
W = exp(-dist.^2/(10*median(dist(:))^2));
W = diag(1./sum(W,1))*W;
[Psi{1},~] = eigs(W,4);
    
for bn = 1:Sensor_num
    dist = squareform(pdist(y((bn-1)*nrel+1:bn*nrel,:).')); 
    W = exp(-dist.^2/(10*median(dist(:))^2));
    W = diag(1./sum(W,1))*W;
    [Psi{bn+1},~] = eigs(W,4);
      
    Psi{bn+1} = Psi{bn+1}(:,2:end);
end

SSD = cell(1,1);
M_mea = cell(1,1);
C_mea = cell(1,1);
for pos2 = 1:snum*2
    M_mea{pos2} = mean(measure_y{pos2});
    t1 = cov(measure_y{pos2});
    t1 = t1 - diag(diag(t1));
    C_mea{pos2} = mean(t1); 
    distSSD = squareform(pdist(measure_y{pos2}'));
    W = exp(-distSSD.^2/(10*median(distSSD(:))^2));
    SSD{pos2} = sum(W)/sum(sum(W)); 
end

M_mea_S2 = squeeze(cat(1,M_mea{:}))';
C_mea_S2 = squeeze(cat(1,C_mea{:}))';
SSD = cat(1,SSD{:});

dist_SSD = squareform(pdist(SSD,'minkowski',1)); 
W = exp(-dist_SSD.^2/(10*median(dist_SSD(:))));
W = diag(1./sum(W,1))*W;
[PsiSSD,~] = eigs(W,4);
%%

inpts = [{PsiSSD},Psi];
titles = {'\pi_{i}','Concatenated','Sensor 1','Sensor 2','Sensor 3','Sensor 4','Sensor 5',...
    'Sensor 6','Sensor 7','Sensor 8','Sensor 9','Sensor 10',...
    'Sensor 11','Sensor 12'};
blue_region = zeros(length(inpts),1);
red_region = zeros(length(inpts),1);
for nin = 1:length(inpts)
    t           = templateSVM('Standardize',1,'KernelFunction','linear');
    mdl_conn    = fitcecoc(inpts{nin},lbls,'Learners',t);
    CVmdl_conn  = crossval(mdl_conn);
    osLoss_conn = kfoldLoss(CVmdl_conn);
    osPred_conn = kfoldPredict(CVmdl_conn);
    conf_conn   = confusionmat(osPred_conn,lbls);
    conf_conn   = bsxfun(@rdivide,conf_conn,sum(conf_conn,1));
    blue_region(nin) = conf_conn(1,1);
    red_region(nin) = conf_conn(2,2);

end
%% plot accuracy
ax = gca;
acc_R = cat(2,blue_region,red_region,0.5*(blue_region+red_region));
figure
bar(acc_R);title(' Classification Accuracy')
legend('Blue Region','Red Region','Mean Acc')
legend('show');xtickangle(45);
xticklabels(titles);%ax.TickLabelInterpreter = 'tex';
set(gca,'FontSize',18)

