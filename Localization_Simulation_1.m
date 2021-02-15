close all; clear; clc ; rng(0102)

snum  = 200;
gdist = 0.2;
nsig  = 0.1;
g1    = [-pi/2+gdist;-gdist]+nsig*randn(2,snum);
g2    = [pi/2+gdist;gdist]+nsig*randn(2,snum);
g3    = [pi/2+gdist;-gdist]+nsig*randn(2,snum);
g4    = [-pi/2+gdist;gdist]+nsig*randn(2,snum);
data  = [g1,g2,g3,g4];
lbls  = [ones(1,snum),2*ones(1,snum),3*ones(1,snum),4*ones(1,snum)];

%% Sensors locations

loc1x = 0;
loc1y = 1/sqrt(2)*1.1;
loc1z = 1/sqrt(2)*1.1;

loc2x = 0;
loc2y = -1/sqrt(2)*1.1;
loc2z = 1/sqrt(2)*1.1;

loc3x = 0;
loc3y = 0;
loc3z = 1.1;

loc4x = 1/sqrt(2)*1.1;
loc4y = 0;
loc4z = 1/sqrt(2)*1.1;

loc5x = -1/sqrt(2)*1.1;
loc5y = 0;
loc5z = 1/sqrt(2)*1.1;

%% Project 2D trajectories on a 3D semi-sphere
r        = 1;
loc(1,:) = r * cos(data(2,:)) .* sin(data(1,:));
loc(2,:) = r * sin(data(2,:)) .* sin(data(1,:));
loc(3,:) = r * cos(data(1,:));

figure
[Xv,Yv,Zv] = sphere;
surf(Xv,Yv,Zv,ones(length(Xv))); set(gca,'FontSize',12)
colormap([0.7,0.7,0.7])
alpha 0.7
hold on
DispRange = 1:size(loc,2);
hPlotHigh1 = scatter3(loc(1,lbls==1),loc(2,lbls==1),loc(3,lbls==1),'b','fill');
hPlotHigh2 = scatter3(loc(1,lbls==2),loc(2,lbls==2),loc(3,lbls==2),'r','fill');
hPlotHigh3 = scatter3(loc(1,lbls==3),loc(2,lbls==3),loc(3,lbls==3),'k','fill');
hPlotHigh4 = scatter3(loc(1,lbls==4),loc(2,lbls==4),loc(3,lbls==4),'y','fill');

hScat1 = scatter3(loc1x,loc1y,loc1z,1000,'g','h','fill'); %text(loc1x,loc1y,loc1z,'1','Fontsize',20);
scatter3(loc2x,loc2y,loc2z,200,'g','h','fill'); text(loc2x,loc2y,loc2z,'2','Fontsize',20);
scatter3(loc3x,loc3y,loc3z,200,'g','h','fill'); text(loc3x,loc3y,loc3z,'3','Fontsize',20);
scatter3(loc4x,loc4y,loc4z,200,'g','h','fill'); text(loc4x,loc4y,loc4z,'4','Fontsize',20);
scatter3(loc5x,loc5y,loc5z,200,'g','h','fill'); text(loc5x,loc5y,loc5z,'5','Fontsize',20);
set(gca,'FontSize',20,'View',[112,8])
axis off;axis equal;

%% Poisson rates
fac1   = 20;
rate1 = fac1*exp(-sqrt((loc(1,:) - loc1x).^2+(loc(2,:) - loc1y).^2+(loc(3,:) - loc1z).^2));
rate2 = fac1*exp(-sqrt((loc(1,:) - loc2x).^2+(loc(2,:) - loc2y).^2+(loc(3,:) - loc2z).^2));
rate3 = fac1*exp(-sqrt((loc(1,:) - loc3x).^2+(loc(2,:) - loc3y).^2+(loc(3,:) - loc3z).^2));
rate4 = fac1*exp(-sqrt((loc(1,:) - loc4x).^2+(loc(2,:) - loc4y).^2+(loc(3,:) - loc4z).^2));
rate5 = fac1*exp(-sqrt((loc(1,:) - loc5x).^2+(loc(2,:) - loc5y).^2+(loc(3,:) - loc5z).^2));

%%  Mmeasurements
nrel    = 100; 
y1 = repmat(rate1,nrel,1).*randn(nrel,length(rate1));
y2 = repmat(rate2,nrel,1).*randn(nrel,length(rate2));
y3 = repmat(rate3,nrel,1).*randn(nrel,length(rate2));
y4 = repmat(rate4,nrel,1).*randn(nrel,length(rate2));
y5 = repmat(rate5,nrel,1).*randn(nrel,length(rate2));

%% Collect trajectories
y = [y1; y2; y3; y4; y5];
figure
yyaxis left
imagesc(y);colormap(jet)
set(gca,'Xtick',[100,300,500,700],'Xticklabel',{'Blue','Red','Black','Yellow'});
set(gca,'Ytick',50:100:500,'Yticklabel',...
    {'Sensor 1','Sensor 2','Sensor 3','Sensor 4','Sensor 5'},'Ycolor','k')
xlabel('Region');
title('Sensor measurement, $$f_{j}(x_{i})$$','interpreter','latex');
yyaxis right
imagesc(y)
ylabel('Dimension');
set(gca,'Ytick',50:100:500,'Yticklabel',{'100','100','100','100','100'},'Ycolor','k')
line([200,200],[0,501],'linewidth',2,'color','k')
line([400,400],[0,501],'linewidth',2,'color','k')
line([600,600],[0,501],'linewidth',2,'color','k')

line([0,801],[100,100],'linewidth',2,'color','k')
line([0,801],[200,200],'linewidth',2,'color','k')
line([0,801],[300,300],'linewidth',2,'color','k')
line([0,801],[400,400],'linewidth',2,'color','k')
c = colorbar;
ylabel(c, 'Measurement value')
set(gca,'Fontsize',20)

%% Histograms

histYN = 0;

if histYN
nbins   = 10;         % # of histogram bins
datAll  = {y1.',y2.',y3.',y4.',y5.'};
datHist = cell(1,length(datAll));

for nd = 1:length(datAll)
    datHist{nd} = zeros(nbins,length(lbls));
    if nd <= 2
        hist_bins = linspace(-20,20,nbins);%min(datAll{nd}(:)), max(datAll{nd}(:)), nbins);%linspace(0,20,nbins);%
    else 
        hist_bins = linspace(0,20,nbins);
    end
    for ii=1:length(lbls)
        datHist{nd}(:,ii) = (hist(datAll{nd}(ii,:), hist_bins))';
    end
end

figure
imagesc(cell2mat(datHist.'))

end
%% Visualization
% Baselines:
colorsch = [0 0 1; ... % blue
            1 0 0; ... % red
            0 0 0; ... % black
            1 1 0];    % yellow

perplexity = 30;
newDims = 2;

if histYN
    datAll   = [{cell2mat(datHist.')}, datHist];
else
    datAll   = {y,y1,y2,y3,y4,y5};
end
Psi      = cell(1,length(datAll));

for bn = 1:length(datAll)
    dist = squareform(pdist(datAll{bn}')); 
    W = exp(-dist.^2/(10*median(dist(:))^2));
    W = diag(1./sum(W,1))*W;
    [Psi{bn},~] = eigs(W,4);
    Psi{bn} = Psi{bn}(:,2:end);
end


%% SSD
SSD = zeros(5,snum);
for sn = 1:size(y,2)
        ytmp = [y1(:,sn),y2(:,sn),y3(:,sn),y4(:,sn),y5(:,sn)];
    distSSD = squareform(pdist(ytmp.'));
    W = exp(-distSSD.^2/(median(distSSD(:))^2));
    W = diag(1./sum(W,1))*W;
    [SSD(:,sn),~] = eigs(W.',1);
    SSD(:,sn) = sign(SSD(1,sn))*SSD(:,sn)/sum(abs(SSD(:,sn)));
end

figure
yyaxis left
imagesc(SSD);title('SSD matrix')
set(gca,'Xtick',[100,300,500,700],'Xticklabel',{'Blue','Red','Black','Yellow'});
set(gca,'Ytick',1:5,'Yticklabel',{'Sensor 1','Sensor 2',...
    'Sensor 3','Sensor 4','Sensor 5'},'Ycolor','k')
xlabel('Region'); 
yyaxis right
imagesc(SSD)
ylabel('Dimension');
set(gca,'Ytick',1:5,'Yticklabel',{'1','1','1','1','1'},'Ycolor','k');
set(gca,'Fontsize',20);c = colorbar;
c.Label.String = 'SSD value';
set(gca,'Fontsize',20)
line([200,200],[0,6],'linewidth',2,'color','k')
line([400,400],[0,6],'linewidth',2,'color','k')
line([600,600],[0,6],'linewidth',2,'color','k')

%%
distSSDall = squareform(pdist(SSD.','minkowski',1)); 
W = exp(-distSSDall.^2/(10*median(distSSDall(:))));
W = diag(1./sum(W,1))*W;
[PsiSSD,~] = eigs(W,4);
figure
scatter3(PsiSSD(:,2),PsiSSD(:,3),PsiSSD(:,4),[],colorsch(lbls,:),'fill'); 
xlabel('\psi_2');ylabel('\psi_3');zlabel('\psi_4');set(gca,'Fontsize',25)
PsiSSD = PsiSSD(:,2:end);
%% Classification

% inpts = {y.',y1.',y2.',y3.',SSD.'};
inpts  = [Psi,{PsiSSD}];
titles = {'All sensors','Sensor 1','Sensor 2','Sensor 3','Sensor 4','Sensor 5','SSD'};
plt_CM = 1;
for nin = 7:length(inpts)
    t           = templateSVM('Standardize',1,'KernelFunction','rbf');
    mdl_conn    = fitcecoc(inpts{nin},lbls,'Learners',t);
    CVmdl_conn  = crossval(mdl_conn);
    osLoss_conn = kfoldLoss(CVmdl_conn);
    osPred_conn = kfoldPredict(CVmdl_conn);
    conf_conn   = confusionmat(osPred_conn,lbls);
    if (plt_CM == 1)
        % Convert this data to a [numClasses x numSample.] matrix
        targets = zeros(4,800);
        outputs = zeros(4,800);
        outputsIdx = sub2ind(size(targets), lbls', [1:800]');
        targetsIdx = sub2ind(size(outputs), osPred_conn,[1:800]');
        targets(targetsIdx) = 1;
        outputs(outputsIdx) = 1;
        % Plot the confusion matrix for a 4-class problem
        figure
        plotconfusion(targets,outputs)
        
        h = gca;
        h.XTickLabel = {'Blue','Red','Black','Yellow',' '};
        h.YTickLabel = {'Blue','Red','Black','Yellow',' '};
        h.YTickLabelRotation = 90;
        h.XTickLabelRotation = 0;
        set(h,'FontSize',18)
        set(findobj(gca,'type','text'),'fontsize',20) 
    end
end


