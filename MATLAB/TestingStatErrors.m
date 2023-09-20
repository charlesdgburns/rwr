%Marek et al. first generate Correlation matrices and then compare to full
%sample correlation pvalues

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATING RANDOM DATA TO TEST %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

rng(2022) %For reproducibility

nSubjects = 1000;
nROIs = 50; %number of regions of interest we are correlating between.
nEdges = nchoosek(nROIs,2);
factor = 15*normrnd(1, 0.5, [nSubjects,1]); %Here we simulate some random factor scored between 0 and 15.
binsize =[25,33,50,70, 100, 135, 200, 265, 375, 525, 725, 1000]; %(so 3 different bin sizes)
iter = 50; %number of iterations


for s = 1:nSubjects
    X = randn(15, nROIs); %the 15 here is arbitrary but represent timepoints of resting state fMRI.
    C= corr(X);    
    rmats(s,:)=C(logical(triu(C,1))); 
end


%We may visualise this in a scatter plot. Given a 'random' population, there should be no covariance here.
randomData = [factor, rmats(:,:)];

figure
scatter(randomData(:,1), randomData(:,2));

%To demonstrate the relation to full sample size,
%We also generate data with 10'000 which we later subsample up to 1000 of;

% 10K random data generation %
nSubjects2 = 10000;
factor2 = 15*normrnd(1, 0.5, [nSubjects2,1]);

for s = 1:nSubjects2
    X = randn(15, nROIs); %the 15 here is arbitrary but represent timepoints of resting state fMRI.
    C= corr(X);    
    rmats2(s,:)=C(logical(triu(C,1))); 
end

rnd2Data = [factor2, rmats2(:,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COPULA 'TRUE' EFFECT SIMULATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Here we simulate data which covaries according to some covariance.
%Since our data is non-normal, we do this using copulas based on generated
%data above.


nTrueEffects = 100; %Number of 'true effects' we want to simulate
rho=0.4; %Size of simulated 'true effects'

nu=50; %This essentially specifies how 'gaussian' our tau distribution is.
N=nSubjects; %sample size of originally generated data.
N2=nSubjects2; 
n=1000; %sample size of COPULA generated data 
n2=10000;

% We comment this out here, but may be used to simulate data:
% tauData = corr(data, 'type', 'Kendall'); %Model data covariation
% RhoData = sin(pi*tauData/2); %Covariance matrix

data = [factor, rmats(:,1:nTrueEffects)];

covMat = repmat(rho, nTrueEffects+1); %Covariance matrix. 
% Note here that technically each edge will now also correlate with each other, b
% However, this is done as a simple way to ensure a symmetric positive
% semi-definite covariance matrix.


for i=1:length(covMat)
covMat(i,i)=1;
end

% Generating data through copulas: 

for v = 1:size(data,2)
EinvCDF{v}=sort(data(:,v));
ns(v)= length(data(:,v));
SinvCDF{v}=ksdensity(data(:,v), EinvCDF{v}, 'function', 'cdf', 'width',.05);
end

nobs = size(rmats,1);

T = mvtrnd(covMat, nu, n); %Generate data modelling covariance using t-distribution.
U = tcdf(T,nu); 

for v=1:size(data,2) %for each variable (edges and factor)
CS(:,v)= [EinvCDF{v}(ceil(ns(1)*U(:,v)))]; %Send the data generated under uniform distribution through our inverse CDF
smCS(:,v) = pchip((1:nobs)/nobs, EinvCDF{v}, U(:,v) ); %smoothed estimates 
end

scatter(smCS(:,1), smCS(:,2))

%The data below takes factors which are correlated with the first
%nTrueEffects edges. 
simData = [smCS(:,:), randomData(:,(nTrueEffects+2):(nEdges+1))];


%% 10k sample true effects %%

data2 = [factor2, rmats2(:,1:nTrueEffects)];

for v = 1:size(data2,2)
EinvCDF2{v}=sort(data2(:,v));
ns2(v)= length(data2(:,v));
SinvCDF2{v}=ksdensity(data2(:,v), EinvCDF2{v}, 'function', 'cdf', 'width',.05);
end

nobs2 = size(rmats2,1);

T = mvtrnd(covMat, nu, n2); %Generate data modelling covariance using t-distribution.
U = tcdf(T,nu); 

for v=1:size(data2,2) %for each variable (edges and factor)
CS2(:,v)= [EinvCDF2{v}(ceil(ns2(1)*U(:,v)))]; %Send the data generated under uniform distribution through our inverse CDF
smCS2(:,v) = pchip((1:nobs2)/nobs2, EinvCDF2{v}, U(:,v) ); %smoothed estimates 
end

scatter(smCS2(:,1), smCS2(:,2))

%The data below takes factors which are correlated with the first
%nTrueEffects edges. 
sim2Data = [smCS2(:,:), rnd2Data(:,(nTrueEffects+2):(nEdges+1))];


%%%%%%%%%%%%%%
%% SCRIPT 1 %%
%%%%%%%%%%%%%%

%We run the random data through bootstrapped resampling and  

%INPUT: the random generated values above


% 1k full sample %
[randomCorrs,randomPvals] = abcd_edgewise_correlation_iterative_reliability_single_factor(randomData(:,2:(nEdges+1)),randomData(:,1),binsize,iter);

[simCorrs,simPvals] = abcd_edgewise_correlation_iterative_reliability_single_factor(simData(:,2:(nEdges+1)),simData(:,1),binsize,iter);


% 10k full sample %
[rnd2Corrs,rnd2Pvals] = abcd_edgewise_correlation_iterative_reliability_single_factor(rnd2Data(:,2:(nEdges+1)),rnd2Data(:,1),binsize,iter);

[sim2Corrs,sim2Pvals] = abcd_edgewise_correlation_iterative_reliability_single_factor(sim2Data(:,2:(nEdges+1)),sim2Data(:,1),binsize,iter);


%The function above is edited only to also output 'TheseBehav' and 'TheseMats', which are the last bootstrap sample iteration, which we use
%as an example of a bootstrap sample at full sample size.

[fscorrs, fspvals] = corr(randomData(:,2:(nEdges+1)),randomData(:,1)); %Simply generating NxP matrix for all subjects (hence all factor scores)

[simfscorrs, simfspvals] = corr(simData(:,2:(nEdges+1)),simData(:,1)); %Simply generating NxP matrix for all subjects (hence all factor scores)

[fs2corrs, fs2pvals] = corr(rnd2Data(:,2:(nEdges+1)),rnd2Data(:,1)); %Simply generating NxP matrix for all subjects (hence all factor scores)

[sim2fscorrs, sim2fspvals] = corr(sim2Data(:,2:(nEdges+1)),sim2Data(:,1)); %Simply generating NxP matrix for all subjects (hence all factor scores)

%%%%%%%%%%%%%%
%% FIGURE 1 %%
%%%%%%%%%%%%%%

%We plot the underlying method, scattering connectivity across factors and
%correlating, for both 'full sample' and bootstrap resample.

figure
% t=tiledlayout(1,2)

% id = find(Corrs(:,iter,size(binsize,2)) == min(Corrs(:,iter,size(binsize,2))))

id = 1; %Specify which edge you wish to plot. For example, we take the first.

ax1=nexttile
scatter(rmats(:,id), factor)
r1=round(fscorrs(id,:),3,'significant')
p1=round(fspvals(id,:),3,'significant')
title({'Random Sample (N=1000)';['r = ' num2str(r1) ' p = ' num2str(p1) ' (Pearson)']})
hold on
lsline
xlabel('Edge 1 connectivity estimates')
ylabel('Factor Score')

%% OBSOLETE - see line 215

%ax2=nexttile
%scatter(randomEdgeSample(:,id),randomFactorSample)
% r2=round(Corrs(id,iter,size(binsize,2)),3,'significant')
% p2=round(Pvals(id,iter,size(binsize,2)),3,'significant')
% title({'Bootstrap Resample (N=1000)';['r = ' num2str(r2) ' p = ' num2str(p2) ' (Pearson)']})
% hold on
% lsline
% xlabel('Edge 1 connectivity estimates')
% ylabel('Factor Score')

% linkaxes([ax1,ax2], 'x','y')

sgtitle('Example edge-behaviour correlation (Edge 1)')

%%%%%%%%%%%%%%
%% FIGURE 2 %%
%%%%%%%%%%%%%%

%We compare Original Sample (OS) and Bootstrap Resample (BR) distributions.

% Corrs = randomCorrs;
% Pvals = randomPvals;
% TheseMats = randomEdgeSample;

%% OBS figure obsolete as simEdgeSample // Thesemats cannot be extracted in parallel loops.
%% This will be re-introduced in R, however.

% Corrs = simCorrs;
% Pvals = simPvals;
% TheseMats = simEdgeSample;
% 
% figure
% t=tiledlayout(3,2)
% 
% %Plot 1: OS rmats
% ax1=nexttile
% hist(rmats)
% title('a')
% xlabel({'Edge-level connectivity estimates'})
% 
% %Plot 2: RS TheseMats
% ax2=nexttile
% hist(TheseMats)
% title('d')
% xlabel({'Edge-level connectivity estimates'})
% 
% %Plot 3: OS corrs
% ax3=nexttile
% hist(fscorrs)
% title({'b'})
% ylabel({'Frequency'})
% xlabel({'Edge-behaviour correlations (Pearson)'})
% 
% %Plot 4: RS corrs
% ax4=nexttile
% hist(Corrs(:,iter,size(binsize,2)))
% title({'e'})
% xlabel({'Edge-behaviour correlations (Pearson)'})
% 
% %Plot 5: OS pvals
% ax5=nexttile
% hist(fspvals)
% title({'e'})
% xlabel({'Edge-behaviour p-values (Pearson)'})
% 
% %Plot 6: BR pvals
% ax6=nexttile
% hist(Pvals(:,iter,size(binsize,2)))
% title({'f'})
% xlabel({'Edge-behaviour p-values (Pearson)'})
% 
% 
% linkaxes([ax1,ax2],'y', 'x')
% linkaxes([ax3, ax4], 'y','x')
% linkaxes([ax5, ax6], 'y','x')
% 
% sgtitle({'Random Sample (N=1000)                Bootstrap resample (N=1000) '})

%%%%%%%%%%%%%%
%% SCRIPT 2 %%
%%%%%%%%%%%%%%

%We run the statistical errors script

[rndtype1,rndtype2,rndtypem,rndtypes,rndtypes_pvals] = abcd_statisticalerrors(randomCorrs,randomPvals,fscorrs,fspvals,iter)

[simtype1,simtype2,simtypem,simtypes,simtypes_pvals] = abcd_statisticalerrors(simCorrs,simPvals,simfscorrs,simfspvals,iter)

[rnd2type1,rnd2type2,rnd2typem,rnd2types,rnd2types_pvals] = abcd_statisticalerrors(rnd2Corrs,rnd2Pvals,fs2corrs,fs2pvals,iter)

[sim2type1,sim2type2,sim2typem,sim2types,sim2types_pvals] = abcd_statisticalerrors(sim2Corrs,sim2Pvals,sim2fscorrs,sim2fspvals,iter)

%Here we can plot as Marek et al. 2020, or generally just check

%%%%%%%%%%%%%%
%% FIGURE 3 %%
%%%%%%%%%%%%%%

%Here we explicitly take a key statistical error estimate: statistical power. We note that resampling bias effects other statistical errors also.

thr = size(simtype2,2);
cmap = [];
cmap(1,:) = [1 100 80]./255;
cmap(2,:) = [103 169 207]./255;
cmap(3,:) = [208 209 230]./255;
[x,y]=meshgrid([1:3],[1:thr]);
cmap = interp2(x([1,round(thr/2),thr],:),y([1,round(thr/2),thr],:),cmap,x,y);

figure;
t= tiledlayout(2,2)

%Plot 1: random effects; full sample 1000; resample 1000
ax1=nexttile

hold on 
for i = 1:thr
    h = plot(1:size(rndtype2,1),rndtype2(:,i),'Color',cmap(i,:),'LineWidth',2);
end

xlim([1 size(rndtype2,1)])
xticks([1:1:size(rndtype2,1)]) 
xticklabels(binsize);
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Type II Error','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 100])
title({'rho = 0, N=1000, r=1000'})

%Plot 2: true effects; full sample 1000; resample 1000

ax2=nexttile 

hold on 
for i = 1:thr
    h = plot(1:size(simtype2,1),simtype2(:,i),'Color',cmap(i,:),'LineWidth',2);
end

xlim([1 size(simtype2,1)])
xticks([1:1:size(simtype2,1)]) 
xticklabels(binsize);
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Type II Error','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 100])
title({'rho = 0.4, N=1000, r=1000'})

%Plot 3: random effects; full sample 10000; resample 1000
ax3=nexttile

hold on 
for i = 1:thr
    h = plot(1:size(rnd2type2,1),rnd2type2(:,i),'Color',cmap(i,:),'LineWidth',2);
end

xlim([1 size(rnd2type2,1)])
xticks([1:1:size(rnd2type2,1)]) 
xticklabels(binsize);
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Type II Error','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 100])
title({'rho = 0, N=10000, r=1000'})

%Plot 4: true effects; full sample 10000; resample 1000

ax4=nexttile 

hold on 
for i = 1:thr
    h = plot(1:size(sim2type2,1),sim2type2(:,i),'Color',cmap(i,:),'LineWidth',2);
end

xlim([1 size(sim2type2,1)])
xticks([1:1:size(sim2type2,1)]) 
xticklabels(binsize);
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Type II Error','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 100])
title({'rho = 0.4, N=10000, r=1000'})


%% Type M plotting %%

% Type M plotting
inflations = [4 9 17];
figure; hold on 
% plot the uncorrected line
%h=plot([1:16],typem_uncorrected_both','-.k*','LineWidth',2);

mtypem = [];    
for f = 1:size(rndtypem,1)
    AllInflatedCorrs = rndtypem{f};
    for x = 1:length(inflations)
       this(:,x) = nanmean(AllInflatedCorrs(:,:,inflations(x)),1)';
    end
    mtypem(:,:,f) = this;
end

mtypem(isnan(mtypem)) = 100;
    
thr = size(mtypem,1) + 5;
cmap = [];
cmap(1,:) = [237 248 251]./255;
cmap(2,:) = [140 150 198]./255;
cmap(3,:) = [110 1 107]./255;
[x,y]=meshgrid([1:3],[1:thr]);
cmap = interp2(x([1,round(thr/2),thr],:),y([1,round(thr/2),thr],:),cmap,x,y);
    
% p < 0.05
% h
% 
% % p < 0.000005
% lineProps.col = {cmap(ti,:)}; % CBCL colors
% lineProps.style = '--';
% mseb(1:size(mtypem,1),mtypem(:,1)',stypem(:,1)',lineProps,1);

figure; hold on 
for p = [1 7]
    for i = 1:length(inflations)
        if p == 1
            h = plot(1:size(mtypem,1),mtypem(:,i,p),'Color',cmap(inflations(i),:),'LineWidth',2,'LineStyle','--');
        else
            h = plot(1:size(mtypem,1),mtypem(:,i,p),'Color',cmap(inflations(i),:),'LineWidth',2);
        end
    end
end

xlim([1 size(mtypem,1)])
xticks([1:1:size(mtypem,1)]) 
xticklabels({'25'; '33' ; '50' ; '70' ; '100' ; '135' ; '200' ; '265' ; '375' ; '519' ; '725' ; '1000' ; '1430' ; '2000' ; '2800'; '4000' ; '10000' ; '20000' ; '30000'});
xtickangle(45)
xlabel('Number of Subjects','FontWeight','bold','FontSize',14)
ylabel('Type M Error','FontWeight','bold','FontSize',14)
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2)
ylim([0 100])

