

%% In Vivo Imaging Analysis Pipeline for Rogers et al., 2024
%
%
%Each part depends on previous parts and will not run without one another
%unless specified.
%
%In the following sections you will:
% 1. Load the data of all animals in a given experimental group from
% storeRogers2024Data.m
% 2. Extract the traces of longitudinally registered cells from all files
% and align to user-defined stimulus times or times of interest
% 3. Plot spatial coordinates of longitudinally registered cells
% 4. Extract and plot average traces in terms of % baseline of all neurons 
% averaged over all trials in each session
% 5. Identify cells that are upregulated or downregulated in response to
% stimuli in a given session according to their average trace
% 6. Measure the number of sessions cells tend to maintain stimulus
% responsiveness
% 7. Calculate freezing encoding in each neuron in each session and compare
% within animals across groups
% 8. Align freezing to trials
% 9. Load & plot TCA models, identifying session-dominant components and
% correlating to freezing
% 10. Calculate and plot normalized trial factors.
% 11. Calculate and plot normalized temporal factors.
% 12. Assess cumulative distribution of neuron weights
% 13. Identify Acq-Dominant, Ext1-Dominant, and Ext3-Dominant ensembles and
% their overlaps
% 14. Use a Fisher linear decoder to predict group based on average
% 15. Calculate changed activity (z-score relative to Acquisition) of each ensemble over time 
% 16. Reconstruct example neurons from TCA
% 17. Stimulus response properties of TCA ensemble neurons.
%
% To run non-shock mice, change for loops that are a=1:25 to a=26:28
%
% For aesthetic purposes, most data was exported from these sections in csv
% and txt files and imported in Prism, where most statistical tests were
% also done
%
%Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania




%% 1. Load the data of all animals in a given experimental group this will take a minute or so


%the following structure contains all the raw data necessary to run the following code.
%the structure is broken down by animals (structures).
%within each animal field are subfields for sessions (structure), 
%cell registration table across sessions (table), 
%cell spatial properties (table), and tca output(table)

%sessions has 5 subfields (6 if broken videos), one for each session
%(habituation, acquisition, extinction 1, extinction 2, extinction 3). in
%each session subfield is calcium, a table of traces; stimulusTimes, a
%vector of stimulus delivery times; and freezing, a table of freezing


load('/Users/Rogers1/Documents/rogers2024CalciumData.mat')
%%
% metadata of experimental design.
trialStructure = table2array(rogers2024.trialStructure);
start = trialStructure(1); % trial start. defined here ashow much time back from stimulus times to begin collecting data. (seconds)
fin = trialStructure(2);   % how much time after stimulus start to stop collecting data. (seconds)
dt = trialStructure(3);    % sampling rate of recording (Hz)
nSesh = trialStructure(4); %number of sessions
trialLength = start+fin+1; %length of trial (seconds)
nTrials = trialStructure(5); %number of total trials
nAnimals = length(rogers2024.animals); %number of animals is the number of substructures
videoOffset = trialStructure(6); %

%experimental groups
groupIDs{1}= rogers2024.groups(1).members; %rapid
groupIDs{2}= rogers2024.groups(2).members; %slow
groupIDs{3}= rogers2024.groups(3).members; %responders
groupIDs{4}= rogers2024.groups(4).members;  %nonresponders
groupNames = {rogers2024.groups(1).name,rogers2024.groups(2).name,rogers2024.groups(3).name,rogers2024.groups(4).name};
nGroups = length(groupIDs); %number of sessions

%define trial event times
bl = rogers2024.stimuli(1).times;  %seconds
stim = rogers2024.stimuli(2).times; %seconds
trace = rogers2024.stimuli(3).times; %seconds
shock = rogers2024.stimuli(4).times; %seconds

stimuliNames = {rogers2024.stimuli(1).name,rogers2024.stimuli(2).name,rogers2024.stimuli(3).name,rogers2024.stimuli(4).name};

%assign trials to sessions
hab = rogers2024.sessions(1).trials; %trials
acq = rogers2024.sessions(2).trials; %trials
ext1 = rogers2024.sessions(3).trials; %trials
ext2 = rogers2024.sessions(4).trials; %trials
ext3 = rogers2024.sessions(5).trials; %trials
sessions = {hab, acq, ext1, ext2, ext3};
sessionNames = {rogers2024.sessions(1).name,rogers2024.sessions(2).name,rogers2024.sessions(3).name,rogers2024.sessions(4).name,rogers2024.sessions(5).name};


% 2. Extract traces of longitudinal cells and align to stimuli

%CUSTOM FUNCTIONS REQUIRED:
%   crossSeshAuto
%   aligntraces

%outputs
nCells = zeros(nAnimals,1);  %initialize array to store number of cells per animal
lr = cell(nAnimals,nSesh);   %initialize cell to store original activity matrices downsampled to 15Hz, for freezing decoding later
tensors = cell(nAnimals,1);  %initialize cell to store same information as poolMat but arranged in tensors of t seconds X c cells X T trials. Export for use in tensortools python kit
spatial = cell(nAnimals,1);  %initialize cell to store spatial coordinates of longitudinally registered cells
poolMat = [];

for n=1:nAnimals
    
    %initialize temp storage
    reRegistered={};
    dA = {};
    timesA = {};
    
    %loop through number of recordings and store raw traces
    nRec= length(rogers2024.animals(n).sessions);
    for m=1:nRec
        dA{m,1} = rogers2024.animals(n).sessions(m).calcium;
    end
    
    %call temporal regristration and spatial registration tables
    L = table2array(rogers2024.animals(n).cellreg);
    locs = rogers2024.animals(n).cellprops;
        
   %longitudinally register and spatially locate cells across 5 sessions
    [longRegistered, coordinates]= crossSeshAuto(L,dA,locs,nRec);
    
   %store spatial data
   spatial{n,1} = coordinates;
   
   %repair broken recordings in 3 animals
    if n == 12
        reRegistered{1,1} = [longRegistered{1,1}];
        reRegistered{2,1} = [longRegistered{2,1}];
        reRegistered{3,1} = [longRegistered{3,1}; longRegistered{4,1}];
        reRegistered{4,1} = [longRegistered{5,1}];
        reRegistered{5,1} = [longRegistered{6,1}];
        longRegistered = reRegistered;
    elseif n == 13
        x=longRegistered{1,1};
        y=longRegistered{2,1};
        reRegistered{1,1} = [x; x(end-16.8*dt/2+1:end,:); y(1:16.8*dt/2,:); y];
        reRegistered{2,1} = [longRegistered{3,1}];
        reRegistered{3,1} = [longRegistered{4,1}];
        reRegistered{4,1} = [longRegistered{5,1}];
        reRegistered{5,1} = [longRegistered{6,1}];
        longRegistered = reRegistered;
    elseif n == 28
        reRegistered{1,1} = [longRegistered{1,1}; longRegistered{2,1}];
        reRegistered{2,1} = [longRegistered{3,1}];
        reRegistered{3,1} = [longRegistered{4,1}];
        reRegistered{4,1} = [longRegistered{5,1}];
        reRegistered{5,1} = [longRegistered{6,1}];
        longRegistered = reRegistered;
    end
    
    %loop through sessions, downsample registered traces to 15Hz and store.
    %also store stimulus delivery times
     for m=1:nSesh
        x=resample(longRegistered{m,1},3,4);
        lr{n,m} = x;
        
        timesA{m,1} = rogers2024.animals(n).sessions(m).stimulusTimes;
     end
    
    %store number of cells
    nCells(n,1) = size(longRegistered{1,1},2);
    
    %align traces to stimulus times and store 1Hz data at the predefined
    %trial windows
    alignedTraces = aligntraces(longRegistered,timesA,nSesh,dt,start,fin);
    
    %concatenate aligned traces
    allAligned = [alignedTraces{1,1}; alignedTraces{2,1}; alignedTraces{3,1}; alignedTraces{4,1}; alignedTraces{5,1}];
    
    %transform into t x c x T tensor and store
    ten = zeros(trialLength, nCells(n), nTrials);
    for t=1:nTrials
        ten(:,:,t) = allAligned(trialLength*(t-1)+1:trialLength*t,:);
    end
    tensors{n,1} = ten;
    
    poolMat= [poolMat allAligned];
end


%% 3.Plot spatial coordinates of longitudinally registered cells. Fig 2B 

%enter animal number of interest; dataset 1 animal 1 is the representative image
animal = 1; 
propM39{1,1} = readtable('mm39hab-props.csv'); propM39{2,1} = readtable('mm39acq-props.csv'); propM39{3,1} = readtable('mm39ext1-props.csv'); propM39{4,1} = readtable('mm39ext2-props.csv'); propM39{5,1} = readtable('mm39ext3-props.csv');
data = spatial{animal};

for m=1:nSesh
    dat = propM39{m,1};
    posX{1,m} = table2array(dat(:,6));
    posY{1,m} = table2array(dat(:,7));
    diameter{1,m} = table2array(dat(:,9));
end


%plot cells at x coordinate by y coordinate at their size in pixels, scaled up to be visible
colors = {[.5 .5 .5], [1 0 0], [0.9290 0.6940 0.1250], [0 1 0], [0 1 1]};
figure
for m=1:nSesh
    hold on
    scatter(posX{1,m},posY{1,m},diameter{1,m}.*(20-3*(m-1)),colors{m},'filled')
end
hold on
scatter(data(:,1),data(:,2),data(:,3).*8,[0.4940 0.1840 0.5560]) 
%% 4. for fig 2D - heatmaps of average traces

%CUSTOM FUNCTIONS
% percBaselineAvg

%initialize storage cell average session activity of each animal by session 
%avPB = cell(nAnimals, nSesh); 
avPBt = cell(nAnimals, nSesh); 
%loop through animals
for n=1:25%nAnimals %can look at specific groups by replacing "1:nAnimals" with saline, responders, nonresponders, or non-shock
     
    %call t x c x T data from tensors and convert to time x cells matrix
    data = tensors{n};
    %zscore to baseline
    data = (data-mean(data(1:10,:,:)))./std(data(1:10,:,:));
    
    %initialize temp storage
    avPB = zeros(trialLength,nCells(n),nSesh);
    
    %loop through sessions
    for m=1:nSesh
        
        trials = sessions{m};
        
        act = median(data(:,:,sessions{m}),3);
        
        %store median zscore in time x cells x session tensor
        avPB(:,:,m) = act;
        
    end
    %store tensor for each animal
    avPBt{n} = avPB;
end

%pool cells from all animals within sessions
poolAvgAct = cell(nSesh,nGroups);

%loop through animals to pool mean trial activity in a session
for g=1:4
    for n = rogers2024.groups(g).members
    
    %call mean % baseline tensor from animal
    data1 = avPBt{n};
    
    %loop through sessions
    for m=1:nSesh
        
        %call a specific session
        data = data1(:,:,m);
        
        %store all average traces in a session from every animal
        pool = [poolAvgAct{m,g}, data];
        poolAvgAct{m,g} = pool;
    end
    
    %store pooled traces by session
    end
end

%sort data to shock
acqAct = poolAvgAct{5,1};
[B,I] = sort(mean(acqAct(shock,:)));
c = palette('scheme',4);
%plot
titles = {'Hab','Acq','Ext1','Ext2','Ext3'};
figure
for g=1:4
    acqAct = poolAvgAct{2,g};
    [B,I] = sort(mean(acqAct(shock,:)));
    length(I)
    for m=1:nSesh
        subplot(4,5,m+5*(g-1))

        data2plot = poolAvgAct{m,g};
        h=heatmap(data2plot(:,I)','Colorlimits',[-2 2]);
        colormap(c)
        if g==1
        title(titles{m})
        end
        ylabel('Cells')
        if g==4
        xlabel('Time')
        end
        h.XDisplayLabels = repmat({''}, 1, size(h.ColorData, 2));
        h.YDisplayLabels = repmat({''}, 1, size(h.ColorData, 1));
        if m<5
        h.ColorbarVisible = 'off';
        end
        grid off
    end
end



%% 5. Identify cells that are upregulated or downregulated in response to stimuli in a given session according to their average trace for fig 2J-M

%CUSTOM FUNCTIONS
%permTest
tone = [11:35];
trace = [36:55];
shock = [56:60];
bl = {[1:10],[1:10],[1:10],[1:10]};
stimuli = {tone, trace, [tone trace], shock};
allCellsEver = cell(nAnimals, 8);
allCellsTrials = cell(nAnimals, 8);
numIt = 1000;

for a=1:nAnimals
    data = poolMat(:,sum(nCells(1:a-1))+1:sum(nCells(1:a)));
    
    disp(a)
    for m=1:nSesh
        trialSet = sessions{m};
        
        %make data tensors
        x = zeros(60,size(data,2),length(trialSet));
        for t=1:length(trialSet)
            x(:,:,t) = data(60*(trialSet(t)-1)+1:60*(trialSet(t)),:);
        end
        %zscore compared to baseline
        x = (x-mean(x(1:10,:,:)))./std(x(1:10,:,:));
        
        %median zscore
        data2 = squeeze(median(x,3));
        
        %for tone, trace, tone+trace, shock
        for k=1:4
            stim = stimuli{k};
            bL = bl{k};
            
            
            responses = isSig(data2,bL,stim);
            
            sigCells{a,m,k} = responses;
        end
    end
end

%calculate fractions of activated neurons
for a=1:nAnimals
    for m=1:nSesh
        for k=1:4
            fracs(a,m,k)=length(sigCells{a,m,k}{1})/nCells(a);
        end
    end
end

%plot
figure
for k=1:4
    subplot(2,2,k)
    for g = 1:nGroups
    errorbar(mean(fracs(rogers2024.groups(g).members,:,k)),std(fracs(rogers2024.groups(g).members,:,k))./sqrt(length(rogers2024.groups(g).members)))
    hold on
    end
end



%% 6. Measure the number of sessions cells tend to maintain stimulus responsiveness. For fig 2N-O

%CUSTOM FUNCTIONS
%permTest

%loop through animals to identify tone, trace, and shock responsive neurons
for a = 1:25%nAnimals
    %call average activity tensor 
    toTake = 0;
    if a>1
        toTake = sum(nCells(1:a-1));
    end
    %loop through stimuli after baseline
    for k = 1:length(rogers2024.stimuli)
            counter = zeros(5,nCells(a));
            for c=1:nCells(a)
                for m=1:8
                    counter(m,c) = ismember(c,sigCells{a,m,k}{1});
                    
                end
            end
            
            for m=1:8
                frac(a,m,k) = length(find(sum(counter)==m))./nCells(a);
            end
    end
         disp(a)
 end
    

%plot
for s=1:3
    labels{s} = strcat('Fraction  ',rogers2024.stimuli(s+1).name,'-responsive');
end

plotPersist(frac,groupIDs,groupNames,labels(1),labels)

%% 7. Measure the freezing encoding of single neurons. For fig 2P-Q, Supp Fig. 2B
%Representative image is animal g=10, session=3

%initialize storage matrix for the fraction of neurons and average freezing
%encoding
freezeCells = zeros(nSesh,nAnimals);
AUCsMean = zeros(nSesh,nAnimals);

delay = 30; %offset between miniscope and freezing recording (seconds)
dt = 15; %sampling rate of freezing (Hz)

%loop through animals
for g = 1:25
    
    %loop through sessions
    for h = 3%1:nSesh
        
        %initialize temporary storage counter and vector
        numCells = 0;
        AUCtot = [];
        
        data = lr{g,h};    % downsampled traces
        ac = zscore(data(delay*15:end,:)); %remove temporal offset of recordings, z-score over the session
        
        freezing = table2array(rogers2024.animals(g).sessions(h).freezing);      %loop through freezing
        freezing(freezing==100) = 1;  %turn into binary
        
        if length(freezing)>length(ac)  %make matrices the same length
            freezing = freezing(1:length(ac(:,1)));
        else
            ac = ac(1:length(freezing),:);
        end
        
        %initialize storage vectors for auROC and cell indices for auROC > .6
        AUCs = [];
        posCells = [];
        
        %loop through neuronsΩç
        for m=1:length(ac(1,:))
                act = ac(:,m);
                
                %create shuffle vector
                randN = randperm(round(length(freezing)));
                
                %pick random 50% of activity freezing values for training
                %set
                trainAct = act(randN(1:round(length(randN)/2)),:);
                trainFreeze = freezing(randN(1:round(length(randN)/2)),:);
                
                %pick random 50% of activity freezing values for test
                %set
                testAct = act(randN(round(length(randN)/2)+1:end),:);
                testFreeze = freezing(randN(round(length(randN)/2)+1:end),1);
                
                %train logistic regressor
                [B,DEV,STATS] = glmfit(trainAct,trainFreeze,'binomial','Link','logit');
                
                %test
                pred = glmval(B,testAct,'logit');
                
                %calculate auROC
                [X,Y,T,AUC] = perfcurve(testFreeze,pred,1, 'XCrit','FPR','YCrit','TPR');
                
                %record auROC
                AUCs = [AUCs AUC];
               
                if STATS.p(2)<0.0001
                    
                    %count cells with auROC > 0.6
                    numCells = numCells+1;
                    
                    %record their index
                    posCells = [posCells m];

                end
            
            freezeEncoding{g,h} = posCells;
        end
        
        %record average auROC of animal in session
        AUCsMean(h,g) = mean(AUCs);
        
        %record fraction of freezing encoding neurons
        freezeCells(h,g) = length(posCells)/nCells(g);
        ac(isnan(ac)) = 0;
        c = palette('scheme',2);
        
       
        [B,I] = sort(AUCs);
        cofI = I(ismember(I,posCells))';
        %[B,I] = sort(mean(zscore(ac(132+1141:194+1141,cofI))));
        
        %Uncomment for representative image
%         figure
%         subplot(211)
%         heatmap(freezing(1241:1440)')
%         grid off
%         subplot(212)
%         if nCells(g)<50
%             f=nCells(g)-1;
%         else
%             f=49;
%         end
%         heatmap(zscore(ac(1241:1440,cofI))','colorlimits',[-2 2])
%         colormap(c)
%         grid off

    end
    

  disp(g)
    
end

%plot
for a = 1:nAnimals
    blPop = freezeEncoding{a,1};
    for m=2:nSesh
        pop = freezeEncoding{a,m};
        pop(ismember(pop,blPop)) = [];
        fearEncoding{a,m} = pop;
        freezeFrac(m-1,a) = length(pop)/nCells(a);
    end
end
ylab1 = {'Mean freezing encoding'};
ylab2 = {'Fraction of Cells'};
ylim1 = [.5 .7];
ylim2 = [0 .7];

plot2metrics(AUCsMean,freezeFrac,sessionNames,groupNames,groupIDs,ylab1,ylab2,ylim1,ylim2)
%% 8. Extract freezing and align to stimuli
%calculate tensors from Williams et al. TCA code (https://github.com/ahwillia/tensortools)
% load tensors 

%get trial by trial freezing
freezingInTrials = {};
dt = 15;
timesF = cell(5,1);
%get standard stimulus times
for m=1:nSesh
    timesF{m,1} = rogers2024.animals(24).sessions(m).stimulusTimes-videoOffset;
end

%create second by second % freezing
for n=1:25%nAnimals
    fr = cell(nSesh,1);
    for m=1:nSesh
        %read each animal's binary freezing data from each session
         fz = table2array(rogers2024.animals(n).sessions(m).freezing);
         fr{m} = fz;
         
         %downsample to 1 second by thresholding (if freezing is greater
         %than 50% in that second
         fzH = zeros(round(length(fz)/dt)-1,1);
        for t=1:round(length(fz)/dt)-1
            
            val = mean(fz((t-1)*dt+1:t*dt));
            if val>=50
                fzH(t) = 1;
            else
                fzH(t) = 0;
            end
            
        end
        
        %get that animal's stimulus times
        stimTimes = rogers2024.animals(n).sessions(m).stimulusTimes;
        offset = stimTimes(1)-videoOffset;
        fzHs = [];
        for s =1:length(stimTimes)
            fzHs = [fzHs; fzH(stimTimes(s)-offset-10:stimTimes(s)-offset+49)];
        end
        fzHz{n,m} = fzHs;
    end

    %fix broken video
    if n==5
        timesf = timesF;
        timesf{3,1} = [120; 120+round(2844/15); 234+round(2844/15); 331+round(2844/15); 441+round(2844/15); 543+round(2844/15)];
    else
        timesf = timesF;
    end
    
    %align freezing to stimulus times
    fIT = aligntraces(fr,timesf,nSesh,dt,start,fin);
    freezingInTrials = [freezingInTrials fIT];
end


%take freezing 10 seconds before each tone to 50 seconds after on each
%trial and concatenate
for n=1:25%nAnimals
    a=0;
    for m=1:nSesh
        fre = freezingInTrials{m,n};

        for t = 1:length(sessions{m})
           a=a+1;
           freTot(a,n) = mean(fre(trialLength*(t-1)+1:trialLength*t));
        end

    end
end

%% 9. Load & plot TCA models, identifying session-dominant components and correlating to freezing for fig. 3A,D-F

nDims = 5;

coefs = zeros(nAnimals,nDims);
%to plot representative TCAs, uncomment the figure file in the for-loop -
%the representative image in the paper is from animal 12

for a=1:25
    
    fr = freTot(:,a);
    if a<22 
    tt = table2array(rogers2024.animals(a).tca);
    timeFactor = tt(2:61,:);
    neuronFactor = tt(62:end-34,:);
    trialFactor = tt(end-33:end,:);
    elseif ismember(a,22:25)
        tt = table2array(rogers2024.animals(a).tca(:,2:end));
        timeFactor = tt(1:60,:);
        neuronFactor = tt(61:end-34,:);
        trialFactor = tt(end-33:end,:);
    end
    
    %store neuron factors - used in Fig. 3G, 4-5
    neuFac{a,1} = neuronFactor;
    timFac{a,1} = timeFactor;
    triFac{a,1} = trialFactor;
    
    %Identify session-dominant factors - used in Fig. 3D-F, 4, 5
    for n=1:nSesh
        facLoads = mean(trialFactor(sessions{n},:));
        domFacs(a,n) = find(facLoads==max(facLoads));
    end
    
    %calculate strength of each factor in each session - used in Fig. D,E
    for n=1:nSesh
        facLoads = mean(trialFactor(sessions{n},:));
        
        %for m=1:nDims
            strength(a,n) = facLoads(domFacs(a,n))/sum(facLoads(domFacs(a,:)));
           
        %end
    end
    
    %calculate correlation of trial factor loadings with trial-by-trial
    %freezing over extinction with multiple linear regression. save stats
    [b, bint, r, rint,stats] = regress(fr, [trialFactor ones(nTrials,1)]);
    coefs(a,:) = b(2:6);
    R2(a,:) = stats(1); 
    fs(a,:) = stats(2); 
    ps(a,:) = stats(3); 

    %uncomment to plot representative figure
    %plotTCAmodel(timeFactor,neuronFactor,trialFactor,fr)
end
%%
r=zeros(5);
p=zeros(5);
for m=1:5
    for n=1:5
        [r(n,m),p(n,m)] = corr(mean(freTot(sessions{n},1:25))',strength(1:25,m));
    end
end
r(p>0.05) = 0;
figure
c = palette('scheme',4);
s=sign(r);
heatmap((r))
colormap(c)

%% 10. Plot normalized trial factors. for Supp Fig. 4A-E

%initialize trial weight storage vectors
normTrialWeights = zeros(nTrials, nAnimals, nSesh);

for n=1:25%nAnimals
    
    
    %call trial weights from a given animal
    tWeights = triFac{n};
    
    %store normalized trial weights
    for m=1:nSesh
        normTrialWeights(:,n,m) = tWeights(:,domFacs(n,m))./max(tWeights(:,domFacs(n,m)));
       
        facLoads = mean(normTrialWeights(sessions{m},n,:));
        
        for m=1:nDims
            strength(n,m) = facLoads(domFacs(n,m))/sum(facLoads(domFacs(n,:)));
           
        end
   
    end
    
end

%plot
ys = [0 1];% {'Hab-Dom','Acq-Dom','Ext1-Dom','Ext2-Dom','Ext3-Dom'};
colors = {[.5 .5 .5],'k','c','r'};
figure
for n=1:nSesh
subplot(5,1,n)
    for m=1:nGroups
        data=normTrialWeights(:,rogers2024.groups(m).members,n);
        plot(1:length(data),mean(data,2),'color',colors{m},'Linewidth',1)
        dat=data;
        hold on
        y_upper = mean(dat,2)+std(dat,[],2)./sqrt(size(dat,2));
        y_lower = mean(dat,2)-std(dat,[],2)./sqrt(size(dat,2));
        fill([1:size(dat, 1), fliplr(1:size(dat, 1))], [y_upper; flipud(y_lower)], colors{m}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on
    end
    ylim([ys])
    %ylabel('Normalized trial weights')
    %xlabel('Time')
    %title(strcat(sessionNames{n},'-Dominant Component'))
    %if n==5
    %    legend(groupNames)
    %end
end


%% 11. Plot normalized temporal factors. for Supp Fig. 3F-J
%initialize trial weight storage vectors
normTrialWeights = zeros(nTrials, nAnimals, nSesh);

for n=1:25%nAnimals
    
    %call trial weights from a given animal
    tWeights = timFac{n};
    
    %store normalized trial weights
    for m=1:nSesh
        normTempWeights(:,n,m) = tWeights(:,domFacs(n,m))./max(tWeights(:,domFacs(n,m)));
    end
    
end

%plot
ys = [0 1];% {'Hab-Dom','Acq-Dom','Ext1-Dom','Ext2-Dom','Ext3-Dom'};
colors = {[.5 .5 .5],'k','c','r'};
figure
for n=1:nSesh
subplot(1,5,n)
    for m=1:nGroups
        data=normTempWeights(:,rogers2024.groups(m).members,n);
        plot(1:length(data),mean(data,2),'color',colors{m},'Linewidth',1)
        hold on
        dat=data;
        y_upper = mean(dat,2)+std(dat,[],2)./sqrt(size(dat,2));
        y_lower = mean(dat,2)-std(dat,[],2)./sqrt(size(dat,2));
        fill([1:size(dat, 1), fliplr(1:size(dat, 1))], [y_upper; flipud(y_lower)], colors{m}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on
    end
    ylim([ys])
%     ylabel('Normalized trial weights')
%     xlabel('Time')
%     title(strcat(rogers2024.sessions(n).name,'-Dominant Component'))
%     if n==5
%         legend(groupNames{1},groupNames{2},groupNames{3},groupNames{4})
%     end
end

%% 12. Assess cumulative distribution of neuron weights. for Fig. 4G,H

% set threshold range
threshold = linspace(0,10,11);

%define storage vectors for fraction of cells included in Acq-, Ext1-, and
%Ext3-dominant ensembles
nComps = 3;
comps = [2 3 5];
cellsInComp = zeros(nAnimals, length(threshold),nComps);

for t = 1:length(threshold)
    aCells = [];
    e1Cells = [];
    e3Cells = [];
    
    for n = 1:25%nAnimals
        
        %call neuron factor weights from animal
        weights = neuFac{n};
        
        for c = 1:nComps
        %calculate fractions of neurons with w > threshold and store
            cellsInComp(n,t,c) = length(find(weights(:,domFacs(n,comps(c)))>threshold(t)))/nCells(n);
        end
        
        %record ensembles at each threshold
        if n == 1
            c = 0;
        else
            c = sum(nCells(1:n-1));
        end
    
        aCells = [aCells; find(weights(:,domFacs(n,2))>threshold(t))+c];
        e1Cells = [e1Cells; find(weights(:,domFacs(n,3))>threshold(t))+c];
        e3Cells = [e3Cells; find(weights(:,domFacs(n,5))>threshold(t))+c];

    end
    
    %collect all the cells in the Acq, Ext1, and Ext3 ensembles at the 
    %predetermined varying thresholds
    collection{t,1} = aCells;
    collection{t,2} = e1Cells;
    collection{t,3} = e3Cells;
end

%plot
colors = {'r','y','c'};
figure
for c=1:nComps
    data = cellsInComp(:,:,c);
    errorbar(threshold,mean(data),std(data)/sqrt(nAnimals-1),colors{c},'Linewidth',3)
    hold on
end
hold on
xline(1)
xlabel('Threshold weight')
ylabel('Fraction of neurons')
legend({'Acq-Dominant','Ext1-Dominant','Ext3-Dominant'})
title('Choosing dominant neurons')


%% 13. Identify Acq-Dominant, Ext1-Dominant, and Ext3-Dominant ensembles and their overlaps for Fig. 5B, 6A

%set threshTest to 0 if reproducing main figs. set threshTest to 1 if
%reproducing Supp Fig 6. note that null thresholds were not determined for
%non-shock mice
threshTest = 0;

%define groups
responders = rogers2024.groups(3).members;
nonresponders = rogers2024.groups(4).members;
rapid = rogers2024.groups(1).members;
slow = rogers2024.groups(2).members;
nonshock = rogers2024.groups(5).members;

if threshTest == 1
    threshes = table2array(readtable('fMeans.csv'));
    threshold = [threshes([[15:21] [1:14]]); [1.55; 1.28; 1.7; 0.975]];
    nAnimals = 21;
end

%create storage vectors for all Acq-, Ext1-, and Ext3-dominant cells
acqDomCells = [];
ext1DomCells = [];
ext3DomCells = [];

%create storage vectors for their overlaps
acqOnly = [];
ext1Only = [];
ext3Only = [];
acqExt1 = [];
acqExt3 = [];
ext1Ext3 = [];
acqExt1Ext3 = [];

%allocate vectors to assign cell indices to for each group
cR = [];
cNR = [];
cSS = [];
cSR = [];
cNS = [];
for a = 1:25%1:21%nAnimals
    if threshTest == 1
        th = threshold(a);
    else
        th=1;
    end
    %number to add to convert within-animal cell index to pooled cell index
    if a == 1
        c = 0;
    else
        c = sum(nCells(1:a-1));
    end
    
    %load cell indices into their group
    if ismember(a,responders)
        cR = [cR [1+c:nCells(a)+c]];
    elseif ismember(a,nonresponders)
        cNR = [cNR [1+c:nCells(a)+c]];
    elseif ismember(a,rapid)
        cSR = [cSR [1+c:nCells(a)+c]];
    elseif ismember(a,slow)
        cSS = [cSS [1+c:nCells(a)+c]];
    elseif ismember(a,nonshock)
        cNS = [cNS [1+c:nCells(a)+c]];
    end
    
    %call neuron factor weights
    weights = neuFac{a};
    
    %dominant cells are defined by those weights > 1 in the component
    %dominating that session


    aCells = find(weights(:,domFacs(a,2))>th);
    e1Cells = find(weights(:,domFacs(a,3))>th);
    e3Cells = find(weights(:,domFacs(a,5))>th);
    
    weightMat{a,1} = weights(aCells,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,2} = weights(e1Cells,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,3} = weights(e3Cells,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    
    overlap = aCells(ismember(aCells,e1Cells));
    
    %cells dominating all three sessions
    ae1e3 = overlap(ismember(overlap,e3Cells));
    
    %cells dominating all acq, ext1
    ae1 = overlap(~ismember(overlap,e3Cells));
    
    %cells dominating all acq, ext3
    overlap = aCells(ismember(aCells,e3Cells));
    
    ae3 = overlap(~ismember(overlap,e1Cells));
    
    %cells dominating all ext1, ext3
    overlap = e1Cells(ismember(e1Cells,e3Cells));
    mat(a,1) = length(overlap)./nCells(a);
    e1e3 = overlap(~ismember(overlap,aCells));
    
    %cells dominating only acq
    aOnly = aCells(~ismember(aCells,ae1e3));
    aOnly = aOnly(~ismember(aOnly,ae1));
    aOnly = aOnly(~ismember(aOnly,ae3));
    
    %cells dominating only ext1
    e1Only = e1Cells(~ismember(e1Cells,ae1e3));
    e1Only = e1Only(~ismember(e1Only,ae1));
    e1Only = e1Only(~ismember(e1Only,e1e3));
    
    %cells dominating only ext3
    e3Only = e3Cells(~ismember(e3Cells,ae1e3));
    e3Only = e3Only(~ismember(e3Only,e1e3));
     e3Only = e3Only(~ismember(e3Only,ae3));
     
    %store fraction of overlap for individual animals
    ensOverlap(a,1,1) = length(aOnly)/length(aCells);
    ensOverlap(a,2,1) = length(ae1)/length(aCells);
    ensOverlap(a,3,1) = length(ae3)/length(aCells);
    ensOverlap(a,4,1) = length(ae1e3)/length(aCells);
    
    ensOverlap(a,1,2) = length(e1Only)/length(e1Cells);
    ensOverlap(a,2,2) = length(ae1)/length(e1Cells);
    ensOverlap(a,3,2) = length(e1e3)/length(e1Cells);
    ensOverlap(a,4,2) = length(ae1e3)/length(e1Cells);
    
    ensOverlap(a,1,3) = length(e3Only)/length(e3Cells);
    ensOverlap(a,2,3) = length(ae3)/length(e3Cells);
    ensOverlap(a,3,3) = length(e1e3)/length(e3Cells);
    ensOverlap(a,4,3) = length(ae1e3)/length(e3Cells);
    
    %assign neurons indices for activity pooled across animals and store in
    %vectors
    weightMat{a,4} = weights(aOnly,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,5} = weights(e1Only,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,6} = weights(e3Only,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,7} = weights(ae1,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,8} = weights(ae3,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,9} = weights(e1e3,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,10} = weights(ae1e3,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    
    ensAns(a,1).acqcells = aCells+c;
    ensAns(a,1).ext1cells = e1Cells+c;
    ensAns(a,1).ext3cells = e3Cells+c;
    ensAns(a,1).a = aOnly+c;
    ensAns(a,1).e1 = e1Only+c;
    ensAns(a,1).e3 = e3Only+c;
    ensAns(a,1).ae1 = ae1+c;
    ensAns(a,1).ae3 = ae3+c;
    ensAns(a,1).e1e3 = e1e3+c;
    ensAns(a,1).ae1e3 = ae1e3+c;
    
    acqDomCells = [acqDomCells; aCells + c];
    ext1DomCells = [ext1DomCells; e1Cells + c];
    ext3DomCells = [ext3DomCells; e3Cells + c];
    
    acqOnly = [acqOnly; aOnly+c];
    ext1Only = [ext1Only; e1Only+c];
    ext3Only = [ext3Only; e3Only+c];
    acqExt1 = [acqExt1; ae1+c];
    acqExt3 = [acqExt3; ae3+c];
    ext1Ext3 = [ext1Ext3; e1e3+c];
    acqExt1Ext3 = [acqExt1Ext3; ae1e3+c];
end
% 
% 
labels = {'Overlap with Acq-Dom Neurons','Overlap with Ext1-Dom Neurons','Overlap with Ext3-Dom Neurons'};
labels2{1} = ['Acq Only', 'Acq/Ext1', 'Acq/Ext3', 'Acq/Ext1/Ext3'];
labels2{2} = ['Ext1 Only', 'Acq/Ext1', 'Ext1/Ext3', 'Acq/Ext1/Ext3'];
labels2{3} = ['Ext3 Only', 'Acq/Ext3', 'Ext1/Ext3', 'Acq/Ext1/Ext3'];

nComps=3;
figure
for n = 1:nComps
    for m = 1:nGroups
        subplot(3,4,m+4*(n-1))
        barWithError([1:4],ensOverlap(rogers2024.groups(m).members,:,n),1)
        ylabel('Fraction of Cells')
        xticks([1:4])
        xticklabels(labels2{n})
        ylabel(labels{n})
        ylim([0 .8])
        if n==1
        title(rogers2024.groups(m).name)
        end
    end
end

%% 14. Use a Fisher linear decoder to predict group based on average activity (z-score relative to Acquisition) of each ensemble over time for fig. 6B

%TO AVOID HAVING TO RE-RUN PREVIOUS SCRIPTS, LOAD THESE FILES WITH
%ENSEMBLES' CELL INDICES
%
% IF RUNNING THRESHTEST SKIP THIS SECTION
% 
% poolMat = [];
% 
% for n=1:nAnimals
%     poolMat = [poolMat reshape(tensors{n},60*34,nCells(n))];
% end

%define time windows
timeHab = [1:trialLength*max(hab)];
timeAcq = [trialLength*max(hab)+1:trialLength*max(acq)];
timeExt1 = [trialLength*max(acq)+1:trialLength*max(ext1)];
timeExt2 = [trialLength*max(ext1)+1:trialLength*max(ext2)];
timeExt3 = [trialLength*max(ext2)+1:trialLength*max(ext3)];

%normalize activity to acquisition
poolMat = (poolMat-mean(poolMat(481:960,:)))./std(poolMat(481:960,:));


%define ensembles to loop through
ens = fieldnames(ensAns(1));



%create matrices of average ensembles acts over time in the session of
%interest
sess = 'Ex1'; %set to Acq, Ex1, or Ex3

if sess == 'Ex1'
    times = timeExt1;
    s=3;
elseif sess == 'Ex3'
    times = timeExt3;
    s=5;
elseif sess == 'Ex2'
    times = timeExt2;
    s=4;
elseif sess == 'Hab'
    times = timeHab;
    s=1;
elseif sess == 'Acq'
    times = timeAcq;
    s=2;
end
timeCell= {timeHab, timeAcq, timeExt1, timeExt2, timeExt3};
%define binary class vector
classes = [zeros(length(times),1); ones(length(times),1)];

ens = fieldnames(ensAns(1));
numEnsembles = length(ens);

sepActs = cell(25,10,2,5);
set1 = zeros(10,360);
set2 = zeros(25,10,360);
set3 = zeros(25,10,360);
for a=1:25
    if a==1
       c=0;
    else
       c=sum(nCells(1:a-1));
    end
    for m=[3 5]
        for l = 1:numEnsembles+1
            if l<11
            pop = ensAns(a).(ens{l});
            
                
            x=0;
            psave = [];
            if ismember(a,rapid)
                if ismember(l,[3,6])
                for p=1:length(pop)
                if mean(poolMat(timeExt3,pop(p)))>8
                    psave = [psave p];
                    disp(psave)
                end
                end
                end
                
                if ~isempty(psave)
                if ~isempty(pop)
                    pop(psave) = [];
                end
                end
            end
            
            
            %pop(ismember(pop,freezeEncoding{a,1}+c)) = [];
            idx=find(fzHz{a,m}==1);
            sepActs{a,l,1,m} = mean(poolMat(timeCell{m}(idx),pop),2);
            idx=find(fzHz{a,m}==0);
            sepActs{a,l,2,m} = mean(poolMat(timeCell{m}(idx),pop),2);
            sepActs{a,l,3,m} = mean(poolMat(timeCell{m},pop),2);
            lens(a,m,l,1) =  size(sepActs{a,l,1,m},1);
            lens(a,m,l,2) =  size(sepActs{a,l,2,m},1);
            lens(a,m,l,3) =  size(sepActs{a,l,3,m},1);
            else
                for ll = 1:10
                    set1(ll,:) = sepActs{a,ll,3,m};
                end
                sepActs{a,l,3,m} = set1;
                lens(a,m,l,1) =  size(sepActs{a,l,1,m},1);
                lens(a,m,l,2) =  size(sepActs{a,l,2,m},1);
                lens(a,m,l,3) =  size(sepActs{a,l,3,m},1);
            end
        end
    end
end

%number of iterations of decoder
numIt = 100;

%set of comparisons - res vs. nonres, res vs. sal, nonres vs. sal
comparisons = [0:2];

%store 100 iterations of the model (rows), for each ensemble plus a round
%where every ensemble is a predictor (columns), for each comparison
%(columns')
%accuracy = zeros(numIt,numEnsembles+1,2);
%accuracyShuff = zeros(numIt,numEnsembles+1,2,3,2);
sOI = [3 5];
for f=3
    for s=1:2
    for l = 1:10
        g1 = rogers2024.groups(3).members;
        g2 = rogers2024.groups(4).members;
        g3 = rogers2024.groups(1).members;
  
        
        if l==11
            ensActs1 = [];
            mL = min(lens(g1,sOI(s),l,f));
            for a=1:length(g1)
                set = sepActs{g1(a),l,f,sOI(s)};
                ensActs1 = [ensActs1 set];
            end
            
            ensActs2 = [];
            mL = min(lens(g2,sOI(s),l,f));
            for a=1:length(g2)
                set = sepActs{g2(a),l,f,sOI(s)};
                ensActs2 = [ensActs2 set];
            end
            
            ensActs3 = [];
            mL = min(lens(g3,sOI(s),l,f));
            for a=1:length(g3)
                set = sepActs{g3(a),l,f,sOI(s)};
                ensActs3 = [ensActs3 set];
            end
            
            mL = min([length(ensActs1) length(ensActs2) length(ensActs3)]);
  
        
            mL = min([length(ensActs1) length(ensActs2) length(ensActs3)]);
  
        
        data = [ensActs1(:,(length(ensActs1)-mL)/2+1:end-(length(ensActs1)-mL)/2)'; ensActs2(:,(length(ensActs2)-mL)/2+1:end-(length(ensActs2)-mL)/2)'; ensActs3(:,(length(ensActs3)-mL)/2+1:end-(length(ensActs3)-mL)/2)']; %compare responders to nonresponders
        data(isnan(data))=0;
        
        else
        
        ensActs1 = [];
        mL = min(lens(g1,sOI(s),l,f));
        for a=1:length(g1)
            set = sepActs{g1(a),l,f,sOI(s)};
            ensActs1 = [ensActs1; set((length(set)-mL)/2+1:end-(length(set)-mL)/2)];
        end
        %ensActs1 = mean(ensActs1)
        ensActs2 = [];
        mL = min(lens(g2,sOI(s),l,f));
        for a=1:length(g2)
            set = sepActs{g2(a),l,f,sOI(s)};
            ensActs2 = [ensActs2; set((length(set)-mL)/2+1:end-(length(set)-mL)/2)];
        end
        
        ensActs3 = [];
        mL = min(lens(g3,sOI(s),l,f));
        for a=1:length(g3)
            set = sepActs{g3(a),l,f,sOI(s)};
            ensActs3 = [ensActs3; set((length(set)-mL)/2+1:end-(length(set)-mL)/2)];
        end
        
        mL = min([length(ensActs1) length(ensActs2) length(ensActs3)]);
  
        
        data = [ensActs1((length(ensActs1)-mL)/2+1:end-(length(ensActs1)-mL)/2); ensActs2((length(ensActs2)-mL)/2+1:end-(length(ensActs2)-mL)/2); ensActs3((length(ensActs3)-mL)/2+1:end-(length(ensActs3)-mL)/2)]; %compare responders to nonresponders
        data(isnan(data))=0;
        
        end
        
        
        lD = length(data)/3;
        for c=1:3
            if c==1
                dat = data(1:2*lD);
                
            elseif c==2
                dat = [data(1:lD); data(1+2*lD:end)];
            else
                dat = data(1+lD:end);
            end
            
            for n=1:100
                classes = [zeros(mL,1); ones(mL,1)];
                [X,Y,x,y] = split_data(dat,classes,.5); %split_data is a custom function, described below that randomly selects test_size % of your data to test on and 1 - test_size to train on
                model = fitcdiscr(x,y); %fit discr fits a linear model to your training data & classes

                predictedLabels  = predict(model, X); %predict applies your model to your test data to generate class predictions

                %count true labels in your test set to normalize confusion
                %matrix
                normVec = zeros(length(unique(y)),1);
                for k=1:length(unique(y))
                    normVec(k,1) = sum(Y==k-1);
                end

                %different sessions have different numbers of unique behaviors 
                accuracy(n,l,f,c,1,s) = mean(predictedLabels == Y);

                shuffle = randperm(length(Y));
                Y = Y(shuffle);
                
                predictedLabels  = predict(model, X); %predict applies your model to your test data to generate class predictions
                accuracy(n,l,f,c,2,s) = mean(predictedLabels == Y); %the accuracy of your model is the number of instances in which your predicted classes matched your a priori classes divided by total number of predictions

            end
        end
    end
    end
end

titles = {'Acq-Dom','Ext1-Dom','Ext3-Dom','Acq Only', 'Ext1 Only', 'Ext3 Only', 'Acq/Ext1', 'Acq/Ext3', 'Ext1/Ext3', 'Acq/Ext1/Ext3','All ensembles'};

tits = {'Decoding between responders and nonresponders during Ext3','Decoding between responders and rapid mice during Ext3','Decoding between nonresponders and rapid mice during Ext3'};
figure; 
for c=1:3
    subplot(3,1,c)
barWithError([1:10],squeeze(accuracy(:,1:10,2,c,1,2)),0.5)
hold on
barWithError([1:10],squeeze(accuracy(:,1:10,1,c,1,2)),0.5)
hold on
barWithError([1:10],squeeze(accuracy(:,1:10,3,c,1,2)),0.5)
hold on
ylim([.5,1])
legend({'Motion','','Freezing'})
xticklabels(titles(1:10))
ylabel('Accuracy')
title(tits{c})
xtickangle(45)
end

tits = {'Decoding between responders and nonresponders during Ext1','Decoding between responders and rapid mice during Ext1','Decoding between nonresponders and rapid mice during Ext1'};
figure; 
for c=1:3
    subplot(3,1,c)
barWithError([1:10],squeeze(accuracy(:,1:10,2,c,1,1)),0.5)
hold on
barWithError([1:10],squeeze(accuracy(:,1:10,1,c,1,1)),0.5)
hold on
barWithError([1:10],squeeze(accuracy(:,1:10,3,c,1,1)),0.5)
hold on
ylim([.5,1])
legend({'Motion','','Freezing'})
xticklabels(titles(1:10))
ylabel('Accuracy')
title(tits{c})
xtickangle(45)
end

%% 15. Extract average activity (z-score relative to Acquisition) of each cell in the ensemble during Extinction 1, Extinction 3. for fig. 5D-F, fig. 6C-I, fig. 7B-D

for d=1:3
for n = 1:25
    if n==1
       c=0;
    else
       c=sum(nCells(1:n-1));
    end
    %loop through all ensembles
    for l=1:10

        cellPop = ensAns(n).(ens{l});
        %cellPop(ismember(cellPop,unique([freezeEncoding{n,1}])+c)) = [];
        %calculate average activity of cells during Extinction 1 and
        %Extinction 3
        %cellPop(ismember(cellPop,freezeEncoding{n,3}+c))= [];
        %(fzHz{n,5}==0)
        %(fzHz{n,3}==0)
         if d==1
             idx1 = (fzHz{n,3}==0);
             idx2 = (fzHz{n,5}==0);
         elseif d==2
             idx1 = (fzHz{n,3}==1);
             idx2 = (fzHz{n,5}==1);
         else
             idx1 = [1:360];
             idx2 = [1:360];
         end
        x =  [mean(poolMat(timeExt1(idx1),cellPop))' mean(poolMat(timeExt3(idx2),cellPop))']; % fzHz{n,5}==0)
        
        
        %store activity matrices
        changedActivity2{l,n,d} = x;
       
    end
    
end
%pool cells

for n=1:4
    group = rogers2024.groups(n).members;
    for l=1:10
        pop=[];
        pgroup=[];
        for g=1:length(group)
            an = group(g);
       
            if l==1
                pop=[pop; changedActivity2{4,an,d}; changedActivity2{7,an,d}; changedActivity2{8,an,d}; changedActivity2{10,an,d}];
            elseif l==2
                pop=[pop; changedActivity2{5,an,d}; changedActivity2{7,an,d}; changedActivity2{9,an,d}; changedActivity2{10,an,d}];
            elseif l==3
                  pop=[pop; changedActivity2{6,an,d}; changedActivity2{5,an,d}; changedActivity2{9,an,d}; changedActivity2{10,an,d}]; %
            else
                pop = [pop; changedActivity2{l,an,d}];
            end
            %pop = [pop; changedActivity2{l,an}];
            pgroup = [pgroup; ensAns(an).(ens{l})];
        end
        
    
    popGroup{l,n} = pgroup;
    pop(isnan(pop)) = 0;
    changedActivity{l,n,d} = pop;    
 
    end

end
sOI=[1:10];
titles = {'Acq','Ext1','Ext3','Acq-Only','Ext1-Only','Ext3-Only','Acq/Ext1','Acq/Ext3','Ext1/Ext3','Acq/Ext1/Ext3'};
col = {[.5 .5 .5],'k', 'c','r'};%,
grOI = [1 2 3 4];
figure
for l=1:10
subplot(1,10,l)
for g=1:4
    if size(changedActivity{sOI(l),grOI(g),d},1)<2
        continue
    else
        errorbar([1:2],mean(changedActivity{sOI(l),grOI(g)}),std(changedActivity{sOI(l),grOI(g)})./(sqrt(length(changedActivity{sOI(l),grOI(g)}))-1),'Linewidth',3,'Color',col{g})
    end
hold on
end
xlim([0 3])
xticks([1:2])
xticklabels({'Ext1','Ext3'})
title(titles{sOI(l)})
ylabel('Change in Activity (zscore)')
end
end

for n=1:25
    pop=ensAns(n).acqcells;
    %pop2 = sigCells{n,3,4}{1};
    mA(n) = mean([changedActivity2{4,n}(:,1); changedActivity2{7,n}(:,1); changedActivity2{8,n}(:,1); changedActivity2{10,n}(:,1)]);
    %pop=ensAns(n).ext3cells;
    pop=mean([changedActivity2{6,n}(:,2); changedActivity2{8,n}(:,2); changedActivity2{9,n}(:,2); changedActivity2{10,n}(:,2)]);
    if ismember(n,rogers2024.groups(1).members)
        pop(pop>8) = [];
    end
        
    mE(n) = mean(pop);
end
mA(isnan(mA)) = 0;
mE(isnan(mE)) = 0;



%% 16. TCA reconstructions using representative tca neurons from animal 18
an=18;
weights = neuFac{an};
pop = ensAns(an).a-sum(nCells(1:an-1));
fac = domFacs(an,2);
maxW = max(weights(pop,fac));

neu(1) = find(weights(:,fac)==maxW);

pop = ensAns(an).e1-sum(nCells(1:an-1));
fac = domFacs(an,3);
maxW = max(weights(pop,fac));

neu(2) = find(weights(:,fac)==maxW);

pop = ensAns(an).e3-sum(nCells(1:an-1));
fac = domFacs(an,5);
maxW = max(weights(pop,fac));

neu(3) = find(weights(:,fac)==maxW);

pop = ensAns(an).ae1-sum(nCells(1:an-1));
fac = domFacs(an,2);
maxW = max(weights(pop,fac));

neu(4) = find(weights(:,fac)==maxW);

pop = ensAns(an).ae3-sum(nCells(1:an-1));
fac = domFacs(an,2);
maxW = max(weights(pop,fac));

neu(5) = find(weights(:,fac)==maxW);

pop = ensAns(an).e1e3-sum(nCells(1:an-1));
fac = domFacs(an,3);
maxW = max(weights(pop,fac));

neu(6) = find(weights(:,fac)==maxW);

pop = ensAns(an).ae1e3-sum(nCells(1:an-1));
fac = domFacs(an,2);
maxW = max(weights(pop,fac));

neu(7) = find(weights(:,fac)==maxW);


dat = poolMat(timeAcq,sum(nCells(1:an-1))+1:sum(nCells(1:an)));
datA = zeros(trialLength,nCells(an),8);
for t=1:8
    datA(:,:,t) = dat(1+trialLength*(t-1):trialLength*t,:);
end
datA = zscore(datA);


dat = poolMat(timeExt1,sum(nCells(1:an-1))+1:sum(nCells(1:an)));
datE1 = zeros(trialLength,nCells(an),6);
for t=1:6
    datE1(:,:,t) = dat(1+trialLength*(t-1):trialLength*t,:);
end
datE1 = zscore(datE1);

dat = poolMat(timeExt3,sum(nCells(1:an-1))+1:sum(nCells(1:an)));
datE3 = zeros(trialLength,nCells(an),6);
for t=1:6
    datE3(:,:,t) = dat(1+trialLength*(t-1):trialLength*t,:);
end
datE3 = zscore(datE3);



for n=1:7
    wMat(n,1) = weights(neu(n),domFacs(an,2));
    wMat(n,2) = weights(neu(n),domFacs(an,3));
    wMat(n,3) = weights(neu(n),domFacs(an,5));
end


titles={'AcqOnly','Ext1Only','Ext3Only','Acq/Ext1','Acq/Ext3','Ext1/Ext3','Acq/Ext1/Ext3'};
sOI = [2 3 5];
for m=1:3
tim = sessions{sOI(m)};
figure
for n=1:7
fH = timFac{an}(:,domFacs(an,1)).*weights(neu(n),domFacs(an,1)).*mean(triFac{an}(tim,domFacs(an,1)));
fA = timFac{an}(:,domFacs(an,2)).*weights(neu(n),domFacs(an,2)).*mean(triFac{an}(tim,domFacs(an,2)));
fE1 = timFac{an}(:,domFacs(an,3)).*weights(neu(n),domFacs(an,3)).*mean(triFac{an}(tim,domFacs(an,3)));
fE2 = timFac{an}(:,domFacs(an,4)).*weights(neu(n),domFacs(an,4)).*mean(triFac{an}(tim,domFacs(an,4)));
fE3 = timFac{an}(:,domFacs(an,5)).*weights(neu(n),domFacs(an,5)).*mean(triFac{an}(tim,domFacs(an,5)));

subplot(2,4,n)
dat=sum([fH fA fE1 fE2 fE3],2);
dat = dat-min(dat);
dat = dat./max(dat);
plot(dat)
hold on
dat1 = traces{n,m}(1:60)-min(traces{n,3}(1:60));
dat1 = dat1./max(dat1);
plot(dat1)
ylim([-.1 1.1])
legend('Reconstructed','Real')
title(titles{n})

[rRec(n,m),pRec(n,m)] = corr(dat,dat1);
pRec(pRec>.05) = 0;
%rRec(pRec==0) = 0;
end
end
rRec=abs(rRec);
figure
subplot(121)
heatmap(rRec,'colorlimits',[0 1])
colormap('parula')
subplot(122)
heatmap(pRec,'colorlimits',[0 .05])
colormap('parula')

%% 17. identify overlap of TCA ensembles with stim-responsive neurons (fig 8, supp 10)
for a=1:25
    if a==1
        c=0;
    else
        c=sum(nCells(1:a-1));
    end
    for l=1:3
        pop1 = ensAns(a).(ens{l})-c;
        for k=1:3
            for m=1:5
                if k<3
                pop2 = [sigCells{a,m,k}{1} sigCells{a,m,k}{2}];
                else
                    pop2 = [sigCells{a,m,4}{1} sigCells{a,m,4}{2}];
                end
                fracEns(a,l,k,m) = mean(ismember(pop1,pop2));
            end
        end
    end
end