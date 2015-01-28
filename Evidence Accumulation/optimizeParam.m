%script to optimize parameters

%output in this case is ratio of neurons per second
%%
%designate parameter ranges
pRange.feedforwardWeight = 1;%0.01:0.02:2.0;
pRange.intraInhibDirectWeight = -0.05:-0.01:-2.0;
pRange.theta = 0:0.01:2;

%%%%%%%%%% network size and resolution 
em.nClusters = 50; %number of neuronal clusters
em.clusterSize = 1; %number of nuerons in each cluster. Total neurons = clusterSize*nClusters + nInhibitory
em.intraInhibition = true; %boolean of whether or not inhibition within a sequence is included
em.intraDirect = true; %flag of whether or not to directly inhibit neurons as external input or use modeled neurons
em.interInhibition = false; %boolean of whether or not inhibition between sequences is included
em.nIntraInhibitory = 1*em.nClusters; %number of inhibitory neurons within sequence
em.nInterInhibitory = 1*em.nClusters; %number of inhibitory neurons across sequences
em.nTimeBins = 200; %1,000 time bins
em.binSize = 5E-3; %.5 ms/bin*10,000 time bins = 5 second simulation
em.nSequences = 1; %number of sequences

%%%%%%%%%% network weights
em.mistuneWeight = 0; %percentage by which to mistune weights
em.recurrentWeight = 0;
em.clusterWeight = 0.04355; %weight of connections within cluster
em.intraInhibitoryWeight = -.2; %weight of connections from inhibitory to excitatory (should be negative)
em.interInhibitoryWeight = -.02; %weight of connections from inhibitory to excitatory (should be negative)
% em.feedforwardIntraInhibitoryWeight = 2*em.feedforwardWeight; %weight of connections from excitatory to inhibitory
% em.feedforwardInterInhibitoryWeight = 0.1*em.feedforwardWeight; %weight of connections from excitatory to inhibitory
em.scaleInhibDist = true; %boolean of whether or not to scale inhibitory weight by distance
em.inhibScaleFac = 5; %alpha for gaussian
em.routWeight = .005; %weight for rout
em.extWeight = 1; %weight of the external input (corresponds to ai in equation 2)
em.evidenceWeight = 1; %weight of the evidence input

%%%%%%%%%%% time constants
em.tau = .01; %time constant of change in seconds
em.gamma = 1; %decay constant to scale rate of decay

%%%%%%%%%%% direct inhibition constraints
em.capFlag = false; %should cap distance at specific value
em.capVal = .5*em.nClusters*em.clusterSize; %value to cap distance at
em.inhibDirectScaleFac = 5; %scale fac for width of gaussian for inhibitory weights direct

%%%%%%%%%%% connection probabilities
em.clusterConnProb = 1; %probability of neurons within the same cluster being connected
em.feedforwardConnProb = 1;%probability of neurons within cluster i being connected to neuron in cluster i+1
em.recurrentConnProb = 0; %probability of neuron connected to itself
em.feedforwardIntraInhibConnProb = 1;%probability of an excitatory neuron connected an inhibitory neuron
em.feedforwardInterInhibConnProb = 1;%probability of an excitatory neuron connected an inhibitory neuron
em.intraInhibitoryConnProb = 1; %probability of inhibitory neuron connected to excitatory neuron within sequence
em.interInhibitoryConnProb = 1; %probability of inhibitory neuron connected to excitatory neuron within sequence
em.evidenceConnProb = 1; %probability of connections from evidence to excitatory neurons
em.evNeurons = [1:em.nClusters*em.clusterSize;em.nClusters*em.clusterSize+1:2*em.nClusters*em.clusterSize]; %vector of neurons which evidence projects to; each row corresponds to each evidence pulse

%%%%%%%%%%%% input pulse to cluster 1 parameters
em.inputPulse = true; %boolean of whether or not input pulse should occur
em.pulseStart = .09; %time at which the pulse starts in seconds
em.pulseDuration = .1; %pulse duration in seconds
em.pulseAmplitude = 10; %amplitude of pulse
% em.inputNeurons = [1:em.clusterSize em.nClusters*em.clusterSize+1:(em.nClusters+1)*em.clusterSize]; %array containing the identities of the neurons which receive external input
% em.inputNeurons = [1:em.nClusters:(em.nSequences-1)*em.nClusters+1];
em.inputNeurons = 1;

%%%%%%%%%%%% timing sequence
em.timingSeq = false; %boolean of whether or not timing sequence should occur
em.numSteps = em.nClusters; %number of steps in the timing sequence
em.timingAmplitude = .1; %amplitude of each input

%%%%%%%%%%%% global excitatory drive
em.globalExc = false; %boolean of whether or not global drive should be present
em.globalExcNeurons = 1:em.nClusters*em.clusterSize*em.nSequences; %indices of neurons which receive global drive
em.globalExcAmp = .3; %amplitude of global excitatory drive

%%%%%%%%%%% evidence pulse parameters
em.evidencePulseStarts = [NaN]; %vector of pulse start times for evidence pulse; each row corresponds to different evidence
em.evidencePulseDurations = [.3; .3]; %scalar or vector of pulse durations matching length of start times
em.evidencePulseAmplitudes = [1; 1]; %scalar or vector of pulse amplitudes matching length of start times

%%%%%%%%%%% display parameters
em.shouldPlot = ~true; %whether or not output should be plot
em.nShow = 5; %number of neurons to show in the figure
em.alsoShow = [em.nClusters*em.clusterSize]; %neurons to show in addition to nShow

%%%%%%%%%%% debug figures
em.debug = ~true; %enables debug figures
em.maxFR = false; %enables max firing rate vs. time plot
em.seqSpeed = false; %enables index of max firing rate vs. time
em.inputPlot = ~true; %enables plot of inputs to a given neuron
em.inputPlotNeurons = [50]; %indices of neurons to plot
em.plotrout = ~true; %enables figure with rout plot
em.shouldWaitBar = false; %enables waitbar

%% run simulation

%open matlab pool
poolSize = matlabpool('size');
if poolSize == 0
    matlabpool open local 5;
end

%calculate total nSim
fields = fieldnames(pRange);
nSim = 1;
for i=1:length(fields)
    nSim = nSim*length(pRange.(fields{i}));
end

% seqSpeed = zeros(nSim,2); %column 1 - speed 1, column 2 - speed 2, column 3 - ratio
simParam = zeros(nSim,length(fields)); 
simOut = zeros(nSim,4); %1 - maxFR STD 2 - seqSpread (tMaxend - tMax0) 3 - minFR after 20% overall 4 - minFR std across neurons


progBar = waitbar(0,'Percent Complete','CreateCancelBtn',...
            'setappdata(progBar,''canceling'',true)');
setappdata(progBar,'canceling',false);

simInd = 1;
start = tic;
for i=1:length(pRange.(fields{1})) %massive for loop
    for j=1:length(pRange.(fields{2}))
        for k=1:length(pRange.(fields{3}))
%             for l = 1:length(pRange.(fields{4}))
                
                if getappdata(progBar,'canceling')
                    disp('optimizeParam aborted');
                    delete(progBar);
                    return;
                end

                %set paramaters
                em.(fields{1}) = pRange.(fields{1})(i);
                em.(fields{2}) = pRange.(fields{2})(j);
                em.(fields{3}) = pRange.(fields{3})(k);
%                 em.(fields{4}) = pRange.(fields{4})(l);
                
                em.feedforwardIntraInhibitoryWeight = 2*em.feedforwardWeight; %weight of connections from excitatory to inhibitory
                em.feedforwardInterInhibitoryWeight = 0.1*em.feedforwardWeight; %weight of connections from excitatory to inhibitory
                
                %store parameters
                simParam(simInd,1) = em.(fields{1});
                simParam(simInd,2) = em.(fields{2});
                simParam(simInd,3) = em.(fields{3});
%                 simParam(simInd,4) = em.(fields{4});
                
                
                %run model
                [frMatrix, ~, ~, ~] =  evidenceAccumulationModel(em);
                
                %save max and min of each neuron
                frMatrix = frMatrix';
                [maxFR, tMax] = max(frMatrix);
                [minFR] = min(frMatrix(:,round(.2*size(frMatrix,1)):end));
                
                %determine std of maxFR
                simOut(simInd,1) = std(maxFR);
                
                %determine range of times
                simOut(simInd,2) = tMax(end) - tMax(1);
                
                %get min FR after 20%
                simOut(simInd,3) = min(minFR);
                
                %get minFR spread
                simOut(simInd,4) = std(minFR);
                
                if simInd == 1
                    loopTime = toc(start);
                    totalTime = loopTime*nSim;
                end
                
                secRemaining = totalTime-(simInd*totalTime/nSim);
                
                waitbar(simInd/nSim,progBar,['Percent Complete: ',num2str(100*simInd/nSim),'%',...
                    'Minutes Remaining: ',num2str(secRemaining/60)]);
                
                %increment simInd
                simInd = simInd + 1;
%             end
        end
    end
end

delete(progBar);

% seqSpeed(:,3) = seqSpeed(:,2)./seqSpeed(:,1);

%close matlab pool
matlabpool close;