%optimizeWeights.m script to find the optimal weights for the goldman
%simple feedforward model

%%%%%%%%%% designate network parameters 
nNeurons = 100;
nTimeBins = 1E3; %1,000 time bins
binSize = 5E-4; %.5 ms/bin*10,000 time bins = 5 second simulation
extWeight = 5;
feedbackWeight = 0;
recurrentWeight = 0;
tau = .01; %time constant of change in seconds
shouldScaleDist = true; %should the network weights scale with distance such that neurons close receive greater weights
lambda = 0.08;

%%%%%%%%%% designate input pulse to neuron 1 parameters
pulseStart = .1; %time at which the pulse starts in seconds
pulseDuration = .02; %pulse duration in seconds
pulseAmplitude = 5; %amplitude of pulse
inputNeurons = 1; %array containing the identities of the neurons which receive external input

%%%%%%%%%% designate display parameters
shouldPlot = false; %whether or not output should be plot
nShow = 5; %number of neurons to show in the figure

feedforwardWeightOptions = linspace(0,1,100);

maxAmpSpread = zeros(1,100);
timeSpread = zeros(1,100);
for i=1:100
    %%%%%%%%%%% run model
    frMatrix = goldmanFeedforward(...
        nNeurons, nTimeBins, binSize, feedforwardWeightOptions(i),...
        feedbackWeight, recurrentWeight, tau, extWeight, shouldScaleDist, lambda,...
        pulseStart, pulseDuration, pulseAmplitude, inputNeurons, shouldPlot,...
        nShow);
    
    [maxAmps, maxInds] = max(frMatrix(5:end,:),[],2);
    maxAmpSpread(i) = max(maxAmps) - min(maxAmps);
    timeSpread(i) = max(maxInds) - min(maxInds);
end