%%%%%%%%%% designate network parameters 
nClusters = 100; %number of neuronal clusters
clusterSize = 20; %number of nuerons in each cluster. Total neurons = clusterSize*nClusters
nTimeBins = 2E3; %1,000 time bins
binSize = 5E-3; %.5 ms/bin*10,000 time bins = 5 second simulation
feedforwardWeight = .00925; %cluster = .03, feedforward = .022 works
mistuneWeight = 0; %percentage by which to mistune weights
recurrentWeight = 0;
clusterWeight = 0.043; %weight of connections within cluster
routWeight = .005; %weight for rout
tau = .01; %time constant of change in seconds
gamma = 1; %decay constant to scale rate of decay
extWeight = .5; %weight of the external input (corresponds to ai in equation 2)
clusterConnProb = 1; %probability of neurons within the same cluster being connected
feedforwardConnProb = 1;%probability of neurons within cluster i being connected to neuron in cluster i+1
recurrentConnProb = 0; %probability of neuron connected to itself

%%%%%%%%%% designate input pulse to neuron 1 parameters
pulseStart = .1; %time at which the pulse starts in seconds
pulseDuration = .1; %pulse duration in seconds
pulseAmplitude = 20; %amplitude of pulse
inputNeurons = [1:20]; %array containing the identities of the neurons which receive external input

%%%%%%%%%% designate display parameters
shouldPlot = true; %whether or not output should be plot
nShow = 5; %number of neurons to show in the figure
alsoShow = [1]; %neurons to show in addition to nShow

%%%%%%%%%% run model
[frMatrix, rWeights, pulseFunc, rout] =  goldmanFunctionalFeedforward(...
    nClusters, clusterSize, nTimeBins, binSize, feedforwardWeight,...
    recurrentWeight, routWeight, tau, gamma, clusterWeight, extWeight,...
    clusterConnProb, feedforwardConnProb, recurrentConnProb,...
    pulseStart, pulseDuration, pulseAmplitude, inputNeurons, shouldPlot,...
    nShow, alsoShow, mistuneWeight);