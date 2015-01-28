%%%%%%%%%% designate network parameters 
nNeurons = 100;
nTimeBins = 2E3; %1,000 time bins
binSize = 5E-3; %.5 ms/bin*10,000 time bins = 5 second simulation
feedforwardWeight = 2.54; %gamma = 1.42/feedforward = .38 works with full = true;
feedbackWeight = 0;
recurrentWeight = 0;
routWeight = .2; %weight for rout
tau = .1; %time constant of change in seconds
gamma = 1; %decay constant to scale rate of decay
extWeight = 1; %weight of the external input (corresponds to ai in equation 2)
fullFeedforward = true; %should j project to all i's greater than j or only i=j+1
shouldScaleDist = true; %should the network weights scale with distance such that neurons close receive greater weights
lambda = 2; %length constant of the exponential decay used in the distance scale

%%%%%%%%%% designate input pulse to neuron 1 parameters
pulseStart = .1; %time at which the pulse starts in seconds
pulseDuration = .1; %pulse duration in seconds
pulseAmplitude = 20; %amplitude of pulse
inputNeurons = [1 2]; %array containing the identities of the neurons which receive external input

%%%%%%%%%% designate display parameters
shouldPlot = true; %whether or not output should be plot
nShow = 5; %number of neurons to show in the figure
alsoShow = []; %neurons to show in addition to nShow

%%%%%%%%%%% run model
[frMatrix, rWeights, pulseFunc, rout] = goldmanFeedforward(...
    nNeurons, nTimeBins, binSize, feedforwardWeight, fullFeedforward,...
    feedbackWeight, recurrentWeight, routWeight, tau, gamma, extWeight, shouldScaleDist, lambda,...
    pulseStart, pulseDuration, pulseAmplitude, inputNeurons, shouldPlot,...
    nShow, alsoShow);