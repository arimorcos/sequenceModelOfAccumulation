%%%%%%%%%% network size and resolution 
em.nClusters = 50; %number of neuronal clusters
em.clusterSize = 1; %number of nuerons in each cluster. Total neurons = clusterSize*nClusters + nInhibitory
em.intraInhibition = true; %boolean of whether or not inhibition within a sequence is included
em.intraDirect = false; %flag of whether or not to directly inhibit neurons as external input or use modeled neurons
em.interInhibition = false; %boolean of whether or not inhibition between sequences is included
em.nIntraInhibitory = 1*em.nClusters; %number of inhibitory neurons within sequence
em.nInterInhibitory = 1*em.nClusters; %number of inhibitory neurons across sequences
em.nTimeBins = 200; %1,000 time bins
em.binSize = 5E-3; %.5 ms/bin*10,000 time bins = 5 second simulation
em.nSequences = 1; %number of sequences

%%%%%%%%%% network weights
em.feedforwardWeight = 1; %cluster = .04355, feedforward = .00925, intra = -.16 works
em.mistuneWeight = 0; %percentage by which to mistune weights
em.recurrentWeight = 0;
em.clusterWeight = 0.04355; %weight of connections within cluster
em.intraInhibDirectWeight = 0; %base weight of direct inhibitory weight (should be negative)
em.intraInhibitoryWeight = -.2; %weight of connections from inhibitory to excitatory (should be negative)
em.interInhibitoryWeight = -.02; %weight of connections from inhibitory to excitatory (should be negative)
em.feedforwardIntraInhibitoryWeight = 2*em.feedforwardWeight; %weight of connections from excitatory to inhibitory
em.feedforwardInterInhibitoryWeight = 0.1*em.feedforwardWeight; %weight of connections from excitatory to inhibitory
em.scaleInhibDist = true; %boolean of whether or not to scale inhibitory weight by distance
em.inhibScaleFac = 5; %alpha for gaussian
em.routWeight = .005; %weight for rout
em.extWeight = 1; %weight of the external input (corresponds to ai in equation 2)
em.evidenceWeight = 1; %weight of the evidence input

%%%%%%%%%%% time constants
em.tau = .01; %time constant of change in seconds
em.gamma = 1; %decay constant to scale rate of decay
em.theta = 30; %spontaneous firing rate factor

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
em.pulseAmplitude = 5; %amplitude of pulse
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
em.shouldPlot = true; %whether or not output should be plot
em.nShow = 5; %number of neurons to show in the figure
em.alsoShow = [em.nClusters*em.clusterSize]; %neurons to show in addition to nShow

%%%%%%%%%%% debug figures
em.debug = true; %enables debug figures
em.maxFR = false; %enables max firing rate vs. time plot
em.seqSpeed = false; %enables index of max firing rate vs. time
em.inputPlot = ~true; %enables plot of inputs to a given neuron
em.inputPlotNeurons = [50]; %indices of neurons to plot
em.plotrout = ~true; %enables figure with rout plot
em.shouldWaitBar = true; %enables waitbar

%%%%%%%%%% run model
[frMatrix, rWeights, pulseFunc, rout, em] =  evidenceAccumulationModel(em);

