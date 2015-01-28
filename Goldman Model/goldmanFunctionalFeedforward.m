function [frMatrix, rWeights, pulseFunc, rout] =  goldmanFunctionalFeedforward(...
    nClusters, clusterSize, nTimeBins, binSize, feedforwardWeight,...
    recurrentWeight, routWeight, tau, gamma, clusterWeight, extWeight,...
    clusterConnProb, feedforwardConnProb, recurrentConnProb,...
    pulseStart, pulseDuration, pulseAmplitude, inputNeurons, shouldPlot,...
    nShow, alsoShow, mistune)
%goldmanFunctionalFeedforward.m Master function to replicate the
%functionally feedforward memory network presented in the 2009 Mark Goldman
%Neuron paper.
%
%OUTPUTS
%frMatrix - nNeurons x nTimeBins array containing the firing rate of each
%           neuron at each time point
%rWeights - nNeurons x nNeurons matrix of weights from neuron j to neuron i
%pulseFunc - 1 x nTimeBins row vector containing value of input pulse
%rout - 1 x nTimeBins row vector containing FR of readout neuron
%
%
%INPUTS
%nClusters - number of reccurrently connected clusters in the network (i.e.
%               number of stages)
%clusterSize = number of neurons in each recurrently connected cluster
%nTimeBins - number of time bins in the network 
%binSize - bin length in seconds. Total simulation length =
%          binSize*nTimeBins
%feedforwardWeight - weight of projections between clusters
%feedbackWeight - base weight of projections from neurons in which j > i
%clusterWeight - weight of recurrent connections within clusters
%recurrentWeight - base weight of projections from neurons in which j=i
%clusterConnProb - connection probability of neurons within the same
%                   cluster
%feedforwardConnProb - connection probability of neurons between clusters 
%recurrentConnProb - connection probaiblity of neurons to themselves
%routWeight - weight of readout
%tau - time constant in seconds
%extWeight - weight of projections from input pulse
%pulseStart - start time in seconds of pulse
%pulseDuration - duration in seconds of pulse
%pulseAmplitude - max firing rate of pulse 
%inputNeurons - neuron ids which pulse projects to
%shouldPlot - boolean of whether or not to plot
%nShow - number of neurons to show on plot
%mistune - percentage by which to mistune feedforward weight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate the number of neurons required 
nNeurons = nClusters*clusterSize; %total neurons in network

%initialize array of firing rates and rout
frMatrix = zeros(nNeurons,nTimeBins);
rout = zeros(1,nTimeBins);

%mistune weights
feedforwardWeight = mistuning(feedforwardWeight,mistune);
% clusterWeight = mistuning(clusterWeight,mistune);

%generate nNeurons x nNeurons weight matrix
rWeights = generateWeightMatrix(clusterWeight, feedforwardWeight,...
    recurrentWeight, nNeurons, nClusters, clusterSize, clusterConnProb,...
    feedforwardConnProb, recurrentConnProb);

%generate rout weight matrix
routWeightMat = routWeight*ones(1,nNeurons);

%generate pulse function
pulseFunc = generatePulseFunc(pulseStart, pulseDuration, pulseAmplitude,...
    nTimeBins, binSize);

%generate waitbar
progBar = waitbar(0,'Percent Complete');

%run the model
for i=2:nTimeBins
    
    frMatrix(:,i) = updateFR(frMatrix(:,i-1), rWeights, tau, gamma,...
        binSize, pulseFunc(i), extWeight, inputNeurons); %update firing rate
    
    rout(i) = updateRout(rout(i-1), frMatrix(:,i-1), routWeightMat, tau, binSize); %update rout
    
    waitbar(i/nTimeBins,progBar,['Percent Complete: ',num2str(100*i/nTimeBins),'%']); %update waitbar
end

close(progBar); %close waitbar

%plot
if shouldPlot
    plotGoldman(binSize, nTimeBins, pulseFunc, nNeurons, nShow, alsoShow, frMatrix, rout);
end

end

function out = mistuning(in,mistune)
%function to mistune weights


mistune = (100 + mistune)/100; %turn into percentage of total
out = mistune*in;

end

function out = rectify(in)
%function to perform half-wave rectification of the firing rate

out = max(in,0);

end

function pulseFunc = generatePulseFunc(start, duration, amplitude, nTimeBins, binSize)
%function to generate an 1 x nTimeBins row vector describing a pulse with
%max amplitude defined by amplitude which begins at start (in seconds) and
%lasts for duration seconds. binSize specifies the length of one bin in
%seconds.

startBin = start/binSize; %convert the start time in seconds to the start bin
durationBins = duration/binSize; %convert the duration in seconds to duration in bins
startBin = max(startBin,1); %set startBin to 1 if 0

pulseFunc = zeros(1,nTimeBins); %initialize row vector 
pulseFunc(1,startBin:startBin+durationBins) = amplitude; %set pulse

end

function connMatrix = generateConnMatrix(nNeurons, nClusters, clusterSize,...
    clusterConnProb, feedforwardConnProb, recurrentConnProb)
%generates boolean connectivity matrix with size nNeurons x nNeurons based
%on connection probabilities between neurons within clusters and between
%clusters

connMatrix = zeros(nNeurons); %initialize matrix

%perform modifications based on each cluster
ind = 1;
for i = 1:nClusters
    
    %set within cluster probability
    tempCluster = reshape(randsample([0 1],clusterSize^2,true,...
        [1-clusterConnProb clusterConnProb]),clusterSize,clusterSize); %generate weighted connectivity and reshape matrix into clusterSize x clusterSize
    
    %set recurrent prob
    tempCluster(logical(eye(size(tempCluster)))) = randsample([0 1],clusterSize,true,...
        [1-recurrentConnProb recurrentConnProb]); %set diagonal to proper recurrent weights
    
    connMatrix(ind:i*clusterSize,ind:i*clusterSize) = tempCluster;
    
    %set forward prob
    tempForCluster = reshape(randsample([0 1],clusterSize^2,true,...
        [1-feedforwardConnProb feedforwardConnProb]),clusterSize,clusterSize); %generate weighted connectivity and reshape matrix into clusterSize x clusterSize
    
    if i ~= nClusters %only set this if not the last cluster
        connMatrix(ind+clusterSize:(i+1)*clusterSize,ind:i*clusterSize) = tempForCluster;
    end
        
    
    ind = ind + clusterSize;
    
end

connMatrix = logical(connMatrix);

end

function rWeights = generateWeightMatrix(cluster, forward, recurrent, nNeurons,...
    nClusters, clusterSize, clusterConnProb, feedforwardConnProb,...
    recurrentConnProb)
%function to generate weight matrix with uniform values for feedforward,
%feedback, and recurrent connections

%generate connectivity matrices
connMatrixCluster = generateConnMatrix(nNeurons, nClusters, clusterSize,... 
    clusterConnProb, 0, 0);

connMatrixFeedforward = generateConnMatrix(nNeurons, nClusters, clusterSize,... 
    0, feedforwardConnProb, 0);

connMatrixRecurrent = generateConnMatrix(nNeurons, nClusters, clusterSize,... 
    0, 0, recurrentConnProb);

%initialize rWeights
rWeights = zeros(size(connMatrixCluster));

%set cluster weights
rWeights(connMatrixCluster) = cluster;

%set feedforward weights
rWeights(connMatrixFeedforward) = forward;

%set recurrent weights
rWeights(connMatrixRecurrent) = recurrent;

end

function newFR = updateFR(oldFR, weights, tau, gamma, deltat,...
    extVal, extWeight, inputNeurons)
%function to convert input to a neuron into a firing rate based on the
%weights contained in the weights matrix

feedforwardInput = oldFR'*weights'; %sum(WijRj)

externalInput = zeros(size(oldFR)); %initialize array same as number of neurons
externalInput(inputNeurons) = extVal*extWeight; %ai*x(t)

taudFRdt = feedforwardInput'-gamma*oldFR + externalInput; %-Ri + sum(Wij*Rj) + ai*x(t)

dFR = (taudFRdt.*deltat)./tau; %dFR = dFR*dt/tau

newFR = oldFR + dFR; %new firing rate is old firing rate plus change in firing rate

newFR = rectify(newFR); %half-wave rectify new firing rate (set any negative values to 0)


end

function plotGoldman(binSize, nTimeBins, pulseFunc,...
    nNeurons, nShow, alsoShow, frMatrix, rout)
%function to plot activity and pulse

figure('Name','Neuronal Activity in Goldman Model');
xTimes = linspace(0,binSize*nTimeBins,nTimeBins);

pulsePlot = subplot(3,1,1);
plot(xTimes,pulseFunc,'b');
ylim([0 2*max(pulseFunc)]);
xlabel('Time (seconds)');
ylabel('Input');
title('External Input to Neuron 1');

neuronPlot = subplot(3,1,2);
hold on;
whichNeurons = round(unique([alsoShow linspace(1,nNeurons,nShow)]));
colors = hsv(length(whichNeurons));
legendStr = cell(1,length(whichNeurons));
for i=1:length(whichNeurons)
    plot(xTimes,frMatrix(whichNeurons(i),:),'Color',colors(i,:));
    legendStr{i} = num2str(whichNeurons(i));
end
xlabel('Time (seconds)');
ylabel('Firing Rate');
title('Neuronal Activity');
legend(legendStr,'Location','East');

routPlot = subplot(3,1,3);
plot(xTimes,rout,'r');
ylim([0 2*max(rout)]);
xlabel('Time (seconds)');
ylabel('r_{out}');
title('r_{out}');

linkaxes([pulsePlot neuronPlot routPlot],'x');

end

function newRout = updateRout(oldRout, networkFR, weights, tau, deltat)
%function to update FR of rout
    
input = -oldRout + sum(networkFR'.*weights); %decay + sum(ri*Wij)

newRout = oldRout + input*deltat/tau;
newRout = rectify(newRout);
end
