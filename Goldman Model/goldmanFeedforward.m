function [frMatrix, rWeights, pulseFunc, rout] =  goldmanFeedforward(...
    nNeurons, nTimeBins, binSize, feedforwardWeight, fullFeedforward,...
    feedbackWeight, recurrentWeight, routWeight, tau, gamma, extWeight,...
    shouldScaleDist, lambda,...
    pulseStart, pulseDuration, pulseAmplitude, inputNeurons, shouldPlot,...
    nShow, alsoShow)
%goldmanFeedforward.m Master function to replicate the feedforward memory network
%presented in the 2009 Mark Goldman Neuron paper. 
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
%nNeurons - network size
%nTimeBins - number of time bins in the network 
%binSize - bin length in seconds. Total simulation length =
%          binSize*nTimeBins
%%feedforwardWeight - base weight of projections from neurons in which j <
%                      i
%feedbackWeight - base weight of projections from neurons in which j > i
%fullFeedforward - should j project to all i's greater than j or only i=j+1
%recurrentWeight - base weight of projections from neurons in which j=i
%routWeight - weight of readout
%tau - time constant in seconds
%extWeight - weight of projections from input pulse
%shouldScaleDist - boolean of whether or not to scale weights by distance
%lambda - length constant of the exponential decay used in the distance
%           scale
%pulseStart - start time in seconds of pulse
%pulseDuration - duration in seconds of pulse
%pulseAmplitude - max firing rate of pulse 
%inputNeurons - neuron ids which pulse projects to
%shouldPlot - boolean of whether or not to plot
%nShow - number of neurons to show on plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize array of firing rates and rout
frMatrix = zeros(nNeurons,nTimeBins);
rout = zeros(1,nTimeBins);

%generate nNeurons x nNeurons weight matrix
rWeights = generateWeightMatrix(feedforwardWeight, fullFeedforward,... 
    feedbackWeight, recurrentWeight, nNeurons, shouldScaleDist, lambda);

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

function rWeights = generateWeightMatrix(forward, full, back, recurrent,...
    nNeurons, shouldScaleDist,lambda)
%function to generate weight matrix with uniform values for feedforward,
%feedback, and recurrent connections

rWeights = eye(nNeurons);
if shouldScaleDist
    distMatrix = ones(nNeurons);
end

for i=1:nNeurons
    eyeIndex = find(rWeights(:,i)==1); %find the diagonal for the ith column
    rWeights(1:eyeIndex-1,i) = back; %all projections in which j>i are set to feedbackWeight
    rWeights(eyeIndex,i) = recurrent; %all projections in which j=i are set to recurrentWeight
    if full
        rWeights(eyeIndex+1:end,i) = forward; %all projections in which i>j are set to feedbackWeight
    else
        if i<nNeurons
            rWeights(eyeIndex+1,i) = forward;
        end
    end
    
    %scale according to distance
    if shouldScaleDist
        for j=1:nNeurons
%             distMatrix(i,j) = 1 - (abs(i-j)/nNeurons);
            distMatrix(i,j) = forward*exp(-(abs(i-j))*lambda);
        end
    end
    
end

%scale by distance
if shouldScaleDist
    rWeights = rWeights.*distMatrix;
end

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

function newFR = updateFRSimple(oldFR, weights, tau, gamma, deltat,...
    extVal, extWeight, inputNeurons)
%function to convert input to a neuron into a firing rate based on the
%weights contained in the weights matrix

newFR = zeros(size(oldFR));
newFR(1) = extVal;
newFR(2:end) = oldFR(1:end-1)-oldFR(2:end);

newFR = (newFR.*deltat)./tau;
newFR = rectify(newFR);

end

function newRout = updateRout(oldRout, networkFR, weights, tau, deltat)
%function to update FR of rout
    
input = -oldRout + sum(networkFR'.*weights); %decay + sum(ri*Wij)

newRout = oldRout + input*deltat/tau;
newRout = rectify(newRout);
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