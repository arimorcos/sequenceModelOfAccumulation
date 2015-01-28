function [frMatrix, rWeights, pulseFunc, rout, em] =  evidenceAccumulationModel(em)
%evidenceAccumulationModel.m Master function to model evidence accumulation
%by sequences of neuronal activation
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
%inhibition - boolean of whether or not to include inhibitory neurons
%nInhibitory - number of inhibitory interneurons
%inhibitoryWeight - base weight for inhibition
%inhibConnProb - connection probability of inhibitory neuron to excitatory
%                   neurons
%feedforwardInhibConnProb - connection probability between neurons of a
%                           cluster and their inhibitory neurons
%scaleInhibDist - boolean of whether or not to scale the inhibitory weights
%                   by distance
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
%evidenceWeight - connection weight of evidence projection
%evidenceConnProb - connection probability of evidence signal to excitatory
%                   neurons
%evidencePulseStarts - vector of pulse start times for evidence
%evidencePulseDurations - vector of evidence pulse durations
%evidencePulseAmplitudes - vector of evidence pulse amplitudes
%evNeurons - vector of neurons which evidence connects to

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if em.interInhibition && ~em.intraInhibition
    error('cannot have inter inhibition without intra inhibition');
end

%calculate the number of neurons required 
em.nNeurons = em.nSequences*em.nClusters*em.clusterSize; %total neurons in network
if em.intraInhibition && ~em.intraDirect %if inhibition incorporated
    em.nNeurons = em.nNeurons + em.nSequences*em.nIntraInhibitory; %add the inhibitory neurons
end
if em.interInhibition %if inhibition between sequences incorporated
    em.nNeurons = em.nNeurons + em.nSequences*em.nInterInhibitory;
end

%initialize array of firing rates and rout
frMatrix = zeros(em.nNeurons,em.nTimeBins);
rout = zeros(em.nSequences,em.nTimeBins);

%if intraInihibition direct, generate distance matrix
distMat = generateDistMatrix(em);

%mistune weights
em.feedforwardWeight = mistuning(em.feedforwardWeight,em.mistuneWeight);
% em.clusterWeight = mistuning(em.clusterWeight,em.mistuneWeight);

%generate nNeurons x nNeurons weight matrix
rWeights = generateWeightMatrix(em);

%generate rout weight matrix
routWeightMat = em.routWeight*ones(1,em.nNeurons);

%generate pulse function
if em.inputPulse
    pulseFunc = generatePulseFunc(em.pulseStart, em.pulseDuration, em.pulseAmplitude,...
        em.nTimeBins, em.binSize);
else
    pulseFunc = zeros(1,em.nTimeBins);
end

%generate timing vals
em.timingThresh = round(linspace(1,em.nTimeBins,em.numSteps)); %generate threshold indices for each timing stage specifiying where input should go
em.clusterIDs = [];
for i = 1:em.nSequences %generate matrix containing indices of neurons in each cluster
    em.clusterIDs = cat(2,em.clusterIDs,reshape((i-1)*em.nClusters*em.clusterSize+1:i*em.nClusters*em.clusterSize,...
        em.clusterSize,em.nClusters)'); 
end
if ~em.timingSeq %if not on, set amplitude to 0
    em.timingAmplitude = 0;
end

%check for global excitatory
if ~em.globalExc
    em.globalExcAmp = 0;
end

%generate evidence function for each piece of evidence
for i=1:size(em.evidencePulseStarts,1)
    evidenceFunc(i,:) = generatePulseFunc(em.evidencePulseStarts(i,:), em.evidencePulseDurations(i,:),...
        em.evidencePulseAmplitudes(i,:), em.nTimeBins, em.binSize);
end

%generate waitbar
if em.shouldWaitBar
    progBar = waitbar(0,'Percent Complete');
end

%initialize inputs 
em.inputs.feedforward = zeros(em.nNeurons,em.nTimeBins);
em.inputs.evidence = zeros(em.nNeurons,em.nTimeBins);
em.inputs.external = zeros(em.nNeurons,em.nTimeBins);
em.inputs.globalExc = zeros(em.nNeurons,em.nTimeBins);
em.inputs.timing = zeros(em.nNeurons,em.nTimeBins);
em.inputs.dirInhib = zeros(em.nNeurons,em.nTimeBins);
em.inputs.theta = em.theta*ones(em.nNeurons,em.nTimeBins);

%run the model
for i=2:em.nTimeBins
    
    %determine current cluster
    currCluster = find(i >= em.timingThresh,1,'last');
    currClusterIDs = em.clusterIDs(currCluster,:);
    
    [frMatrix(:,i), em.inputs] = updateFR(frMatrix(:,i-1), rWeights, em.tau, em.gamma,...
        em.binSize, pulseFunc(i), em.extWeight, em.inputNeurons, evidenceFunc(:,i),...
        em.evidenceWeight, em.evNeurons, currClusterIDs, em.timingAmplitude,...
        em.globalExcNeurons, em.globalExcAmp, em.inputs, i, distMat); %update firing rate
    
    for j=1:em.nSequences
        seqInds = (j-1)*em.clusterSize*em.nClusters+1:j*em.clusterSize*em.nClusters;
        rout(j,i) = updateRout(rout(j,i-1), frMatrix(seqInds,i-1),...
            routWeightMat(seqInds), em.tau, em.binSize); %update rout
    end
    
    if em.shouldWaitBar
        waitbar(i/em.nTimeBins,progBar,['Percent Complete: ',num2str(100*i/em.nTimeBins),'%']); %update waitbar
    end
end

if em.shouldWaitBar
    close(progBar); %close waitbar
end

%plot
if em.shouldPlot
    plotModel(em.binSize, em.nTimeBins, pulseFunc, em.nNeurons, em.nShow, em.alsoShow,...
        frMatrix, rout, evidenceFunc, em.plotrout);
end

if em.debug
    plotDebug(em,frMatrix);
end

end

function distMat = generateDistMatrix(em)
%generateDistMatrix.m generates distance matrix 

neuronArray = 1:em.nSequences*em.nClusters*em.clusterSize;

[j, i] = meshgrid(neuronArray,neuronArray);
distMat = abs(i-j);

if em.capFlag 
    distMat(distMat > capFlag) = em.capVal;
end

%scale distMat
distMat = distMat.*em.intraInhibDirectWeight;

if ~em.intraDirect
    distMat = [];
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

function connMatrix = generateConnMatrix(nNeurons, nClusters, clusterSize,...
    nSeq, clusterConnProb, feedforwardConnProb, recurrentConnProb,...
    intraInhibConnProb, intraInhibition, intraDirect, nIntraInhibitory, feedforwardInhibConnProb,...
    interInhibConnProb, interInhibition, nInterInhibitory)
%generates boolean connectivity matrix with size nNeurons x nNeurons based
%on connection probabilities between neurons within clusters and between
%clusters

connMatrix = zeros(nNeurons); %initialize matrix

for i=1:nSeq %for each sequence
    %generate connectivity from intra inhibitory neurons
    nExc = clusterSize*nClusters; %number of excitatory neurons
    if intraInhibition && ~intraDirect %if inhibitory neurons are included

        %set output connectivity
        connMatrix(max(nExc*(i-1)+1,1):i*nExc, nSeq*nExc+(i-1)*nIntraInhibitory+1:nSeq*nExc+i*nIntraInhibitory) =...
            reshape(randsample([0 1],nExc*nIntraInhibitory,...
            true, [1 - intraInhibConnProb intraInhibConnProb]), nExc, nIntraInhibitory);
    end
    
    if interInhibition
        rowInd = 1:nSeq*nExc;
        rowInd = rowInd(~ismember(rowInd,max(nExc*(i-1)+1,1):i*nExc));
        colInd = nSeq*(nExc+nIntraInhibitory)+(i-1)*nInterInhibitory+1:...
            nSeq*(nExc+nIntraInhibitory)+i*nInterInhibitory;
        connMatrix(rowInd,colInd) =...
            reshape(randsample([0 1],(nSeq-1)*nExc*nInterInhibitory,...
            true, [1 - interInhibConnProb interInhibConnProb]), (nSeq-1)*nExc, nInterInhibitory);
    end
    
    intraInhibPerCluster = nIntraInhibitory/nClusters;
    interInhibPerCluster = nInterInhibitory/nClusters;
end


%perform modifications based on each cluster
for i=1:nSeq
    ind = (i-1)*nExc+1;
    for j = 1:nClusters
        
        %set within cluster probability
        tempCluster = reshape(randsample([0 1],clusterSize^2,true,...
            [1-clusterConnProb clusterConnProb]),clusterSize,clusterSize); %generate weighted connectivity and reshape matrix into clusterSize x clusterSize
        
        %set recurrent prob
        tempCluster(logical(eye(size(tempCluster)))) = randsample([0 1],clusterSize,true,...
            [1-recurrentConnProb recurrentConnProb]); %set diagonal to proper recurrent weights
        
        connMatrix(ind:(i-1)*nExc+j*clusterSize,ind:(i-1)*nExc+j*clusterSize) = tempCluster;
        
        %set forward prob
        tempForCluster = reshape(randsample([0 1],clusterSize^2,true,...
            [1-feedforwardConnProb feedforwardConnProb]),clusterSize,clusterSize); %generate weighted connectivity and reshape matrix into clusterSize x clusterSize
        
        if j ~= nClusters %only set this if not the last cluster
            connMatrix(ind+clusterSize:(i-1)*nExc+(j+1)*clusterSize,ind:(i-1)*nExc+j*clusterSize) = tempForCluster;
        end
        
        %set inhibitory feedforward projection
        if intraInhibition && ~intraDirect
            %set input connectivity
            connMatrix(nSeq*nExc+(i-1)*nIntraInhibitory+1+(j-1)*intraInhibPerCluster:...
                nSeq*nExc+(i-1)*nIntraInhibitory+1+j*intraInhibPerCluster-1, ind:(i-1)*nExc+j*clusterSize)...
                = reshape(randsample([0 1],intraInhibPerCluster*clusterSize,true,...
                [1 - feedforwardInhibConnProb feedforwardInhibConnProb]),...
                intraInhibPerCluster, clusterSize);
        end
        
        if interInhibition
            %set input connectivity
            connMatrix(nSeq*(nExc+nIntraInhibitory)+(i-1)*nInterInhibitory+1+(j-1)*interInhibPerCluster:...
                nSeq*(nExc+nIntraInhibitory)+(i-1)*nInterInhibitory+1+j*interInhibPerCluster-1, ind:(i-1)*nExc+j*clusterSize)...
                = reshape(randsample([0 1],interInhibPerCluster*clusterSize,true,...
                [1 - feedforwardInhibConnProb feedforwardInhibConnProb]),...
                interInhibPerCluster, clusterSize);
        end
        
        
        ind = ind + clusterSize;
        
    end
end

connMatrix = logical(connMatrix);

end

function gaussian = generateGaussian(nExc,clusterID,clusterSize,alpha)
%function to generate a gaussian offset by cluster location

%find dimension id
nClusters = nExc/clusterSize;
dimID = nClusters - clusterID + 1;

%make gaussian window with twice number of neurons
tempGauss = gausswin(2*nClusters,alpha);

%chop ouot portion corresponding to longest dimension
tempGauss = tempGauss(dimID:dimID+nClusters-1);

%expand gaussian so that all neurons in a given cluster receive the same
%input
gaussian = repmat(tempGauss',clusterSize,1);
gaussian = reshape(gaussian,1,nExc)';

end

function rWeights = generateWeightMatrix(em)
%function to generate weight matrix with uniform values for feedforward,
%feedback, and recurrent connections

%generate connectivity matrices
connMatrixCluster = generateConnMatrix(em.nNeurons, em.nClusters, em.clusterSize,... 
    em.nSequences, em.clusterConnProb, 0, 0, 0, em.intraInhibition, em.intraDirect,...
    em.nIntraInhibitory, 0, 0, em.interInhibition, em.nInterInhibitory);

connMatrixFeedforward = generateConnMatrix(em.nNeurons, em.nClusters, em.clusterSize,... 
    em.nSequences,0, em.feedforwardConnProb, 0, 0, em.intraInhibition, em.intraDirect,...
    em.nIntraInhibitory, 0, 0, em.interInhibition, em.nInterInhibitory);

connMatrixRecurrent = generateConnMatrix(em.nNeurons, em.nClusters, em.clusterSize,... 
    em.nSequences, 0, 0, em.recurrentConnProb, 0, em.intraInhibition, em.intraDirect,...
    em.nIntraInhibitory, 0, 0, em.interInhibition, em.nInterInhibitory);

if em.intraInhibition
    connMatrixFeedforwardIntraInhibitory = generateConnMatrix(em.nNeurons, em.nClusters, em.clusterSize,...
        em.nSequences, 0, 0, 0, 0, em.intraInhibition,  em.intraDirect, em.nIntraInhibitory,...
        em.feedforwardIntraInhibConnProb, 0, false, em.nInterInhibitory);
    
    connMatrixIntraInhibitory = generateConnMatrix(em.nNeurons, em.nClusters,...
        em.clusterSize,em.nSequences, 0, 0, 0,...
        em.intraInhibitoryConnProb,em.intraInhibition,  em.intraDirect, em.nIntraInhibitory,...
        0, 0, em.interInhibition, em.nInterInhibitory);
end

if em.interInhibition
    connMatrixFeedforwardInterInhibitory = generateConnMatrix(em.nNeurons, em.nClusters, em.clusterSize,...
        em.nSequences, 0, 0, 0, 0, false,  em.intraDirect, em.nIntraInhibitory, ...
        em.feedforwardInterInhibConnProb, 0, em.interInhibition, em.nInterInhibitory);
    
    connMatrixInterInhibitory = generateConnMatrix(em.nNeurons, em.nClusters,...
        em.clusterSize,em.nSequences, 0, 0, 0,...
        em.intraInhibitoryConnProb,false,  em.intraDirect, em.nIntraInhibitory,...
        0, em.interInhibitoryConnProb, em.interInhibition, em.nInterInhibitory);
end

if em.intraInhibition && em.scaleInhibDist && ~em.intraDirect %if inhibitory weights scale with distance
    nExc = em.clusterSize*em.nClusters; %number of excitatory neurons
    intraInhibPerCluster = em.nIntraInhibitory/em.nClusters;
    connMatrixIntraInhibitory = double(connMatrixIntraInhibitory);
    if em.interInhibition
        interInhibPerCluster = em.nInterInhibitory/em.nClusters;
        connMatrixInterInhibitory = double(connMatrixInterInhibitory);
    end
    
    for i=1:em.nSequences %for each sequence
        for j=1:em.nClusters %for each cluster
            
            %intra
            intraGaussian = 1-generateGaussian(nExc,j,em.clusterSize,em.inhibScaleFac); %generate gaussian centered on proper cluster and invert (stronger weights farther away)
            tempWeightedIntraInhib = connMatrixIntraInhibitory(nExc*(i-1)+1:i*nExc,...
                em.nSequences*nExc+(i-1)*em.nIntraInhibitory+1+(j-1)*intraInhibPerCluster:...
                em.nSequences*nExc+(i-1)*em.nIntraInhibitory+j*intraInhibPerCluster); %move to temp array and convert to double
            
            for k=1:size(tempWeightedIntraInhib,2) %for each inhibitory neuron in cluster i
                tempWeightedIntraInhib(:,k) = tempWeightedIntraInhib(:,k).*intraGaussian; %multiply connectivity matrix by gaussian
            end

            connMatrixIntraInhibitory(nExc*(i-1)+1:i*nExc,...
                em.nSequences*nExc+(i-1)*em.nIntraInhibitory+1+(j-1)*intraInhibPerCluster:...
                em.nSequences*nExc+(i-1)*em.nIntraInhibitory+j*intraInhibPerCluster) =...
                tempWeightedIntraInhib;  %put back in main conn matrix       
            
            %inter
            if em.interInhibition
                interGaussian = 1-intraGaussian;

                rowInd = 1:em.nSequences*nExc;
                rowInd = rowInd(~ismember(rowInd,max(nExc*(i-1)+1,1):i*nExc));
                colInd = em.nSequences*(nExc+em.nIntraInhibitory)+(i-1)*em.nInterInhibitory...
                    +1+(j-1)*interInhibPerCluster:em.nSequences*(nExc+em.nIntraInhibitory)+...
                    (i-1)*em.nInterInhibitory+j*interInhibPerCluster;

                tempWeightedInterInhib = connMatrixInterInhibitory(rowInd,colInd); %move to temp array and convert to double

                for k=1:size(tempWeightedInterInhib,2) %for each inhibitory neuron in cluster i
                    tempWeightedInterInhib(:,k) = tempWeightedInterInhib(:,k).*interGaussian; %multiply connectivity matrix by gaussian
                end

                connMatrixInterInhibitory(rowInd,colInd) =...
                    tempWeightedInterInhib;  %put back in main conn matrix   
            end
        end
    end
end

%initialize rWeights
rWeights = zeros(size(connMatrixCluster));

%set cluster weights
rWeights(connMatrixCluster) = em.clusterWeight;

%set feedforward weights
rWeights(connMatrixFeedforward) = em.feedforwardWeight;

%set recurrent weights
rWeights(connMatrixRecurrent) = em.recurrentWeight;

if em.intraInhibition && ~em.intraDirect
    %set feedforward inhibitory weights
    rWeights(connMatrixFeedforwardIntraInhibitory) = em.feedforwardIntraInhibitoryWeight;
    
    %set inhibitory weights
    tempInhibMat = connMatrixIntraInhibitory*em.intraInhibitoryWeight;
    rWeights(connMatrixIntraInhibitory~=0) = tempInhibMat(tempInhibMat~=0);
end

if em.interInhibition
    %set feedforward inhibitory weights
    rWeights(connMatrixFeedforwardInterInhibitory) = em.feedforwardInterInhibitoryWeight;
    
    %set inhibitory weights
    tempInhibMat = connMatrixInterInhibitory*em.interInhibitoryWeight;
    rWeights(connMatrixInterInhibitory~=0) = tempInhibMat(tempInhibMat~=0);
end

end

function [newFR, inputs] = updateFR(oldFR, weights, tau, gamma, deltat,...
    extVal, extWeight, inputNeurons, evVal, evidenceWeight, evNeurons,...
    currClusterIDs, timingAmp, globalNeurons, globalAmp, inputs, i,...
    distMat)
%functidon to convert input to a neuron into a firing rate based on the
%weights contained in the weights matrix

inputs.feedforward(:,i) = (oldFR'*weights')'; %sum(WijRj)

inputs.external(inputNeurons,i) = extVal*extWeight; %ai*x(t)

inputs.timing(currClusterIDs,i) = timingAmp;   

if ~isempty(distMat)
    inputs.dirInhib(:,i) = oldFR'*distMat;
else
    inputs.dirInhib(:,i) = zeros(size(oldFR));
end

for j=1:length(evVal)
    inputs.evidence(evNeurons(j,:),i) = evVal(j)*evidenceWeight;%ai*x(t) for evidence
end

inputs.globalExc(globalNeurons,i) = globalAmp;

fx = inputs.feedforward(:,i) +...
    inputs.external(:,i) + inputs.evidence(:,i) +...
    inputs.timing(:,i) + inputs.globalExc(:,i) +...
    inputs.dirInhib(:,i) + inputs.theta(:,i); %f(x) = sum(Wij*Rj) + ai*x(t) + theta

fx = rectify(fx); %half-wave rectify new firing rate (set any negative values to 0)

taudFRdt = fx-gamma*oldFR; %-Ri + f(x)

dFR = (taudFRdt.*deltat)./tau; %dFR = dFR*dt/tau

newFR = oldFR + dFR; %new firing rate is old firing rate plus change in firing rate


end

function plotModel(binSize, nTimeBins, pulseFunc,...
    nNeurons, nShow, alsoShow, frMatrix, rout, evidenceFunc, plotrout)
%function to plot activity and pulse

xTimes = linspace(0,binSize*nTimeBins,nTimeBins);

if plotrout
    figure('Name','Neuronal Activity in Goldman Model');

    pulsePlot = subplot(3,1,1);
    plot(xTimes,pulseFunc,'b');
    hold on;
    colors = distinguishable_colors(size(evidenceFunc,1));
    for i=1:size(evidenceFunc,1)
        plot(xTimes,evidenceFunc(i,:),'Color',colors(i,:));
    end
    ylim([0 max(1,2*max(pulseFunc))]);
    legend('Input to cluster 1','Evidence')
    xlabel('Time (seconds)');
    ylabel('Input');
    title('External Input');

    neuronPlot = subplot(3,1,2);
    hold on;
    whichNeurons = round(unique([alsoShow linspace(1,nNeurons,nShow)]));
    colors = distinguishable_colors(length(whichNeurons));
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
    hold on;
    for i=1:size(rout,1)
        plot(xTimes,rout(i,:),'Color',colors(i,:));
    end
    ylim([0 max(1,2*max(rout(:)))]);
    xlabel('Time (seconds)');
    ylabel('r_{out}');
    title('r_{out}');

    linkaxes([pulsePlot neuronPlot routPlot],'x');
end

figure('Name','Neuronal Activity in Goldman Model');
warning('off','MATLAB:linkaxes:RequireDataAxes');
pulsePlot2 = subplot(5,1,1);
plot(xTimes,pulseFunc,'b');
hold on;
colors = distinguishable_colors(size(evidenceFunc,1));
for i=1:size(evidenceFunc,1)
    plot(xTimes,evidenceFunc(i,:),'Color',colors(i,:));
end
ylim([0 max(1,2*max(pulseFunc))]);
legend('Input to cluster 1','Evidence')
xlabel('Time (seconds)');
ylabel('Input');
title('External Input');

seqPlot = subplot(5,1,2:5);
imagesc(frMatrix);
axisXLabel = num2str(str2num(get(gca,'xticklabel'))*binSize);
set(gca,'xticklabel',axisXLabel);
xlabel('Time (seconds)');
ylabel('Neuron #');
seqPos = get(seqPlot,'Position');
c=colorbar;
cpos = get(c,'Position');
set(c,'Position',[seqPos(1)+seqPos(3)+ 0.01 cpos(2:4)])
ylabel(c,'Firing Rate (Hz)');
end

function plotDebug(em,frMatrix)

if em.maxFR
    figure('Name','Max FR Plot');
    xTimes = linspace(0,em.binSize*em.nTimeBins,em.nTimeBins);

    colors = distinguishable_colors(em.nSequences);

    for i=1:em.nSequences
        hold on;
        plot(xTimes,max(frMatrix((i-1)*em.nClusters+1:i*em.nClusters,:)),'color',colors(i,:));

    end
    legend(cellfun(@(x) ['Sequence ',num2str(x)],num2cell(1:em.nSequences),'UniformOutput',false));
    xlabel('Time (sec)');
    ylabel('Firing Rate (Hz)');
    title('Max FR vs. Time');
end

if em.seqSpeed
    figure('Name','Sequence Speed Plot');
    xTimes = linspace(0,em.binSize*em.nTimeBins,em.nTimeBins);

    colors = distinguishable_colors(em.nSequences);

    for i=1:em.nSequences
        hold on;
        [~,ind] = max(frMatrix((i-1)*em.nClusters+1:i*em.nClusters,:));
        plot(xTimes,ind,'color',colors(i,:));
    end
    
    ylim([0 em.nClusters]);
    legend(cellfun(@(x) ['Sequence ',num2str(x)],num2cell(1:em.nSequences),'UniformOutput',false));
    xlabel('Time (sec)');
    ylabel('Neuron Index of Max FR');
end

if em.inputPlot
    for i=1:length(em.inputPlotNeurons)
        figure('Name',['Neuron ',num2str(em.inputPlotNeurons(i)),' Inputs Over Time']);
        
        xTimes = linspace(0,em.binSize*em.nTimeBins,em.nTimeBins);
        
        inputNames = fieldnames(em.inputs);
        
        colors = distinguishable_colors(length(inputNames));
        
        for j=1:length(inputNames)
            hold on;
            plot(xTimes,em.inputs.(inputNames{j})(em.inputPlotNeurons(i),:),'Color',colors(j,:));
        end
        
        legend(inputNames);
        
        xlabel('Time (sec)');
        ylabel('Input Amplitude');
        title(['Neuron ',num2str(em.inputPlotNeurons(i)),' Inputs Over Time']);
    end
end
end

function newRout = updateRout(oldRout, networkFR, weights, tau, deltat)
%function to update FR of rout
    
input = -oldRout + sum(networkFR'.*weights); %decay + sum(ri*Wij)

newRout = oldRout + input*deltat/tau;
newRout = rectify(newRout);
end
