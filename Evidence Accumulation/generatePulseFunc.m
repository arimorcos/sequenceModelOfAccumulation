function pulseFunc = generatePulseFunc(start, duration, amplitude, nTimeBins, binSize)
%function to generate an 1 x nTimeBins row vector describing a pulse with
%max amplitude defined by amplitude which begins at start (in seconds) and
%lasts for duration seconds. binSize specifies the length of one bin in
%seconds.
%
%Can create multiple pulses if start is a vector of start times
%Duration and amplitude can also be vectors to specify individual durations
%and amplitudes for each pulse. If duration and amplitude are scalar and
%start is a vector, the same pulse duration/amplitude will be used for
%every pulse. If start and duration/amplitude are vectors of different
%sizes, an error will be thrown

if (length(duration) > 1 || length(amplitude) > 1) &&... %throw error if duration/amp contains multiple values but not the same amount as start
        ((length(start) ~= length(duration) || length(start) ~= length(amplitude)) ||...
        length(duration) ~= length(amplitude))
    error('start and duration/amplitude must have same lengths');
end

if length(start) > 1 && length(duration) == 1 %if duration is only one value, replicate that value into vector matching size of start
    duration = duration*ones(size(start));
end

if length(start) > 1 && length(amplitude) == 1 %if amplitude is only one value, replicate that value into vector matching size of start
    amplitude = amplitude*ones(size(start));
end

startBin = start/binSize; %convert the start time in seconds to the start bin
durationBins = duration/binSize; %convert the duration in seconds to duration in bins
startBin = max(startBin(~isnan(startBin)),1); %set startBin to 1 if 0

pulseFunc = zeros(1,nTimeBins); %initialize row vector 


warning('off','MATLAB:colon:nonIntegerIndex');
for i=1:length(startBin) %for each pulse
    if isnan(startBin(i)) %skip if nan
        continue;
    end
    pulseFunc(1,startBin(i):startBin(i)+durationBins(i)) = amplitude(i); %set pulse
end

end