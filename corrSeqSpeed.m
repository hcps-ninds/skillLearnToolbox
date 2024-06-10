function varargout = corrSeqSpeed(varargin)
%
% [t,seqSpeed[,tranTime]] = corrSeqSpeed(timeStamps,typedSeq,targSeq[,tSpan,dt,interpMeth,circShift,tranOut])
%
% Required Inputs:
%   timeStamps - numeric vector of length M with keypress timestamps in seconds or ms [NOTE:
%                timeStamps, tSpan and dt must all be specified in the same units]
%   typedSeq - numeric or string vector of length M with keypress IDs (i.e.- '42314423144231442314423')
%   targSeq - numeric or string vector of length N with target sequence (i.e. - '42314')
%
% Optional Inputs
%   tSpan - bounding limits for t in seconds or ms (i.e. - [0 10000]). ALL timestamps MUST fall within these limits. DEFAULT is [timeStamp(1)
%           timeStamp(end)] [NOTE: timeStamps, tSpan and dt must all be specified in the same units]
%   dt -    delta t in seconds or ms (determines sampling rate of output; DEFAULT = 1e-4 of
%           tSpan) [NOTE: timeStamps, tSpan and dt must all be specified in the same units]
%   interpMeth - interpolation method (type "help interp1" for all options. DEFAULT is 'makima' cubic hermite spline interpolation)
%   circShift - BOOLEAN switch.  DEFAULT is TRUE (looks for all circular shifts of targSeq in typedSeq). 
%               Not recommended if long pauses between sequence repetitions are observed.
%   tranOut - BOOLEAN switch.  DEFAULT is TRUE (outputs keypress transition speed for all individual transitions in correct sequences).
%
% Required Outputs
%   t - time vector (in ms)
%   seqSpeed - 1 x t vector of correct sequence speed (keypresses/s)
%   
% Optional Outputs
%   tranTime - N x t vector of individual transition times observed
%   during a correct sequence (in ms)
%
%--------------------------------------------------------------------------
% v1.0 21 August 2019 by Ethan R Buch
%
% v1.1 23 August 2019 by Ethan R Buch
%   Changed extrapolation behavior.  Instead of appending NaNs to all samples before first and after last completed correct sequences,
%   the speed value for first complete correct sequence is replicated across all prior samples, and the speed value for last complete  
%   correct sequence is replicated across all subsequent samples.
%
% v1.2 4 September 2019 by Ethan R Buch
%   Added check to return NaNs if no correct sequences found
%
% v1.3 16 July 2021 by Ethan R Buch
%   Improved measure of instantaneous correct sequence speed. Removed
%   cyclical bias caused by different mix of transitions for different circular
%   shifts
%
%v1.4 20 December 2022 by Ethan R Buch
%   Updated instantaneous correct sequence speed. Now calculated with
%   wrap-around transition so that all circular shifts have same mix of
%   transitions.
%
% v1.5 16 March 2023 by Ethan R Buch
%   Updated sequence speed measure to improve generalization to different
%   length sequences. Changed sequence speed units from correct sequences/s
%   to keypresses/s. Added optional individual transition-speed output.
%
% v1.6 27 May 2024 by Ethan R Buch
%   Updated handling of input time data variables (timestamps, tspan and dt)
%   to accomodate different units (ms or s)
%
% v1.7 6 June 2024 by Ethan R Buch
%   Updated handling of typed/target sequence comparison to be compatible with cases where there is only a single iteration of the target sequence 
%
% v1.8 7 June 2024 by Ethan R Buch
%   1) Added check to tspan input parameter to make sure that it brackets all
%      timestamp data. If it doesn't, tspan is reset to default.
%   2) Switched to using "strfind" instead of "regexp" so that speed profile is updated correctly when consecutive
%      correct sequences occur (which are now overlapping after adding wrap-around
%      transition in v1.4 update)

vFun = 'v1.8';
%disp(['Running ' mfilename('fullpath') ' ' vFun]);

%% Input parsing and I/O checks

if nargin < 3
    error('A minimum of three inputs (timeStamps, typedSeq and targSeq) are required.');
end
if nargin >= 3
    timeStamps = varargin{1}(:)';
    if ~isnumeric(timeStamps)
        error('timeStamps input must be a numeric array.');
    end
    %Make sure timeStamps array is in ascending order
    [timeStamps,iA] = sort(timeStamps,'ascend');
    typedSeq = varargin{2}(:)';
    if isnumeric(typedSeq) %Convert numeric array to string
        typedSeq = num2str(typedSeq(:))';
    end
    targSeq = varargin{3}(:)';
    if isnumeric(targSeq) %Convert numeric array to string
        targSeq = num2str(targSeq(:))';
    end
end
if length(timeStamps) ~= length(typedSeq)
    error('timeStamps and typedSeq input vectors must be the same length.');
end
typedSeq = typedSeq(iA); %Make sure typed seq array matches ascending ordered timestamp array
    
if nargin > 3
    tSpan = varargin{4};
    if tSpan(1)>timeStamps(1)
        tSpan(1) = timeStamps(1);
        disp('tSpan not inclusive of all timestamp inputs. Resetting lower limit to first timestamp.')
    end
    if tSpan(end)<timeStamps(end)
        tSpan(end) = timeStamps(end);
        disp('tSpan not inclusive of all timestamp inputs. Resetting upper limit to last timestamp.')
    end
else
    tSpan = [timeStamps(1) timeStamps(end)]; %DEFAULT
end
if nargin > 4
    dt = varargin{5};
else
    dt = diff(tSpan)*1e-4; %DEFAULT: 10,000 samples over each 1s of timespan (covers Presentation timestamp resolution)
end
if nargin > 5
    interpMeth = varargin{6};
else
    interpMeth = 'makima'; %DEFAULT
end
if nargin > 6
    circShift = varargin{7};
else
    circShift = true; %DEFAULT
end
if nargout < 3
    if nargin > 7
        tranOut = varargin{8};
    else
        tranOut = false;
    end
else
    tranOut = true; %DEFAULT
end

%% Compute time-resolved correct sequence speed estimate (keypresses/s)
%Store typed and target sequence lengths
M = length(typedSeq);
N = length(targSeq);
if M <= N
    N = N - 1; %Correct transition limit if typed sequence is not longer than one iteration of the target sequence 
end

%Determine time units
if dt < 1e-1
    tUnits = 's';
    tGain = 1;
else
    tUnits = 'ms';
    tGain = 1e3;
end
%disp(['Units of input time variables were determined to be "' tUnits '".']);

% Initialize variables
speedMat = NaN(1,M);
if tranOut, allTT = NaN(N,M); end
t = tSpan(1):dt:tSpan(end); %generate time vector for use in interpolation step and output
seqSpeed = NaN(size(t)); %initialize "seqSpeed" output
if tranOut
    allTT = NaN(N,M);
    tranTime = NaN(N,length(t));
end
if circShift
    allShifts = 0:N-1;
else
    allShifts = 0;
end

anyCorrSeqs = false;
for jS = allShifts %loop through all circular shifts of target sequence
    checkSeq = circshift(targSeq,-jS);
    if length(typedSeq)>length(targSeq)
        checkSeq = [checkSeq checkSeq(1)]; %If more than one sequence iteration then add wrap-around transition so that each shift has same mix of transitions
    end
    iSeqStart = strfind(typedSeq,checkSeq);
    if ~isempty(iSeqStart)
        anyCorrSeqs = true;
    end
    for curSS = iSeqStart %loop through each correct sequence onset index
        speedMat(end+1,:) = NaN; %#ok<AGROW>
        speedMat(end,curSS:curSS+N) = N./(timeStamps(curSS+N) - timeStamps(curSS)).*tGain; %Replicate speed estimate over full length of correct sequence (keypresses/s)
        if tranOut
            allTT(:,:,end+1) = NaN;
            curOrd = circshift(1:N,-jS);
            for k = 1:N
                allTT(curOrd(k),curSS+k,end) = timeStamps(curSS+k) - timeStamps(curSS+k-1);
            end; clear k
        end
            
    end
end
mnSpeed = mean(speedMat,1,'omitnan'); %take mean speed over all shifts for each keypress timepoint
if tranOut, mnAllTT = mean(allTT,3,'omitnan'); end
iObs = isfinite(mnSpeed) & isfinite(timeStamps); %index finite speed values for interpolation step

if anyCorrSeqs && sum(iObs)>0
    seqSpeed = interp1(timeStamps(iObs),mnSpeed(iObs),t,interpMeth,NaN);  %interpolate over full timespan using method "interpMeth" and setting all extrapulated points to NaN
    seqSpeed(t<timeStamps(find(iObs,1,'first'))) = mnSpeed(find(iObs,1,'first')); %Fill-in beginning part with first observed value
    seqSpeed(t>timeStamps(find(iObs,1,'last'))) = mnSpeed(find(iObs,1,'last')); %Fill-in end part with last observed value
    if tranOut %Calculate instantaneous speeds for individual transitions
        for k = 1:N
            jTT = isfinite(mnAllTT(k,:));
            if sum(jTT)>0
                if sum(jTT)>1
                    tranTime(k,:) = interp1(timeStamps(jTT),mnAllTT(k,jTT),t,interpMeth,NaN); %interpolate over full timespan using method "interpMeth" and setting all extrapulated points to NaN
                end
                tranTime(k,t<=timeStamps(find(jTT,1,'first'))) = mnAllTT(k,find(jTT,1,'first')); %Fill-in beginning part with first observed value
                tranTime(k,t>=timeStamps(find(jTT,1,'last'))) = mnAllTT(k,find(jTT,1,'last')); %Fill-in end part with last observed value
            end
        end; clear k
    end
else
    warning('No correct sequences found. Returning all NaNs for sequence speed output.');
end

%% Define outputs
if nargout < 2
  error('At least two outputs (t and seqSpeed) are required.');
else
    varargout(1) = {t};
    varargout(2) = {seqSpeed};
end
if nargout > 2
    varargout(3) = {tranTime};
end
if nargout > 3
    warning('There is a maximum of three outputs (t, seqSpeed and tranTime). All additional specified ouputs will return empty.')
    for iOut = 4:nargout
        varargout(iOut) = cell(0);
    end
end