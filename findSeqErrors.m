function varargout = findSeqErrors(varargin)
%
% [isSeqError, seqErrorType] = findSeqErrors(typedSeq,targSeq,[printToScreen,useRegExp])
%
% Required Inputs:
%   typedSeq - numeric or string vector of length M with keypress IDs (i.e.- '42314423144231442314423')
%   targSeq - numeric or string vector of length N with target sequence (i.e. - '42314')
%
% Optional Inputs:
%   printToScreen - BOOLEAN switch. If TRUE, information about where the errors occur in the keypress vector is printed to the screen. This is
%                   useful for trouble-shooting purposes. (Default is FALSE)
%   useRegExp - BOOLEAN switch. If TRUE, forces use of method based on regular expressions instead of default method based on Smith-Waterman algorithm. (Default is FALSE).
%
% Required Outputs
%   isSeqError - BOOLEAN vector indicating which keypresses are sequences errors or not (i.e. - correct keypresses). By default, this function
%                will determine errors using the Smith-Waterman algorithm (swalign; requires Bioinformatics toolbox). If swalign
%                is not available (or if optional input useRegExp is set to TRUE), an alternate method based on regular expressions is used.
%
% Optional Outputs [NOT IMPLEMENTED YET]
%   seqErrorType - vector containing sequence error type. Coded as follows:
%                   0 = No Error
%                   1 = Insertion Error (i.e. - an extra keypress is added in the middle of a sequence. duplication error is a special case of this)
%                   2 = Deletion Error (i.e. - a key is ommitted)
%                   3 = Substitution Error (i.e. - the wrong key is pressed)
%                   4 = Transposition Error (i.e. - order of two keys is swapped)
%--------------------------------------------------------------------------
% v1.0 23 Feb 2023 by Ethan R Buch
%   Keypress observation padding added at both ends to eliminate edge effects.
%
% v2.0 27 May 2024 by Ethan R Buch
%   Added option of using Smith-Waterman algorithm for local alignment of two sequences
%   (swalign). This is now the default algorithm. Also added seqErrorType
%   output when using swalign.

vFun = 'v2.0';
%disp(['Running ' mfilename('fullpath') ' ' vFun]);

%% I/O house-keeping
if nargin < 2
    error('This function requires at least 2 inputs: typedSeq and targSeq.')
else
    typedSeq = varargin{1};
    targSeq = varargin{2};
    if nargin > 2
        printToScreen = varargin{3};
    else
        printToScreen = false;
    end
    if nargin > 3
        useRegExp = varargin{4};
    else
        useRegExp = false;
    end
end
if isnumeric(typedSeq)
    if exist('swalign')==0 || useRegExp %Regular Expression method requires input to be string row vector
        typedSeq = num2str(typedSeq(:))'; %Convert numeric array to string row vector
    end
else
    if exist('swalign')~=0 && ~useRegExp %S-W Align method requires input to be numeric row vector
        typedSeq = str2num(typedSeq(:))'; %Convert string array to numeric row vector
    end
end
typedSeq = typedSeq(:)'; %Make sure input array is row vector
if isnumeric(targSeq)
    if exist('swalign')==0 || useRegExp %Regular Expression method requires input to be string row vector
        targSeq = num2str(targSeq(:))'; %Convert numeric array to string row vector
    end
else
    if exist('swalign')~=0 && ~useRegExp  %S-W Align method requires input to be numeric row vector
        targSeq = str2num(targSeq(:))'; %Convert string array to numeric row vector
    end
end
targSeq = targSeq(:)'; %Make sure input array is row vector

%% Perform target sequence error check
nKP = length(typedSeq); %Number of total keypresses
nTargSeq = length(targSeq); %Length of a single target sequence iteration

if exist('swalign')==0 || useRegExp %Use method based on Regular Expressions instead of SWALIGN
    seqErrorType = cell(size(typedSeq));
    prePad = targSeq; %Create pre-pad (should just be target sequence)
    postPad = targSeq; %Create post-pad (initialize as target sequence). Run loop to compare similarity of all circular shifts against last 5 keypresses to get best estimate. This should handle situations where insertion/deletion errors could cause a shift in the beginning of correct sequence iterations
    if nKP>=nTargSeq %Check to make sure there are at least enough keypresses to complete a full sequence before proceeding
        hamSim = 0; %Initialize Hamming similarity (i.e. - number of matched ordinal keypresses)
        for jS = 0:nTargSeq-1 %Loop through all circular shifts of target sequence
            curSim = sum(typedSeq(end-nTargSeq+1:end)==circshift(targSeq,-jS));
            if curSim > hamSim
                hamSim = curSim;
                postPad = circshift(targSeq,-jS);
            end
        end
    end
    testID = [prePad typedSeq postPad]; %Pad keypress vector to eliminate edge effects by adding correct sequence to beginning and end (rotated appropriately given the length)
    checkKPmat = repmat(testID,nTargSeq,1); %Replicate KP sequence into nTargSeq x nKP matrix so that all circular rotations of the target sequence can be compared against keypress vector
    for jS = 0:nTargSeq-1 %Loop through all possible circular rotations of target sequence
        iCS = regexp(testID,circshift(targSeq,-jS)); %Locate all instances of complete (rotated) target sequence in keypress vector
        for kCS = iCS
            checkKPmat(jS+1,kCS:kCS+nTargSeq-1) = '0'; %Replace keypress ID in matrix with '0' for all keypressed associated with a complete correct sequence
        end
    end
    iCorrSeqKP = sum(checkKPmat(:,nTargSeq+1:end-nTargSeq)=='0',1)>0; %Remove padding from checkKPmat and convert to correct/incorrect BOOLEAN matrix. If KP is correct for any rotation than it is considered correct (this does not mark any keypresses associated with a deletion error).
    checkKP = typedSeq; %Initialize vector used to display where keypress locations are
    checkKP(iCorrSeqKP) = '0';
    isSeqError = checkKP~='0';
    varargout(1) = {isSeqError};

    %%  ADD section later that classifies ERROR TYPE

    %% Add sequence error type to output if requested
    if nargout == 2
        seqErrorType(~isSeqError) = {'No Error'};
        seqErrorType(isSeqError) = {'Unknown'};
        disp('Sequence Error Type has not yet been implemented for the Regular Expressions method. Returning "Unknown" for all errors.')
        varargout(2) = {seqErrorType};
    end

else
    [~, seqAlign, ~] = swalign(...
        int2aa(typedSeq), ...
        int2aa(repmat(targSeq,1,ceil(length(typedSeq)/length(targSeq)))), ...
        'ScoringMatrix', 2.*eye(max(targSeq), max(targSeq))-1, 'GapOpen', 1);
    alignTyped = seqAlign(1,:);
    seqMatch = seqAlign(2,:);
    alignTarg = seqAlign(3,:);
    nAlign = length(alignTyped);
    isSeqError = seqMatch == ' ';

    seqErrorType = cell(size(seqMatch));
    seqErrorType(~isSeqError) = {'No Error'};
    seqErrorType(isSeqError) = {'Substitution'}; %Start with default that errors are simple substitution errors, then go back and correct


    %Next assign deletion/insertion simple error types
    isDeletion = alignTyped == '-';
    iRem = find(isDeletion);
    isSeqError(iRem+1) = true; %For "deletions" we mark the next KP as an error since we can't mark KPs that didn't happen (although we could revisit this decision)
    seqErrorType([iRem iRem+1]) = {'Deletion'};
    isInsertion = alignTarg == '-';
    seqErrorType(isInsertion) = {'Insertion'};

    %Next deal with compound errors (i.e. - transpositions or "combination" errors)
    if sum(isDeletion)>0 && sum(isInsertion)>0 && nAlign>nKP
        for jD = iRem(:)' %Find all sample indices in typed sequence where deletions holders are marked
            if jD>=3 && seqMatch(jD-1)=='|' && alignTarg(jD-2)=='-'
                if alignTarg(jD) == alignTyped(jD-2)
                    seqErrorType(jD-2:jD) = {'Transposition'};
                end 
            end
            if jD<=nAlign-2 && seqMatch(jD+1)=='|' && alignTarg(jD+2)=='-'
                if alignTarg(jD) == alignTyped(jD+2)
                    seqErrorType(jD:jD+2) = {'Transposition'};
                end
            end
        end
    end

    %Now remove "deletion" holders so that output lengths match typed KP sequence input array
    isSeqError(iRem) = [];
    varargout(1) = {isSeqError};

    %Add sequence error type to output if requested
    if nargout == 2
        seqErrorType(iRem) = [];
        varargout(2) = {seqErrorType};
    end
end

%% Print results to screen if requested
if exist('swalign')==0 || useRegExp && printToScreen
    disp('Target sequence:');
    disp(targSeq);
    disp('Pre-pad:')
    disp(prePad);
    disp('Post-pad:')
    disp(postPad);
    disp('Keypress observations:')
    disp(typedSeq);
    disp('Target sequence ERROR locations:');
    disp(checkKP);
    disp([num2str(sum(isSeqError)) ' of ' num2str(nKP) ' key presses for this trial were determined to be sequence errors.']);
    if nKP==sum(isSeqError)
        disp('No correct keypresses were found for this trial.');
    end
elseif exist('swalign')~=0 && ~useRegExp && printToScreen
    disp('Target sequence:');
    disp(num2str(targSeq')');
    disp('Keypress observations:')
    disp(num2str(typedSeq')');
    checkKP = typedSeq; checkKP(~isSeqError) = 0;
    disp('Target sequence ERROR locations:');
    disp(num2str(checkKP')');
    disp([num2str(sum(isSeqError)) ' of ' num2str(nKP) ' key presses for this trial were determined to be sequence errors.']);
    if nKP==sum(isSeqError)
        disp('No correct keypresses were found for this trial.');
    end
end

















