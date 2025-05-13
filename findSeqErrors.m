function varargout = findSeqErrors(varargin)
%
% [isSeqError[,seqErrorType,align]] = findSeqErrors(typedSeq,targSeq[,method,printToScreen,gapCost])
%
% Required Inputs:
%   typedSeq - numeric or string vector of length M with keypress IDs (e.g.- '4132441324413244132441324413')
%   targSeq - numeric or string vector of length N with target sequence (e.g. - '41324')
%
% Optional Inputs:
%   method - Determines which method is used for error-labeling. Options are:
%                   'regexp' - pattern matching of the target sequence to the typed sequence
%                   'regexpcs' [DEFAULT if Bioinformatics toolbox is UNAVAILABLE] - pattern matching of all possible circular shifts of the target sequence to the typed sequence
%                   'nwalign' [DEFAULT if Bioinformatics toolbox is AVAILABLE] - Needleman-Wunsch global sequence alignment algorithm (nwalign; requires Bioinformatics toolbox)
%                   'swalign' - Smith-Waterman partial sequence alignment algorithm (swalign; requires Bioinformatics toolbox)
%   printToScreen - BOOLEAN switch. If TRUE, information about where the errors occur in the typed sequence and what type of errors they are is printed to the screen. (Default is TRUE)
%   gapOpen - This parameter determines how likelihood of
%             DELETION/INSERTION error versus SUBSITUTION error classifications for "nwalign" or "swalign" methods.
%             (Default is 4)
%
% Required Outputs
%   isSeqError - BOOLEAN vector indicating which keypresses are sequences errors or not (i.e. - correct keypresses).
%
% Optional Outputs [NOT IMPLEMENTED YET]
%   seqErrorType - vector containing sequence error type label. Labels are one of the following:
%                   1) No Error (-)
%                   2) Insertion Error (I; i.e. - an extra keypress is added in the middle of a sequence. DUPLICATION error [R] is a special case of this)
%                   3) Deletion Error (D; i.e. - a keypress is skipped or ommitted)
%                   4) Substitution Error (S; i.e. - the wrong key is pressed. A TRANSPOSITION error [T; i.e. - the order of two consecutive keypresses is reversed] is a special case of this)
%   align - outputs same text strings as printed to screen when
%           printToScreen set to TRUE (automatic if this output is specificed)
%
%--------------------------------------------------------------------------

%% I/O house-keeping
if nargin < 2
    error('This function requires at least 2 inputs: typedSeq and targSeq.')
else
    %Required Input 1 - typedSeq
    typedSeq = varargin{1};
    if isempty(typedSeq)
        error('The "typedSeq" input array is empty.')
    end
    %Required Input 2 - targSeq
    targSeq = varargin{2};
    if isempty(targSeq)
        error('The "targSeq" input array is empty.')
    end
    %Optional Input 3 - method
    if nargin > 2
        method = varargin{3};
        if isempty(method) && exist('nwalign')==0
            method = 'regexpcs';
            disp('Setting method to "regexpcs"');
        elseif isempty(method) && exist('nwalign')~=0
            method = 'nwalign';
            disp('Setting method to "nwalign"');
        end
    elseif nargin<=2 && exist('nwalign')==0
        method = 'regexpcs';
        disp('Setting method to "regexpcs"');
    else
        method = 'nwalign';
        disp('Setting method to "nwalign"');
    end
    %Optional Input 4 - printToScreen
    if nargin > 3
        printToScreen = varargin{4};
        if isempty(printToScreen)
            printToScreen = true;
        end
    else
        printToScreen = true;
    end
    if nargout == 3 && ~printToScreen
        printToScreen = true;
        warning('Switching "printToScreen" to TRUE since "align" output is requested.')
    end
    %Optional Input 5 - gapCost
    if nargin > 4
        gapCost = varargin{5};
    else
        gapCost = 4;
    end
end

%Check for necessary "typedSeq" conversions
if isnumeric(typedSeq)
    if strcmp(method,'regexp') || strcmp(method,'regexpcs') %Regular Expression methods require input to be string row vector
        typedSeq = num2str(typedSeq(:))'; %Convert numeric array to string row vector
    end
else
    if strcmp(method,'nwalign') || strcmp(method,'swalign') %"nwalign" and "swalign" methods require input to be numeric row vector
        typedSeq = str2num(typedSeq(:))'; %Convert string array to numeric row vector
    end
end
typedSeq = typedSeq(:)'; %Make sure input array is row vector

%Check for necessary "targSeq" conversions
if isnumeric(targSeq)
    if strcmp(method,'regexp') || strcmp(method,'regexpcs') %Regular Expression methods require input to be string row vector
        targSeq = num2str(targSeq(:))'; %Convert numeric array to string row vector
    end
else
    if strcmp(method,'nwalign') || strcmp(method,'swalign') %"nwalign" and "swalign" methods require input to be numeric row vector
        targSeq = str2num(targSeq(:))'; %Convert string array to numeric row vector
    end
end
targSeq = targSeq(:)'; %Make sure input array is row vector
seqLen = length(targSeq);

%% Perform target sequence error check
nKP = length(typedSeq); %Number of total keypresses
nTargSeq = length(targSeq); %Length of a single target sequence iteration
perfectSeq = repmat(targSeq,1,ceil(nKP/nTargSeq));
perfectSeq = perfectSeq(1:nKP);


%% Check to make sure there are at least enough keypresses to evaluate one full iteration of the sequence.
if nKP < nTargSeq
    warning('The number of total keypresses is LESS THAN the sequence length. Returning all NaNs.')
    iSeqError = NaN(1,nKP);
    varargout(1) = {iSeqError};
    if nargout == 2
        seqErrorType = repmat({'Unknown'},1,nKP);
        varargout(2) = {seqErrorType};
    end
    if nargout == 3
        if isnumeric(targSeq)
            align.targetSeq = num2str(targSeq(:))';
        else align.targetSeq = targSeq(:)';
        end
        if isnumeric(typedSeq)
            align.typedSeq = num2str(typedSeq(:))';
        else align.typedSeq = typedSeq(:)';
        end
        align.compare = '';
        varargout(3) = {align};
    end
    return %Exit function
end

%% Perform error detection
if strcmp(method,'regexp') %Regular Expression method (no circular shift)
    seqErrorType = cell(size(typedSeq));
    checkKP = typedSeq;
    errorGap = zeros(1,length(typedSeq));
    iCS = regexp(typedSeq,targSeq); %Locate all instances of complete (unrotated) target sequence in keypress vector
    prevCorrSeqStop = 0;
    for kCS = iCS
        corrSeqStart = kCS;
        corrSeqStop = kCS+nTargSeq-1;
        checkKP(corrSeqStart:corrSeqStop) = '-'; %Replace keypress ID in matrix with '0' for all keypressed associated with a complete correct sequence
        if corrSeqStart > prevCorrSeqStop + 1
            if corrSeqStart-1 == nKP
                errorGap(prevCorrSeqStop+1:corrSeqStart-1) = Inf; %If errors happen at the end then we don't really know what the gap is
            else
                errorGap(prevCorrSeqStop+1:corrSeqStart-1) = corrSeqStart - (prevCorrSeqStop + 1); %Keep track of error gap lengths between correct sequences
            end
        end
    end
    iCorrSeqKP = checkKP=='-'; %Convert to correct/incorrect BOOLEAN matrix.
    iSeqError = ~iCorrSeqKP;
    checkKP(iSeqError) = 'E';

    % Add sequence error type to output if requested
    if nargout >= 2
        seqErrorType(iCorrSeqKP) = {'No Error'};
        seqErrorType(iSeqError) = {'Unknown Error'};
        for k = find(iSeqError)
            if errorGap(k) == nTargSeq - 1 %Possible Deletion
                seqErrorType(k) = {'Deletion'};
                checkKP(k) = 'D';
            elseif errorGap(k) == nTargSeq %Possible Substitution
                seqErrorType(k) = {'Substitution'};
                checkKP(k) = 'S';
            elseif errorGap(k) == nTargSeq + 1 %Possible Insertion
                seqErrorType(k) = {'Insertion'};
                checkKP(k) = 'I';
            end
        end
    end
elseif strcmp(method,'regexpcs') % Regular Expression with all circular shifts of target sequence method
    prePad = targSeq; %Create pre-pad (should just be target sequence)
    postPad = targSeq; %Create post-pad (initialize as target sequence). Run loop to compare similarity of all circular shifts against last 5 keypresses to get best estimate. This should handle situations where insertion/deletion errors could cause a shift in the beginning of correct sequence iterations
    nAppendErrors(1:nTargSeq) = 2*nTargSeq; %Initialize errors for testAppend (i.e. - start with all errors so we can minimize)
    for jS = 0:nTargSeq-1 %Loop through all circular shifts of target sequence
        curPostPad = circshift(targSeq,-jS);
        testPostPad = [typedSeq(end-nTargSeq+1:end) curPostPad];
        checkPPmat = repmat(testPostPad,nTargSeq,1);
        for jSS = 0:nTargSeq-1 %Loop through all possible circular rotations of target sequence
            iCS = regexp(testPostPad,circshift(targSeq,-jSS)); %Locate all instances of complete (rotated) target sequence in keypress vector
            for kCS = iCS
                checkPPmat(jSS+1,kCS:kCS+nTargSeq-1) = '-'; %Replace keypress ID in matrix with '-' for all keypressed associated with a complete correct sequence
            end
        end; clear jSS iCS kCS 
        nMatches = sum(checkPPmat=='-',1);
        % First check for obvious correct KPs (i.e. - belongs to at least two overlapping pattern matches)
        iCorrAppend = nMatches >= 2;
        % Now check for non-overlapping pattern matches that indicate a DELETION error (i.e. - "11" patterns in nMatches)
        iDE = regexp(num2str(nMatches(:))','11');
        iCorrAppend(iDE) = true;
        iCorrAppend(iDE+1) = false;
        % Now check for non-overlapping pattern matches that indicate a SUBSTITUTION error (i.e. - "101" patterns in nMatches)
        iSE = regexp(num2str(nMatches(:))','101');
        iCorrAppend([iSE iSE+2]) = true;
        iCorrAppend(iSE+1) = false;
        nAppendErrors(jS+1) = sum(~iCorrAppend);
    end; clear iCorrAppend iDE iSE checkPPmat testPostPad curPostPad jS % Housekeeping
    jPP = find(nAppendErrors==min(nAppendErrors),1,'first') - 1; clear nAppendErrors
    postPad = circshift(targSeq,-jPP); clear jPP
    testID = [prePad typedSeq postPad]; %Pad keypress vector to eliminate edge effects by adding correct sequence to beginning and end (rotated appropriately given the length)
    targetID = [prePad perfectSeq postPad];

    %First pattern match padded typed sequence (testID) to target sequence
    %with no shifts to get error type info from error gap lengths
    seqErrorType = cell(size(testID));
    checkKP = testID;
    errorGap = zeros(1,length(testID));
    iCS = regexp(testID,targSeq); %Locate all instances of complete (unrotated) target sequence in keypress vector
    prevCorrSeqStop = 0;
    for kCS = iCS
        corrSeqStart = kCS;
        corrSeqStop = kCS+nTargSeq-1;
        checkKP(corrSeqStart:corrSeqStop) = '-'; %Replace keypress ID in matrix with '0' for all keypressed associated with a complete correct sequence
        if corrSeqStart > prevCorrSeqStop + 1
            if corrSeqStart-1 == nTargSeq+nKP %Length of pre-pad + typed sequence
                errorGap(prevCorrSeqStop+1:corrSeqStart-1) = Inf; %If errors happen at the end then we don't really know what the gap is
            else
                errorGap(prevCorrSeqStop+1:corrSeqStart-1) = corrSeqStart - (prevCorrSeqStop + 1); %Keep track of error gap lengths between correct sequences
            end
        end
        prevCorrSeqStop = corrSeqStop;
    end
    iCorrSeqKP = checkKP=='-'; %Convert to correct/incorrect BOOLEAN matrix.
    iSeqError = ~iCorrSeqKP;
    checkKP(iSeqError) = 'E';

    % Add sequence error type to output if requested
    if nargout >= 2
        seqErrorType(iCorrSeqKP) = {'No Error'};
        seqErrorType(iSeqError) = {'Unknown Error'};
        for k = find(iSeqError)
            if errorGap(k) == nTargSeq - 1 %Possible Deletion
                seqErrorType(k) = {'Deletion'};
                checkKP(k) = 'D';
            elseif errorGap(k) == nTargSeq %Possible Substitution
                seqErrorType(k) = {'Substitution'};
                checkKP(k) = 'S';
            elseif errorGap(k) == nTargSeq + 1 %Possible Insertion
                seqErrorType(k) = {'Insertion'};
                checkKP(k) = 'I';
            end
        end
    end

    checkKPmat = repmat(testID,nTargSeq,1); %Replicate KP sequence into nTargSeq x nKP matrix so that all circular rotations of the target sequence can be compared against keypress vector
    for jS = 0:nTargSeq-1 %Loop through all possible circular rotations of target sequence
        iCS = regexp(testID,circshift(targSeq,-jS)); %Locate all instances of complete (rotated) target sequence in keypress vector
        for kCS = iCS
            checkKPmat(jS+1,kCS:kCS+nTargSeq-1) = '-'; %Replace keypress ID in matrix with '-' for all keypressed associated with a complete correct sequence
        end
    end
    nMatches = sum(checkKPmat=='-',1);
    % First check for obvious correct KPs (i.e. - belongs to at least two overlapping pattern matches)
    iCorrSeqKP = nMatches >= 2;
    % Now check for non-overlapping pattern matches that indicate a DELETION or DUPLICATION error (i.e. - "11" patterns in nMatches)
    iE = regexp(num2str(nMatches(:))','11');
    for cE = iE(:)'
        iCorrSeqKP(cE:cE+1) = [true false];
        if errorGap(cE) == nTargSeq-1
            checkKP(cE:cE+1) = '-D';
            seqErrorType(cE:cE+1) = {'No Error','Deletion'};
        elseif errorGap(cE) == nTargSeq+1
            checkKP(cE:cE+1) = '-R';
            seqErrorType(cE:cE+1) = {'No Error','Duplication'};
        else
            checkKP(cE:cE+1) = '-E';
            seqErrorType(cE:cE+1) = {'No Error','Unknown Error'};
        end
    end
    % Now check for non-overlapping pattern matches that indicate a SUBSTITUTION or INSERTION error (i.e. - "101" patterns in nMatches)
    iE = regexp(num2str(nMatches(:))','101');
    for cE = iE(:)'
        iCorrSeqKP(cE:cE+2) = [true false true];
        if errorGap(cE+1) == nTargSeq
            checkKP(cE:cE+2) = '-S-';
            seqErrorType(cE:cE+2) = {'No Error','Substitution','No Error'};
        elseif errorGap(cE+1) == nTargSeq+1
            checkKP(cE:cE+2) = '-I-';
            seqErrorType(cE:cE+2) = {'No Error','Insertion','No Error'};
        else
            checkKP(cE:cE+2) = '-E-';
            seqErrorType(cE:cE+2) = {'No Error','Unknown Error','No Error'};
        end
    end
    % Now check for consecutive SUBSTITUTION or possible TRANSPOSITION errors (i.e. - "1001" patterns in nMatches)
    iE = regexp(num2str(nMatches(:))','1001');
    for cE = iE(:)'
        iCorrSeqKP(cE:cE+3) = [true false false true];
        if targetID(cE+1)==testID(cE+2) && targetID(cE+2)==testID(cE+1)
            checkKP(cE:cE+3) = '-TT-';
            seqErrorType(cE:cE+3) = {'No Error','Transposition','Transposition','No Error'};
        else
            checkKP(cE:cE+3) = '-SS-';
            seqErrorType(cE:cE+3) = {'No Error','Substitution','Substitution','No Error'};
        end
    end
    iC = regexp(num2str(nMatches(:))','21');
    for cC = iC(:)'
        iCorrSeqKP(cC:cC+1) = [true true];
        checkKP(cC:cC+1) = '--';
        seqErrorType(cC:cC+1) = {'No Error','No Error'};
    end    
    %Now perform final update of iSeqError and checkKP
    iSeqError = ~iCorrSeqKP;
    checkKP(iCorrSeqKP) = '-';
    checkKP(strcmp(checkKP,'E')) = 'S'; %Convert all remaining "Unknown Errors" to substitutions
    seqErrorType(iCorrSeqKP) = {'No Error'};
    seqErrorType(strcmp(seqErrorType,'Unknown Error')) = {'Substitution'}; %Convert all remaining "Unknown Errors" to substitutions

    %Now remove padding from outputs    
    iCorrSeqKP = iCorrSeqKP(nTargSeq+1:end-nTargSeq);
    iSeqError = iSeqError(nTargSeq+1:end-nTargSeq);
    checkKP = checkKP(nTargSeq+1:end-nTargSeq);
    seqErrorType = seqErrorType(nTargSeq+1:end-nTargSeq);

elseif strcmp(method,'nwalign') || strcmp(method,'swalign') %Needleman-Wunsch/Smith-Waterman methods
    scoringMat = -abs(unique(targSeq) - unique(targSeq)') + eye(max(unique(targSeq)),max(unique(targSeq))); % scoringMat = 2.*eye(max(targSeq), max(targSeq))-1;

    % 1) Perform global align to handle DELETION and INSERTION errors
    seqObs = int2aa(typedSeq); %Convert typed integer sequence to amino acid sequence
    seqComp = int2aa(perfectSeq); %Convert comparison sequence to amino acid sequence of same length as typed sequence
    if strcmp(method,'nwalign')
        [~, seqAlign, ~] = nwalign(... %1st Global align
            seqObs, ...
            seqComp, ...
            'ScoringMatrix', scoringMat, 'GapOpen', 4);
    elseif strcmp(method,'swalign')
        [~, seqAlign, ~] = swalign(... %1st Global align
            seqObs, ...
            seqComp, ...
            'ScoringMatrix', scoringMat, 'GapOpen', 4);
    end
    alignTyped = seqAlign(1,:); %Row 1 is the observed/typed sequence aligned to the target sequence
    alignTarg = seqAlign(3,:); %Row 3 is the target sequence aligned to the observed/typed sequence
    nAlign = length(alignTyped); %The number of alignments might differ from the number of keypresses if insertions/deletions are present

    % Check for erroneous DELETION error added at the end. This can happen if INSERTION error is marked earlier in sequence
    % but deletions are not possible after the last typed KP
    while alignTyped(end)=='-'
        alignTyped(end) = [];
        alignTarg(end) = [];
        nAlign = length(alignTyped);
    end

    % Check for erroneous INSERTION error added at the end. This can happen if DELETION error is marked earlier in sequence
    if alignTarg(end)=='-' && strcmp(alignTyped(end-nTargSeq+1:end),int2aa(targSeq))
        iStartInsert = length(alignTarg);
        while iStartInsert > 1 && alignTarg(iStartInsert) == '-'
            iStartInsert = iStartInsert - 1;
        end
        if alignTarg(iStartInsert)~='-', iStartInsert = iStartInsert+1; end
        alignTarg(iStartInsert:end) = alignTyped(iStartInsert:end);
    end

    % Assign initial DELETION error types
    iDeletion = find(alignTyped == '-'); %Any deletions (i.e. - sequence elements that were skipped over) will be marked as "-" in alignTyped (i.e. - row 1 of seqAlign)
    iOmit = [];
    % Loop through each DELETION error and make adjustments if needed
    for curDel =  iDeletion(:)'
        if ismember(curDel,iOmit), continue; end %Check against running list of errors already corrected
        %Deal with special case of repeated keypresses (i.e. - double-tap or
        %triple-tap). Make sure deletions are marked for later KPs.
        %Eventually we can incorporate KTT anomoly scores to mark most
        %likely deletion.
        curKPid = alignTarg(curDel);
        iBack = 0;
        while curDel+iBack > 1 && alignTarg(curDel+iBack) == curKPid
            iBack = iBack - 1;
        end
        if alignTarg(curDel+iBack) ~= curKPid, iBack = iBack + 1; end %Move back to the right since it will overshoot by 1
        iFwd = 0;
        while curDel+iFwd<nAlign && alignTarg(curDel+iFwd) == curKPid
            iFwd = iFwd + 1;
        end
        if alignTarg(curDel+iFwd) ~= curKPid, iFwd = iFwd - 1; end %Move back to the left since it will overshoot by 1
        nSeg = length(curDel+iBack:curDel+iFwd);
        nDel = sum(alignTyped(curDel+iBack:curDel+iFwd) == '-');
        segCorr = [repmat(curKPid,1,nSeg-nDel) repmat('-',1,nDel)];
        alignTyped(curDel+iBack:curDel+iFwd) = segCorr;
        iOmit = [iOmit curDel+iBack:curDel+iFwd];
    end % Adjustments finished

    % Reassign DELETION error types after adjustments
    % Check for erroneous DELETION error added at the end. This can happen if INSERTION error is marked earlier in sequence
    % but deletions are not possible after the last typed KP
    while alignTyped(end)=='-'
        alignTyped(end) = [];
        alignTarg(end) = [];
        nAlign = length(alignTyped);
    end

    seqObsCorr = alignTyped; %Initialize corrected observed sequence array
    seqCompCorr = alignTarg; %Initialize corrected target/comparison sequence array

    %Correct typed sequence (seqObsCorr variable) by adding DELETIONs (i.e. - missing KPs)
    iDeletion = find(alignTyped == '-'); %Reset array now that adjustments have been made.
    seqObsCorr(iDeletion) = seqCompCorr(iDeletion); %Use aligned comparison sequence to add in missing keypress IDs

    seqErrorType = repmat({'No Error'},1,nAlign); %Initialize sequence error type and size to number of alignments
    iSeqError = false(1,nAlign); %Initialize sequence error boolean vector and size to number of alignments

    iSeqError(iDeletion+1) = true; %For "deletions" we mark the next KP as an error since we can't mark KPs that didn't happen (although we could revisit this decision)
    seqErrorType(iDeletion+1) = {'Deletion'};
    seqErrorType(iDeletion) = {'Remove From List'};

    % Assign INSERTION error types and correct
    iInsertion = find(alignTarg == '-');
    iSeqError(iInsertion) = true; %Mark as an error

    %Correct target sequence (seqCompCorr variable) by adding INSERTIONS (i.e. - extra KPs)
    seqCompCorr(iInsertion) = seqObsCorr(iInsertion); %Add extra KPs to comparison sequence so that insertions match

    %Loop through each insertion error and make adjustments if needed
    if ~isempty(iInsertion)
        for jA = iInsertion(:)'
            isDup = false;
            if jA < nAlign && alignTyped(jA) == alignTyped(jA+1) %Check for DUPLICATIONS (special case of insertion error) and make sure 2nd one is marked as an error if found
                iSeqError(jA:jA+1) = [false true]; %Reassign error to 2nd KP in duplication
                seqErrorType(jA:jA+1) = {'No Error', 'Duplication'};
                isDup = true;
            end
            if jA > 1 && alignTyped(jA-1) == alignTyped(jA)
                iSeqError(jA-1:jA) = [false true]; %Reassign error to 2nd KP in duplication
                seqErrorType(jA-1:jA) = {'No Error', 'Duplication'};
                isDup = true;
            end
            if ~isDup
                seqErrorType(jA) = {'Insertion'};
            end
        end
    end

    % 2) Now compare corrected observed and target keypress lists for SUBSTITUTION errors
    % perfectSeq = repmat(targSeq,1,ceil(nKP/nTargSeq));
    % perfectSeq = perfectSeq(1:nKP);
    iSubstitution = find(seqObsCorr~=seqCompCorr);
    seqErrorType(iSubstitution) = {'Substitution'}; %Label errors as SUBSTITUTIONS, then go back and check for special cases
    iSeqError(iSubstitution) = true;

    % Look for consecutive substitution errors and label special case of TRANSPOSITIONS (reversing order of two keypresses) if found
    if length(iSubstitution)>1
        for jS = 1:length(iSubstitution)-1
            if iSubstitution(jS+1)-iSubstitution(jS)==1 && strcmp(seqCompCorr(iSubstitution(jS)),seqObsCorr(iSubstitution(jS+1))) && strcmp(seqCompCorr(iSubstitution(jS+1)),seqObsCorr(iSubstitution(jS)))
                seqErrorType(iSubstitution(jS:jS+1)) = {'Transposition'};
            end
        end
    end
end

%% Print results to screen if requested
if printToScreen
    if strcmp(method,'nwalign') || strcmp(method,'swalign')
        if isnumeric(alignTarg)
            targetSeqPrint = num2str(alignTarg')';
        else
            targetSeqPrint = alignTarg;
            targetSeqPrint(targetSeqPrint=='A') = num2str(aa2int('A'));
            targetSeqPrint(targetSeqPrint=='D') = num2str(aa2int('D'));
            targetSeqPrint(targetSeqPrint=='N') = num2str(aa2int('N'));
            targetSeqPrint(targetSeqPrint=='R') = num2str(aa2int('R'));
        end
        if isnumeric(alignTyped)
            typedSeqPrint = num2str(alignTyped')';
        else
            typedSeqPrint = alignTyped;
            typedSeqPrint(typedSeqPrint=='A') = num2str(aa2int('A'));
            typedSeqPrint(typedSeqPrint=='D') = num2str(aa2int('D'));
            typedSeqPrint(typedSeqPrint=='N') = num2str(aa2int('N'));
            typedSeqPrint(typedSeqPrint=='R') = num2str(aa2int('R'));
        end
        checkKP = repmat('-',1,nAlign);
        checkKP(strcmp(seqErrorType,'Substitution')) = 'S';
        checkKP(strcmp(seqErrorType,'Transposition')) = 'T';
        checkKP(strcmp(seqErrorType,'Remove From List')) = 'D';
        checkKP(strcmp(seqErrorType,'Deletion')) = '|'; %This is where the deletion error will be marked in the keypress list
        checkKP(strcmp(seqErrorType,'Insertion')) = 'I';
        checkKP(strcmp(seqErrorType,'Duplication')) = 'R';
        disp([targetSeqPrint ' (Target sequence)']);
        disp([typedSeqPrint ' (Typed sequence)']);
        disp([checkKP ' (Error types and locations)']);
        disp([...
            num2str(sum(iSeqError(~strcmp(seqErrorType,'Remove From List')))) ...
            ' of ' ...
            num2str(nKP) ...
            ' key presses for this trial are associated with sequence errors.']);
    elseif strcmp(method,'regexpcs') % || strcmp(method,'regexp')
        %Correct typed sequence for DELETION errors before printing to
        %screen
        corrTypedSeq = typedSeq;
        corrCheckKP = checkKP;
        iDE = find(checkKP=='D');
        for jD = 1:length(iDE)
            curD = iDE(jD);
            corrTypedSeq = [corrTypedSeq(1:curD-1) '-' corrTypedSeq(curD:end)];
            corrCheckKP = [corrCheckKP(1:curD) '|' corrCheckKP(curD+1:end)];
            if jD<length(iDE)
                iDE(jD+1:end) = iDE(jD+1:end) + 1;
            end
        end
        typedSeqPrint = corrTypedSeq;
        checkKP = corrCheckKP;

        %Correct perfect target sequence for INSERTION errors before printing to
        %screen
        corrPerfectSeq = perfectSeq;
        iIE = find(checkKP=='I' | checkKP=='R');
        for jI = 1:length(iIE)
            curI = iIE(jI);
            corrPerfectSeq = [corrPerfectSeq(1:curI-1) '-' corrPerfectSeq(curI:end)];
            if jI<length(iIE)
                iIE(jI+1:end) = iIE(jI+1:end) + 1;
            end
            if curI>1 && (corrTypedSeq(curI) == corrTypedSeq(curI-1))
                checkKP(curI) = 'R';
                seqErrorType(curI) = {'Duplication'};
            end
        end
        targetSeqPrint = corrPerfectSeq;
        %% ================================================================
        %Now attempt to label any remaining Unknown errors that can
        %misalign arrays (NOTE: THIS ROUTINE CAN BE IMPROVED TO ACCOMODATE
        %MORE CASES)
        iUE = find(corrCheckKP=='E');
        if ~isempty(iUE)
            if length(corrCheckKP)==length(typedSeqPrint) && length(corrCheckKP)>length(targetSeqPrint) %Unknown Error is an INSERTION error
                corrCheckKP(iUE) = 'I';
                targetSeqPrint = [targetSeqPrint(1:iUE-1) '-' targetSeqPrint(iUE:end)];
            elseif length(corrCheckKP)==length(targetSeqPrint) && length(corrCheckKP)>length(typedSeqPrint) %Unknown Error is an DELETION error
                corrCheckKP(iUE:iUE+1) = 'D|';
                typedSeqPrint = [targetSeqPrint(1:iUE-1) '-' targetSeqPrint(iUE:end)];
            elseif length(corrCheckKP)==length(targetSeqPrint) && length(corrCheckKP)==length(typedSeqPrint)
                corrCheckKP(iUE) = 'S';
            end
        end
        if length(typedSeqPrint) > length(targetSeqPrint)
            targetSeqPrint = [targetSeqPrint typedSeqPrint(length(targetSeqPrint)+1:end)];
        elseif length(typedSeqPrint) < length(targetSeqPrint)
            targetSeqPrint = targetSeqPrint(1:length(typedSeqPrint));
        end
        
        %% ================================================================
        disp([targetSeqPrint ' (Target sequence)']);
        disp([typedSeqPrint ' (Typed sequence)']);
        disp([corrCheckKP ' (Error types and locations)']);
    elseif strcmp(method,'regexp')
        typedSeqPrint = typedSeq;
        targetSeqPrint = perfectSeq;
    end
    if nKP==sum(iSeqError)
        disp('No correct keypresses were found for this trial.');
    end
    if  nargout == 3
            align.targetSeq = targetSeqPrint;
            align.typedSeq = typedSeqPrint;
            align.compare = checkKP;
    end
end

%Now remove "deletion" holders so that output variables match with
%typed KP sequence input array instead of alignment array
iSeqError(strcmp(seqErrorType,'Remove From List')) = [];
seqErrorType(strcmp(seqErrorType,'Remove From List')) = [];

%% CHECK OUTPUT LENGTH TO MAKE SURE IT MATCHES LENGTH OF ORIGINAL TYPED SEQUENCE INPUT
if length(iSeqError) ~= length(typedSeq) || length(seqErrorType) ~= length(typedSeq)
    warning('Mismatch between size of INPUT & OUTPUT arrays. Please investigate.')
end

%% Write Ouput
varargout(1) = {iSeqError};
if nargout >= 2
    varargout(2) = {seqErrorType};
end
if nargout == 3
    varargout(3) = {align};
end