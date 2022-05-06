function [y, x, S] = RTOReadBin(filename, acquisitions, xInterval, nOutputSamplesMax)
% RTOREADBIN   
% ****************************************************************************
% FUNCTION:      [y, x, S] = RTOReadBin(filename)
%                [y, x, S] = RTOReadBin(filename, acquisitions, xInterval,...
%                            nOutputSamplesMax)
%
% SPECIFICATION: The R&S RTO series oscilloscopes use a proprietary .bin format
%                as the standard export format for waveform data. 
%
%                The format consists of two files. The waveform samples are
%                stored in the file <filename>.Wfm.bin following the IEEE 754
%                standard for floating point arithmetic. A companion file called
%                <filename>.bin contains the waveform properties in a
%                XML-compatible format.
%   
%                The function RTOReadBin reads binary files exported with 
%                .bin format from the RTO and produces the sample array y, 
%                the horizontal vector x, and the structure S with waveform 
%                properties.
%
% PARAMETERS:    filename - Name of the file. File extension is not
%                required. Path name muss be included if the file is not 
%                located in the current working directory.
%
% OPTIONAL PAR.: acquisitions - There are three possibilities for this input
%                       parameter:
%                       1) If omitted or empty: all available acquisitions 
%                           are read.
%                       2) Scalar value: the corresponding acquisition 
%                           (only one) is read using MATLAB Indexing s.
%                           below
%                       3) Acquistion interval [acqMin, acqMax]:
%                           acquisitions from acqMin until acqMax are read.
%                           i.e. [2 4] means acquisition 2, 3 and 4 are
%                           read.
%
%                           Equivalence between RTO and MATLAB index system
%                           Matlab indexing is used in this function:
%                           RTO Indexing      -6  -5  -4  -3  -2  -1   0
%                                             ---|---|---|---|---|---|---
%                           MATLAB Indexing    1   2   3   4   5   6   7
%                           In the figure above the last waveform acquired
%                           in the RTO corresponds to the index 7
%
%                xInterval - User-defined time/frequency interval [xMin, xMax].
%                       Only samples within this interval are read from the
%                       .bin file. The whole waveform is read when this
%                       parameter is omitted or empty.
%                nOutputSamplesMax - Maximum number of output samples per 
%                       acquistion or xInterval. When the number of samples 
%                       within one acquisition (respectively xInterval)
%                       exceeds nOutputSamplesMax, then the samples
%                       are decimated by the smallest integer that yields a 
%                       number of output samples below nOutputSamplesMax.

%
% RETURN VALUES: y - Matrix of waveform samples (double format or for raw 
%                    ADC values - 8 bit integer). 
%                    Each column corresponds to one acquisition.
%                    Last column corresponds to the last acquisition.
%                    In case of multi channel export the variable y is a 
%                    3-dimensional array.
%                    y[NofSamples, NofAcq, NofChs]
%                    i.e: y[5000,1,4] indicates a multi channel export (all
%                    4 analog channels), one acquisition pro channel with
%                    5000 values.
%                    For Envelope the NofChs is double because the resulting 
%                    diagram consists of two waveforms: the minimums
%                    (floor) and maximus (roof). i.e. y[1000,1,2] for one
%                    channel export indicates y(1000,1,1) is the waveform
%                    floor and y(1000,1,2) is the waveform roof.
%                x - Vector of time/frequency instants. Remains unchanged for
%                    all acquisitions.
%                S - Structure with waveform properties.
%
% LIMITATIONS:   Multi-Channel-Export:
%                   nNofChannels must coincide with the number of channels 
%                   exported in the file, which means for four-channel 
%                   export it is not possible to read two of them. For this 
%                   case those channels can be excluded outside this
%                   function.
%                Export using Envelope:
%                   For exports using trace arithmetic envelope, this script 
%                   supports only the case when all channels (in case of multi
%                   channel export) are exported with envelope. 
%                   i.e. two channel export, Ch1 without waveform arithmetic 
%                   and Ch2  with envelope is not supported in this script.
%
% EXCEPTIONS:    Not used.
% ****************************************************************************

% ****************************************************************************
%
% Copyright: (c) 2013 Rohde & Schwarz. All rights reserved.
%                Muehldorfstr. 15, D-81671 Munich, Germany
%
% TERMS OF USE:  Rohde & Schwarz products and services are supplied to
%                customers subject to certain contractual terms and 
%                conditions. In addition, there are some requirements that 
%                apply especially to certain products, customers or 
%                circumstances.
%                see http://www.termsofuse.rohde-schwarz.com/
%
% PROJECT:       Matlab tools.
%
% COMPILER:      Matlab 8.2 (R2013b).
%
% LANGUAGE:      Matlab Interpreter.
%
% AUTHOR:        Ruben Villarino
%                Johannes Ganzert
%                Rafael Ruiz
%                Andreas Ritter
%
% PREMISES:      None.
%
% REMARKS:       None.
%
% HISTORY:       $Log: $
%
% COMMENTS:      RTOReadBin V. 2.3 supports vertical offsets due to position scaling
%                RTOReadBin V. 2.2 supports export of history including
%                time stamp for each waveform. S.vecTimeStamps contains the 
%                relative time stamp values. This also works for the
%                ultra segmentation mode
%                RTOReadBin V. 2.1 supports digital waveforms (MSO export)
%                RTOReadBin V. 2.0 supports multiple acquisition and 
%                multiple channel export 
%
% VERSION:       2.3
%
% ****************************************************************************

%% --- Check input parameters ---

% Check number of arguments
narginchk(1, 4);   

%   Check existence of input file
[pathstr,name,ext] = fileparts(filename);
if isempty(pathstr)
    pathstr = '.';
end
if isempty(ext)
    filenameHeader = fullfile(pathstr, [name '.bin']); 
    filenameSamples = fullfile(pathstr, [name '.Wfm.bin']);
elseif strcmpi(ext, '.bin')
    % Check if the user has selected the waveform data "Wfm.bin"
    %   If so, eliminate the ".Wfm" part
    if strcmpi(name(max(1,numel(name) - 3):end), '.Wfm'),
        name = name(1:end-4);
    end
    filenameHeader = fullfile(pathstr, [name '.bin']); 
    filenameSamples = fullfile(pathstr, [name '.Wfm.bin']);
else
    filenameHeader = fullfile(pathstr, [name ext '.bin']); 
    filenameSamples = fullfile(pathstr, [name ext '.Wfm.bin']);
end

if ~exist(filenameSamples,'file')
    error('File <%s> not found', filenameSamples)
end

if ~exist(filenameHeader,'file')
    warning([mfilename ':HeaderMissing'],...
        'Parameter file <%s> not found, reverting to standard parameters',filenameHeader);
    filenameHeader = [];
end


%   Check validity of horizontal interval
if ~exist('xInterval', 'var')
    xInterval = [];
end
if ~isempty(xInterval)
    if length(xInterval) == 1
        error('Please specify the parameter xInterval as [xMin, xMax]')
    elseif xInterval(1) > xInterval(2) || xInterval(1) == xInterval(2)
        error('Please specify the parameter xInterval as [xMin, xMax]')
    end
end

%   Check validity of maximum number of output samples
if ~exist('nOutputSamplesMax', 'var')
    nOutputSamplesMax = [];
end
if ~isempty(nOutputSamplesMax)
    if length(nOutputSamplesMax) > 1
        error('The parameter nOutputSamplesMax must be a scalar value')
    end
    if nOutputSamplesMax < 1 || mod(nOutputSamplesMax, 1)
        error('The parameter nOutputSamplesMax must be a positive integer')
    end
end


%% --- Load waveform parameters from file ---

% Use default parameters if the previous file check has failed
if isempty(filenameHeader)
    % Set up parameter structure
    S.RecordLength = [];
    S.XStart = [];
    S.XStop = [];
else
    % Import parameters from companion file
    % We use Java to retrieve data from a DOM representation of the XML
    % input file
    
    % Import waveform parameters
    xDoc = xmlread(filenameHeader);
    
    % Get property nodes
    properties = xDoc.getElementsByTagName('Prop');
    
    % Get name and value of all properties
    S = [];     % Initilize output structure
    for k = 0 : properties.getLength-1
        % Get all attributes of current property node
        attributes = properties.item(k).getAttributes;
        
        % Get property name (one of the node attributes)
        if isempty(attributes.getNamedItem('Name'))
            continue;
        else
            propertyNameStr = (attributes.getNamedItem('Name').getNodeValue.toCharArray)';
        end
        
        % Get property value (one of the node attributes)
        if ~isempty(attributes.getNamedItem('Value'))
            propertyValueStr = (attributes.getNamedItem('Value').getNodeValue.toCharArray)';
            % Check if it is numeric
            propertyValueDouble = str2double(propertyValueStr);
            if isnan(propertyValueDouble)
                propertyValue = propertyValueStr;
            else
                propertyValue = propertyValueDouble;
            end
            
        elseif ~isempty(attributes.getNamedItem('I_0'))
            %Properties of MultiChannel Export
            i=0;
            if ~isempty(attributes.getNamedItem('Size'))
                size=str2double(attributes.getNamedItem('Size').getNodeValue.toCharArray);  
                propertyValue=cell(1,size);
                while ~isempty(attributes.getNamedItem(['I_' int2str(i)]))
                    propertyValueStr = (attributes.getNamedItem(['I_' int2str(i)]).getNodeValue.toCharArray)';
                    % Check if it is numeric
                    propertyValueDouble = str2double(propertyValueStr);
                    if isnan(propertyValueDouble)
                        propertyValue{i+1} = propertyValueStr;
                    else
                        propertyValue{i+1} = propertyValueDouble;
                    end
                    i=i+1;
                end
            else
            propertyValue = '';
            end
        else
                propertyValue = '';
        end
        
        % Store properties in structure
        S.(propertyNameStr) = propertyValue;
    end
end

% Order the fields in S alphabetically
S = orderfields(S);

% Get Number of Channels from Header File
nNofChannels=1;
if isfield(S, 'MultiChannelExport')
    if strcmp(S.MultiChannelExport, 'eRS_ONOFF_ON') 
        nNofChannels=0;
        for(count=1:length(S.MultiChannelExportState))
            if(strcmp(S.MultiChannelExportState{count},'eRS_ONOFF_ON'))
                nNofChannels=nNofChannels+1;
            end
        end    
    end 
end
% Get header information from waveform file and finalize the initialization 
% of the parameter structure
%   Open the file and read the first eight bytes
[fid, message] = fopen(filenameSamples, 'r');
if ~isempty(message)
    error('Error while opening waveform file: %s', message)
end

% Read header segment information. The second value is the number of samples
% in one acquisition. For files with multiple channels or
% acquisitions, this in not an indication of the file length.
fileInfo = fread(fid, 2, 'uint32');

% Close file
fclose(fid);

% Initialization depends on existence of parameter file
if isempty(filenameHeader)
    % Initialice the fields that were left empty
    S.SignalHardwareRecordLength      = fileInfo(2);
    S.SignalRecordLength              = fileInfo(2);
    S.SignalResolution                = 1;
    S.XStart                          = 0;
    S.XStop                           = S.RecordLength - 1;
    S.HardwareXStart                  = S.XStart;
    S.HardwareXStop                   = S.XStop;
    nNofSamplesPerAcqPerCh            = fileInfo(2) / nNofChannels;
else
    % For compatibility reasons with FW versions older than 2.0 we need to check 
    % that the number of samples in waveform and parameter files match.
    % With the release of FW version 2.0 the firmware version is available 
    % in the header file. According to this we decide for one of them for the 
    % following computations in case the parameters differ from each other.
    if S.SignalHardwareRecordLength ~= fileInfo(2);
        warning('Number of samples in waveform file and header file do not match')
        if ~isfield(S, 'FirmwareVersion')
            nNofSamplesPerAcqPerCh = fileInfo(2) / nNofChannels;   % FW < 2.0
        else
            nNofSamplesPerAcqPerCh = S.SignalHardwareRecordLength; % FW >= 2.0
        end
    else
        nNofSamplesPerAcqPerCh = S.SignalHardwareRecordLength;
    end  
end   


%% --- Check for waveform arithmetic ENVELOPE ---
% Using waveform arithmetic Envelope two envelope waveforms are exported,
% because of that the number of channels is double internally 

if isfield(S, 'TraceType')
    if strfind(S.TraceType, 'ENVELOPE')   
        nNofChannels = nNofChannels * 2; 
    end
end

%% --- Check for digital waveforms ---
bIsDigitalSource = false;
if isfield(S,'SourceType');
    if strfind(S.SourceType,'SOURCE_TYPE_DIGITAL');
        bIsDigitalSource = true;
    end
end

%% --- Check time stamp information ---
bTimeStampState = false;
if isfield(S,'TimestampState');
    if strfind(S.TimestampState,'eRS_ONOFF_ON');
        bTimeStampState = true;
    end
end


%% --- Following parameters must be calculated ---

%   a) Preamble length
%   b) Sample interval
%   c) Sample decimation

%% --- a) Compute preamble length ---

%   The sample file may have some preamble samples. This samples are needed to
%   reproduce measurements performed on the waveform
if isfield(S,'LeadingSettlingSamples')
    nPreambleSamples = S.LeadingSettlingSamples;
else
    % Estimate number of preamble samples (backwards compatibility)
    nPreambleSamples = round((S.XStart-S.HardwareXStart) / S.SignalResolution);
end

if nPreambleSamples < 0
    warning('Error while computing number of preamble samples')
    nPreambleSamples = 0;
end

%% ---  b) Compute sample interval ---

%   The first valid sample has the index 1 (Matlab indexing)
%   Map the user defined sample interval into a sample index range

tolerance = 1e-15;  % Used for comparing floating values (equals 1 fs for time domain)

if isempty(xInterval)
    % Default range
    idxRange = [1, S.SignalRecordLength];
else
    if ( S.XStart-xInterval(1) > tolerance ) && ( xInterval(2)-S.XStop > tolerance )
        warning('User-defined range lies outside waveform horizontal range')
    elseif (S.XStart-xInterval(1)) > tolerance
        warning('User-defined lower range lies outside waveform horizontal range')
    elseif (xInterval(2) - S.XStop) > tolerance
        warning('User-defined upper range lies outside waveform horizontal range')
    end
           
    % Compute lower and upper index range
    idxRange(1) = max(1, round((xInterval(1) - S.XStart) / S.SignalResolution));
    idxRange(2) = min(S.SignalRecordLength, round((xInterval(2) - S.XStart) / S.SignalResolution));
end

%% --- c) Compute sample decimation ---

%   Read a decimated version of the waveform if the user constrains the number
%   of imported samples
if isempty(nOutputSamplesMax)
    % Read all data
    decimationFactor = 1;
    nSamples = idxRange(2) - idxRange(1) + 1;
elseif diff(idxRange) <= nOutputSamplesMax    
    % Number of samples in range is less than upper limit
    % No decimation required, read all data
    decimationFactor = 1;
    nSamples = idxRange(2) - idxRange(1) + 1;
else
    % Decimation required
    decimationFactor = ceil((idxRange(2) - idxRange(1) + 1) / nOutputSamplesMax);
    nSamples = floor((idxRange(2) - idxRange(1)) / decimationFactor) + 1;
end    

%% --- Open file and determine size and precision per sample  ---

% Open file
[fid, message] = fopen(filenameSamples, 'r');
if ~isempty(message)
    error('Error while opening waveform file: %s',message)
end

% Get number of bytes and precision per sample
dataType = fileInfo(1);
sampleSizeX = 0; % samplSizeX is different from zero only for the X/Y interleaved export
switch dataType
    case 0              % samples are int8
        sampleSizeY = 1;
        precisionY = '*int8';
    case 1              % samples are int16
        sampleSizeY = 2;
        precisionY = '*int16';
    case 2              % samples are int24
        sampleSizeY = 3;
        precisionY = '*int24';
    case 3              % samples are int32
        sampleSizeY = 4;
        precisionY = '*int32';
    case 4              % samples are single (float)
        sampleSizeY = 4;
        precisionY = '*single';
    case 5              % samples are double 
        sampleSizeY = 8;
        precisionY = '*double';
    case 6              % samples are XY interleaved - X: double, Y:float
        sampleSizeX = 8;
		precisionX = '*double';
		sampleSizeY = 4;
        precisionY = '*single';
    otherwise
        error('Error unkown sample data type in waveform file')
end

% Define size and precision for digital waveforms
if bIsDigitalSource       
    sampleSizeY = 1;
    precisionY = '*ubit1';
end

%% --- Check acquisitions ---

% Define constant
headerSegmentInBytes = 8;

% Get size (in Bytes) of samples file (<filename>.Wfm.bin)
dirInfo = dir(pathstr);
[~, name, ~] = fileparts(filenameSamples);

index = strcmpi({dirInfo.name}, [name '.bin']);
fileSizeSamples = dirInfo(index).bytes;

% Compute number of available acquisitions in file
if bIsDigitalSource
    % MSO signals are saved bitwise. Each byte on file corresponds to 
    % 8 sequential digital values 
    nNofAvailableAcq = floor((fileSizeSamples - headerSegmentInBytes)*8 / (S.SignalRecordLength * sampleSizeY * nNofChannels)); % Interleaved x/y is not available for MSO channels
else
    if bTimeStampState
        nNofAvailableAcq = floor((fileSizeSamples - headerSegmentInBytes) / (nNofSamplesPerAcqPerCh * (sampleSizeX + sampleSizeY * nNofChannels)+8));
    else
        nNofAvailableAcq = floor((fileSizeSamples - headerSegmentInBytes) / (nNofSamplesPerAcqPerCh * (sampleSizeX + sampleSizeY * nNofChannels)));
    end
end

if ~exist('acquisitions', 'var') || isempty(acquisitions)    % all acquisitions will be read
    
    % No precedent acquisitions to be skiped
    NofPreAcq = 0;
    
    % Set number of acquisitions to the number of available acq. in file
    nNofAcquisitions = nNofAvailableAcq;
  
elseif length (acquisitions) == 1   % One acquisition will be read specified as an AcqIndex
       
    %Verify acquisition index       
    if acquisitions < 1 || mod(acquisitions, 1)
        error('Acquistion index must be positive integer')
    elseif acquisitions > nNofAvailableAcq
        error('Acquisition index exceeds %s which is the biggest acquisition index available',  num2str(nNofAvailableAcq))
    end
    
    % Determine number of precedent acquisitions
    NofPreAcq = acquisitions - 1;

    % Just one acquisition will be read.
    nNofAcquisitions = 1;
    
elseif length (acquisitions) == 2   % Acquisitions betwenn [acqMin, acqMax] will be read
        
    %Verify acquisition interval       
    if acquisitions(1) < 1 || mod(acquisitions(1), 1) || acquisitions(2) < 1 || mod(acquisitions(2), 1)
        error('Values for acquisition interval must be positive integers')
    elseif acquisitions(1) > acquisitions(2)
        error('Please specify acquisition interval in the form [acqMin, acqMax]')
    elseif acquisitions(2) > nNofAvailableAcq
        error('Acquisition interval is out of range. For the current file the biggest acquisition index available is %s',  num2str(nNofAvailableAcq))
    end
    
    % Determine number of precedent acquisitions
    NofPreAcq = acquisitions(1) - 1;
    
    % Compute number of acquistions to be read
    nNofAcquisitions = acquisitions(2) - acquisitions(1) + 1;    
else
    error('Please specify the parameter acquisition either as a scalar value or as an interval [acqMin, acqMax]')    
end


%% --- Verify that the required number of bytes or bits (MSO data) are available in file ---
if bIsDigitalSource
    % Number of samples given in Bit for all channels
    nNofBit_AllCHs = nSamples * (sampleSizeX + sampleSizeY * nNofChannels);
       
    % Number of samples given in Bit for all acquisitions and all channels
    nNofBit_AllCHs_AllAcqs = nNofBit_AllCHs * nNofAcquisitions;
       
    % Check if the number of required Bit are available in file
    if nNofBit_AllCHs_AllAcqs > (fileSizeSamples - headerSegmentInBytes) * 8
        error('Number of required samples exceeds number of samples available in file')
    end
else
    % Number of samples given in byte for all channels
    nNofBytes_AllCHs = nSamples * (sampleSizeX + sampleSizeY * nNofChannels);

    % Number of samples given in byte for all acquisitions and all channels
    nNofBytes_AllCHs_AllAcqs = nNofBytes_AllCHs * nNofAcquisitions;

    % Check if the number of required bytes are available in file
    if nNofBytes_AllCHs_AllAcqs > (fileSizeSamples - headerSegmentInBytes)
        error('Number of required samples exceeds number of samples available in file')
    end
end

%% --- Skip 8 byte header ---
if fseek(fid, headerSegmentInBytes, 'bof') ~= 0
    error('Error while seeking header')
end

%% --- Skip precedent acquisitions (if necessary) ---
if fseek(fid, (nNofSamplesPerAcqPerCh * NofPreAcq) * (sampleSizeX + sampleSizeY * nNofChannels), 'cof') ~= 0
    error('Error while seeking precedent acquisitions')
end

% skip previous timestamps if necessary
if bTimeStampState
    if fseek(fid, 8*NofPreAcq, 'cof')~=0
        error('Error while skipping precedent Timestamps')
    end
end


%% --- Read waveform samples ---

% Initialize matrix or three dimensional array of waveform samples
y = zeros(nSamples, nNofAcquisitions, nNofChannels);


vecTimeStamps = zeros(nNofAcquisitions,1);

% Define constants for skipping preamble and postamble
nPreambleInBytes = (nPreambleSamples + idxRange(1)-1) * (sampleSizeX + sampleSizeY * nNofChannels);
nPostambleInBytes = (nNofSamplesPerAcqPerCh - nSamples*decimationFactor - nPreambleSamples - (idxRange(1)-1)) * (sampleSizeX + sampleSizeY * nNofChannels);
                      
if dataType ~= 6        % samples are not XY
    for nAcqIndex = 1 : nNofAcquisitions  
        if bTimeStampState
            % Read time stamp. Notice the time stamp is a double value and
            % it precedes all vertical values of each acquisition
            vecTimeStamps(nAcqIndex) = fread(fid, 1, '*double', 0);
        end
        % Skip preamble
        if fseek(fid, nPreambleInBytes, 'cof') ~= 0
            error('Error while seeking postamble')
        end   
        
        % Get vertical values
        %   Decimate data by jumping decimationFactor samples
        yAux = fread(fid, nSamples * nNofChannels, precisionY, (decimationFactor-1) * sampleSizeY * nNofChannels);
        
        % Use reshape to separate channels (undo the interleaving)
        yCHs = reshape(yAux, nNofChannels, []); %Each row vector corresponds to one channel
        
        % Assign values to the three dimensional array in case of multiple
        % acquisitions and multiple channels
        y(:, nAcqIndex, :) = yCHs';
        
        if dataType==0 ||dataType==1 %Convert Data from RAW to Voltages
            VertOffsetByPosition = S.VerticalScale*S.VerticalPosition; %this was missing in the formula in the user guide; Jazz# 470158; change 31-Jul-2018, AR
            if isfield(S, 'MultiChannelExport')
                if strcmp(S.MultiChannelExport, 'eRS_ONOFF_ON') 
                    %MultiChannel Export 
                    ChanIndex=1;
                    for(ChanNum=1:length(S.MultiChannelExportState))
                        if strcmp(S.MultiChannelExportState{ChanNum}, 'eRS_ONOFF_ON')
                            VertOffsetByPosition = S.MultiChannelVerticalScale{ChanNum}*S.MultiChannelVerticalPosition{ChanNum}; %calculate vertical offset for each channel
                            ConversionFactor=(1/S.NofQuantisationLevels)*S.MultiChannelVerticalScale{ChanNum}*S.VerticalDivisionCount;
                            y(:,nAcqIndex,ChanIndex)=y(:,nAcqIndex,ChanIndex)*ConversionFactor+S.MultiChannelVerticalOffset{ChanNum}-VertOffsetByPosition;  %Jazz# 470158; change 31-Jul-2018, AR
                            ChanIndex=ChanIndex+1;
                        end     
                    end 
                else
                    %single Channel Export
                    ConversionFactor=(1/S.NofQuantisationLevels)*S.VerticalScale*S.VerticalDivisionCount;
                    y(:,nAcqIndex,:)=y(:,nAcqIndex,:)*ConversionFactor+S.VerticalOffset-VertOffsetByPosition;    %Jazz# 470158; change 31-Jul-2018, AR
                end
            else
                %single Channel Export
                ConversionFactor=(1/S.NofQuantisationLevels)*S.VerticalScale*S.VerticalDivisionCount;
                y(:,nAcqIndex,:)=y(:,nAcqIndex,:)*ConversionFactor+S.VerticalOffset-VertOffsetByPosition;   %Jazz# 470158; change 31-Jul-2018, AR
            end
        end
            
            
        % Skip postamble if this was not the last acquisition
        if nAcqIndex < nNofAcquisitions
            if fseek(fid, nPostambleInBytes, 'cof') ~= 0
                error('Error while seeking postamble')
            end          
        end
    end
    
    % Reconstruct horizontal vector
    x = (idxRange(1)-1 : decimationFactor : idxRange(2)-1) * S.SignalResolution + S.XStart;

else                    % interleaved samples (XY)
    
    x = zeros(nSamples, nNofAcquisitions);
    for nAcqIndex = 1 : nNofAcquisitions
        
        if bTimeStampState
            % Read time stamp. Notice the time stamp is a double value and
            % it precedes all vertical values of each acquisition
            vecTimeStamps(nAcqIndex) = fread(fid, 1, '*double', 0);
        end
        
        % Skip preamble
        if fseek(fid, nPreambleInBytes, 'cof') ~= 0
            error('Error while seeking postamble')
        end   
        
        % Get horizontal and vertical values:  <X,Y> pair
        for indexSample = 1 : nSamples
            x(indexSample, nAcqIndex) = fread(fid, 1, precisionX);
            for indexCH = 1 : nNofChannels
                y(indexSample, nAcqIndex, indexCH) = fread(fid, 1, precisionY);
            end
            if fseek(fid, (decimationFactor-1) * (sampleSizeX + sampleSizeY * nNofChannels), 'cof') ~= 0
                error('Error while decimating interleaved samples')
            end
        end 
     
        % Skip postamble if this was not the last acquisition
        if nAcqIndex < nNofAcquisitions
            if fseek(fid, nPostambleInBytes, 'cof') ~= 0
                error('Error while seeking postamble')
            end 
        end
    end
end

    
if bTimeStampState
    % Assign vector with time stamps to the output structure
    S.('vecTimeStamps') = vecTimeStamps;
end

%% --- Close file ---
fclose(fid);

%% --- Check horizontal length of y and x ---
if (length(y) ~= length(x))
    warning([mfilename, ':NumelMismatch'], ...
        'Sample vector y and horizontal vector x have different lenghts')
end
