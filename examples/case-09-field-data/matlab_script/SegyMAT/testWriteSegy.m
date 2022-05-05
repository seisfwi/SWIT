% testWriteSegy : Script to test WriteSegy and WriteSegyStructure
%
n=200;
seisdata=peaks(n);
imagesc(seisdata);colorbar;
xlabel('Traces');
ylabel('Time')
title('The data that should be written to disk')

% use WriteSegy to write an IEEE Revision 1 formatted file
%
WriteSegy('test.segy',seisdata);

% Read the file we just created with all of its header values
%
[Data,SegyTraceHeaders,SegyHeader,HeaderInfo]=ReadSegy('test.segy');

% NOW WE CAN WRITE ALL THE FORMATS SUPPORTED BY SegyMAT :


% WRITE IBM FLOATING POINT REV 0
% THIS IS THE ONLY IMPLEMENTED FORMAT FOR REV 0
SegyHeader.SegyFormatRevisionNumber=0;
SegyHeader.DataSampleFormat=1;
WriteSegyStructure('data_IBM_REV0.segy',SegyHeader,SegyTraceHeaders,Data);

% WRITE IBM FLOATING POINT REV 1
SegyHeader.SegyFormatRevisionNumber=100;
SegyHeader.DataSampleFormat=1;
WriteSegyStructure('data_IBM_REV1.segy',SegyHeader,SegyTraceHeaders,Data);


% WRITE 4BYTE INT POINT REV 1
SegyHeader.SegyFormatRevisionNumber=100;
SegyHeader.DataSampleFormat=2;
WriteSegyStructure('data_4byteINT.segy',SegyHeader,SegyTraceHeaders,Data);

% WRITE 2BYTE INT POINT REV 1
SegyHeader.SegyFormatRevisionNumber=100;
SegyHeader.DataSampleFormat=3;
WriteSegyStructure('data_2byteINT.segy',SegyHeader,SegyTraceHeaders,Data);

% WRITE 1BYTE INT POINT REV 1
SegyHeader.SegyFormatRevisionNumber=100;
SegyHeader.DataSampleFormat=8;
WriteSegyStructure('data_1byteINT.segy',SegyHeader,SegyTraceHeaders,Data);


% WRITE IEEE POINT REV 1
SegyHeader.SegyFormatRevisionNumber=100;
SegyHeader.DataSampleFormat=5;
WriteSegyStructure('data_IEEE.segy',SegyHeader,SegyTraceHeaders,Data);

