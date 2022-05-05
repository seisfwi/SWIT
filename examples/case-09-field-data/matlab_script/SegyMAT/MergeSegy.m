% MergeSegy : Merge multiple SEGY files
%
% Example :
%    MergeSegy('*.sgy','merge.sgy')
%
%    f{1}='file1.sgy';
%    f{2}='file2.sgy';
%    f{3}='file3.sgy';
%    MergeSegy(f,'merge.sgy')
%
%
% Note: All imput segy files must have the same constant trace length
%       The SEGY header of the merged SEGY file will be the SEGY header
%       form the first input SEGY file.
%
%

function [file_out,D,STH,SegyHeader]=MergeSegy(files,file_out)

if nargin<2
    file_out='segymerge.sgy';
end

if isstr(files)
    file_names=dir(files);
end

if iscell(files);
    for i=1:length(files);
        file_names(i).name=files{i};
    end
end

for i=1:length(file_names);
    [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(file_names(i).name);

    if i==1;
        ss=size(Data);
    else
        %check size
        if (length(find(ss==size(Data)))~=2)
           SegymatVerbose(sprintf('%s : Data sizes differs',mfilename),-1); 
        end
    end

    if i==1;
        D=Data;
        STH=SegyTraceHeaders;
    else
        n_sth=length(STH);
        for j=1:length(SegyTraceHeaders)
            STH(j+n_sth)=SegyTraceHeaders(j);
            STH(j+n_sth).TraceNumber=j;
        end
        D=[D Data];
    end
    
end

WriteSegyStructure(file_out,SegyHeader,STH,D);
