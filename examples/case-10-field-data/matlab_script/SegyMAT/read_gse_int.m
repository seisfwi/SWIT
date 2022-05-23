function GSE=read_gse_int(filename);

fid=fopen(filename,'r');

i=0;
while ~feof(fid);
    i=i+1;
    SegymatVerbose(sprintf('Reading trace #%02d',i))

    line1=fgetl(fid);
    
    try,GSE{i}.type=line1(45:47);end
       
    try,GSE{i}.year=str2num(line1(6:9));,end
    try,GSE{i}.month=str2num(line1(11:12));,end
    try,GSE{i}.day=str2num(line1(14:15));,end
    
    try,GSE{i}.hour=str2num(line1(17:18));,end
    try,GSE{i}.minute=str2num(line1(20:21));,end
    try,GSE{i}.second=str2num(line1(23:24));,end
    try,GSE{i}.mills=str2num(line1(26:28));,end
    try,GSE{i}.nspersec=str2num(line1(57:68));end
    try,GSE{i}.dt=1./GSE{i}.nspersec;end
    try;GSE{i}.ns=str2num(line1(52:56));end
    
    line2=fgetl(fid);
    line3=fgetl(fid);
    
    for j=1:GSE{i}.ns
        l=fgetl(fid);
        GSE{i}.data(j)=sscanf(l,'%f');
    end
    
    %GSE{i}.data=fscanf(fid,'%f',GSE{i}.ns)
    
    line4=fgetl(fid);
    GSE{i}.CHK2=str2num(line4(6:length(line4)));
        
end
fclose(fid);