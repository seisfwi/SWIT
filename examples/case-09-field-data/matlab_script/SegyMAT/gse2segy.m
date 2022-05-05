function file_out=gse2segy(file_in,file_out);

if nargin<2 
    [p,f]=fileparts(file_in);
    file_out=[f,'.sgy'];
end

GSE = read_gse_int(file_in);

for i=1:length(GSE)
    YearDataRecorded(i)=GSE{i}.year;
    DayOfYear(i)=GSE{i}.day;
    HourOfDay(i)=GSE{i}.hour;
    MinuteOfOur(i)=GSE{i}.minute;
    SecondOfMinute(i)=GSE{i}.second;
    TimeBaseCode(i) = 1; %Local(1), GMT(2); other(3), UTC(4)
    data(:,i)=GSE{i}.data(:);
end
dt=GSE{1}.dt;
WriteSegy(file_out,data,    'dt',dt,'YearDataRecorded',YearDataRecorded,...
    'DayOfYear',DayOfYear,...
    'HourOfDay',HourOfDay,...
    'MinuteOfOur',MinuteOfOur,...
    'SecondOfMinute',SecondOfMinute,...
    'TimeBaseCode',TimeBaseCode);

