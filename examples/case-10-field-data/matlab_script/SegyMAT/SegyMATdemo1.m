% SegyMATdemo1 : Creates, Reads and plots a Segy File;

% load test data set
load Clown;

% choose a file name;
file='segy_test.sgy';

% Write segy file, using default setting
WriteSegy(file,X);

% Read the segy file with Segy Header and Segy Trace Header
[Data,STH,SH]=ReadSegy(file);


x=[STH.TraceNumber];
t=SH.time;


subplot(2,2,1)
wiggle(x,t,Data) % variable Area
title('Wiggle')
subplot(2,2,2)
title('Wiggle')
wiggle(x,t,Data,'VA') % variable Area
title('Variable Area')

subplot(2,1,2)
scale=1; % MAX Data Value
showmax=300; % Max number of traces to display using wiggle
plimage=1; % plot image plot
wiggle(x,t,Data,'wiggle',[],showmax,plimage);
title('Wiggle + image')