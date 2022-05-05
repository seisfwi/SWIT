% progress_txt : console based progress bar
%
% Ex1 : 
%   for i=1:10000;
%     progress_txt(i,10000,'Ciao');
%   end
%
% Ex1 :
%
%   for i=1:10;
%   for j=1:10;
%   for k=1:10;
%     progress_txt([i j k],[10 100 1000],'i','j','k');
%   end
%   end
%   end
%
% TMH/2005, thomas@cultpenguin.com
%
function progress_txt(i,max,varargin);
  
  if nargin==0
    help progress_txt
    return;
  end
  
  ncols=length(i);
  
  %
  nchar=45;  
  
  % 
  pc=i./max;
  
  % clear command window
  clc; 

  for m=1:ncols
    
    try
      txt=varargin{m};
    catch
      txt='';
    end
    
    char_prog='';
    for j=1:nchar
      if j<=(pc(m)*nchar);
        char_prog=[char_prog,'#'];
      else
        char_prog=[char_prog,'_'];
      end
    end
    disp(sprintf('%10s %s %3.1f%% %d/%d',txt,char_prog,100*pc(m),i(m),max(m)))
    
  end