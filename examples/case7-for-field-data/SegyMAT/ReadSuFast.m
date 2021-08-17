% ReadSuFast
%
% PURPOSE : reads a SEISMIC section i  SU format in big endian format, 
%           strips the headers and returns the field in the matrix seis.
%           If nx==0 and nt<>0, nx will be computed
%           If nt==0 and nx<>0, nt will be computed           
%
% Call : function seis=ReadSuFast(fileid,nt,nx,'byteorder');
%           byteorder : 'l' for little or 'b' for big endian (Default : Native )
%
% BY : TMH 1/8 1997
% Updated by Thomas Mejer Hansen : 22-03-1999
%
function [seis,nt,nx]=ReadSuFast(fileid,nt,nx,byte);
  
  if nargin==4 | nargin==2,
    if nargin==2, byte=nt; nt=0; end
    if byte=='b'
      b_order='ieee-be';
    else
      b_order='ieee-le';
    end
    SegymatVerbose(['SU2MAT : Byteorder : ',b_order],2)
  else
    % USE DEFAULT BYTEORDER
    SegymatVerbose(['SU2MAT : Byteorder : DEFAULT'],2)
  end
  
  
  if nargin==1 | nargin==2,
    if exist('b_order')==1, fid=fopen(fileid,'r',b_order);
    else, fid=fopen(fileid,'r'); end
    
    S0=fread(fid,60,'int16');nt=S0(58);
    if nt<0,
      % THIS CANNOT BE, MAYBE WRONG BYTE ORDER
      SegymatVerbose([mfilename,' : Wrong NT, maybe wrong byte order'])
      seis=[];nt=nt;nx=[];
      return
    end
    SegymatVerbose(['nt=',num2str(nt)],2)
    nx=0;
    fclose(fid);
  end
  
  if exist('b_order')==1, fid=fopen(fileid,'r',b_order);
  else, fid=fopen(fileid,'r'); end
  S=fread(fid,inf,'float');
  
  l=length(S);
  
  
  if nt==0, nt=l/nx-60; SegymatVerbose(['nt=',num2str(nt)],2);end
  if nx==0, nx=l/(nt+60); SegymatVerbose(['nx=',num2str(nx)],2); end,
  
  S=reshape(S,nt+60,nx);
  
  seis=S(61:nt+60,:);
  
