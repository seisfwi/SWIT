function SegymatHelp(file);
  if nargin==0
    file='index';
  end
  
  [SegyMAT_root]=fileparts(which('ReadSegy'));
  
  docpath=[SegyMAT_root, filesep , 'html'];
  
  if isdir(docpath)
    helpfile=[docpath, filesep, file, '.html'];
    if exist(helpfile)
      web(['file://',helpfile])
    else
      SegymatVerbose(sprintf('Could not open help file : "%s"',helpfile))
    end
  else
    SegymatVerbose(sprintf('"%s" does not exist',docpath))
    SegymatVerbose('Cannot show local help')
    SegymatVerbose(sprintf('Trying to open the online help',docpath))
    web('http://segymat.sourceforge.net/doc/')
  end