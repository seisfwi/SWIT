% isoctave : checks of octave
function r=isoctave
  v=version;
  if (str2num(v(1)))>4
    r=0;
  else
    r=1;
  end
