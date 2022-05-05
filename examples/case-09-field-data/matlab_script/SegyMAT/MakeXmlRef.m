f{1}=dir('*.m')
f{2}=dir('visim/*.m')

fid=fopen('segymat-functions.xml','w');

fprintf(fid,'<sect1 id=\"%s\"><title>%s</title>\n','reference','M-file Reference');
for ff=1:length(f)
     
     for i=1:length(f{ff})
     [p,name,ext]=fileparts(f{ff}(i).name);
     disp(name)
     h=help(name);
     
     fprintf(fid,'<sect2 id=\"%s\"><title>%s</title>\n',['_',name],name);
     fprintf(fid,'<para><programlisting><![CDATA[%s]]></programlisting></para>\n',h);
     fprintf(fid,'</sect2>\n\n');
     
     end
     
end
fprintf(fid,'</sect1>');

fclose(fid);
  
