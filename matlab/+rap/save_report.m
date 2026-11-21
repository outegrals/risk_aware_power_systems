function save_report(folder,name,value)
if ~isfolder(folder), mkdir(folder); end
f=fopen(fullfile(folder,[name '.json']),'w'); cleanup=onCleanup(@()fclose(f));
fprintf(f,'%s\n',jsonencode(value,PrettyPrint=true));
end
