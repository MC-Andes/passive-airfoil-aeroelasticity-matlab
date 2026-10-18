function write_json(path,value)
%WRITE_JSON Checked JSON output with guaranteed handle cleanup.
fid=fopen(path,'w'); assert(fid~=-1,'cad:Report','Cannot write %s.',path);
c=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
