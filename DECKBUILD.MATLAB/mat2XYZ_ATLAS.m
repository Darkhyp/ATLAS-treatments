function outXYZ = mat2XYZ_ATLAS(filename,isSave)
% filename without file extension
% X,J_SC,V_OC,CE,FF are stored in (filename).m

if ~exist('isSave','var')
    isSave = false;
end

load(filename);

outXYZ = [];
str = '';
for n=length(X):-1:1
    outXYZ = [outXYZ,X{n}(:)];
    str = [str,vars{length(X)-n+1,1},'	'];
end
outXYZ = [outXYZ,J_SC(:),V_OC(:)*1e3,CE(:),FF(:)];
str = [str,'J_SC	V_OC	CE	FF'];

ind = find(CE(:)==max(CE(:)));
fprintf('\nMaximum efficiency (#%i): ',ind)
fprintf('%g; ',outXYZ(ind,:))
fprintf('\t%s\t',filename)
fprintf('\n')

if isSave
    if strcmpi(filename(end-1:end),'.m')
        filename = filename(1:end-2);
    end
    fout = fopen([filename,'XYZ.dat'],'w');
    fprintf(fout,'%s\n',str);
    fclose(fout);
    save([filename,'XYZ.dat'],'outXYZ','-ascii','-double','-tabs','-append')
end
