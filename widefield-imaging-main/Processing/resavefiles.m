function resavefiles(root,oldfiletype,newfiletype)

%root = 'C:\NHP data summary\ORISFODexpts';
%oldfiletype = '.fig'
%newfiletype = '.eps'

files = dir([root,'\*' oldfiletype]);

for i = 1:length(files)
    open([root '\' files(i).name])

    if strcmp('.eps',newfiletype)
        print('-depsc2',[root '\' files(i).name(1:end-4) '.eps'])
    elseif  strcmp('.jpeg',newfiletype)
        print('-djpeg',[root '\' files(i).name(1:end-4) newfiletype])
    end
    
    close
end