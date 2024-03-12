function TfileID =  getTextureID(imid,noiseflag)

%gets a file identifier for texture stimuli from Corey

if ~noiseflag
    fid(1:5)  = {'im13-t*','im48-t*','im56-t*','im60-t*','im71-t*'};
else
    fid(1:5) = {'im13-n*','im48-n*','im56-n*','im60-n*','im71-n*'};
end

TfileID = fid(imid);
TfileID = TfileID{1};