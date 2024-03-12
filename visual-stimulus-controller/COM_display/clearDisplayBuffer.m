function clearDisplayBuffer   

global DcomState 

comhandle = DcomState.serialPortHandle;

%Clear the buffer
n = get(comhandle,'BytesAvailable')
if n > 0
    fread(comhandle,n); %clear the buffer
end
