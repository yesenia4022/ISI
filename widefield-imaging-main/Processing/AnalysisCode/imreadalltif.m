function A = imreadalltif(filename, framestart,framestop)
% IMREADALLTIFF  Utility function to read all the frames from a TIFF stack.
% It avoids most of the overhead of IMREAD when trying to access and read
% all the data from a TIFF stack by assuming that every frame is the same
% type and size as the first, and so not checking this.
% 
%   [A, FRAMESREAD] = IMREADALLTIFF(FILENAME, NFRAMES) reads NFRAMES 
%   from the TIFF file represented by the string FILENAME.
% 
%   Input:
%   NFRAMES should be a scaler representing the number of frames to read.
% 
%   Output:
%   A is a 3-dimensional array containing the image data.
%   
% 
%%%12/1/07: modified code to make it compatible with previous matlab
%%%releases


% read the first frame
tf = imformats('tif');
I = feval(tf.read, filename, 1);

% preallocate output
A = zeros(size(I,1), size(I,2), framestop-framestart+1);

% copy in first image
A(:,:,1) = I;

% try to read all frames
for frame=framestart:framestop
    A(:,:,frame-framestart+1) = feval(tf.read, filename, frame);
end

