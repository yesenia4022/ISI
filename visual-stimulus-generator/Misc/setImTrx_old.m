function H = setImTrx(InputPts,OutputPts,gun)

%If you measure 'up' as positive Y (second column), and 'right' is positive X,
%this is the correct transform. 3 reasons: 
%1) Monitor is rotated 90 degrees (flipping of dimensions), 
%2) Rear projection (image is flipped around vertical (horizontal since rotated) axis)  
%3) My code (and Matlab) treats the y-axis as going from top to bottom.

global imcolortrx

if isempty(gun)
    
    imcolortrx = cell(1); %Reset it.
    imcolortrx{1} = eye(4,2); %If 'make' sees one element it will not trx any gun images
    
else
    
    OutputPts(:,1) = -OutputPts(:,1);
    InputPts(:,1) = -InputPts(:,1);
    OutputPts = fliplr(OutputPts); InputPts = fliplr(InputPts); %monitor is rotated
    
    %% Affine warping

    InputPts = [InputPts InputPts(:,1).*InputPts(:,2)]; 
    
    InputPts = [InputPts ones(size(InputPts,1),1)];  %[x y x*y 1]
    H = inv(InputPts'*InputPts)*InputPts'*OutputPts;
    
    Outhat = InputPts*H;
    
    Np = size(InputPts,2);
    
    imcolortrx{1} = eye(Np,2);    
    imcolortrx{2} = eye(Np,2);
    imcolortrx{3} = eye(Np,2);
    
    imcolortrx{gun} = H; %Trx the given gun; R=1;G=2;B=3
    
end
