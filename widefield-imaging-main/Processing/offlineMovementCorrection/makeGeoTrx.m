function CHtrx = makeGeoTrx(CH,Px,Py)

%Make geometric transformation for the given movie based on x and y
%position waveforms

global ACQinfo


M = ACQinfo.linesPerFrame;
N = ACQinfo.pixelsPerLine;

m = 1:M;
n = 1:N;

CHtrx = zeros(size(CH));
for f = 1:length(CH(1,1,:))
  
    idframe = ((f-1)*M+1):f*M;
    
    y = m - Py(idframe);    
    y = y'*ones(1,N);
    x = ones(M,1)*n - Px(idframe)'*ones(1,N);

    %CHtrx(:,:,f) = griddata(x,y,double(CH(:,:,f)),np,mp,'cubic');
    %ZI = interp2(X,Y,Z,XI,YI)
    
    warning('off', 'all') 
    
    CHdum = CH(:,:,f);
    for col = 1:N        
        [CHdum(:,col)] = interp1(y(:,col),CHdum(:,col),m','spline');        
    end
    
    for row = 1:M        
        [CHdum(row,:)] = interp1(x(row,:),CHdum(row,:),n,'spline');        
    end
    
    warning('on', 'all') 
    
    CHtrx(:,:,f) = CHdum;
    
end