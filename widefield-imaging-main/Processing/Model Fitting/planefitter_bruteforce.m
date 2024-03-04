function err = planefitter_bruteforce(xpts,ypts,RF)

%Ian Nauhaus

xslope_dom = linspace(-60,60,200);
yslope_dom = linspace(-60,60,200);
base_dom = linspace(-360,360,100);

for i = 1:length(xslope_dom)
    for j = 1:length(yslope_dom)
        for k = 1:length(base_dom)
            
            img = xpts*xslope_dom(i) + ypts*yslope_dom(j) + base_dom(k);
            
            img = angle(exp(1i*img*pi/180))*180/pi; %wrap it
            id = find(img(:)<0);
            img(id) = img(id)+360;
            
        end
    end
end

img = xpts*param(1) + ypts*param(2) + param(3);

img = angle(exp(1i*img*pi/180))*180/pi; %wrap it
id = find(img(:)<0);
img(id) = img(id)+360;

err = nansum(abs(oridiff(img(:)*pi/180,RF(:)*pi/180)));


function dist = oridiff(angle1,angle2)

w1 = exp(1i*angle1);
w2 = exp(1i*angle2);
dist = angle(w1 ./ w2);

