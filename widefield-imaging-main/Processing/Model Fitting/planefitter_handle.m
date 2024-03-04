function err = planefitter_handle(param)

%Ian Nauhaus

global RF xpts ypts;

img = xpts*param(1) + ypts*param(2) + param(3);

img = angle(exp(1i*img*pi/180))*180/pi; %wrap it
id = find(img(:)<0);
img(id) = img(id)+360;

err = nansum(abs(oridiff(img(:)*pi/180,RF(:)*pi/180)));


function dist = oridiff(angle1,angle2)

w1 = exp(1i*angle1);
w2 = exp(1i*angle2);
dist = angle(w1 ./ w2);

