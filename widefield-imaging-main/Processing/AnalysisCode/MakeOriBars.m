function MakeOriBars(im)

N = length(im(:,1));

hyp = 7; %Hypotenuse

cidx = hsv;

oris = 0:20:180;
Nori = length(oris);
ydom = N*(0:Nori-1)/Nori;
ydom = ydom+(ydom(2)-ydom(1))/2;
for i = 1:Nori
    
    center = [-8 ydom(i)]; 
    
    x = hyp*cos(oris(i)*pi/180);
    y = hyp*sin(oris(i)*pi/180);
    
    x = [center(1)-x center(1)+x];
    y = [center(2)-y center(2)+y];
    
    idx = 1+round(63*oris(i)/180);
    
    line(x,y,'Clipping','off','linewidth',3,'color',cidx(idx,:))
    
end