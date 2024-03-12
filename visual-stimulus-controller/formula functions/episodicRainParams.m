function a = episodicRainParams(ori,barWidth,x_size)

y_size = x_size;
mask_radius = x_size/2; 
s_freq = 1/x_size; 
%y_zoom = abs(round(10*(cos(ori*pi/180)+.2)));  
%x_zoom = abs(round(10*(sin(ori*pi/180)+.2)));
x_zoom = 5;
y_zoom = 5;
s_duty = barWidth/x_size;

a(1)=y_size; 
a(2)=mask_radius; 
a(3)=s_freq; 
a(4)=y_zoom; 
a(5)=x_zoom; 
a(6)=s_duty;