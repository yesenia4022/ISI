%uiopen('C:\Mandi\xq8_vertical_ret_20.fig',1)
ax1 = gca; % get current axis.
ax1.Children; % returns ax1 as an array. 
get(ax1.Children,'CData');% gets the RGB value of the pixels. 
n1 = get(ax1.Children,'CData'); %saving the RGB values of the pixels as a variable.
%n2 = imresize(n1,0.2) % downsizing the data to remove some of the noise.
figure
imagesc(diff(n1)) 
imagesc(diff(n1,[],1)) % finds the difference between the pixels in the x direction.
imagesc(diff(n1,[],2)) % finds the difference between the pixels in the y direction.
x = diff(n1,[],1); % saving as a variable. 
y = diff(n1,[],2); % saving as a variable. 
m1 = roipoly; % saving region of interest that I draw as a variable. 
x(m1); % finding the difference between the pixels in x direction of my ROI.
y(m1); % finding the difference between the pixels in y direction of my ROI. 
x1 = x(m1); % saving as a variable. 
y1 = y(m1); % saving as a variable. 
%figure
%plot(sort(x1)) this shows the distribution. 
plot(x1)
%figure
plot([x1 y1])
%figure
plot([x1' y1']) % Plots the pixel differences in x and y transposed. 
xy1 = [x1' y1']; % saving as a variable.
%figure
plot(sort(xy1)) % Shows the distribution.
signal = mean(xy1)
noise = std(xy1) % Finds the standard deviation of the distribution.
%noise = var(xy1)
snr = signal./noise

figure
imagesc(n1)
bw = roipoly;
12:15
ddmap = diff(diff(n1,1),2);
12:16
ddvector = ddmap(find(bw));
12:16

mn=mean(ddvector);
md=median(ddvector);
stdv=std(ddvector);

figure, hist(ddvector)
mn
md
stdv

