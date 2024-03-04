function dist = oridiff(angle1,angle2)

%pepOriDiff        Returns the difference between the angles in angle1
%                  and angle2 in the orientation domain (that is they
%                  wrap around at pi radians!!!  The angles should be in rad.

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2;

%Another way to do it that I can put in an equation:
%dist = mod(angle1-angle2,pi);
%dist = mod(dist+pi/2,pi)-pi/2;

%Yet another (probably the easiest way to define it in text)
% dist = angle1-angle2;
% id1 = find(dist>pi/2);
% id2 = find(dist<-pi/2);
% dist(id1) = dist(id1)-pi;
% dist(id2) = dist(id2)+pi/2;
