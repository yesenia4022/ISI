function G = Travelerfitguess

%Double check Initial guesses

global RF;

f = RF;


[ma_y ma_x] = find(f == max(f(:)));

G(1) = ma_x(1);  %Spatial center
G(4) = ma_y(1);  %time to peak of center
G(2) = length(f(1,:))/8;  %spatial sigma
G(5) = 2;   %temporal sigma
G(3) = 1.5;   %slope (wave velocity)

% G(6) = max(f(:))-min(f(:));  %amplitude
% G(7) = min(f(:));  %baseline
