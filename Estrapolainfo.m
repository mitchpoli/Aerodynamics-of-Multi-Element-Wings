clc;
clear; 
close; %% 1) Carica immagine
I = imread('ala10.png');   % metti il nome giusto del file
figure; imshow(I); axis on; hold on;

%% 2) Calibrazione assi
disp('Click 2 punti sull''asse X: sinistra (x=-1) e destra (x=2)');
[x1,y1] = ginput(1);   % corrisponde a x = -1
[x2,y2] = ginput(1);   % corrisponde a x =  2

disp('Click 2 punti sull''asse Y: basso (y=-1) e alto (y=1)');
[x3,y3] = ginput(1);   % y = -1
[x4,y4] = ginput(1);   % y =  1

% ipotizziamo assi ortogonali (buona approssimazione)
% relazioni lineari pixel -> coord. fisiche
mX = (2 - (-1)) / (x2 - x1);
qX = -1 - mX * x1;

mY = (1 - (-1)) / (y4 - y3);
qY = -1 - mY * y3;

px2x = @(px) mX*px + qX;
py2y = @(py) mY*py + qY;

%% 3) Definisci la riga da cui estrarre il profilo
disp('Click 2 punti sull''immagine lungo la riga da cui vuoi il profilo');
[xa,ya] = ginput(1);
[xb,yb] = ginput(1);

N = 200;  % numero di punti del profilo
xp = linspace(xa, xb, N);
yp = linspace(ya, yb, N);

% coord. fisiche corrispondenti
x_phys = px2x(xp);
y_phys = py2y(yp);

%% 4) Estrazione colore lungo la riga
I_double = im2double(I);
R = I_double(:,:,1);
G = I_double(:,:,2);
B = I_double(:,:,3);

% interpola i canali ai punti xp,yp
[Xgrid,Ygrid] = meshgrid(1:size(I,2), 1:size(I,1));
Rp = interp2(Xgrid, Ygrid, R, xp, yp);
Gp = interp2(Xgrid, Ygrid, G, xp, yp);
Bp = interp2(Xgrid, Ygrid, B, xp, yp);

rgb_profile = [Rp(:), Gp(:), Bp(:)];

%% 5) Mappa colore -> valore di velocità usando la stessa colormap
% assumiamo che la figura usi 'jet' da 0 a 1.5 (come nella barra colori)
cmap = jet(256);
val  = linspace(0, 1.5, 256);

vel_profile = zeros(N,1);
for k = 1:N
    diffRGB = cmap - rgb_profile(k,:);
    d2 = sum(diffRGB.^2,2);
    [~,idx] = min(d2);
    vel_profile(k) = val(idx);
end

%% 6) Plot del profilo estratto
figure; 
plot(x_phys, vel_profile, '-o');
grid on;
xlabel('x/c (approssimato)');
ylabel('\sqrt{U^2+V^2}/U_\infty (da colori)');
title('Profilo velocità estratto dalla PIV a \alpha = 10^\circ');
