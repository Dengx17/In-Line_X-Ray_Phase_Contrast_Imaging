%% Betatron辐射光源同轴相衬成像数值模拟 %%
clc;clear;
tic()
%% Betatron辐射光源（光子能量为1-100keV）及成像系统参数
E = linspace(1,100,100);
lambda = linspace(1,100,100);
k = linspace(1,100,100);
for m = 1:100
    lambda(m) = (1.241E-6)/E(m); % mm 波长=hc/E
    k(m) = 2*pi/lambda(m); % /mm 波矢
end
r1 = 1000; % mm 物距
r2 = 1000; % mm 像距
M = (r1+r2)/r1; % 放大倍数
z = r2/M; % 点光源等效于平行光时的等效物距

%% 物平面
D = 0.1; % mm 圆柱半径
h = 1; % mm 圆柱高度
lx = 2; % mm 取样长度
ly = 0.2; % mm 取样宽度
N = 500; % 取样点数
% 以圆柱体中心为原点，构建直角坐标系
x0 = linspace(-lx/2,lx/2,N);
y0 = linspace(-ly/2,ly/2,N);
th = linspace(0,0,N);
for m = 1:N
    if(abs(y0(m))<= D/2)
        th(m) = 2*sqrt((D/2).^2-(y0(m)).^2); % 厚度与纵坐标的关系
    end
end
th = th';
P = zeros(N);
for n = 1:N
    if(abs(x0(n)) <= h/2)
        P(:,n) = th; % 绘制物平面
    end
end
%figure;imagesc(x0,y0,P)
%plot(y0,th)

% Index of Refraction = (1-delta)-i(beta)
% C7H5N3O6 Density=1.654 @PE = 10keV
delta = 3.51828226E-06;
beta = 6.33161923E-09; 
for m = 1:100
    mu(:,:,m) = (k(m)*beta).*P; % 吸收系数
    phi(:,:,m) = (-k(m)*delta).*P; % 相移系数
end
u0 = exp(-mu+i*phi); % 物平面的复振幅
I0 = zeros(N);
for m = 1:100
    U_0(:,:,m) = fftshift(fft2(u0(:,:,m))); % Fourier变换后的频域上的物平面
    I0 = I0 + abs(u0(:,:,m)).^2; % 物平面光强
end
figure('NumberTitle','off','Name','物平面光强');imagesc(x0,y0,I0)
figure('NumberTitle','off');plot(y0,I0(:,N/2));grid on;title('物平面沿y轴光强变化曲线')

%% 系统传输矩阵
% 频域坐标
x1 = linspace(-N/2/lx,N/2/lx,N);
y1 = linspace(-N/2/ly,N/2/ly,N);
y1 = y1'; % or [x y] = meshgrid(x1,y1)
for m = 1:100
    H(:,:,m) = exp(i*k(m)*z-i*pi*lambda(m)*z*(x1.^2+y1.^2));
end

%% 像平面
U_1 = U_0.*H; % Fourier变换后的频域上的像平面
I1 = zeros(N);
for m = 1:100
    u1(:,:,m) = ifft2(U_1(:,:,m)); % 像平面的复振幅
    I1 = I1 + abs(u1(:,:,m)).^2; % 像平面光强
end
figure('NumberTitle','off','Name','物像光强对比');imagesc(x0*M,y0*M,I1)
figure('NumberTitle','off');plot(y0*M,I1(:,N/2));grid on;title('像平面沿y轴光强变化曲线')
toc()
%{
figure('NumberTitle','off','Name','物像光强对比')
subplot(2,2,1)
imagesc(x0,y0,I0)
subplot(2,2,2)
plot(y0,I0(:,N/2))
subplot(2,2,3)
imagesc(x0*M,y0*M,I1)
subplot(2,2,4)
plot(y0*M,I1(:,N/2))
%}