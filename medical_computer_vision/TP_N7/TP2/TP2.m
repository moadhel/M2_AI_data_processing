load("scatters.mat")
colormap("gray")
figure()
imagesc(scatters)
y=abs(hilbert(scatters));
figure()
hold on
plot(scatters(:,60) )
plot(abs(hilbert(scatters(:,60))))

figure()
hold on
plot(abs(hilbert(scatters)))
hold off
figure()
hist(y(:),64);
figure()
hist(log(y(:)),64)
%%
y1=abs(hilbert(us));
denoise=imgaussfilt(us);
y1_denoise=abs(hilbert(denoise));
load("us.mat")
figure()
colormap("gray")
imagesc(us)
figure()
hold on
plot(us(:,10) )
plot(abs(hilbert(us(:,10))))
figure()
hold on
plot(abs(hilbert(us)))
hold off
figure()
hist(y1(:),64);
figure()
hist(log(5+6*y1(:)),64)
figure()
colormap("gray")
y1_denoise1=imgaussfilt(log(5+6*y1(:)));
imagesc(y1_denoise1)

