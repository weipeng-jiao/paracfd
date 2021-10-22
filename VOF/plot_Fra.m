function [] = plot_Fra( F,istep,xm,ym,imax,imin,jmax,jmin,dt)

fig = figure('Visible','off'); 
pcolor(xm(imin:imax),ym(jmin:jmax),F(imin:imax,jmin:jmax)') 
frame = getframe(fig); 
img = frame2im(frame); 
imwrite(img,['./image/','ʱ�䲽����',num2str(istep,'%06d\n'),'.png']); "
fig2 = figure('Visible','off'); 
contourf(xm(imin:imax),ym(jmin:jmax),F(imin:imax,jmin:jmax)',1);
frame2 = getframe(fig2); 
img = frame2im(frame2); 
imwrite(img,['./result/','ʱ�䲽����',num2str(istep,'%06d\n'),'.png']); 
fprintf('����ͼƬ:%s\n \n',num2str(istep));
end
