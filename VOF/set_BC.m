function [u,v,F] = set_BC(u,v,F,imax,imin,jmax,jmin)
v(:,jmin-1)=v(:,jmin);
v(:,jmax+1)=v(:,jmax);
u(imin-1,:)=u(imin,:);
u(imax+1,:)=u(imax,:);
%为F设置不可穿透条件
F(imin-1,:)=F(imin,:);%left
F(imax+1,:)=F(imax,:);%right
F(:,jmin-1)=F(:,jmin);%bottom;
F(:,jmax+1)=F(:,jmax);%top
end
