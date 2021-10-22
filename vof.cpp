#include<iostream>
#include<cmath>
using namespace std;



double** Matrix(int m,int n)
{
	double** matrix = new double*[m];
	for(int i=0; i<m; ++i)
	{
		matrix[i] = new double[n];
	}
	return matrix;
}


int main()
{
//输入网格参数
int nx, ny, Lx, Ly;
double rho_water, rho_air, nu_water, nu_air, gx, gy;
double x1, y1, x2, y2;
double  dt;
nx=32;          //x方向网格点数
ny=32;          //y方向网格点数
Lx=1;           //区域长度
Ly=1;           //区域高度
//输入物理参数
rho_water=1;    //水的密度
rho_air=0.001;    //空气密度
nu_water=0.01;  //运动粘性系数
nu_air=0.005;   //运动粘性系数
gx=0;
gy=-1;
//初始体积分数区域（方形区域）
x1=1/3;
x2=2/3;
y1=0.5;
y2=1;
//求解参数
dt=0.001;
//======  ===================生成网格================================================
int imax, imin, jmax, jmin, dxi, dyi;
imin=2; imax=imin+nx-1;
jmin=2; jmax=jmin+ny-1;

x(imin : imax+1)=linspace (0 ,Lx, nx+1);
y(jmin : jmax+1)=linspace (0 ,Ly, ny+1);

xm(imin : imax)=0.5*(x(imin : imax)+x(imin+1:imax+1));
ym(jmin : jmax)=0.5*(y(jmin : jmax)+y(jmin+1:jmax+1));

dx=x(imin+1)-x(imin );
dy=y(jmin+1)-y(jmin );

dxi=1/dx ;
dyi=1/dy ;
//参数设置

//初始条件
u=zeros(imax+1,jmax+1);
v=zeros(imax+1,jmax+1);
p=zeros(imax,jmax);
F=zeros(imax+1,jmax+1);
//设置初始体积分数
for j=1:1:jmax
    for i=1:1:imax
        if(xm(i)>=x1)&&(xm(i)<=x2)&&(ym(j)>=y1)&&(ym(j)<=y2)
            F(i,j)=1;
        end
    end
end
//初始化rho，mu
rho=zeros(imax+1,jmax+1);
mu=zeros(imax+1,jmax+1);
//=================================================
// 创建Laplace算子
L=Laplace_operator(nx,ny,dxi,dyi);
//=====================================================
istep=0;
istep_max=2000;
Residual=zeros(istep_max/100,3);
check_mass=zeros(istep_max/100,1);//检查质量
R_limit=15;
count=0;
while (istep<istep_max)
    //设置边界条件
    [u,v,F]=set_BC(u,v,F);
    //由F更新rho，mu
    [rho,mu] = cal_mu_rho( F,rho,mu,imin,imax,jmin,jmax);
    //投影法求解压力possion方程
    [un,vn,pn] = M_Possion(L,u,v,p,mu,rho);
    //显示，差分法求解F
    [Fn] = solve_F( un,vn,F);
    for j=1:1:jmax
        for i=1:1:imax
            Fn(i,j)=var(0,1,Fn(i,j));
        end
    end
    istep=istep+1;//时间步+1
    u_R=norm(un-u)/(nx*ny);
    v_R=norm(vn-v)/(nx*ny);
    p_R=norm(pn-p)/(nx*ny);
    u=un;
    v=vn;
    p=pn;
    F=Fn;
    [u,v,Fn]=set_BC(u,v,Fn);
    if (u_R>R_limit)||(v_R>R_limit)
        fprintf('残差过大不收敛');
        break;
    end
   
   
   return 0;
}