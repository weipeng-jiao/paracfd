//二维流动 tony jiao

#include <math.h>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <time.h>
#include <stdlib.h>
using namespace std;




//int threadnums=omp_get_max_threads();
int threadnums;


//常数
double Pi = 3.14159265359;
double ad = 96*Pi/1199;
double c0 = 30.0;

//系数声明
int i, j, k, psi, q, qi;
//粒子参数声明
int n, nb, d, ni, ne, stencil_size, bt;
double  mass, mu, Rho, V, V1, V2, h, space;
//网格变量声明
int r, c, rb, cb, n1, n2, row_fact, col_fact;
double Dx, Dy, dx, dy, L, H, l1, l2, h1, h2, gx, gy;
//时间声明
double dt, t_sim, t_total, steps;


//输出一维矩阵
template <typename T>
void print1D(T a, int cols) 
{

    for (int jj=0; jj<cols; jj++) 
    {
        cout << a[jj] << " ";
    }
    cout << endl;
}

template <typename T>
void setVal(T a, int cols, double val) 
{
    for (int jj=0; jj<cols; jj++)
     {
        a[jj] = val;
    }
}

//输出二维矩阵
template <typename T>
void print2D(T a, int rows, int cols)
 {
    for (int ii=0; ii<rows; ii++) 
    {
        for (int jj=0; jj<cols; jj++) 
        {
            cout << a[ii*cols+jj] << " ";
        }
        cout << endl;
    }
}

template <typename T>
auto slice3D(T a, int* size, int r1, int r2, int c1, int c2, int d1, int d2) 
{
    int rows = size[0]; int cols = size[1]; int dep = size[2];
    double* res = new double[(r2-r1)*(c2-c1)*(d2-d1)]{0};
    int cnt = 0;
    for (int ii=r1; ii<r2; ii++) 
    {
        for (int jj=c1; jj<c2; jj++) 
        {
            for (int kk=d1; kk<d2; kk++) 
            {
                res[cnt] = a[ (ii*cols+jj)*dep + kk];
                cnt++;
            }
        }
    }
    return res;
}


//打印指向文件的指针
void print2file(double* a, int rows, int cols, string path)
 {
    ofstream myfile (path, ios_base::app);
    for (int ii=0; ii<rows; ii++) 
    {
        for (int jj=0; jj<cols; jj++) 
        {
            myfile << a[ii*cols+jj] << " ";
        }
        myfile << endl;
    }
    myfile.close();
}
//检查计算域边界
bool stencilBoundaryCheck(int rows, int cols)
 {
    if (rows < 0 || cols < 0)
        return 1;

    else if (rows >= rb || cols >= cb)
        return 1;
    return 0;
}

void indexFixSoft(int& ii, int& jj) 
{
    
    if (ii == -1) {ii = 0;}
}
//检查粒子是否在域之外
bool coordBoundaryCheck(double xx, double yy) 
{ 
    if (yy<0 || yy>=Dy*rb)
        return 1;

    else if (xx<0 || xx>=Dx*cb)
        return 1;

    return 0;
}
//检查粒子是否在域之外
bool boundaryCheckArray(double* a, double* b, int size, bool print) 
{
    bool res;
    for (int jj=0; jj<size; jj++) {
        res = coordBoundaryCheck(a[jj], b[jj]);
        if (res == 1) {
            if (print == 1) {
                cout << "粒子超出边界" << endl;
                exit (EXIT_FAILURE);
            }
            return res;
        }
    }
    return res;
}
//定义核函数
double kernel(double dr) 
{
   double o = abs(dr)/h;

	if (0 <= o && o <= 0.5) 
	{
		return (ad/h)*(pow(5/2-o, 4) - 5*pow(3/2-o, 4) + 10*pow(0.5-o, 4));
	}
	else if (0.5 < o && o <= 1.5)
	{
		return (ad/h)*(pow(5/2-o, 4) - 5*pow(3/2-o, 4));
	}
	else if (1.5 < o && o <= 2.5) 
	{
		return (ad/h)*(pow(5/2-o, 4));
	}
    return 0;
}
//核函数导数
double gradKernel(double dr) {
 
  double o = abs(dr)/h;

	if (0 <= o && o <= 0.5) {
		return (ad/h)*(-4/h*pow(5/2-o, 3) + 20/h*pow(3/2-o, 3) - 40/h*pow(0.5-o, 3));
	}
	else if (0.5 < o && o <= 1.5) {
		return (ad/h)*(-4/h*pow(5/2-o, 3) + 20/h*pow(3/2-o, 3));
	}
	else if (1.5 < o && o <= 2.5) {
		return (ad/h)*(-4/h*pow(5/2-o, 3));
	}
    return 0;
}

double randZeroToOne() {
    return rand() / (RAND_MAX + 1.);
}

//初始化流动粒子
void initialise(double* a, double* b, double l1, double l2) {


    bool randomIni = 1;
    space = row_fact*Dx;

    if (randomIni == 1) 
    {

        int jj = 0; int ii = 0;
        for (int kk=0; kk<n1; kk++) 
        {
            if (jj*dx > l1) {ii++; jj = 0;}
            a[kk] = space + jj*dx;
            b[kk] = space + ii*dx;
            jj++;
        }
        jj = 0; ii = 0;
        for (int kk=n1; kk<n; kk++) 
        {
            if (jj*dx > l2) {ii++; jj = 0;}
            a[kk] = space + jj*dx;
            b[kk] = space + dx/2 + h1 + dx + ii*dx;
            jj++;
        }
    }
    else {
         //其他粒子
        for (int kk = 0; kk<n; kk++) {
            a[kk] = space + randZeroToOne()*L;
            b[kk] = space + randZeroToOne()*H;
        }
    }
}
  // 边界粒子初始化
void initialiseBoundaries(double* x, double* y, int nbh, int nbv, int nbh_side, int nbv_side) {
 
    double space2 = row_fact*Dx-(bt-0)*dx;

    int jj = 0; int ii = 0;
    for (int kk=n; kk<n+nbh/2; kk++) 
    {
        if (jj*dx > L+2*bt*dx) {ii++; jj = 0;}
        x[kk] = space2 + jj*dx;
        y[kk] = space2 + ii*dx;
        jj++;
    }
    jj = 0; ii = 0;
    for (int kk=n+nbh/2; kk<n+nbh; kk++) 
    {
        if (jj*dx > L+2*bt*dx) {ii++; jj = 0;}
        x[kk] = space2 + jj*dx;
        y[kk] = space2 + dx + bt*dx + H + ii*dx;  
        jj++;
    }
    jj = 0; ii = 0;
    for (int kk=n+nbh; kk<n+nbh+nbv/2; kk++) 
    {
        if (jj*dx > (bt-1)*dx) {ii++; jj = 0;}
        x[kk] = space2 + jj*dx;
        y[kk] = space2 + bt*dx + ii*dx;
        jj++;
    }
    jj = 0; ii = 0;
    for (int kk=n+nbh+nbv/2; kk<n+nbv+nbh; kk++) 
    {
        if (jj*dx >= 2*bt*dx) {ii++; jj = 0;}
        x[kk] = space2 + dx + bt*dx + L + jj*dx;  
        y[kk] = space2 + bt*dx + ii*dx;
        jj++;
        if (jj%nbv_side == 0) {ii++; jj = 0;}
    }
}
//初始化文件
void fileIni(string path, int size) {

    int no = 12;
    ofstream ofs;
    ofs.open(path, ofstream::out | ofstream::trunc);
    ofstream myfile (path, ios_base::app);
    myfile << dt << " " << steps << " " << t_total << " ";
    myfile << Dx << " " << Dy << " " << dx << " ";
    myfile << bt << " " << L << " " << H << " " << n << " " << nb << " " << n1 << " ";
    for (int jj=0; jj<size-no; jj++) {myfile << 0 << " ";}
    myfile << endl;
    ofs.close();
}

//检查点是否在半径范围内
bool circleCheck(double a, double b, double radius)
{

    a = abs(a); b = abs(b);
    if (a + b <= radius)
        return 1;
    if (a > radius)
        return 0;
    if (b > radius)
        return 0;
    if (a*a + b*b <= radius*radius)
        return 1;
    else
        return 0;
}

//计算粘滞张量
double viscTensor(double dvel, double difr, double mui, double muj, double rhoi, double rhoj, double deltaw) 
{

    double eh = 0.01*h*h;

    return -16*(mui+muj)/(pow(difr, 2)*(rhoi+rhoj) + eh)*(dvel*difr)*deltaw;
}

//计算压力
void setPressure(double* a, double* rho, int cols) 
{
    double gamma = 7.0;

    for (int jj=0; jj<cols; jj++) {
        a[jj] = Rho*c0*c0/gamma*(pow(rho[jj]/Rho, gamma) - 1);
    }

}

//计算压力
void setPressure2(double* a, double rrho, int cols)
{
    double gamma = 7.0;
    a[cols] = Rho*c0*c0/gamma*(pow(rrho/Rho, gamma) - 1);
}

//边界处理
void mirrorBC(double* xx, double* yy, double* uu, double* vv, int col) 
{
	if (xx[col] < Dx) {
		xx[col] = -xx[col] / 2.0;
		uu[col] = -uu[col];
	}
	if (yy[col] < Dx) {
		yy[col] = -yy[col] / 2.0;
		vv[col] = -vv[col];
	}
	if (xx[col] >= Dx+L) {
		xx[col] = Dx+L - xx[col] + Dx+L;
		uu[col] = -uu[col];
	}
	if (yy[col] >= Dy+H) {
		yy[col] = Dy+H - yy[col] + Dy+H;
		vv[col] = -vv[col];
	}
}


int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);


cout<<"请输入线程数量： " <<endl;

cin>>threadnums;
cout<<endl;

    //帧数
    double fps = 30;
    //输出次数
    double print_count = 0;

//*******************粒子数量******************
  
  cout<<"请输入粒子数量： " <<endl;

cin>>n;
cout<<endl;
  
    //粒子数
   // n = 500;
    //密度
    Rho = 1000;


 //******************时间步*********************
    dt = 0.0001;
    steps = 40000;
    t_total = dt*steps;
    t_sim = 0.0;


 //******************计算域设置********************
    //边界厚度
    bt = 3;
    //尺寸
    L = 20; H = 10;

    //几何初始化参数
    h1 = 0.2*H; l1 = 0.15*L;
    h2 = 0.3*H; l2 = 0.15*L;
    //加速度设置
    gx = 0;
    gy = -9.81;
    
 //********************初始化*********************  
    
    V1 = h1*l1; V2 = h2*l2;
    V = V1 + V2;
    double ratio = V1/V;
    n1 = ratio*n; n2 = n-n1;
    double factor = L*H/V;
    int nx = sqrt(L/H*n*factor + pow((L-H), 2)/(4*H*L)) - (L-H)/(2*H);
    int ny = n*factor/nx;
    dx = L/(nx-1);
    int nbh_side = (L/dx+1+2*bt);
    int nbh = bt*nbh_side*2;
    int nbv_side = bt;
    int nbv = (H/dx+1)*nbv_side*2;
     nb = nbh+nbv;

    //指针初始化
    double* ax = new double[n+nb]{0};
    double* ay = new double[n+nb]{0};
    double* fx = new double[n+nb]{0};
    double* fy = new double[n+nb]{0};
    double* u = new double[n+nb]{0};
    double* v = new double[n+nb]{0};
    double* u_old = new double[n+nb]{0};
    double* v_old = new double[n+nb]{0};
    double* x = new double[n+nb]{0};
    double* y = new double[n+nb]{0};
    double* x_old = new double[n+nb]{0};
    double* y_old = new double[n+nb]{0};
    double* rho = new double[n+nb]{0};
    double* rho_old = new double[n+nb]{0};
    double* drho = new double[n+nb]{0};
    double* p = new double[n+nb]{0};
    setVal(rho, n+nb, Rho);
    setVal(rho_old, n+nb, Rho);
    setPressure(p, rho, n+nb);

    //粒子核的网格参数
    h = sqrt(2)*dx;
    Dx = 2.5*h; Dy = Dx;
    d = (Dx/dx+1)*(Dx/dx+1); 
    r = H/Dy+0.5; c = L/Dx+0.5;
    row_fact = Dy / (dx*(bt-1)) + 0;
    col_fact = Dx / (dx*(bt-1)) + 0;
    rb = r+2*row_fact; cb = c+2*col_fact;
    int* jaret = new int[rb*cb*d]{-1};
    int* jaret_shape = new int[3]{rb, cb, d};
    int* ndeg = new int[rb*cb]{0};
    stencil_size = 9;
    int* stencil_x = new int[stencil_size]{-1, -1, -1, 0, 0, 0, 1, 1, 1};
    int* stencil_y = new int[stencil_size]{-1, 0, 1, -1, 0, 1, -1, 0, 1};
    long long iter_cnt = 0;

    double dt_max = 0.1*h/c0;


    //最终初始化
    initialise(x, y, l1, l2);
    initialise(x_old, y_old, l1, l2);
    initialiseBoundaries(x, y, nbh, nbv, nbh_side, nbv_side);

    //流体属性
    mass = Rho * dx * dx;
    mu = 0.001;

    //文件初始化
    string pathx = "x.txt";
    string pathy = "y.txt";
    fileIni(pathx, n+nb); 
    fileIni(pathy, n+nb);


    
    cout<<"初始化完成，开始计算"<<endl;
//****************************************计算************************* 
    clock_t startime=clock();
        
    
    // 时间步循环
    while (t_sim < t_total)
    {
 
        //计算每个粒子的周围信息（每4步更新一次）
        if (iter_cnt % 1 == 0) {
            setVal(ndeg, rb*cb, 0);  // 初始化
            for (ni=0; ni<n+nb; ni++) {

                //跳过泄漏粒子
                if (x[ni] == 0 && y[ni] == 0 && ni < n) {continue;}

                //计算i，j
                i = (rb-1) - (int)(y[ni]/Dy); j = x[ni]/Dx;
                indexFixSoft(i, j);
                //计算粒子指数
                q = (i*cb + j)*d + ndeg[i*cb+j];
                //将粒子存储在适当的网格单元中
                jaret[q] = ni;
                ndeg[i*cb+j]++;
            }
        }
        
       
        clock_t startime1=clock();
        // 粒子循环
        //实验结果：大循环进行多线程向导式调度性能更好些
     #pragma omp parallel for num_threads(threadnums) schedule(guided)
        for (ni=0; ni<n; ni++)
        {
           
            if (x[ni] == 0 && y[ni] == 0) 
            {
              continue;
            }

            int i0, j0, q, qi, i, j;
            double weight0 = 0; 
            double weight1 = 0;
            double difx, dify, weightx, weighty, visc_termx, visc_termy, du, dv, dr0, dwr;

     
            i0 = (rb-1) - (int)(y[ni]/Dy); 
            j0 = x[ni]/Dx;

            //计算周围信息
            for (int psi=0; psi<stencil_size; psi++) 
            {
                i = i0 + stencil_x[psi];
                j = j0 + stencil_y[psi];

                bool out_of_bounds = stencilBoundaryCheck(i, j);
                if (out_of_bounds == 1) 
                {
                  continue;
                }

              
                for (int k=0; k<ndeg[i*cb+j]; k++) 
                {

                    // 当前粒子系数
                    q = (i*cb + j)*d + k; 

                    // 周围粒子系数
                    qi = jaret[q];

                    // 跳过相同粒子
                    if (q == qi) 
                      {
                        continue;
                      }

                    // 计算位移
                    difx = x[qi] - x[ni];
                    dify = y[qi] - y[ni];

                    //检查粒子是否在范围内
                    bool in_circle = circleCheck(difx, dify, Dx);

                    if (in_circle == 1) 
                    {
                        
                        //速度差
                        du = u[qi] - u[ni];
                        dv = v[qi] - v[ni];
                        //距离
                        dr0 = sqrt(difx*difx + dify*dify);
                        //核的梯度
                        dwr = gradKernel(dr0);
                        //权重
                        weightx = dwr*difx/(dr0+10e-20);
                        weighty = dwr*dify/(dr0+10e-20);

                        // 粘性
                        double mui = mu;
                        double muj = mu;
                        visc_termx = viscTensor(du, dr0, mui, muj, rho[ni], rho[qi], dwr);
                        visc_termy = viscTensor(dv, dr0, mui, muj, rho[ni], rho[qi], dwr);

                        //计算总和
                        double sum_termx = mass*((p[ni]/rho[ni]/rho[ni] + p[qi]/rho[qi]/rho[qi])*weightx + visc_termx);
                        double sum_termy = mass*((p[ni]/rho[ni]/rho[ni] + p[qi]/rho[qi]/rho[qi])*weighty + visc_termy);
                        
                        // 外力计算
                        #pragma omp atomic
                        fx[ni] += sum_termx;
                        #pragma omp atomic
                        fy[ni] += sum_termy;

                        //更新密度
                        #pragma omp atomic
                        drho[ni] += mass*(du*weightx + dv*weighty);
                        //更新边界密度
                        if (qi > n) 
                        {
                            #pragma omp atomic
                            rho[qi] += mass*(du*weightx + dv*weighty)*dt;
                        }

                        //光滑处理
                        if (iter_cnt % 10 == 0) 
                        {
                            weight0 += kernel(dr0);
                            weight1 += kernel(dr0)/rho[qi];
                        }

                    }
                } 
            } 

            //光滑处理
            if (iter_cnt % 10 == 0) 
            {
              rho[ni] =  weight0/weight1;
            }

        } 

        double maxrho = -1.0;
        
      #pragma omp parallel for num_threads(threadnums) schedule(guided)
        //时间积分更新
        for (ni=0; ni<n+nb; ni++)
        {

            if (ni < n) 
            {
                //粒子密度
                rho[ni] = rho_old[ni] + drho[ni]*dt;
                double perc = abs(rho[ni]-rho_old[ni])/rho_old[ni]*100;
                if (perc > maxrho) {maxrho = perc;}
                rho_old[ni] = rho[ni];
                drho[ni] = 0;

                //粒子加速
                ax[ni] = fx[ni] + gx;
                ay[ni] = fy[ni] + gy;
                fx[ni] = 0; fy[ni] = 0;

               //粒子速度
                u[ni] = u_old[ni] + ax[ni]*dt;
                v[ni] = v_old[ni] + ay[ni]*dt;

                //粒子位置
                x[ni] = x_old[ni] + u[ni]*dt;
                y[ni] = y_old[ni] + v[ni]*dt;

                //旧值更新
                u_old[ni] = u[ni];
                v_old[ni] = v[ni];
                x_old[ni] = x[ni];
                y_old[ni] = y[ni];

                //检查粒子越界
                if (coordBoundaryCheck(x[ni], y[ni]) == 1) 
                {
                
                    x[ni] = 0; 
                    y[ni] = 0;  
                    
                }
            }

            setPressure2(p, rho[ni], ni);

        } //粒子求解结束
          clock_t endtime1=clock();
         

         //数据导出
        //将x、y打印到文件
        if (t_sim/(1/fps) > print_count) 
        {
            print2file(x, 1, n+nb, pathx);
            print2file(y, 1, n+nb, pathy);
            print_count = print_count + 1;
            cout << "已完成（%）： " << t_sim/t_total*100 << "  时间："<<(endtime1-startime1)<< "ms" <<endl;
        }

        //更新运行时间
        t_sim = t_sim + dt;
        iter_cnt++; 
        

    } 
    //时间循环结束
    clock_t endtime=clock();
    cout<<"总的时间："<<(endtime-startime)/1000<< "s" <<endl;
    cout<<"迭代次数："<<iter_cnt<<endl;
    //释放动态内存
    delete[] jaret, jaret_shape, stencil_x, stencil_y, ndeg;
    delete[] ax, ay, fx, fy, u, v, u_old, v_old, x, y, x_old, y_old, rho, rho_old, drho, p;


system("pause");
    return 0;
}
