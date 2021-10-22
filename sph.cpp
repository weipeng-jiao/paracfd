#include <math.h>
#include <fstream>
#include <iostream>
#include <omp.h>
using namespace std;


#define CHUNKSIZE 40  //��̬������������ʹ�õ�ÿ���̵߳Ŀ��С
#define THREAD_NUM 16


// ȫ�ֲ���
int i, j, k, psi, q, qi;
double M_pi = 3.14159265359;
double ad = 96*M_pi/1199;
double c0 = 30.0;

// ���Ӳ���
int n, nb, d, ni, ne, stencil_size, bt;
double  mass, mu, Rho, V, V1, V2, h, space;

// ����
int r, c, rb, cb, n1, n2, row_fact, col_fact;
double Dx, Dy, dx, dy, L, H, l1, l2, h1, h2, gx, gy;

// ʱ��
double dt, t_sim, t_total, steps;

//�����������

//һά����
template <typename T>
void print1D(T a, int cols) {

    for (int jj=0; jj<cols; jj++) {
        cout << a[jj] << " ";
    }
    cout << endl;
}

template <typename T>
void setVal(T a, int cols, double val) {
    for (int jj=0; jj<cols; jj++) {
        a[jj] = val;
    }
}

//��ά����
template <typename T>
void print2D(T a, int rows, int cols) {


    for (int ii=0; ii<rows; ii++) {
        for (int jj=0; jj<cols; jj++) {
            cout << a[ii*cols+jj] << " ";
        }
        cout << endl;
    }
}

template <typename T>
auto slice3D(T a, int* size, int r1, int r2, int c1, int c2, int d1, int d2) {

    int rows = size[0]; int cols = size[1]; int dep = size[2];
    double* res = new double[(r2-r1)*(c2-c1)*(d2-d1)]{0};
    int cnt = 0;
    for (int ii=r1; ii<r2; ii++) {
        for (int jj=c1; jj<c2; jj++) {
            for (int kk=d1; kk<d2; kk++) {
                res[cnt] = a[ (ii*cols+jj)*dep + kk];
                cnt++;
            }
        }
    }
    return res;
}


//��ӡָ���ļ���ָ��
void print2file(double* a, int rows, int cols, string path) {

    ofstream myfile (path, ios_base::app);
    for (int ii=0; ii<rows; ii++) {
        for (int jj=0; jj<cols; jj++) {
            myfile << a[ii*cols+jj] << " ";
        }
        myfile << endl;
    }
    myfile.close();
}
// ��������߽�
bool stencilBoundaryCheck(int rows, int cols) {

    if (rows < 0 || cols < 0)
        return 1;

    else if (rows >= rb || cols >= cb)
        return 1;
    return 0;
}

void indexFixSoft(int& ii, int& jj) {
    
    if (ii == -1) {ii = 0;}
}
//��������Ƿ�����֮��
bool coordBoundaryCheck(double xx, double yy) {
   
    if (yy<0 || yy>=Dy*rb)
        return 1;

    else if (xx<0 || xx>=Dx*cb)
        return 1;

    return 0;
}
//��������Ƿ�����֮��
bool boundaryCheckArray(double* a, double* b, int size, bool print) {
 

    bool res;
    for (int jj=0; jj<size; jj++) {
        res = coordBoundaryCheck(a[jj], b[jj]);
        if (res == 1) {
            if (print == 1) {
                cout << "Particle out of bounds." << endl;
                exit (EXIT_FAILURE);
            }
            return res;
        }
    }
    return res;
}

double kernel(double dr) {
    // Kernel function definition

    double o = abs(dr)/h;

	if (0 <= o && o <= 0.5) {
		return (ad/h)*(pow(5/2-o, 4) - 5*pow(3/2-o, 4) + 10*pow(0.5-o, 4));
	}
	else if (0.5 < o && o <= 1.5) {
		return (ad/h)*(pow(5/2-o, 4) - 5*pow(3/2-o, 4));
	}
	else if (1.5 < o && o <= 2.5) {
		return (ad/h)*(pow(5/2-o, 4));
	}
    return 0;
}
//�˺����ĵ���
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

//��ʼ����������
void initialise(double* a, double* b, double l1, double l2) {


    bool randomIni = 1;
    space = row_fact*Dx;

    if (randomIni == 1) {
        // ����

        int jj = 0; int ii = 0;
        for (int kk=0; kk<n1; kk++) {
            if (jj*dx > l1) {ii++; jj = 0;}
            a[kk] = space + jj*dx;
            b[kk] = space + ii*dx;
            jj++;
        }
        jj = 0; ii = 0;
        for (int kk=n1; kk<n; kk++) {
            if (jj*dx > l2) {ii++; jj = 0;}
            a[kk] = space + jj*dx;
            b[kk] = space + dx/2 + h1 + dx + ii*dx;
            jj++;
        }
    }
    else {
        //�������
        for (int kk = 0; kk<n; kk++) {
            a[kk] = space + randZeroToOne()*L;
            b[kk] = space + randZeroToOne()*H;
        }
    }
}
  // �߽����ӳ�ʼ��
void initialiseBoundaries(double* x, double* y, int nbh, int nbv, int nbh_side, int nbv_side) {
 
    double space2 = row_fact*Dx-(bt-0)*dx;

    int jj = 0; int ii = 0;
    for (int kk=n; kk<n+nbh/2; kk++) {
        if (jj*dx > L+2*bt*dx) {ii++; jj = 0;}
        x[kk] = space2 + jj*dx;
        y[kk] = space2 + ii*dx;
        jj++;
    }
    jj = 0; ii = 0;
    for (int kk=n+nbh/2; kk<n+nbh; kk++) {
        if (jj*dx > L+2*bt*dx) {ii++; jj = 0;}
        x[kk] = space2 + jj*dx;
        y[kk] = space2 + dx + bt*dx + H + ii*dx;  // dx + ??
        jj++;
    }
    jj = 0; ii = 0;
    for (int kk=n+nbh; kk<n+nbh+nbv/2; kk++) {
        if (jj*dx > (bt-1)*dx) {ii++; jj = 0;}
        x[kk] = space2 + jj*dx;
        y[kk] = space2 + bt*dx + ii*dx;
        jj++;
    }
    jj = 0; ii = 0;
    for (int kk=n+nbh+nbv/2; kk<n+nbv+nbh; kk++) {
        if (jj*dx >= 2*bt*dx) {ii++; jj = 0;}
        x[kk] = space2 + dx + bt*dx + L + jj*dx;  // dx + ??
        y[kk] = space2 + bt*dx + ii*dx;
        jj++;
        if (jj%nbv_side == 0) {ii++; jj = 0;}
    }
}
//��ʼ���ļ�
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

//�����Ƿ��ڰ뾶��Χ��
bool circleCheck(double a, double b, double radius) {

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

//����ճ������
double viscTensor(double dvel, double difr, double mui, double muj, double rhoi, double rhoj, double deltaw) {

    double eh = 0.01*h*h;

    return -16*(mui+muj)/(pow(difr, 2)*(rhoi+rhoj) + eh)*(dvel*difr)*deltaw;
}

//����ѹ��
void setPressure(double* a, double* rho, int cols) {

    
    double gamma = 7.0;

    for (int jj=0; jj<cols; jj++) {
        a[jj] = Rho*c0*c0/gamma*(pow(rho[jj]/Rho, gamma) - 1);
    }

}

//����ѹ��
void setPressure2(double* a, double rrho, int cols) {
 
    
    double gamma = 7.0;

    a[cols] = Rho*c0*c0/gamma*(pow(rrho/Rho, gamma) - 1);

}

//�߽紦��
void mirrorBC(double* xx, double* yy, double* uu, double* vv, int col) {

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


int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);

    //֡��
    double fps = 30;
    double print_count = 0;

    //������
    n = 5000;
    //�ܶ�
    Rho = 1000;

    //ʱ�䲽
    dt = 0.0001;
    steps = 40000;
    t_total = dt*steps;
    t_sim = 0.0;

    //�߽���
    bt = 3;
    //�ߴ�
    L = 20; H = 10;

    //���γ�ʼ������
    h1 = 0.2*H; l1 = L;
    h2 = 0.3*H; l2 = 0.15*L;
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

    //ָ���ʼ��
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

    //���Ӻ˵��������
    h = sqrt(2)*dx;
    Dx = 2.5*h; Dy = Dx;
    d = (Dx/dx+1)*(Dx/dx+1); // max particle inside cell
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

    //��ӡ���
    double dt_max = 0.1*h/c0;
    cout << "dt: " << dt << endl;
    cout << "dt_max " << dt_max << endl;
    cout << "dx: " << dx << endl;
    cout << "Dx: " << Dx << endl;
    cout << "Max particles in cell: " << d << endl;

    // ���ճ�ʼ��
    initialise(x, y, l1, l2);
    initialise(x_old, y_old, l1, l2);
    initialiseBoundaries(x, y, nbh, nbv, nbh_side, nbv_side);

    //��������
    mass = Rho * dx * dx;
    mu = 0.001;
    gx = 0; gy = -9.81;

    //�ļ���ʼ��
    string pathx = "x.txt";
    string pathy = "y.txt";
    fileIni(pathx, n+nb); fileIni(pathy, n+nb);

    
    // ʱ�䲽ѭ��
    while (t_sim < t_total) {

        //����ÿ�����ӵ���Χ��Ϣ��ÿ4������һ�Σ�
        if (iter_cnt % 1 == 0) {
            setVal(ndeg, rb*cb, 0);  // ��ʼ��
            for (ni=0; ni<n+nb; ni++) {

                //����й©����
                if (x[ni] == 0 && y[ni] == 0 && ni < n) {continue;}

                //����i��j
                i = (rb-1) - (int)(y[ni]/Dy); j = x[ni]/Dx;
                indexFixSoft(i, j);
                //��������ָ��
                q = (i*cb + j)*d + ndeg[i*cb+j];
                //�����Ӵ洢���ʵ��ĵ�Ԫ��
                jaret[q] = ni;
                ndeg[i*cb+j]++;
            }
        }
        
        // ����ѭ��
        #pragma omp parallel for num_threads(THREAD_NUM) schedule(guided, CHUNKSIZE)
        for (ni=0; ni<n; ni++) {

           //����й©����
            if (x[ni] == 0 && y[ni] == 0) {continue;}

            int i0, j0, q, qi, i, j;
            double weight0 = 0; double weight1 = 0;
            double difx, dify, weightx, weighty, visc_termx, visc_termy, du, dv, dr0, dwr;

            //i j
            i0 = (rb-1) - (int)(y[ni]/Dy); j0 = x[ni]/Dx;

            //������Χ��Ϣ
            for (int psi=0; psi<stencil_size; psi++) {

                // ��Χ��i, j
                i = i0 + stencil_x[psi];
                j = j0 + stencil_y[psi];

                bool out_of_bounds = stencilBoundaryCheck(i, j);
                if (out_of_bounds == 1) {continue;}

                // Find all particle neighbours @ stencil cell
                for (int k=0; k<ndeg[i*cb+j]; k++) {

                    // ��ǰ����ϵ��
                    q = (i*cb + j)*d + k; 

                    // ��Χ����ϵ��
                    qi = jaret[q];

                    // ������ͬ����
                    if (q == qi) {continue;}

                    // ����λ��
                    difx = x[qi] - x[ni];
                    dify = y[qi] - y[ni];

                    //����ھ��Ƿ��ڷ�Χ��
                    bool in_circle = circleCheck(difx, dify, Dx);

                    if (in_circle == 1) {
                        
                        //�ٶȲ�
                        du = u[qi] - u[ni];
                        dv = v[qi] - v[ni];
                        //����
                        dr0 = sqrt(difx*difx + dify*dify);
                        //�˵��ݶ�
                        dwr = gradKernel(dr0);
                        //����
                        weightx = dwr*difx/(dr0+10e-20);
                        weighty = dwr*dify/(dr0+10e-20);

                        // ճ��
                        double mui = mu;
                        double muj = mu;
                        visc_termx = viscTensor(du, dr0, mui, muj, rho[ni], rho[qi], dwr);
                        visc_termy = viscTensor(dv, dr0, mui, muj, rho[ni], rho[qi], dwr);

                        //�����ܺ�
                        double sum_termx = mass*((p[ni]/rho[ni]/rho[ni] + p[qi]/rho[qi]/rho[qi])*weightx + visc_termx);
                        double sum_termy = mass*((p[ni]/rho[ni]/rho[ni] + p[qi]/rho[qi]/rho[qi])*weighty + visc_termy);
                        
                        // ��������
                        // #pragma omp atomic
                        fx[ni] += sum_termx;
                        // #pragma omp atomic
                        fy[ni] += sum_termy;

                        //�����ܶ�
                        // #pragma omp atomic
                        drho[ni] += mass*(du*weightx + dv*weighty);
                        //���±߽��ܶ�
                        if (qi > n) {
                            // #pragma omp atomic
                            rho[qi] += mass*(du*weightx + dv*weighty)*dt;
                        }

                        //�⻬����
                        if (iter_cnt % 10 == 0) {
                            weight0 += kernel(dr0);
                            weight1 += kernel(dr0)/rho[qi];
                        }

                    }
                } 
            } 

            //�⻬����
            if (iter_cnt % 10 == 0) {rho[ni] =  weight0/weight1;}

        } 

        double maxrho = -1.0;
        //ʱ����ָ���
        for (ni=0; ni<n+nb; ni++) {

            if (ni < n) {
                //�����ܶ�
                rho[ni] = rho_old[ni] + drho[ni]*dt;
                double perc = abs(rho[ni]-rho_old[ni])/rho_old[ni]*100;
                if (perc > maxrho) {maxrho = perc;}
                rho_old[ni] = rho[ni];
                drho[ni] = 0;

                //���Ӽ���
                ax[ni] = fx[ni] + gx;
                ay[ni] = fy[ni] + gy;
                fx[ni] = 0; fy[ni] = 0;

               //�����ٶ�
                u[ni] = u_old[ni] + ax[ni]*dt;
                v[ni] = v_old[ni] + ay[ni]*dt;

                //����λ��
                x[ni] = x_old[ni] + u[ni]*dt;
                y[ni] = y_old[ni] + v[ni]*dt;

                // ��ֵ����
                u_old[ni] = u[ni];
                v_old[ni] = v[ni];
                x_old[ni] = x[ni];
                y_old[ni] = y[ni];

                //����Ƿ���й©����
                if (coordBoundaryCheck(x[ni], y[ni]) == 1) {
                
                    x[ni] = 0; 
                    y[ni] = 0;  
                    
                }
            }

            setPressure2(p, rho[ni], ni);

        } //����������


         //���ݵ���
        //��x��y��ӡ���ļ�
        if (t_sim/(1/fps) > print_count) {
            print2file(x, 1, n+nb, pathx);
            print2file(y, 1, n+nb, pathy);
            print_count = print_count + 1;
            cout << "Time elapsed (%): " << t_sim/t_total*100 << endl;
            cout << "Rho change: " << maxrho << endl;
        }

        //��������ʱ��
        t_sim = t_sim + dt;
        iter_cnt++;  //����������

    } 
    //ʱ��ѭ������
    
    
    //�ͷŶ�̬�ڴ�
    delete[] jaret, jaret_shape, stencil_x, stencil_y, ndeg;
    delete[] ax, ay, fx, fy, u, v, u_old, v_old, x, y, x_old, y_old, rho, rho_old, drho, p;

    return 0;
}
