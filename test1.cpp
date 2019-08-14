#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#define M 101													//�ڵ���
#define N 100													//������
double phi[M][M], phi1[M][M], Ni[M][M], E[M][M], Ex[M][M], Ey[M][M], Num[N][N], D[N][N];		//�ڵ����ֵphi,�ڵ����ܶ�Ni,�糡E,������Num,�������ܶ�D
double vx[N * N * 10], vy[N * N * 10];                       			//�����ٶ�v
double Px[N * N * 10], Py[N * N * 10], W[4];								//����λ��P,���Ӷ���������ڵ��Ȩ��ϵ��W
double kB = 1.38e-23, T = 500.0, mi = 2.18e-25;							//������������kB,�����¶�T,��������M
double ep0 = 8.854e-12, e = 1.6e-19, ni = 10e15;						//��ս�糣��ep0��Ԫ���e���������ܶ�ni
double vth, pi = 3.1415926;										//��������ѧ�ٶ�vth,Բ����pi
double dx = 1e-7, dt = 1e-11;										//������dx�����沽��dt
int k = 0, l = 0;
void initial()													//��ʼ������λ�ú��ٶ�
{
	int i, j, k, n = 0;
	double R1, R2;
	//vth = pow(kB * T / mi, 0.5);
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			for (k = 0; k < 10; k++)
			{
				R1 = (rand() % 1000) / 999.0;
				R2 = (rand() % 1000) / 999.0;
				Px[n] = double(j) * dx + R1 * dx;
				Py[n] = double(i) * dx + R2 * dx;
				vx[n] = 0;
				vy[n] = 0;
				n++;
			}
		}
	}
}
void weighting()												//����Ȩ��ϵ��������ܶȷ��䵽����
{
	int i, j, x, y;
	double Si, Si1, Sj, Sj1;
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < M; j++)
		{
			Ni[i][j] = 0;
		}
	}
	for (i = 0; i < N * N * 10; i++)
	{
		x = int(Px[i] / dx);
		y = int(Py[i] / dx);
		Si = (Py[i] - y * dx) / dx;
		Si1 = ((y + 1.0) * dx - Py[i]) / dx;
		Sj = (Px[i] - x * dx) / dx;
		Sj1 = ((x + 1.0) * dx - Px[i]) / dx;
		W[0] = Si * Sj;
		W[1] = Si1 * Sj;
		W[2] = Si * Sj1;
		W[3] = Si1 * Sj1;
		Ni[x][y] += W[0] * e / ep0 * ni;
		Ni[x + 1][y] += W[1] * e / ep0 * ni;
		Ni[x][y + 1] += W[2] * e / ep0 * ni;
		Ni[x + 1][y + 1] += W[3] * e / ep0 * ni;
	}
}
void E_field()                                                  //����ڵ�糡
{
	int i, j;
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < M; j++)
		{
			if (i == 0)                                            //���ñ߽�ڵ�糡�������
			{
				if (j == 0)
				{
					Ex[i][j] = -(phi[i][j + 1] - phi[i][j]) / (2 * dx);
					Ey[i][j] = -(phi[i + 1][j] - phi[i][j]) / (2 * dx);
				}
				else if (j == M - 1)
				{
					Ex[i][j] = -(phi[i][j] - phi[i][j - 1]) / (2 * dx);
					Ey[i][j] = -(phi[i + 1][j] - phi[i][j]) / (2 * dx);
				}
				else
				{
					Ex[i][j] = -(phi[i][j + 1] - phi[i][j - 1]) / (2 * dx);
					Ey[i][j] = -(phi[i + 1][j] - phi[i][j]) / (2 * dx);
				}
			}
			else if (i == M - 1)
			{
				if (j == 0)
				{
					Ex[i][j] = -(phi[i][j + 1] - phi[i][j]) / (2 * dx);
					Ey[i][j] = -(phi[i][j] - phi[i - 1][j]) / (2 * dx);
				}
				else if (j == M - 1)
				{
					Ex[i][j] = -(phi[i][j] - phi[i][j - 1]) / (2 * dx);
					Ey[i][j] = -(phi[i][j] - phi[i - 1][j]) / (2 * dx);
				}
				else
				{
					Ex[i][j] = -(phi[i][j + 1] - phi[i][j - 1]) / (2 * dx);
					Ey[i][j] = -(phi[i][j] - phi[i - 1][j]) / (2 * dx);
				}
			}
			else if (j == 0)
			{
				Ex[i][j] = -(phi[i][j + 1] - phi[i][j]) / (2 * dx);
				Ey[i][j] = -(phi[i + 1][j] - phi[i - 1][j]) / (2 * dx);
			}
			else if (j == M - 1)
			{
				Ex[i][j] = -(phi[i][j] - phi[i][j - 1]) / (2 * dx);
				Ey[i][j] = -(phi[i + 1][j] - phi[i - 1][j]) / (2 * dx);
			}
			else
			{
				Ex[i][j] = -(phi[i][j + 1] - phi[i][j - 1]) / (2 * dx);
				Ey[i][j] = -(phi[i + 1][j] - phi[i - 1][j]) / (2 * dx);
			}
		}
	}
}
void movement()								//�����˶�
{
	int i, x, y;
	double Si, Si1, Sj, Sj1;
	for (i = 0; i < N * N * 10; i++)
	{

		if (Px[i] <= 0 || Px[i] >= 100 * dx || Py[i] <= 0 || Py[i] >= 100 * dx)				//�߽���������
		{
			vx[i] = 0;
			vy[i] = 0;
			continue;
		}
		x = int(Px[i] / dx);
		y = int(Py[i] / dx);
		Si = (Py[i] - y * dx) / dx;
		Si1 = ((y + 1.0) * dx - Py[i]) / dx;
		Sj = (Px[i] - x * dx) / dx;
		Sj1 = ((x + 1.0) * dx - Px[i]) / dx;
		W[0] = Si * Sj;
		W[1] = Si1 * Sj;
		W[2] = Si * Sj1;
		W[3] = Si1 * Sj1;
		vx[i] += e / mi * (Ex[x][y] * W[0] + Ex[x + 1][y] * W[1] + Ex[x][y + 1] * W[2] + Ex[x + 1][y + 1] * W[3]) * dt;
		vy[i] += e / mi * (Ey[x][y] * W[0] + Ey[x + 1][y] * W[1] + Ey[x][y + 1] * W[2] + Ey[x + 1][y + 1] * W[3]) * dt;
		Px[i] += vx[i] * dt;
		Py[i] += vy[i] * dt;
	}
}
void density()										//�����������ܶ�
{
	int i, j, x, y;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			Num[i][j] = 0;
		}
	}
	for (i = 0; i < N * N * 10; i++)
	{
		x = int(Px[i] / dx);
		y = int(Py[i] / dx);
		Num[x][y]++;
	}
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			D[i][j] = Num[i][j] / (dx * dx);
		}
	}
}
int main()
{
	int i, j;
	double w;
	using namespace std;
	ofstream fp("position.dat", ofstream::out);
	ofstream fp1("density.dat", ofstream::out);
	ofstream fp2("potential.dat", ofstream::out);
	ofstream fp3("electric field.dat", ofstream::out);
	ofstream fp4("speed.dat", ofstream::out);
	initial();
	w = 2 / (1 + sqrt(1 - pow((cos(pi / M) + cos(pi / N)) / 2, 2)));
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < M; j++)
		{
			phi[i][j] = 0;
			phi1[i][j] = 0;
		}
	}
	phi[50][50] = -30;
	phi1[50][50] = -30;
	for (k = 0; k < 100; k++)
	{
		weighting();
		double t, tmax = 1;
		while (tmax > 0.001)
		{
			l++;
			tmax = 0;
			for (i = 1; i < M - 1; i++)
			{
				for (j = 1; j < M - 1; j++)
				{
					if (i == 50 && j == 50)
						continue;
					phi[i][j] = phi1[i][j] + (w / 4) * (phi[i - 1][j] + phi1[i + 1][j] + phi[i][j - 1] + phi1[i][j + 1] - 4 * phi1[i][j]) - pow(dx, 2) / 4 * Ni[i][j];
					t = fabs(phi[i][j] - phi1[i][j]);
					if (t > tmax)
						tmax = t;
				}
			}
			for (i = 1; i < M - 1; i++)
			{
				for (j = 1; j < M - 1; j++)
				{
					phi1[i][j] = phi[i][j];
				}
			}
		}
		E_field();
		movement();
	}
	density();
	fp1 << "title=\"Density\"" << endl << "variables = \"x\",\"y\",\"num\"" << endl << "zone i=100,j=100,f=point" << endl;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			fp1 << j + 1 << " " << i + 1 << " " << D[i][j] << endl;
		}
	}
	fp2 << "title=\"Potential\"" << endl << "variables = \"x\",\"y\",\"V\"" << endl << "zone i=101,j=101,f=point" << endl;
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < M; j++)
		{
			fp2 << j + 1 << " " << i + 1 << " " << phi[i][j] << endl;
		}
	}
	fp3 << "title=\"Electric Field\"" << endl << "variables = \"x\",\"y\",\"Ex\",\"Ey\",\"E\"" << endl << "zone i=101,j=101,f=point" << endl;
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < M; j++)
		{
			E[i][j] = sqrt(pow(Ex[i][j], 2) + pow(Ey[i][j], 2));
			fp3 << j + 1 << " " << i + 1 << " " << Ex[i][j] << " " << Ey[i][j] << " " << E[i][j] << endl;
		}
	}
	fp4 << "title=\"Speed\"" << endl << "variables = \"x\",\"y\",\"vx\",\"vy\"" << endl << "zone i=100000,f=point" << endl;
	for (i = 0; i < N * N * 10; i++)
	{
		fp4 << Py[i] << " " << Px[i] << " " << vx[i] << " " << vy[i] << endl;
	}
	printf("k=%d l=%d\n", k, l);
	fp << "title=\"Position\"" << endl << "variables=\"x\",\"y\"" << " " << "zone i=" << N * N * 10 << endl;
	for (i = 0; i < N * N * 10; i++)
	{
		fp << Px[i] << " " << Py[i] << endl;
	}
	fp.close();
	fp1.close();
	fp2.close();
	fp3.close();
	fp4.close();
}