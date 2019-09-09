#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#include<iomanip>
#include<iostream>
#define M 101															//�ڵ���
#define N 100															//������
#define Nrel 1e4
double phi[M][M], Ni[M][M], Ne[M][M], Ex[M][M], Ey[M][M];				//�ڵ����ֵphi,�ڵ��������ܶ�Ni,�ڵ㸺�����ܶ�Ne,�糡E
double Num[N][N], Nume[N][N], D[N][N], Di[N][N], De[N][N];				//��������Num,��������Nume,�������ܶ�D,���������ܶ�Di,���������ܶ�De
double phi1[M][M], E[M][M];
double vx[N * N * 10], vy[N * N * 10];                       			//�������ٶ�v
double vex[N * N * 10], vey[N * N * 10];								//�����ٶ�ve
double Px[N * N * 10], Py[N * N * 10];									//������λ��P
double Pex[N * N * 10], Pey[N * N * 10];								//����λ��Pe
double kB = 1.38e-23, T = 500.0;										//������������kB,�����¶�T
double mi = 2.18e-25, me = 2.18e-25;									//����������mi,����������me
double ep0 = 8.854e-12, e = 1.6e-19, ni = 10e15, ne = 10e15;			//��ս�糣��ep0,Ԫ���e,���������ܶ�ni,���������ܶ�ne
double vth, pi = 3.1415926;												//��������ѧ�ٶ�vth,Բ����pi
double dx = 1e-5, dt = 1e-9;											//������dx,���沽��dt
int K = 1e5, L=1e5;
int k = 0, l = 0;
void initial()															//��ʼ������λ�ú��ٶ�
{
	int i, j, k, n = 0;
	double R1, R2, R3, R4;
	//vth = pow(kB * T / mi, 0.5);
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			for (k = 0; k < 10; k++)
			{
				R1 = (rand() % 1000) / 999.0;
				R2 = (rand() % 1000) / 999.0;
				R3 = (rand() % 1000) / 999.0;
				R4 = (rand() % 1000) / 999.0;
				Px[n] = double(j) * dx + R1 * dx;
				Py[n] = double(i) * dx + R2 * dx;
				Pex[n] = double(j) * dx + R3 * dx;
				Pey[n] = double(i) * dx + R4 * dx;
				vx[n] = 0;
				vy[n] = 0;
				vex[n] = 0;
				vey[n] = 0;
				n++;
			}
		}
	}
}
void weighting()												//����Ȩ��ϵ��������ܶȷ��䵽����
{
	int i, j, x, y;
	double Si, Si1, Sj, Sj1, W[4];
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < M; j++)
		{
			Ni[i][j] = 0;
			Ne[i][j] = 0;
		}
	}
	for (i = 0; i < K; i++)									//������Ȩ��ϵ������͵���ۻ�
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
		Ni[y][x] += W[3] * e / ep0 * (Nrel / pow(dx, 2));
		Ni[y + 1][x] += W[2] * e / ep0 * (Nrel / pow(dx, 2));
		Ni[y][x + 1] += W[1] * e / ep0 * (Nrel / pow(dx, 2));
		Ni[y + 1][x + 1] += W[0] * e / ep0 * (Nrel / pow(dx, 2));
	}
	for (i = 0; i < L; i++)									//������Ȩ��ϵ������͵���ۻ�
	{
		x = int(Pex[i] / dx);
		y = int(Pey[i] / dx);
		Si = (Pey[i] - y * dx) / dx;
		Si1 = ((y + 1.0) * dx - Pey[i]) / dx;
		Sj = (Pex[i] - x * dx) / dx;
		Sj1 = ((x + 1.0) * dx - Pex[i]) / dx;
		W[0] = Si * Sj;
		W[1] = Si1 * Sj;
		W[2] = Si * Sj1;
		W[3] = Si1 * Sj1;
		Ne[y][x] += W[3] * e / ep0 * (Nrel / pow(dx, 2));
		Ne[y + 1][x] += W[2] * e / ep0 * (Nrel / pow(dx, 2));
		Ne[y][x + 1] += W[1] * e / ep0 * (Nrel / pow(dx, 2));
		Ne[y + 1][x + 1] += W[0] * e / ep0 * (Nrel / pow(dx, 2));
	}
}
void E_field()													//����ڵ�糡
{
	int i, j;
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < M; j++)
		{
			if (i == 0)											//���ñ߽�ڵ�糡�������
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
void movement()												//�����˶�
{
	int i = 0, j = 0, x, y;
	double Si, Si1, Sj, Sj1, W[4];
	while (i < K)									//�������˶�
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
		vx[i] += e / mi * (Ex[y][x] * W[3] + Ex[y + 1][x] * W[2] + Ex[y][x + 1] * W[1] + Ex[y + 1][x + 1] * W[0]) * dt;
		vy[i] += e / mi * (Ey[y][x] * W[3] + Ey[y + 1][x] * W[2] + Ey[y][x + 1] * W[1] + Ey[y + 1][x + 1] * W[0]) * dt;
		Px[i] += vx[i] * dt;
		Py[i] += vy[i] * dt;
		if (Px[i] <= 0 || Px[i] >= 100 * dx || Py[i] <= 0 || Py[i] >= 100 * dx)				//�߽�����������
		{
			K--;
			vx[i] = vx[K];
			vy[i] = vy[K];
			Px[i] = Px[K];
			Py[i] = Py[K];
		}
		i++;
	}
	while (j < L)									//�������˶�
	{
		x = int(Pex[j] / dx);
		y = int(Pey[j] / dx);
		Si = (Pey[j] - y * dx) / dx;
		Si1 = ((y + 1.0) * dx - Pey[j]) / dx;
		Sj = (Pex[j] - x * dx) / dx;
		Sj1 = ((x + 1.0) * dx - Pex[j]) / dx;
		W[0] = Si * Sj;
		W[1] = Si1 * Sj;
		W[2] = Si * Sj1;
		W[3] = Si1 * Sj1;
		vex[j] += -e / mi * (Ex[y][x] * W[3] + Ex[y + 1][x] * W[2] + Ex[y][x + 1] * W[1] + Ex[y + 1][x + 1] * W[0]) * dt;
		vey[j] += -e / mi * (Ey[y][x] * W[3] + Ey[y + 1][x] * W[2] + Ey[y][x + 1] * W[1] + Ey[y + 1][x + 1] * W[0]) * dt;
		Pex[j] += vex[j] * dt;
		Pey[j] += vey[j] * dt;
		if (Pex[j] <= 0 || Pex[j] >= 100 * dx || Pey[j] <= 0 || Pey[j] >= 100 * dx)				//�߽�����������
		{
			L--;
			vex[j] = vex[L];
			vey[j] = vey[L];
			Pex[j] = Pex[L];
			Pey[j] = Pey[L];
		}
		j++;
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
			Nume[i][j] = 0;
		}
	}
	for (i = 0; i < K; i++)
	{
		x = int(Px[i] / dx);
		y = int(Py[i] / dx);
		Num[y][x]++;
	}
	for (i = 0; i < L; i++)
	{
		x = int(Pex[i] / dx);
		y = int(Pey[i] / dx);
		Nume[y][x]++;
	}
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			Di[i][j]= Num[i][j] * Nrel / (dx * dx);
			De[i][j] = Nume[i][j] * Nrel / (dx * dx);
			D[i][j] = (Num[i][j] + Nume[i][j]) * Nrel / (dx * dx);
		}
	}
}
int main()
{
	int i, j;
	double w;
	using namespace std;
	ofstream fp("position.dat", ofstream::out);
	ofstream fp1("idensity.dat", ofstream::out);
	ofstream fp2("potential.dat", ofstream::out);
	ofstream fp3("electric field.dat", ofstream::out);
	ofstream fp4("speed.dat", ofstream::out);
	ofstream fp5("espeed.dat", ofstream::out);
	ofstream fp6("eposition.dat", ofstream::out);
	ofstream fp7("edensity.dat", ofstream::out);
	initial();
	w = 2 / (1 + sqrt(1 - pow((cos(pi / M) + cos(pi / M)) / 2, 2)));
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
	for(k = 0; k < 20000; k++)
	{
		weighting();
		double t, tmax = 1;
		while (tmax > 1e-3)
		{
			l++;
			tmax = 0;
			for (i = 1; i < M - 1; i++)
			{
				for (j = 1; j < M - 1; j++)
				{
					if (i == 50 && j == 50)
						continue;
					phi[i][j] = phi1[i][j] + (w / 4) * (phi[i - 1][j] + phi1[i + 1][j] + phi[i][j - 1] + phi1[i][j + 1] - 4 * phi1[i][j] + pow(dx, 2) / 4 * (Ni[i][j] - Ne[i][j]));
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
		cout << k << endl;
	}
	density();
	cout << "K=" << K << " " << "L=" << L << endl;
	fp << "title=\"Position\"" << endl << "variables=\"x\",\"y\"" << " " << "zone i=" << K << endl;
	for (i = 0; i < K; i++)
	{
		fp << Px[i] << " " << Py[i] << endl;
	}
	fp1 << "title=\"iDensity\"" << endl << "variables = \"x\",\"y\",\"num\"" << endl << "zone i=100,j=100,f=point" << endl;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			fp1 << j + 1 << " " << i + 1 << " " << Di[i][j] << endl;
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
	fp4 << "title=\"Speed\"" << endl << "variables = \"x\",\"y\",\"vx\",\"vy\"" << endl << "zone i=" << K << endl;
	for (i = 0; i < K; i++)
	{
		fp4 << Px[i] << " " << Py[i] << " " << vx[i] << " " << vy[i] << endl;
	}
	fp5 << "title=\"eSpeed\"" << endl << "variables = \"x\",\"y\",\"vex\",\"vey\"" << endl << "zone i=" << L << endl;
	for (i = 0; i < L; i++)
	{
		fp5 << Pex[i] << " " << Pey[i] << " " << vex[i] << " " << vey[i] << endl;
	}
	fp6 << "title=\"ePosition\"" << endl << "variables=\"x\",\"y\"" << " " << "zone i=" << L << endl;
	for (i = 0; i < L; i++)
	{
		fp6 << Pex[i] << " " << Pey[i] << endl;
	}
	fp7 << "title=\"eDensity\"" << endl << "variables = \"x\",\"y\",\"num\"" << endl << "zone i=100,j=100,f=point" << endl;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			fp7 << j + 1 << " " << i + 1 << " " << De[i][j] << endl;
		}
	}
	fp.close();
	fp1.close();
	fp2.close();
	fp3.close();
	fp4.close();
	fp5.close();
	fp6.close();
	fp7.close();
}