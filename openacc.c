#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <openacc.h>

#define deltaX 0.1
#define deltaT 0.1
#define alpha 0.01
#define itmax 3

/*

1
5 5
10.0 10.0 10.0 10.0 10.0
10.0  0.0  0.0  0.0 10.0
10.0  0.0  0.0  0.0 10.0
10.0  0.0  0.0  0.0 10.0
10.0 10.0 10.0 10.0 10.0

*/


//--------------------------------------------------------
int geramatriz(float** dominio, float** matA, float* vetb, float* vetb_ext, int m, int n)
{
	int k, i, j;
	float C, X;

	k = 0;
	C = (float)1 + (4 * alpha * deltaT) / (deltaX * deltaX);
	X = (float)-(alpha * deltaT) / (deltaX * deltaX);

	for (i = 1; i < m - 1; i++)
	{
		for (j = 1; j < n - 1; j++)
		{
			k = (i - 1) * (m - 2) + (j - 1);

			//valor da 1a diagonal
			if (i == 1)
			{
				vetb_ext[k] += -X * dominio[i - 1][j];
			}
			else
				matA[k][0] = X;

			//valor da 2a diagonal
			if (j == 1)
			{
				vetb_ext[k] += -X * dominio[i][j - 1];
			}
			else
			{
				matA[k][1] = X;
			}

			//valor da diagonal central
			matA[k][2] = C;
			vetb[k] = dominio[i][j];

			//valor da 4a diagonal
			if (j == n - 2)
			{
				vetb_ext[k] += -X * dominio[i][j + 1];
			}
			else
				matA[k][3] = X;

			//valor da 5a diagonal
			if (i == n - 2)
			{
				vetb_ext[k] += -X * dominio[i + 1][j];
			}
			else
				matA[k][4] = X;
		}
	}

	return (m - 2) * (n - 2);
}

//---------------------------------------------------------

int main()
{
	float** dominio;
	float** matA; int* offset;
	float* vetb;
	float* vetb_ext;
	float* vetx;
	clock_t start, end;

	unsigned int m, n, size;
	unsigned int i, j;
	unsigned int ciclos;

	// leitura e alocação dimensão do domínio
	scanf("%d", &ciclos);
	scanf("%d", &m);
	scanf("%d", &n);

	offset = (int*)calloc(5, sizeof(int));
	offset[0] = -(n - 2);
	offset[1] = -1;
	offset[2] = 0;
	offset[3] = 1;
	offset[4] = n - 2;

	dominio = (float**)calloc(m, sizeof(float*));
	for (i = 0; i < m; i++)
		dominio[i] = (float*)calloc(n, sizeof(float));

	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			scanf("%f", &dominio[i][j]);

	//alocação estrutura de dados da pentadiagonal
	//(n-2) e (m-2) por causa da borda
	matA = (float**)calloc(((n - 2) * (m - 2)), sizeof(float*));
	for (i = 0; i < ((n - 2) * (m - 2)); i++)
		matA[i] = (float*)calloc(5, sizeof(float));

	//alocação do vetor dos termos independentes b e b_ext e de resposta x
	vetb = (float*)calloc(((n - 2) * (m - 2)), sizeof(float));
	vetb_ext = (float*)calloc(((n - 2) * (m - 2)), sizeof(float));
	vetx = (float*)calloc(((n - 2) * (m - 2)), sizeof(float));

	//geração dos sistemas
	size = geramatriz(dominio, matA, vetb, vetb_ext, m, n);
	// float valor = dominio[0][0];
	unsigned int alloc = size * sizeof(float);

	//solução iterativa - ciclos
	start = clock();

	float* r = (float*)malloc(alloc);
	float* d = (float*)malloc(alloc);
	float* q = (float*)malloc(alloc);
	float* aux = (float*)malloc(alloc);
#pragma acc data copyin(size,matA[0:size][0:5],offset[0:5],vetb[0:size],vetb_ext[0:size],vetx[0:size],r[0:size],d[0:size],q[0:size],aux[0:size]) copyout(vetx[0:size])
	{
		for (j = 0;j < ciclos;j++)
		{
			float sigman = 0, den = 0, alfa = 0, sigmav = 0, beta = 0;
#pragma acc enter data create(sigman,den,alfa,sigmav,beta)

#pragma acc parallel default(none)  loop
			for (int i = 0; i < size; i++)
				vetb[i] = vetb[i] + vetb_ext[i];

#pragma acc parallel default(none)  loop
			for (int i = 0; i < size; i++)
			{
				float tmp = 0.00;
				if (((i + offset[0]) >= 0) && (matA[i][0] != 0))
					tmp += matA[i][0] * vetx[i + offset[0]];
				if (((i + offset[1]) >= 0) && (matA[i][1] != 0))
					tmp += matA[i][1] * vetx[i + offset[1]];
				if (((i + offset[2]) >= 0) && (matA[i][2] != 0))
					tmp += matA[i][2] * vetx[i + offset[2]];
				if (((i + offset[3]) >= 0) && (matA[i][3] != 0))
					tmp += matA[i][3] * vetx[i + offset[3]];
				if (((i + offset[4]) >= 0) && (matA[i][4] != 0))
					tmp += matA[i][4] * vetx[i + offset[4]];

				r[i] = tmp;
			}

#pragma acc parallel default(none)  loop
			for (int i = 0; i < size; i++)
			{
				r[i] = vetb[i] - r[i];
			}

#pragma acc parallel default(none)  loop
			for (int i = 0; i < size; i++)
			{
				d[i] = r[i];
			}




			//ZONA CRITTICA 
#pragma acc parallel default(none)  loop reduction(+:sigman)
			for (int i = 0; i < size; i++)
			{
				sigman += r[i] * r[i];
			}

#pragma acc update host(sigman)

			int it = 0;
			do
			{
#pragma acc parallel default(none)  loop
				for (int i = 0; i < size; i++)
				{
					float tmp = 0.00;
					if (((i + offset[0]) >= 0) && (matA[i][0] != 0))
						tmp += matA[i][0] * d[i + offset[0]];
					if (((i + offset[1]) >= 0) && (matA[i][1] != 0))
						tmp += matA[i][1] * d[i + offset[1]];
					if (((i + offset[2]) >= 0) && (matA[i][2] != 0))
						tmp += matA[i][2] * d[i + offset[2]];
					if (((i + offset[3]) >= 0) && (matA[i][3] != 0))
						tmp += matA[i][3] * d[i + offset[3]];
					if (((i + offset[4]) >= 0) && (matA[i][4] != 0))
						tmp += matA[i][4] * d[i + offset[4]];

					q[i] = tmp;
				}




				//ZONA CRITTICA 
				den = 0;
				// #pragma acc enter data create(den)
#pragma acc update device(den)
#pragma acc parallel default(none)  loop reduction(+:den)
				for (int i = 0; i < size; i++)
				{
					den += d[i] * q[i];
				}
#pragma acc update host(den)




				//ZONA CRITTICA 
#pragma acc parallel default(none) 
				{
					alfa = sigman / den;
					sigmav = sigman;
#pragma acc exit data copyout(alfa,sigmav)
				}




#pragma acc parallel default(none)  loop
				for (int i = 0; i < size; i++)
				{
					aux[i] = alfa * d[i];
				}



#pragma acc parallel default(none)  loop
				for (int i = 0; i < size; i++)
				{
					vetx[i] = vetx[i] + aux[i];
				}



#pragma acc parallel default(none)  loop
				for (int i = 0; i < size; i++)
				{
					aux[i] = alfa * q[i];
				}




#pragma acc parallel default(none)  loop
				for (int i = 0; i < size; i++)
				{
					r[i] = r[i] - aux[i];
				}



				//ZONA CRITTICA 
				sigman = 0;
				// #pragma acc enter data create(sigman)
#pragma acc update device(sigman)
#pragma acc parallel default(none)  loop reduction(+:sigman)
				for (int i = 0; i < size; i++)
					sigman += r[i] * r[i];
#pragma acc update host(sigman)




				//ZONA CRITTICA 
#pragma acc parallel default(none) 
				{
					beta = sigman / sigmav;
#pragma acc exit data copyout(beta)
				}


#pragma acc parallel default(none)  loop
				for (int i = 0; i < size; i++)
				{
					aux[i] = beta * d[i];
				}


#pragma acc parallel default(none)  loop
				for (int i = 0; i < size; i++)
				{
					d[i] = r[i] + aux[i];
				}


				it++;
			} while (it < itmax);

#pragma acc parallel default(none)  loop
			for (int i = 0; i < size; i++)
			{
				vetb[i] = vetx[i];
			}

		}
	}

	end = clock();

	//saida
	for (i = 0; i < size; i++)
	{
		printf("%4.1f\n", vetx[i]);
	}
	printf("\n ");

	printf("%lf \n ", (double)(end - start) / CLOCKS_PER_SEC);

	return(0);
}


