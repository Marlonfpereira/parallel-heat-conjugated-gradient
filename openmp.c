#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

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
	float start, end;

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
	unsigned int alloc = size * sizeof(float);

	//solução iterativa - ciclos
	start = omp_get_wtime();
	float* r = (float*)malloc(alloc);
	float* d = (float*)malloc(alloc);
	float* q = (float*)malloc(alloc);
	float* aux = (float*)malloc(alloc);

	for (j = 0;j < ciclos;j++)
	{
		float sigman = 0, den = 0, alfa = 0, sigmav = 0, beta = 0;

#pragma omp parallel num_threads(8) 
		//small 16= 24.332031  10= 4.610840  8= 2.894531  4= 3.468750  2= 5.449219 
		//big   16= 80.180054  10= 36.992065  8= 32.939453 4= 35.902954 2= 46.782227
		{
#pragma omp for
			for (int i = 0; i < size; i++)
				vetb[i] = vetb[i] + vetb_ext[i];

#pragma omp for
			for (int i = 0; i < size; i++)
			{
				float tmp = 0.00;
				for (int j = 0; j < 5; j++)
				{
					if (((i + offset[j]) >= 0) && (matA[i][j] != 0))
						tmp += matA[i][j] * vetx[i + offset[j]];
				}
				r[i] = tmp;
			}

#pragma omp for
			for (int i = 0; i < size; i++)
				r[i] = vetb[i] - r[i];

#pragma omp for
			for (int i = 0; i < size; i++)
				d[i] = r[i];

#pragma omp single
			sigman = 0;
#pragma omp for reduction(+:sigman)
			for (int i = 0; i < size; i++)
				sigman += r[i] * r[i];

			int it = 0;
			do
			{
#pragma omp for
				for (int i = 0; i < size; i++)
				{
					float tmp = 0.00;
					for (int j = 0; j < 5; j++)
					{
						if (((i + offset[j]) >= 0) && (matA[i][j] != 0))
							tmp += matA[i][j] * d[i + offset[j]];
					}
					q[i] = tmp;
				}

#pragma omp single
				den = 0;
#pragma omp for reduction(+:den)
				for (int i = 0; i < size; i++)
					den += d[i] * q[i];

#pragma omp single
				{
					alfa = sigman / den;
					sigmav = sigman;
					sigman = 0;
				}

#pragma omp for
				for (int i = 0; i < size; i++)
					aux[i] = alfa * d[i];

#pragma omp for
				for (int i = 0; i < size; i++)
					vetx[i] = vetx[i] + aux[i];

#pragma omp for
				for (int i = 0; i < size; i++)
					aux[i] = alfa * q[i];

#pragma omp for
				for (int i = 0; i < size; i++)
					r[i] = r[i] - aux[i];

				// float sigmav = sigman;

// #pragma omp critical
// 				sigman = 0;
#pragma omp for reduction(+:sigman)
				for (int i = 0; i < size; i++)
					sigman += r[i] * r[i];

#pragma omp single
				beta = sigman / sigmav;

#pragma omp for
				for (int i = 0; i < size; i++)
					aux[i] = beta * d[i];

#pragma omp for
				for (int i = 0; i < size; i++)
					d[i] = r[i] + aux[i];

#pragma omp atomic
				it++;
			} while (it < itmax);

#pragma omp for
			for (int i = 0; i < size; i++)
				vetb[i] = vetx[i];

		} // end parallel
	}
	free(r);
	free(d);
	free(q);
	free(aux);
	end = omp_get_wtime();

	//saida
	// for (i = 0; i < size; i++)
	// {
	// 	printf("%4.1f\n", vetx[i]);
	// }
	// printf("\n ");

	printf("%lf \n ", (end - start));

	return(0);
}


