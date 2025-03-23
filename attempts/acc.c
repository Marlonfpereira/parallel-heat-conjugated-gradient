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

3
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
	float valor = dominio[0][0];
	unsigned int alloc = size * sizeof(float);

	//solução iterativa - ciclos
	start = omp_get_wtime();

#pragma acc enter data copyin(size,matA[0:size][0:5],vetb[0:size],vetb_ext[0:size],vetx[0:size])

	for (j = 0;j < ciclos;j++)
	{
		printf("iteracao %d\n", j);
#pragma acc update host(vetx[0:size])
		if (vetx[size / 2] == valor)
			break;
		float* r = (float*)malloc(alloc);
		float* d = (float*)malloc(alloc);
		float* q = (float*)malloc(alloc);
		float* aux = (float*)malloc(alloc);
		float sigman = 0, den = 0, alfa = 0, sigmav = 0, beta = 0;
		// #pragma acc exit data if(acc_deviceptr(r)) delete(r[0:size],d[0:size],q[0:size],aux[0:size],sigman,den,alfa,sigmav,beta)
#pragma acc enter data copyin(r[0:size],d[0:size],q[0:size],aux[0:size],sigman,den,alfa,sigmav,beta)

		// #pragma omp for
#pragma acc parallel num_gangs(1)  loop
		for (int i = 0; i < size; i++)
			vetb[i] = vetb[i] + vetb_ext[i];

		// #pragma omp for
#pragma acc parallel num_gangs(1)  loop
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

		// #pragma omp for
#pragma acc parallel num_gangs(1)  loop
		for (int i = 0; i < size; i++)
			r[i] = vetb[i] - r[i];

		// #pragma omp for
#pragma acc parallel num_gangs(1)  loop
		for (int i = 0; i < size; i++)
			d[i] = r[i];


		// #pragma omp single
		sigman = 0;
		// #pragma omp for reduction(+:sigman)
#pragma acc parallel num_gangs(1)  loop reduction(+:sigman)
		for (int i = 0; i < size; i++)
			sigman += r[i] * r[i];

#pragma acc update host(sigman)
#pragma acc update device(sigman)

		int it = 0;
		do
		{
			// #pragma omp for
#pragma acc parallel num_gangs(1)  loop
			for (int i = 0; i < size; i++)
			{
				// printf("d inical = %f\n", d[i]);

				float tmp = 0.00;
				// Unrolled loop
				if (((i + offset[0]) >= 0) && (matA[i][0] != 0))
					tmp += matA[i][0] * d[i + offset[0]];
				// printf("%d, %f, %f", i, )
				if (((i + offset[1]) >= 0) && (matA[i][1] != 0))
					tmp += matA[i][1] * d[i + offset[1]];
				if (((i + offset[2]) >= 0) && (matA[i][2] != 0))
					tmp += matA[i][2] * d[i + offset[2]];
				if (((i + offset[3]) >= 0) && (matA[i][3] != 0))
					tmp += matA[i][3] * d[i + offset[3]];
				if (((i + offset[4]) >= 0) && (matA[i][4] != 0))
					tmp += matA[i][4] * d[i + offset[4]];
				// printf("%f, %d, %d\n",d[i + offset[j]], i, offset[j]);
				q[i] = tmp;
				// printf("q = %f\n", q[i]);
			}

			// #pragma omp single
			den = 0;
#pragma acc update device(den)

			// #pragma omp for reduction(+:den)
#pragma acc parallel num_gangs(1)  loop reduction(+:den)
			for (int i = 0; i < size; i++)
			{
				den += d[i] * q[i];
			}

#pragma acc update host(den)
#pragma acc update device(den)
			// #pragma omp single
#pragma acc parallel num_gangs(1) 
			{

				alfa = sigman / den;
				sigmav = sigman;
#pragma acc exit data copyout(alfa,sigmav)
			}


			// alfs, den, sigmav: 0.806452 24.799999 20.000000
			// alfs, den, sigmav: 0.761404 0.382666 0.291363
			// alfs, den, sigmav : 0.618758 0.004457 0.002758


			// #pragma omp for
#pragma acc parallel num_gangs(1)  loop
			for (int i = 0; i < size; i++)
			{

				aux[i] = alfa * d[i];

			}
#pragma acc update host(aux)
#pragma acc update device(aux)

			// #pragma omp for
#pragma acc parallel num_gangs(1)  loop
			for (int i = 0; i < size; i++)
			{
				vetx[i] = vetx[i] + aux[i];
			}
#pragma acc update host(vetx)
#pragma acc update device(vetx)

			// #pragma omp for
#pragma acc parallel num_gangs(1)  loop
			for (int i = 0; i < size; i++)
			{
				aux[i] = alfa * q[i];
			}
#pragma acc update host(aux)
#pragma acc update device(aux)

			// #pragma omp for
#pragma acc parallel num_gangs(1)  loop
			for (int i = 0; i < size; i++)
			{
				r[i] = r[i] - aux[i];
			}
#pragma acc update host(r)
#pragma acc update device(r)

			// float sigmav = sigman;

	// #pragma omp critical
			sigman = 0;
#pragma acc update device(sigman)
			// #pragma omp for reduction(+:sigman)
#pragma acc parallel num_gangs(1)  loop reduction(+:sigman)
			for (int i = 0; i < size; i++)
				sigman += r[i] * r[i];

#pragma acc update host(sigman)

			// #pragma omp single
#pragma acc parallel num_gangs(1) 
			{
				beta = sigman / sigmav;
#pragma acc exit data copyout(beta)
			}

			// #pragma acc update host(beta)
			// #pragma acc update device(beta)

						// #pragma omp for
#pragma acc parallel num_gangs(1)  loop
			for (int i = 0; i < size; i++)
			{
				aux[i] = beta * d[i];
			}

			// #pragma acc update host(aux)
			// #pragma acc update device(aux)

						// #pragma omp for
#pragma acc parallel num_gangs(1)  loop
			for (int i = 0; i < size; i++)
			{
				d[i] = r[i] + aux[i];
			}

#pragma acc update host(d)
#pragma acc update device(d)

			// #pragma omp atomic
			it++;
		} while (it < itmax);

		// #pragma omp for
#pragma acc parallel num_gangs(1)  loop
		for (int i = 0; i < size; i++)
		{
			vetb[i] = vetx[i];
		}

	}
#pragma acc exit data delete(size,matA[0:size][0:5],vetb[0:size],vetb_ext[0:size]) copyout(vetx[0:size])

	end = omp_get_wtime();

	//saida
	for (i = 0; i < size; i++)
	{
		printf("%4.1f\n", vetx[i]);
	}
	printf("\n ");

	printf("%lf \n ", (end - start));

	return(0);
}


