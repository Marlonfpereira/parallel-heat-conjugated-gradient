#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// #include <omp.h>
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

//---------------------------------------------------------
void copia_vetor(float* dest, float* orig, int tdom)
{
	int i;

	for (i = 0;i < tdom;i++)
	{
		dest[i] = orig[i];
	}
}

void mult_mat_vet(float** matriz, int* offset, float* vetor, float* result, int tdom)
{
	float tmp;
	int i, j;

	// #pragma acc parallel loop copyin(matriz[0:tdom][0:5], offset[0:5], vetor[0:tdom]) copyout(result[0:tdom]) private(tmp,j)
	for (i = 0;i < tdom;i++)
	{
		tmp = 0.00;
		for (j = 0;j < 5;j++)
		{
			if (((i + offset[j]) >= 0) && (matriz[i][j] != 0))
				tmp += matriz[i][j] * vetor[i + offset[j]];
		}
		result[i] = tmp;
	}
}

float produto_escalar(float* vetor1, float* vetor2, int tdom)
{
	int i;
	float resposta = 0;

	for (i = 0;i < tdom;i++)
		resposta += vetor1[i] * vetor2[i];

	return (resposta);
}


void escalar_vetor(float* vetor, float escalar, float* resposta, int tdom)
{
	int i;

	for (i = 0;i < tdom;i++)
		resposta[i] = escalar * vetor[i];
}

void soma_vetor(float* vetor1, float* vetor2, float* resposta, int tdom)
{
	int i;

	for (i = 0;i < tdom;i++)
		resposta[i] = vetor1[i] + vetor2[i];
}

void sub_vetor(float* vetor1, float* vetor2, float* resposta, int tdom)
{
	int i;

	for (i = 0;i < tdom;i++)
		resposta[i] = vetor1[i] - vetor2[i];
}


//////////////////////////////////////////////////////
//              Gradiente Conjugado                 //
//////////////////////////////////////////////////////
void GC(float** matriz,
	int* offset,
	float* vet1,
	float* resp,
	int tdom)
{
	int it;
	float* r, * d, * q, * aux, * aux2;
	float alfa, beta, sigman, sigmav, den;

	r = (float*)calloc(tdom, sizeof(float));
	d = (float*)calloc(tdom, sizeof(float));
	q = (float*)calloc(tdom, sizeof(float));
	aux = (float*)calloc(tdom, sizeof(float));
	aux2 = (float*)calloc(tdom, sizeof(float));

	it = 0;
	// mult_mat_vet(matriz, offset, resp, r, tdom);
	// 
	// float tmp;
	int i, j;

	for (i = 0;i < tdom;i++)
	{
		float tmp = 0.00;
		for (j = 0;j < 5;j++)
		{
			if (((i + offset[j]) >= 0) && (matriz[i][j] != 0))
				tmp += matriz[i][j] * resp[i + offset[j]];
		}
		r[i] = tmp;
	}
	// 

	sub_vetor(vet1, r, r, tdom);
	copia_vetor(d, r, tdom);
	sigman = produto_escalar(r, r, tdom);
	copia_vetor(aux2, resp, tdom);

#pragma acc data copyin(matriz[0:tdom][0:5], offset[0:5], tdom)
	do
	{
#pragma acc parallel
		{
			//# mult_mat_vet(matriz, offset, d, q, tdom);
			// float tmp;
			// int i, j;
#pragma acc data copyin(d[0:tdom]) create(q[0:tdom])
#pragma acc loop private(j)
			for (int i = 0;i < tdom;i++)
			{

				float tmp = 0.00;
				for (int j = 0;j < 5;j++)
				{
					if (((i + offset[j]) >= 0) && (matriz[i][j] != 0))
						tmp += matriz[i][j] * d[i + offset[j]];
				}
				q[i] = tmp;
			}
			//#

			//# den = produto_escalar(d, q, tdom);
			// int i;
			den = 0;
#pragma acc data create(den)
#pragma acc loop
			for (int i = 0;i < tdom;i++)
				den += d[i] * q[i];
			// #pragma data copyout(den)
						//# 

#pragma acc data create(alfa) copyin(sigman)
			alfa = sigman / den;

			//# escalar_vetor(d, alfa, aux, tdom);
			// int i;
#pragma acc data create(aux[0:tdom])
#pragma acc loop 
			for (int i = 0;i < tdom;i++)
				aux[i] = alfa * d[i];
			// #pragma acc data copyout(aux[0:tdom])
					//#

					//# soma_vetor(resp, aux, resp, tdom);
					// int i;
#pragma acc data copyin(resp[0:tdom])
#pragma acc loop
			for (int i = 0;i < tdom;i++)
				resp[i] = resp[i] + aux[i];
#pragma acc data copyout(resp[0:tdom])
			//# 

			//# escalar_vetor(q, alfa, aux, tdom);
			// int i;
#pragma acc loop
			for (int i = 0;i < tdom;i++)
				aux[i] = alfa * q[i];
			// #pragma acc data copyout(aux[0:tdom])
					//#

					//# sub_vetor(r, aux, r, tdom);
					// int i;
#pragma acc data copyin(r[0:tdom])
#pragma acc loop
			for (int i = 0;i < tdom;i++)
				r[i] = r[i] - aux[i];
#pragma acc data copyout(r[0:tdom])
			//#

#pragma acc data create(sigmav)
			sigmav = sigman;

			//# sigman = produto_escalar(r, r, tdom);
			// int i;
			sigman = 0;
#pragma acc loop
			for (int i = 0;i < tdom;i++)
				sigman += r[i] * r[i];
#pragma acc data copyout(sigman)
			//# 

#pragma acc data create(beta)
			beta = sigman / sigmav;

			//# escalar_vetor(d, beta, aux, tdom);
			// int i;
#pragma acc loop
			for (int i = 0;i < tdom;i++)
				aux[i] = beta * d[i];
			// #pragma acc data copyout(aux[0:tdom])
					//# 

					//# soma_vetor(r, aux, d, tdom);
					// int i;
#pragma acc loop
			for (int i = 0;i < tdom;i++)
				d[i] = r[i] + aux[i];
#pragma acc data copyout(d[0:tdom])
			//#

		}
		it++;
	} while (it < itmax);
}

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

	//solução iterativa - ciclos
	// timea = omp_get_wtime();
	start = clock();

	for (j = 0;j < ciclos;j++)
	{
		soma_vetor(vetb, vetb_ext, vetb, size);
		GC(matA, offset, vetb, vetx, size);
		copia_vetor(vetb, vetx, size);
	}
	// timeb = omp_get_wtime();
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


