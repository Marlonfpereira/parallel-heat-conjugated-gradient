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
// void copia_vetor(float* dest, float* orig, int tdom)
// {
// 	int i;

// 	for (i = 0;i < tdom;i++)
// 	{
// 		dest[i] = orig[i];
// 	}
// }

// void mult_mat_vet(float** matriz, int* offset, float* vetor, float* result, int tdom)
// {
// 	float tmp;
// 	int i, j;

// 	for (i = 0;i < tdom;i++)
// 	{
// 		tmp = 0.00;
// 		for (j = 0;j < 5;j++)
// 		{
// 			if (((i + offset[j]) >= 0) && (matriz[i][j] != 0))
// 				tmp += matriz[i][j] * vetor[i + offset[j]];
// 		}
// 		result[i] = tmp;
// 	}
// }

// float produto_escalar(float* vetor1, float* vetor2, int tdom)
// {
// 	int i;
// 	float resposta = 0;

// 	for (i = 0;i < tdom;i++)
// 		resposta += vetor1[i] * vetor2[i];

// 	return (resposta);
// }


// void escalar_vetor(float* vetor, float escalar, float* resposta, int tdom)
// {
// 	int i;

// 	for (i = 0;i < tdom;i++)
// 		resposta[i] = escalar * vetor[i];
// }

// void soma_vetor(float* vetor1, float* vetor2, float* resposta, int tdom)
// {
// 	int i;

// 	for (i = 0;i < tdom;i++)
// 		resposta[i] = vetor1[i] + vetor2[i];
// }

// void sub_vetor(float* vetor1, float* vetor2, float* resposta, int tdom)
// {
// 	int i;

// 	for (i = 0;i < tdom;i++)
// 		resposta[i] = vetor1[i] - vetor2[i];
// }


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
	start = clock();

	for (j = 0;j < ciclos;j++)
	{
		#pragma acc parallel copyin(matA[0:size][0:5], offset[0:5], vetb[0:size], vetb_ext[0:size]) copyout(vetx[0:size])
		{
			for (int i = 0;i < size;i++)
				vetb[i] = vetb[i] + vetb_ext[i];

			///////////////////////////////////////////////
			float *r = (float*)calloc(size, sizeof(float));
			//# mult_mat_vet(matriz, offset, vetx, r, size);
			for (int i = 0;i < size;i++)
			{
				float tmp = 0.00;
				for (int j = 0;j < 5;j++)
				{
					if (((i + offset[j]) >= 0) && (matA[i][j] != 0))
						tmp += matA[i][j] * vetx[i + offset[j]];
				}
				r[i] = tmp;
			}
			//#

			//# sub_vetor(vetb, r, r, size);
			for (int i = 0;i < size;i++)
				r[i] = vetb[i] - r[i];
			//#

			float* d = (float*)calloc(size, sizeof(float));
			//# copia_vetor(d, r, size);
			for (i = 0;i < size;i++)
				d[i] = r[i];
			//#

			//# sigman = produto_escalar(r, r, size);
			float sigman = 0;
			for (int i = 0;i < size;i++)
				sigman += r[i] * r[i];
			//#

			int it = 0;
			do
			{
				float* q = (float*)calloc(size, sizeof(float));
				//# mult_mat_vet(matA, offset, d, q, size);
				for (int i = 0;i < size;i++)
				{
					float tmp = 0.00;
					for (int j = 0;j < 5;j++)
					{
						if (((i + offset[j]) >= 0) && (matA[i][j] != 0))
							tmp += matA[i][j] * d[i + offset[j]];
					}
					q[i] = tmp;
				}
				//#

				//# den = produto_escalar(d, q, size);
				float den = 0;
				for (int i = 0;i < size;i++)
					den += d[i] * q[i];
				//# 

				float alfa;
				alfa = sigman / den;

				float* aux = (float*)calloc(size, sizeof(float));
				//# escalar_vetor(d, alfa, aux, size);
				for (int i = 0;i < size;i++)
					aux[i] = alfa * d[i];
				//#

				//# soma_vetor(vetx, aux, vetx, size);
				for (int i = 0;i < size;i++)
					vetx[i] = vetx[i] + aux[i];
				//# 

				//# escalar_vetor(q, alfa, aux, size);
				for (int i = 0;i < size;i++)
					aux[i] = alfa * q[i];
				//#

				//# sub_vetor(r, aux, r, size);
				for (int i = 0;i < size;i++)
					r[i] = r[i] - aux[i];
				//#

				float sigmav = sigman;

				//# sigman = produto_escalar(r, r, size);
				sigman = 0;
				for (int i = 0;i < size;i++)
					sigman += r[i] * r[i];
				//# 

				float beta;
				beta = sigman / sigmav;

				//# escalar_vetor(d, beta, aux, size);
				for (int i = 0;i < size;i++)
					aux[i] = beta * d[i];
				//# 

				//# soma_vetor(r, aux, d, size);
				for (int i = 0;i < size;i++)
					d[i] = r[i] + aux[i];
				//#

				it++;
			} while (it < itmax);
			///////////////////////////////////////////////

			//# copia_vetor(vetb, vetx, size);
			for (int i = 0;i < size;i++)
				vetb[i] = vetx[i];
			//# 
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


