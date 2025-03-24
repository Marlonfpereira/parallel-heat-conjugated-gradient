#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

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

float timea, timeb;

//---------------------------------------------------------
void copia_vetor(float *dest, float *orig, int tdom)
{
    int i;

    for(i=0;i<tdom;i++)
    {
        dest[i]=orig[i];
    }
}

void mult_mat_vet(float **matriz, int *offset, float *vetor, float *result, int tdom)
{
	float tmp;
	int i,j;

	for(i=0;i<tdom;i++)
	{
	      tmp=0.00;
	      for(j=0;j<5;j++)
	      {
		  if(((i+offset[j])>=0)   &&    (matriz[i][j]!=0))
		     tmp+=matriz[i][j]*vetor[i+offset[j]];
	      }
	      result[i]=tmp;
	}
}

float produto_escalar (float *vetor1, float *vetor2, int tdom)
{
	int i;
	float resposta=0;

	for (i=0;i<tdom;i++)
		resposta+= vetor1[i] * vetor2[i];

	return (resposta);
}


void escalar_vetor (float *vetor, float escalar, float *resposta, int tdom)
{
	int i;

	for (i=0;i<tdom;i++)
		resposta[i]=escalar * vetor[i];
}

void soma_vetor (float *vetor1, float *vetor2, float *resposta, int tdom)
{
	int i;

	for (i=0;i<tdom;i++)
		resposta[i]=vetor1[i]+vetor2[i];
}

void sub_vetor (float *vetor1, float *vetor2, float *resposta, int tdom)
{
	int i;

	for (i=0;i<tdom;i++)
		resposta[i]=vetor1[i]-vetor2[i];
}


//////////////////////////////////////////////////////
//              Gradiente Conjugado                 //
//////////////////////////////////////////////////////
void GC(float **matriz,
	int *offset,
        float *vet1,
        float *resp,
        int tdom)
{
	// GC(matA, offset, vetb, vetx,  size);
	int it;
	float *r,*d,*q,*aux,*aux2;
	float alfa, beta, sigman, sigmav,den;

	r=(float *)calloc(tdom,sizeof(float));
	d=(float *)calloc(tdom,sizeof(float));
	q=(float *)calloc(tdom,sizeof(float));
	aux=(float *)calloc(tdom,sizeof(float));
	aux2=(float *)calloc(tdom,sizeof(float));

	it=0;
	mult_mat_vet(matriz, offset, resp, r, tdom);
	// usa: matriz, offset, resp
	// altera: r
	sub_vetor(vet1, r, r, tdom);
	// usa: vet1, r
	// altera: r
	copia_vetor(d, r, tdom);
	// usa: r
	// altera: d
	sigman = produto_escalar(r, r, tdom);
	// usa: r
	// altera: sigman
	copia_vetor(aux2, resp, tdom);

    do
	{
		mult_mat_vet(matriz, offset, d, q, tdom);
		// usa: matriz, offset, d
		// altera: q
		den = produto_escalar(d, q, tdom);
		// usa: d, q
		// altera: den
		alfa = sigman / den;
		// usa: sigman, den
		// altera: alfa
		escalar_vetor(d, alfa, aux, tdom);
		// usa: alfa, d
		// altera: aux
		soma_vetor(resp, aux, resp, tdom);
		// usa: resp, aux
		// altera: resp
		escalar_vetor(q, alfa, aux, tdom);
		// usa: alfa, q
		// altera: aux
		sub_vetor(r, aux, r, tdom);
		// usa: r, aux
		// altera: r
		sigmav = sigman;
		// usa: sigman
		// altera: sigmav
		sigman = produto_escalar(r, r, tdom);
		// usa: r
		// altera: sigman
		beta = sigman / sigmav;
		// usa: sigman, sigmav
		// altera: beta
		escalar_vetor(d, beta, aux, tdom);
		// usa: beta, d
		// altera: aux
		soma_vetor(r, aux, d, tdom);
		// usa: r, aux
		// altera: d
		it++;
	}
	while (it < itmax);

	free(r);
	free(d);
	free(q);
	free(aux);
	free(aux2);
}

//--------------------------------------------------------
int geramatriz(float **dominio, float **matA, float *vetb, float *vetb_ext, int m, int n)
{
    int k, i, j;
    float C, X;

    k=0;
    C=(float) 1+(4*alpha*deltaT)/(deltaX*deltaX);
    X=(float) -(alpha*deltaT)/(deltaX*deltaX);

    for(i = 1; i < m - 1; i++)
    {
        for(j = 1; j < n - 1; j++)
    	{
	        k=(i-1)*(m-2)+(j-1);

	        //valor da 1a diagonal
	        if(i==1)
	        {
	          vetb_ext[k]+=-X*dominio[i-1][j];
	        }
	        else
	          matA[k][0]=X;

	        //valor da 2a diagonal
	        if(j==1)
	        {
	          vetb_ext[k]+=-X*dominio[i][j-1];
	        }
	        else
	        {
	          matA[k][1]=X;
	        }

	        //valor da diagonal central
	        matA[k][2]=C;
	        vetb[k]=dominio[i][j];

	        //valor da 4a diagonal
	        if(j==n-2)
	        {
	          vetb_ext[k]+=-X*dominio[i][j+1];
	        }
	        else
	          matA[k][3]=X;

	        //valor da 5a diagonal
	        if(i==n-2)
	        {
	          vetb_ext[k]+=-X*dominio[i+1][j];
	        }
	        else
	          matA[k][4]=X;
	    }
    }

    return (m-2)*(n-2);
}

//---------------------------------------------------------

int main()
{
    float **dominio;
    float **matA; int *offset;
    float *vetb;
    float *vetb_ext;
    float *vetx;
    float timea, timeb;

    unsigned int m,n,size;
    unsigned int i,j;
    unsigned int ciclos;


// leitura e alocação dimensão do domínio
    scanf("%d", &ciclos);
    scanf("%d", &m);
    scanf("%d", &n);

    offset = (int*) calloc(5, sizeof(int));
    offset[0]=-(n-2);
    offset[1]=-1;
    offset[2]=0;
    offset[3]=1;
    offset[4]=n-2;

    dominio = (float**) calloc(m, sizeof(float *));
    for(i = 0; i < m; i++)
       dominio[i] = (float*) calloc(n, sizeof(float));

    for(i = 0; i < m; i++)
       for(j = 0; j < n; j++)
          scanf("%f", &dominio[i][j]);

//alocação estrutura de dados da pentadiagonal
//(n-2) e (m-2) por causa da borda
    matA = (float**) calloc(((n-2)*(m-2)), sizeof(float *));
    for(i = 0; i < ((n-2)*(m-2)); i++)
	matA[i] = (float*) calloc(5, sizeof(float));

//alocação do vetor dos termos independentes b e b_ext e de resposta x
    vetb = (float*) calloc(((n-2)*(m-2)), sizeof(float));
    vetb_ext = (float*) calloc(((n-2)*(m-2)), sizeof(float));
    vetx = (float*) calloc(((n-2)*(m-2)), sizeof(float));


//geração dos sistemas
    size=geramatriz(dominio, matA, vetb, vetb_ext, m, n);

//solução iterativa - ciclos
timea=omp_get_wtime();
    for(j=0;j<ciclos;j++)
    {
	    soma_vetor (vetb,vetb_ext,vetb,size);
	    GC(matA, offset, vetb, vetx,  size);
	    copia_vetor(vetb, vetx,size);
    }
timeb=omp_get_wtime();

//saida
 	for(i = 0; i < size; i++)
 	{
 	  printf("%4.1f\n", vetx[i]);
 	}
 	printf("\n ");

printf("%f \n ",timeb-timea);

return(0);

}


