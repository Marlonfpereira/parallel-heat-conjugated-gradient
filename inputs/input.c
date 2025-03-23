#include<stdio.h>

#define borda 10.0
#define interior 1.0

int main()
{
    unsigned int m,n,i,j,t;
    float b=borda;
    float w=interior;

    printf("Passos de tempo: ");
    scanf("%u",&t);
    printf("Largura: ");
    scanf("%u",&m);
    printf("Altura: ");
    scanf("%u",&n);
    
    printf("%d\n",t);
    printf("%d %d\n",m,n);
    
    for(i=0;i<m;i++)
        printf("%4.1f ",b);
    printf("\n");    


    for(i=0;i<n-2;i++)
    {
       printf("%4.1f ",b);
       for(j=0;j<n-2;j++)
          printf("%4.1f ",w);  
       printf("%4.1f\n",b);   
    }


    for(i=0;i<m;i++)
        printf("%4.1f ",b);
    printf("\n");    

    return 0;
}


