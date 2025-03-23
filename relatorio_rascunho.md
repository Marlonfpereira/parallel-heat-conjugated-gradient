METODOLOGIA

## Como foi paralizado
Rodar o profiler para identificar os hotspot
para -O0 e -O3

Tranformar todo o codigo para inline -> o -O3 executa tudo em inline

Foi então realizado a paralelização dos loops e reduções

Identificar os pontos de transferencia de dados

Colocar as clausulas de paralelização para loops e reduções

## Analise quantitativa

Coletou tempos cronologicos do serial -O0 e -O3

-O0
    small: 29.532715
    big  : 288.633057 << picos bizarros de swap????>>

-O3
    small: 10.424561 
    big  : 87.129395

// 16= 24.332031  10= 4.610840  8= 2.894531 4= 3.468750 2= 5.449219 1= 10.573242
Tempos do openMP -O3 small
    - 2 theads   = 5.449
    - 4 theads   = 3.42345
    - 8 threads  = 2.89453 
    - 10 threads = 4.6108
    - 16 threads = 24.4563

Tempos do openMP -O3 big
    - 2 theads   = 46.782227
    - 4 theads   = 35.902954
    - 8 threads  = 32.339453
    - 10 threads = 36.992065 
    - 16 threads = 80.180054 


Tempos do openACC com -fast
    -Ofast 
        small: 20.974804
        big  : 100.792057

    -fast 
        small: 21.828314
        big  : 101.411133 





