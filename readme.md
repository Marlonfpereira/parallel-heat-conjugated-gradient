# Transferência de Calor + Gradiente Conjugado
### Aplicação das técnicas de paralelização em OpenMP e OpenAcc.
Discentes: Marlon Fabichacki Pereira e Ronaldo Drecksler

## Compilacao

Para compilar os programas, utilize as regras de makefile `serial`, `openmp` e `openacc`. Ex:

```
make serial openmp openacc
```


## Execucao

Para executar os programas, utilize no comando make a técnica desejada (`xserial`, `xopenmp` ou `xopenacc`), seguida do arquivo de entrada (`small`, `big` ou `input FILE={path}`. "path" toma por ./ a pasta inputs) com `run` ao final. Ex:

```
make xserial small run

make xopenmp input FILE=arquivo_dentro_de_inputs.dat run
```
