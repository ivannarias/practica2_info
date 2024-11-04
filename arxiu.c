#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 512

float Mat[N][N];
float MatDD[N][N];
float V1[N], V2[N], V3[N], V4[N], V5[N];

void InitData(){
        int i,j;
        srand(334411);

        for( i = 0; i < N; i++ ){
            for( j = 0; j < N; j++ ){
                Mat[i][j]=(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
                if ( (abs(i - j) <= 3) && (i != j)){
                    MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0*RAND_MAX));
                } else if ( i == j ){
                    MatDD[i][j]=(((i*j)%3)?-1:1)*(10000.0*(rand()/(1.0*RAND_MAX)));
                } else { MatDD[i][j] = 0.0;
                }
            }
        }
        for( i = 0; i < N; i++ ){
            V1[i]=(i<N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
            V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
            V3[i]=(((i*j)%5)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
            }
        }

//Punt 1
void PrintVect( float vect[N], int from, int numel ){
        for( int i = from; i < (from+numel+1); i++ ){
            printf( "%.2f, ", vect[i] );
        }
}

//Punt 2        row --> fila         numel --> nombre d'elements      from --> posicó
void PrintRow(float Mat[N][N], int row, int from, int numel) {
	for (int i=0;i<numel; i++){
            printf("%.2f, ",Mat[row][from+i]); // Li diem que imprimeixi els valors de la fila "row" fins a "from"
	}
	printf("\n");
}

//Punt 3
void MultEscalar( float vect[N], float vectres[N], float alfa ){
        for( int i = 0; i < N; i++ ){
            vectres[i] = vect[i]*alfa;
        }
}

//Punt 4
float Scalar(float vect1[N], float vect2[N] ) {
    float producte_escalar = 0;
    for (int i=0;i<N;i++){
        producte_escalar += vect1[i]*vect2[i];
    }
    return producte_escalar;
}

//Punt 5
float Magnitude( float vect[N] ){
        float suma = 0;
        for( int i = 0; i < N; i++ ){
            suma += vect[i] * vect[i];
        }
        return sqrt(suma);
}

//Punt 7
void Projection( float vect1[N], float vect2[N], float vectres[N] ){
        float divisor = Scalar(vect1,vect2);
        float magV1 = Magnitude(vect2);
        float escalar = (divisor/magV1);
        vectres[0] = escalar*vect2[0];
        vectres[1] = escalar*vect2[1];
        vectres[2] = escalar*vect2[2];
}

//Punt 8
float Infininorm( float M[N][N] ){
        float maxSuma = 0;
        for( int i = 0; i < N; i++ ){
            float sumaFila = 0;
            for( int j = 0; j < N; j++ ){
                sumaFila += fabs(M[i][j]);
            }
            if( sumaFila > maxSuma ){
                maxSuma = sumaFila;
            }

        }
        return maxSuma;
}

//Punt 9
float Onenorm( float M[N][N] ){
	float maxSuma = 0;
	for( int j = 0; j < N; j++ ){
	    float sumaColumna = 0;
	    for( int i = 0; i < N; i++ ){
	    	sumaColumna += fabs(M[i][j]);
	    }
	    if( sumaColumna > maxSuma ){
	    	maxSuma = sumaColumna;
	    }

	}
	return maxSuma;
}

//Punt 11
int DiagonalDom( float M[N][N] ){
	for (int i = 0; i < N; i++) {
        float sumaNoDiagonal = 0; // Inicialitzem una variable que sumarà aquells valors que no estàn a la diagonal.
		for (int j = 0; j < N; j++) {
        	    if (i != j) {
                	sumaNoDiagonal += fabs(M[i][j]);  // Sumem només els valors que no estàn a la diagonal
            	}
        	}
		if (fabs(M[i][i]) < sumaNoDiagonal) { // Comparem la suma dels valors que no estàn a la diagonal amb el valor de la diagonal de la linia "i".
        	    return 0;  // No es diagonal dominant
        	}
    	} // I axò ho fem per CADA linia. En cas que no s'hagi complert a cap linia...
    return 1;  // És diagonal dominant
}


int main (){
    InitData();
    printf( "El vector 1 és: \n" );
    PrintVect(V1, 0, 9);
    printf( "\n" );

    printf( "\nLa multiplicació escalar és:" );
    MultEscalar(V1, V4, 2);
    printf( "\n" );
    for( int i = 0; i < N; i++ ){
            printf( "%.2f, ", V4[i] );
    }

    printf( "\n" );

    float magnitud = Magnitude(V1);
    printf("\nLa magnitud del vector és: %.2f\n", magnitud);

    Projection(V2,V3,V5);
    printf("\nL'escalar dels vectors demanats és: (%.2f, %.2f, %.2f)\n", V5[0],V5[1],V5[2]);

    float InfNorm = Infininorm(Mat);
    float UNorm = Onenorm(Mat);
    printf("\nLa Infinorma de la matriu és: %.2f\n",InfNorm);
    printf("\nLa Norma Ú de la matriu és: %.2f\n",UNorm);

    int dom = DiagonalDom(Mat);
    if( dom == 1 ){
	printf("\nLa matriu és Diagonal Dominant");
    } else if( dom == 0 ){
	printf("\nLa matriu no és Diagonal Dominant");
    }
}
