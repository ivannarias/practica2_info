#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 512

float Mat[N][N];
float MatDD[N][N];
float V1[N], V2[N], V3[N], V4[N], V5[N], Vaux[N];

void InitData(){ // Apliquem el codi que ens dona l'exercici per generar matrius i vectors amb 512 elements pseudo-aleatoris
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
        for( int i = from; i < (from+numel+1); i++ ){ // Imprimim "numel" elements des de l'índex "from".
            printf( "%.2f, ", vect[i] );
        }
}

//Punt 2        row --> fila         numel --> nombre d'elements      from --> posicó inicial
void PrintRow(float Mat[N][N], int row, int from, int numel) {
	for (int i=0;i<numel; i++){
            printf("%.2f, ",Mat[row][from+i]); // Li diem que imprimeixi "numel" valors de la fila "row" des de "from".
	}
	printf("\n");
}

//Punt 3
void MultEscalar( float vect[N], float vectres[N], float alfa ){
        for( int i = 0; i < N; i++ ){
            vectres[i] = vect[i]*alfa; // Fem la multiplicació escalar multiplicant els elements per alfa.
        }
}

//Punt 4
float Scalar(float vect1[N], float vect2[N] ) {
    float producte_escalar = 0;
    for (int i=0;i<N;i++){
        producte_escalar += vect1[i]*vect2[i]; // Fem la multiplicació escalar seguint la seva definició, multiplicant els vectors entre si i sumant els resultats.
    }
    return producte_escalar;
}

//Punt 5
float Magnitude( float vect[N] ){
        float suma = 0;
        for( int i = 0; i < N; i++ ){
            suma += vect[i] * vect[i]; // Fem la suma dels valors del vector al cuadrat
        }
        return sqrt(suma); // Fem l'arrel de la suma per obtenir la magnitud del vector.
}

//punt 6
int Ortogonal(float vect1[N], float vect2[N]) {
    float orto;
    orto = Scalar(vect1,vect2); // Assignem la variable auxiliar orto al producte escalar dels dos vectors.
    if (orto == 0) { // En cas que el producte escalar sigui 0, els vectors són perpendiculars (ortogonals).
        return 1;
    } else {
        return 0;
    }
}

//Punt 7
void Projection( float vect1[N], float vect2[N], float vectres[N] ){
        float divisor = Scalar(vect1,vect2); // Fem servir la nostra funció per obtenir el producte escalar dels dos vectors.
        float magV1 = Magnitude(vect2); // Fem servir la nostra funció per obtenir la magnitud del segón vector.
        float escalar = (divisor/magV1);
        vectres[0] = escalar*vect2[0]; // Apliquem la definció de la projecció i assignem els valors de cada índex del vector resultat.
        vectres[1] = escalar*vect2[1];
        vectres[2] = escalar*vect2[2];
}

//Punt 8
float Infininorm( float M[N][N] ){
        float maxSuma = 0; // Iniciem la variable maxSuma per posteriorment comparar-la.
        for( int i = 0; i < N; i++ ){
            float sumaFila = 0; // Inicialitzem a CADA ITERACIÓ sobre i la variable que sumarà els valors de cada línia.
            for( int j = 0; j < N; j++ ){
                sumaFila += fabs(M[i][j]); // Es suma el valor absolut flotant de cada valor de la fila.
            }
            if( sumaFila > maxSuma ){ // Es comparen els dos valors i actualitzem maxSuma al valor de la suma de la fila en cas que aquest superi a l'anterior registre maxSuma (anterior fila).
                maxSuma = sumaFila;
            }

        }
        return maxSuma; // Retornem la suma de la fila amb maxSuma.
}

//Punt 9
float Onenorm( float M[N][N] ){
	float maxSuma = 0; // Iniciem la variable maxSuma per posteriorment comparar-la.
	for( int j = 0; j < N; j++ ){
	    float sumaColumna = 0; // Inicialitzem a CADA ITERACIÓ sobre j la variable que sumarà els valors de cada columna.
	    for( int i = 0; i < N; i++ ){
	    	sumaColumna += fabs(M[i][j]); // Es suma el valor absolut flotant de cada valor de la columna.
	    }
	    if( sumaColumna > maxSuma ){ // Es comparen els dos valors i actualitzem maxSuma al valor de la suma de la columna en cas que aquest superi a l'anterior registre maxSuma (anterior columna).
	    	maxSuma = sumaColumna;
	    }

	}
	return maxSuma; // Retornem la suma de la columna amb maxSuma
}

//punt 10
float NormFrobenius(float M[N][N]) {
    float suma = 0; // Inicialitzem la variable per seguir la suma.
    for (int i=0;i<N; i++){
        for (int j=0;j<N; j++){
            suma += ((M[i][j])*(M[i][j])); // Sumem els valors de la matriu a "suma".
        }
    }
    float Frobenius;
    Frobenius = sqrt(suma); // Retornem l'arrel de la suma.
    return Frobenius;
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

//punt 12
void Matriu_x_Vector(float M[N][N], float vect[N], float vectres[N]) {
    float suma_linia=0; // Iniciem una variable per comptar la suma de cada linia.
    for (int i=0;i<N; i++){
        for (int j=0;j<N; j++){
            suma_linia += M[i][j] * vect[j]; // Sumem a la variable "suma" els valors de la matriu multiplicats pel component [j] del vector.
        }
        vectres[i] = suma_linia; // Assignem el resultat de la suma al vector resultat, i posteriorment reinicialitzem la variable suma per sumar la pròxima fila.
        suma_linia=0;
    }
}

//Punt 13
int Jacobi( float M[N][N] , float vect[N], float vectres[N], unsigned iter ){
	float tolerancia = 1e-6;
	int dom = DiagonalDom (M);
	if (dom == 0) { // Abans de començar, veiem si el vector NO és diagonal dominant, per així retornar 0 i veure que no podem fer servir Jacobi.
	    return 0;
	}

	float v_nou[N];  // Vector per a la nova iteració

	for (int i = 0; i < N; i++) {
            vectres[i] = 0.0; // Inicialitzem el vector resultat.
    	}

	for (unsigned k = 0; k < iter; k++) {
            int conver = 1;  // Variable per a comprovar la convergència

	    for (int i = 0; i < N; i++) {
                float suma = 0.0; // Inicialitzem una variable per fer la suma dels valors que no són a la diagonal.
                for (int j = 0; j < N; j++) {
                    if (j != i) {
                        suma += M[i][j] * vectres[j]; // Sumem només els que NO són a la diagonal.
                    }
                }
                v_nou[i] = (vect[i] - suma) / M[i][i]; // Calculem l'aproximació del valor del vector nou.

		if (fabs(v_nou[i] - vectres[i]) > tolerancia) {
                conver = 0;
 	    }
	}

	for (int i = 0; i < N; i++) {
            vectres[i] = v_nou[i]; // Assignem els valors del nostre vector auxiliar al vector resultat.
        }

	if (conver) {
            return 1;  // Hem trobat una solució que compleix la tolerància
        }
    }
}

int main (){
    InitData();
    printf( "El vector 1 és: \n" );
    PrintVect(V1, 0, 9);
    printf( "\n" );

    PrintRow(Mat,0,0,10);
    PrintRow(Mat,100,0,10);
    printf("\n");

    printf( "\nLa multiplicació escalar és:" );
    MultEscalar(V1, V4, 2);
    printf( "\n" );
    for( int i = 0; i < N; i++ ){
            printf( "%.2f, ", V4[i] );
    }

    printf( "\n" );

    float producte_escalar;
    producte_escalar=Scalar(V1,V2);
    printf("%f, ", producte_escalar);
    producte_escalar=Scalar(V1,V3);
    printf("%f, ", producte_escalar);
    producte_escalar=Scalar(V2,V3);
    printf("%f\n", producte_escalar);
    printf("\n");

    float magnitud = Magnitude(V1);
    printf("La magnitud del vector és: %.2f\n", magnitud);

    printf("\n");
    int ortogonal;
    ortogonal = Ortogonal(V1,V2);
    if (ortogonal==1){
        printf("V1 i V2 són ortogonals.\n");
    }
    ortogonal = Ortogonal(V1,V3);
    if (ortogonal==1){
        printf("V1 i V3 són ortogonals.\n");
    }
    ortogonal = Ortogonal(V3,V2);
    if (ortogonal==1){
        printf("V2 i V3 són ortogonals.\n");
    }

    Projection(V2,V3,V5);
    printf("\nL'escalar dels vectors demanats és: (%.2f, %.2f, %.2f)\n", V5[0],V5[1],V5[2]);

    float InfNorm = Infininorm(Mat);
    float UNorm = Onenorm(Mat);
    printf("\nLa Infinorma de la matriu és: %.2f\n",InfNorm);
    printf("\nLa Norma Ú de la matriu és: %.2f\n",UNorm);

    printf("\nNorma de Frobenius per Mat = %.3f\n",NormFrobenius(Mat));
    printf("\nNorma de Frobenius per Mat = %.3f\n",NormFrobenius(MatDD));

    int dom = DiagonalDom(Mat);
    if( dom == 1 ){
	printf("\nLa matriu és Diagonal Dominant\n");
    } else if( dom == 0 ){
	printf("\nLa matriu no és Diagonal Dominant\n");
    }

    // Primer càlcul amb 1 iteració amb la matriu MatDD
    printf("\nEls elements 0 a 9 de la solució (1 iteració) del sistema d'equacions són:\n");
    if (Jacobi(MatDD, V3, Vaux, 1)) {
        for (int i = 0; i < 10; i++) {
            printf("%.6f ", Vaux[i]);
        }
        printf("\n");
    } else {
        printf("\nLa matriu no és diagonal dominant, no es pot aplicar Jacobi.\n");
    }

    // Segon càlcul amb 1000 iteracions amb la matriu MatDD
    printf("\nEls elements 0 a 9 de la solució (1000 iteracions) del sistema d'equacions són:\n");
    if (Jacobi(MatDD, V3, Vaux, 1000)) {
        for (int i = 0; i < 10; i++) {
            printf("%.6f ", Vaux[i]);
        }
        printf("\n");
    } else {
        printf("\nLa matriu no és diagonal dominant, no es pot aplicar Jacobi.\n");
    }

    // Càlcul amb la matriu no diagonal dominant Mat
    printf("\nIntent de resoldre el sistema d'equacions amb una matriu no diagonal dominant:\n"); // Donarà que no es pot fer, però ho comprovem igualment.
    if (Jacobi(Mat, V3, Vaux, 1000)) {
        for (int i = 0; i < 10; i++) {
            printf("%.6f ", Vaux[i]);
        }
        printf("\n");
    } else {
        printf("La matriu M no és diagonal dominant, no es pot aplicar Jacobi.\n");
    }


}
