//#include <stdio.h>  //sprintf, FILE
//#include <string.h> //strcpy, strcat
//#include <math.h>   //pow, floor
#include "maxiutils.h"  //tiene los cuantificadores
#include <time.h>   //time
#include <random>
#include <conio.h> //getch, solo para MS-DOS

using namespace std; // Agrego el espacio de nombres del estandar de C++, zafo de poner std:: adelante de cada l�nea.

int main()
{
    int Bases[] = {2}; // Vector con las bases que quiero probar
        int NBases = (sizeof(Bases)/sizeof(Bases[0]));
    double Margins[2] = {0, 1}; //Los m�rgenes de la PDF-Val
    unsigned long int NInitialConditions = 100; // Es la cantidad de condiciones iniciales diferentes de los que se larga el atractor.
    unsigned long int NIter = 1e7; // Es el largo de cada atractor

    unsigned long int Bins = 1024; // Cantidad de bines del histograma
    unsigned long int DimEmb = 6; // Dimensi�n de embedding para MP, BP y BPW

    double* Map; //Declare the pointer
        Map = (double*) malloc (sizeof(double) * (NIter + 1)); //Creates the array. It has one more postition at first for the length
        Map[0] = (double)NIter;

    double Hval, Qval, Cval, Hbp, Qbp, Cbp, Hbpw, Qbpw, Cbpw, MP, Period; // Vectores en donde van guardados los cuantificadores

    FILE *Results = fopen("Flotante.dat","w"); //Abre archivo de resultados
    fprintf(Results, "Map\tHval\tQval\tCval\tHbp\tQbp\tCbp\tHbpw\tQbpw\tCbpw\tMP\tPeriod\n"); //Escribe encabezado en archivo

    double InitialConditions[NInitialConditions]; // Declaro el vector de condiciones iniciales
    random_device rd;
    mt19937_64 gen(rd()); // Creo un objeto del tipo mt19937 en presici�n double con semilla en 0
    uniform_real_distribution <double> dis(0, 1); // Defino el tipo de distribuci�n
    for (unsigned long int iInitialCondition = 0; iInitialCondition < NInitialConditions; iInitialCondition++) // Lleno un vector de condiciones iniciales
    {
        InitialConditions[iInitialCondition] = dis(gen); // Sortea el primer valor del mapa, es una variable uniformemente distribu�da entre 0 y 1
    }
    printf("Generadas %ld condiciones iniciales\n\n", NInitialConditions);

    /*********Log�stico**************/
    printf("Log�stico\n");
    for (unsigned int iInitialCondition = 0; iInitialCondition < NInitialConditions; iInitialCondition++) // Va sorteando condiciones iniciales
    {

        Map[1] = InitialConditions[iInitialCondition]; // floorl sirve para long double, como son mapas positivos puedo usar floor en vez de floor

        printf("\tCondicion inicial %d/%d = %.32f\n", (int)iInitialCondition+1, (int)NInitialConditions, Map[1]); // Para debuguear

        for (unsigned long int iMap = 1; iMap < NIter; iMap++) // Va riterando el mapa log�stico
        {
            Map[iMap+1] =  4*Map[iMap]*(1-Map[iMap]); // Mapa log�stico, x[n] = r*x[n-1]*(1-x[n-1]), ca�tico con r=4. Ni la resta ni la multiplicaci�n por un entero generan fraccionarios
        } // Ac� ya tengo el atractor guardado en el vector Map

        double* PDFval = PDF_val(Map, Bins, Margins, "normalyzed"); // Genera el histograma de patrones de �rden
        Hval = entropy(PDFval, "normalyzed"); // Le calcula la entrop�a y la suma para el promedio
        Qval = disequilibrum(PDFval, "normalyzed"); // Le calcula el desequilibrio
        Cval = Hval*Qval; // Le calcula la complejidad
        free(PDFval); // Libera el vector que contiene al histograma

        double* PDFbp = PDF_BP(Map, DimEmb, "overlapped", "normalyzed"); // Genera el histograma de patrones de �rden
        Hbp = entropy(PDFbp, "normalyzed"); // Le calcula la entrop�a
        Qbp = disequilibrum(PDFbp, "normalyzed"); // Le calcula el desequilibrio
        Cbp = Hbp*Qbp; // Le calcula la complejidad
        MP = missing_patterns(PDFbp);
        free(PDFbp); // Libera el vector que contiene al histograma

        double* PDFbpw = PDF_BPW(Map, DimEmb, "overlapped"); // Genera el histograma de patrones de �rden
        Hbpw = entropy(PDFbpw, "normalyzed"); // Le calcula la entrop�a
        Qbpw = disequilibrum(PDFbpw, "normalyzed"); // Le calcula el desequilibrio
        Cbpw = Hbpw*Qbpw; // Le calcula la complejidad
        free(PDFbpw); // Libera el vector que contiene al histograma

        Period = find_period(Map, 1); //Para mapas unidimensionales, para los switch voy a tener que usar dimensi�n 2.

        fprintf(Results,"Logistico\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", Hval, Qval, Cval, Hbp, Qbp, Cbp, Hbpw, Qbpw, Cbpw, MP, Period); //Guarda los valores en el archivo de salida, escribo la condici�n inicial para evaluar el comportamiento del rand()
    }

    /*********Tent**************/
    printf("Tent\n");
    for (unsigned int iInitialCondition = 0; iInitialCondition < NInitialConditions; iInitialCondition++) // Va sorteando condiciones iniciales
    {

        Map[1] = InitialConditions[iInitialCondition]; // floorl sirve para long double, como son mapas positivos puedo usar floor en vez de floor

        printf("\tCondicion inicial %d/%d = %.32f\n", (int)iInitialCondition+1, (int)NInitialConditions, Map[1]); // Para debuguear

        for (unsigned long int iMap = 1; iMap < NIter; iMap++) // Va riterando el mapa log�stico
        {
            if (Map[iMap] < 0.5)
            {
                Map[iMap+1] = 2*Map[iMap];
            }
            else
            {
                Map[iMap+1] = 2*(1 - Map[iMap]);
            }
        } // Ac� ya tengo el atractor guardado en el vector Map

        double* PDFval = PDF_val(Map, Bins, Margins, "normalyzed"); // Genera el histograma de patrones de �rden
        Hval = entropy(PDFval, "normalyzed"); // Le calcula la entrop�a y la suma para el promedio
        Qval = disequilibrum(PDFval, "normalyzed"); // Le calcula el desequilibrio
        Cval = Hval*Qval; // Le calcula la complejidad
        free(PDFval); // Libera el vector que contiene al histograma

        double* PDFbp = PDF_BP(Map, DimEmb, "overlapped", "normalyzed"); // Genera el histograma de patrones de �rden
        Hbp = entropy(PDFbp, "normalyzed"); // Le calcula la entrop�a
        Qbp = disequilibrum(PDFbp, "normalyzed"); // Le calcula el desequilibrio
        Cbp = Hbp*Qbp; // Le calcula la complejidad
        MP = missing_patterns(PDFbp);
        free(PDFbp); // Libera el vector que contiene al histograma

        double* PDFbpw = PDF_BPW(Map, DimEmb, "overlapped"); // Genera el histograma de patrones de �rden
        Hbpw = entropy(PDFbpw, "normalyzed"); // Le calcula la entrop�a
        Qbpw = disequilibrum(PDFbpw, "normalyzed"); // Le calcula el desequilibrio
        Cbpw = Hbpw*Qbpw; // Le calcula la complejidad
        free(PDFbpw); // Libera el vector que contiene al histograma

        Period = find_period(Map, 1); //Para mapas unidimensionales, para los switch voy a tener que usar dimensi�n 2.

        fprintf(Results,"Tent\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", Hval, Qval, Cval, Hbp, Qbp, Cbp, Hbpw, Qbpw, Cbpw, MP, Period); //Guarda los valores en el archivo de salida, escribo la condici�n inicial para evaluar el comportamiento del rand()
    }

    /*********SkewTent**************/
    printf("SkewTent\n");

    double gamma = .51; //Es el par�metro del skew tent, para gamma = .5 el mapa es el tent del otro c�digo

    for (unsigned int iInitialCondition = 0; iInitialCondition < NInitialConditions; iInitialCondition++) // Va sorteando condiciones iniciales
    {

        Map[1] = InitialConditions[iInitialCondition]; // floorl sirve para long double, como son mapas positivos puedo usar floor en vez de floor

        printf("\tCondicion inicial %d/%d = %.32f\n", (int)iInitialCondition+1, (int)NInitialConditions, Map[1]); // Para debuguear

        for (unsigned long int iMap = 1; iMap < NIter; iMap++) // Va riterando el mapa log�stico
        {
            if (Map[iMap] < gamma)
            {
                Map[iMap+1] = Map[iMap]/gamma;
            }
            else
            {
                Map[iMap+1] = (1 - Map[iMap])/(1 - gamma);
            }
        } // Ac� ya tengo el atractor guardado en el vector Map

        double* PDFval = PDF_val(Map, Bins, Margins, "normalyzed"); // Genera el histograma de patrones de �rden
        Hval = entropy(PDFval, "normalyzed"); // Le calcula la entrop�a y la suma para el promedio
        Qval = disequilibrum(PDFval, "normalyzed"); // Le calcula el desequilibrio
        Cval = Hval*Qval; // Le calcula la complejidad
        free(PDFval); // Libera el vector que contiene al histograma

        double* PDFbp = PDF_BP(Map, DimEmb, "overlapped", "normalyzed"); // Genera el histograma de patrones de �rden
        Hbp = entropy(PDFbp, "normalyzed"); // Le calcula la entrop�a
        Qbp = disequilibrum(PDFbp, "normalyzed"); // Le calcula el desequilibrio
        Cbp = Hbp*Qbp; // Le calcula la complejidad
        MP = missing_patterns(PDFbp);
        free(PDFbp); // Libera el vector que contiene al histograma

        double* PDFbpw = PDF_BPW(Map, DimEmb, "overlapped"); // Genera el histograma de patrones de �rden
        Hbpw = entropy(PDFbpw, "normalyzed"); // Le calcula la entrop�a
        Qbpw = disequilibrum(PDFbpw, "normalyzed"); // Le calcula el desequilibrio
        Cbpw = Hbpw*Qbpw; // Le calcula la complejidad
        free(PDFbpw); // Libera el vector que contiene al histograma

        Period = find_period(Map, 1); //Para mapas unidimensionales, para los switch voy a tener que usar dimensi�n 2.

        fprintf(Results,"SkewTent\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", Hval, Qval, Cval, Hbp, Qbp, Cbp, Hbpw, Qbpw, Cbpw, MP, Period); //Guarda los valores en el archivo de salida, escribo la condici�n inicial para evaluar el comportamiento del rand()
    }

    /*********Switch**************/
    printf("Switch\n");
    for (unsigned int iInitialCondition = 0; iInitialCondition < NInitialConditions; iInitialCondition++) // Va sorteando condiciones iniciales
    {

        Map[1] = InitialConditions[iInitialCondition]; // floorl sirve para long double, como son mapas positivos puedo usar floor en vez de floor

        printf("\tCondicion inicial %d/%d = %.32f\n", (int)iInitialCondition+1, (int)NInitialConditions, Map[1]); // Para debuguear

        for (unsigned long int iMap = 1; iMap < NIter; iMap = iMap + 2) // Va iterando switch
        {
            //Itero tent
            if (Map[iMap] < 0.5)
            {
                Map[iMap+1] = 2*Map[iMap];
            }
            else
            {
                Map[iMap+1] = 2*(1 - Map[iMap]);
            }

            //Itero logistico
            Map[iMap+2] =  4*Map[iMap+1]*(1-Map[iMap+1]);

        } // Ac� ya tengo el atractor guardado en el vector Map

        double* PDFval = PDF_val(Map, Bins, Margins, "normalyzed"); // Genera el histograma de patrones de �rden
        Hval = entropy(PDFval, "normalyzed"); // Le calcula la entrop�a y la suma para el promedio
        Qval = disequilibrum(PDFval, "normalyzed"); // Le calcula el desequilibrio
        Cval = Hval*Qval; // Le calcula la complejidad
        free(PDFval); // Libera el vector que contiene al histograma

        double* PDFbp = PDF_BP(Map, DimEmb, "overlapped", "normalyzed"); // Genera el histograma de patrones de �rden
        Hbp = entropy(PDFbp, "normalyzed"); // Le calcula la entrop�a
        Qbp = disequilibrum(PDFbp, "normalyzed"); // Le calcula el desequilibrio
        Cbp = Hbp*Qbp; // Le calcula la complejidad
        MP = missing_patterns(PDFbp);
        free(PDFbp); // Libera el vector que contiene al histograma

        double* PDFbpw = PDF_BPW(Map, DimEmb, "overlapped"); // Genera el histograma de patrones de �rden
        Hbpw = entropy(PDFbpw, "normalyzed"); // Le calcula la entrop�a
        Qbpw = disequilibrum(PDFbpw, "normalyzed"); // Le calcula el desequilibrio
        Cbpw = Hbpw*Qbpw; // Le calcula la complejidad
        free(PDFbpw); // Libera el vector que contiene al histograma

        Period = find_period(Map, 1); //Para mapas unidimensionales, para los switch voy a tener que usar dimensi�n 2.

        fprintf(Results,"Switch\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", Hval, Qval, Cval, Hbp, Qbp, Cbp, Hbpw, Qbpw, Cbpw, MP, Period); //Guarda los valores en el archivo de salida, escribo la condici�n inicial para evaluar el comportamiento del rand()
    }

    /*********SwitchEven**************/
    printf("SwitchEven\n");
    for (unsigned int iInitialCondition = 0; iInitialCondition < NInitialConditions; iInitialCondition++) // Va sorteando condiciones iniciales
    {

        Map[1] = InitialConditions[iInitialCondition]; // floorl sirve para long double, como son mapas positivos puedo usar floor en vez de floor

        printf("\tCondicion inicial %d/%d = %.32f\n", (int)iInitialCondition+1, (int)NInitialConditions, Map[1]); // Para debuguear

        for (unsigned long int iMap = 1; iMap < NIter; iMap++) // Va iterando switch
        {
            //Itero tent
            if (Map[iMap] < 0.5)
            {
                Map[iMap+1] = 2*Map[iMap];
            }
            else
            {
                Map[iMap+1] = 2*(1 - Map[iMap]);
            }

            //Itero logistico
            Map[iMap+1] =  4*Map[iMap+1]*(1-Map[iMap+1]);

        } // Ac� ya tengo el atractor guardado en el vector Map

        double* PDFval = PDF_val(Map, Bins, Margins, "normalyzed"); // Genera el histograma de patrones de �rden
        Hval = entropy(PDFval, "normalyzed"); // Le calcula la entrop�a y la suma para el promedio
        Qval = disequilibrum(PDFval, "normalyzed"); // Le calcula el desequilibrio
        Cval = Hval*Qval; // Le calcula la complejidad
        free(PDFval); // Libera el vector que contiene al histograma

        double* PDFbp = PDF_BP(Map, DimEmb, "overlapped", "normalyzed"); // Genera el histograma de patrones de �rden
        Hbp = entropy(PDFbp, "normalyzed"); // Le calcula la entrop�a
        Qbp = disequilibrum(PDFbp, "normalyzed"); // Le calcula el desequilibrio
        Cbp = Hbp*Qbp; // Le calcula la complejidad
        MP = missing_patterns(PDFbp);
        free(PDFbp); // Libera el vector que contiene al histograma

        double* PDFbpw = PDF_BPW(Map, DimEmb, "overlapped"); // Genera el histograma de patrones de �rden
        Hbpw = entropy(PDFbpw, "normalyzed"); // Le calcula la entrop�a
        Qbpw = disequilibrum(PDFbpw, "normalyzed"); // Le calcula el desequilibrio
        Cbpw = Hbpw*Qbpw; // Le calcula la complejidad
        free(PDFbpw); // Libera el vector que contiene al histograma

        Period = find_period(Map, 1); //Para mapas unidimensionales, para los switch voy a tener que usar dimensi�n 2.

        fprintf(Results,"SwitchEven\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", Hval, Qval, Cval, Hbp, Qbp, Cbp, Hbpw, Qbpw, Cbpw, MP, Period); //Guarda los valores en el archivo de salida, escribo la condici�n inicial para evaluar el comportamiento del rand()
    }

    /*********SwitchOdd**************/
    printf("SwitchOdd\n");
    for (unsigned int iInitialCondition = 0; iInitialCondition < NInitialConditions; iInitialCondition++) // Va sorteando condiciones iniciales
    {

        Map[1] = InitialConditions[iInitialCondition]; // floorl sirve para long double, como son mapas positivos puedo usar floor en vez de floor

        printf("\tCondicion inicial %d/%d = %.32f\n", (int)iInitialCondition+1, (int)NInitialConditions, Map[1]); // Para debuguear

        for (unsigned long int iMap = 1; iMap < NIter; iMap++) // Va iterando switch
        {
            //Itero logistico
            Map[iMap+1] =  4*Map[iMap]*(1-Map[iMap]);

            //Itero tent
            if (Map[iMap+1] < 0.5)
            {
                Map[iMap+1] = 2*Map[iMap+1];
            }
            else
            {
                Map[iMap+1] = 2*(1 - Map[iMap+1]);
            }

        } // Ac� ya tengo el atractor guardado en el vector Map

        double* PDFval = PDF_val(Map, Bins, Margins, "normalyzed"); // Genera el histograma de patrones de �rden
        Hval = entropy(PDFval, "normalyzed"); // Le calcula la entrop�a y la suma para el promedio
        Qval = disequilibrum(PDFval, "normalyzed"); // Le calcula el desequilibrio
        Cval = Hval*Qval; // Le calcula la complejidad
        free(PDFval); // Libera el vector que contiene al histograma

        double* PDFbp = PDF_BP(Map, DimEmb, "overlapped", "normalyzed"); // Genera el histograma de patrones de �rden
        Hbp = entropy(PDFbp, "normalyzed"); // Le calcula la entrop�a
        Qbp = disequilibrum(PDFbp, "normalyzed"); // Le calcula el desequilibrio
        Cbp = Hbp*Qbp; // Le calcula la complejidad
        MP = missing_patterns(PDFbp);
        free(PDFbp); // Libera el vector que contiene al histograma

        double* PDFbpw = PDF_BPW(Map, DimEmb, "overlapped"); // Genera el histograma de patrones de �rden
        Hbpw = entropy(PDFbpw, "normalyzed"); // Le calcula la entrop�a
        Qbpw = disequilibrum(PDFbpw, "normalyzed"); // Le calcula el desequilibrio
        Cbpw = Hbpw*Qbpw; // Le calcula la complejidad
        free(PDFbpw); // Libera el vector que contiene al histograma

        Period = find_period(Map, 1); //Para mapas unidimensionales, para los switch voy a tener que usar dimensi�n 2.

        fprintf(Results,"SwitchOdd\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n", Hval, Qval, Cval, Hbp, Qbp, Cbp, Hbpw, Qbpw, Cbpw, MP, Period); //Guarda los valores en el archivo de salida, escribo la condici�n inicial para evaluar el comportamiento del rand()
    }

    fclose(Results); // Cierra el archivo de salida, queda un archivo por base para este oscilador
    printf("\nPresionar una tecla");
    getch(); //Avisa que termin� y espera una tecla
}
