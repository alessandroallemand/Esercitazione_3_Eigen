#include <iostream>
#include <iomanip>

#include "Eigen/Eigen"
using namespace std;
using namespace Eigen;

// funzione per la soluzione del sistema lineare attraverso il metodo palu.
VectorXd palumatrix(const MatrixXd &A, VectorXd &b){
        MatrixXd X = A.lu().solve(b);
        cout << "the solution of the linear sistem is X = " << X << endl;
    return X;
}
// funzione per il calcolo dell'errore relativo 
double ErrorSolution(const VectorXd &xg, const VectorXd &xa){
    double err = (xg-xa).norm()/xg.norm();
    return err;
}

// funzione per la soluzione del sistema lineare attraverso il metodo QR
VectorXd qrmatrix(const MatrixXd &A, VectorXd &b){
    VectorXd X = A.colPivHouseholderQr().solve(b);
    return X;
}

int main()
{
    cout << fixed << setprecision(1) << scientific;

    VectorXd xcorretta(2);  
    xcorretta << -1.0e+00, -1.0e+00;


    //PUNTO 1  

	//Definisco la matrice A1
	MatrixXd A1(2,2);
	A1 << 5.547001962252291e-01, -3.770900990025203e-02, 8.320502943378437e-01, -9.992887623566787e-01;
	
	//Definisco il vettore b1
	VectorXd b1(2);
	b1 << -5.169911863249772e-01, 1.672384680188350e-01;

    // se la matrice è non singolare applico la fattorizzazione PA=LU e QR
    double r1 = A1.determinant();
    if(abs(r1) <1.0e-12 )
		cout << "A1 is singular" << endl;
	else 
		cout << "A1 is not singular" << endl;

        //applico la funzione per la soluzione del sistema lineare ottenuto tramite la fattorizzazione QR
        VectorXd X1palu = palumatrix(A1,b1);
        cout<< "questo è l'output di palumatrix nel main, nonchè soluzione del sistema lineare con fattorizzazione PALU \n"<< X1palu << endl;
        
          
        
        //applico la funzione per il calcolo dell'errore relativo,in questo caso, al sistema 1, risolto con Fattorizzazione PALU
        double err1palu= ErrorSolution(xcorretta, X1palu);
        cout << "l'errore relativo associato al primo risultato con palu e corrisponde a: err1palu = "<< err1palu << endl;
        
        //applico la funzione per la soluzione del sistema lineare ottenuto tramite la fattorizzazione QR
        VectorXd X1qr = qrmatrix(A1,b1);
        cout<< "questo è l'output di qrmatrix nel main, nonchè soluzione del sistema lineare con fattorizzazione qr \n"<< X1qr << endl;
        //applico la funzione per il calcolo dell'errore relativo,in questo caso, al sistema 1, risolto con Fattorizzazione QR
        double err1qr = ErrorSolution(xcorretta, X1qr);
        cout << "l'errore relativo associato al primo risultato con qr e corrisponde a: err1qr = "<< err1qr << endl;
        

    // PUNTO 2
    //Definisco la matrice A2
	MatrixXd A2(2,2);
	A2 << 5.547001962252291e-01, -5.540607316466765e-01, 8.320502943378437e-01, -8.324762492991313e-01;
	
	//Definisco il vettore b2
	VectorXd b2(2);
	b2 << -6.394645785530173e-04, 4.259549612877223e-04;

    // se la matrice è non singolare applico la fattorizzazione PA=LU e QR
    double r2 = A2.determinant();
    if(abs(r2) <1.0e-12 )
		cout << "A2 is singular" << endl;
	else 
		cout << "A2 is not singular" << endl;

        //applico la funzione per la soluzione del sistema lineare ottenuto tramite la fattorizzazione QR
        VectorXd X2palu = palumatrix(A2,b2);
        cout<< "questo è l'output di palumatrix nel main, nonchè soluzione del secondo sistema lineare con fattorizzazione PALU \n"<< X2palu << endl;
        
          
        
        //applico la funzione per il calcolo dell'errore relativo,in questo caso, al sistema 2, risolto con Fattorizzazione PALU
        double err2palu= ErrorSolution(xcorretta, X2palu);
        cout << "l'errore relativo associato al secondo sistema risolto con palu e corrisponde a: err2palu = "<< err2palu << endl;
        
        //applico la funzione per la soluzione del sistema lineare ottenuto tramite la fattorizzazione QR
        VectorXd X2qr = qrmatrix(A2,b2);
        cout<< "questo è l'output di qrmatrix nel main, nonchè soluzione del secondo sistema lineare con fattorizzazione qr \n"<< X2qr << endl;
        //applico la funzione per il calcolo dell'errore relativo,in questo caso, al sistema 2, risolto con Fattorizzazione QR
        double err2qr = ErrorSolution(xcorretta, X2qr);
        cout << "l'errore relativo associato al secondo sistema risolto con qr e corrisponde a: err2qr = "<< err2qr << endl;

    //PUNTO3

    //Definisco la matrice A3
	MatrixXd A3(2,2);
	A3 << 5.547001962252291e-01, -5.547001955851905e-01, 8.320502943378437e-01, -8.320502947645361e-01;
	
	//Definisco il vettore b3
	VectorXd b3(2);
	b3 << -6.400391328043042e-10, 4.266924591433963e-10;

    // se la matrice è non singolare applico la fattorizzazione PA=LU e QR
    double r3 = A3.determinant();
    if(abs(r3) <1.0e-12 )
		cout << "A3 is singular" << endl;
	else 
		cout << "A3 is not singular" << endl;

        //applico la funzione per la soluzione del sistema lineare ottenuto tramite la fattorizzazione QR
        VectorXd X3palu = palumatrix(A3,b3);
        cout<< "questo è l'output di palumatrix nel main, nonchè soluzione del secondo sistema lineare con fattorizzazione PALU \n"<< X3palu << endl;
        
          
        
        //applico la funzione per il calcolo dell'errore relativo,in questo caso, al sistema 3, risolto con Fattorizzazione PALU
        double err3palu= ErrorSolution(xcorretta, X3palu);
        cout << "l'errore relativo associato al terzo sistema risolto con palu e corrisponde a: err3palu = "<< err3palu << endl;
        
        //applico la funzione per la soluzione del sistema lineare ottenuto tramite la fattorizzazione QR
        VectorXd X3qr = qrmatrix(A3,b3);
        cout<< "questo è l'output di qrmatrix nel main, nonchè soluzione del terzo sistema lineare con fattorizzazione qr \n"<< X3qr << endl;
        //applico la funzione per il calcolo dell'errore relativo,in questo caso, al sistema 3, risolto con Fattorizzazione QR
        double err3qr = ErrorSolution(xcorretta, X3qr);
        cout << "l'errore relativo associato al terzo sistema risolto con qr e corrisponde a: err3qr = "<< err3qr << endl;


        


    return 0;
}
