#include<iostream>
#include<iomanip>
#include<fstream>
using namespace std;

void Explicit_Euler(int n);
void Implicit_Euler(int n);


int main() {
	int n = 3;
	double* U0 = new double[n];
	U0[0] = 1;	 U0[1] = 1;	 U0[2] = 1;  //задача Коши
	
	Explicit_Euler(n);
	Implicit_Euler(n);
	
	delete[] U0;
}