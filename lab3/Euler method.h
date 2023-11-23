#pragma once
#include<fstream>

using namespace std;

double foo1(double**, double*, double*, int);

double foo2(double**, double*, double*, int);

double foo3(double**, double*, double*, int);

void readFromFile(double***, istream&, int);

void readFromFile(double**, istream&, int);

void filling(double**, double*, double***, double**, double*, int);

template<typename T>
void print(T***, int);

template<typename T>
void print(T**, int);

template<typename T>
void print(T*, int);

void implictEuler(double**, double**, double*, double (*) (double**, double*, double*, int),
	double (*) (double**, double*, double*, int),
	double (*) (double**, double*, double*, int),
	double, double, double, double, int);

void explictEuler(double*, double**, double*, double(double**, double*, double*, int), double(double**, double*, double*, int)
	, double(double**, double*, double*, int), double, double, double, int);