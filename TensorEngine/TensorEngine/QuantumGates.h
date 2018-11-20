#ifndef QUANTUMGATES_H
#define QUANTUMGATES_H
#include "TensorNode.h"
#include <complex>
using std::complex;

void H(QuantumProgMap & prog_map, qsize_t qubit);

void X(QuantumProgMap & prog_map, qsize_t qubit);

void Y(QuantumProgMap & prog_map, qsize_t qubit);   

void Z(QuantumProgMap & prog_map, qsize_t qubit);

void X1(QuantumProgMap & prog_map, qsize_t qubit);

void Y1(QuantumProgMap & prog_map, qsize_t qubit);

void Z1(QuantumProgMap & prog_map, qsize_t qubit);

void RX(QuantumProgMap & prog_map, qsize_t qubit, double angle);
void RY(QuantumProgMap & prog_map, qsize_t qubit, double angle);
void RZ(QuantumProgMap & prog_map, qsize_t qubit, double angle);

void CZ(QuantumProgMap & prog_map,
    qsize_t qubit1,
    qsize_t qubit2);

void CNOT(QuantumProgMap &prog_map,
    qsize_t qubit1,
    qsize_t qubit2);

vector<pair<qsize_t, double>> PMeasure
(QuantumProgMap &prog_map,
    vector<qsize_t> & qubit_vector);


/*

class H_Gate
{
public:
void addVerticeAndEdge(QuantumProgMap &,qsize_t qubit);
};

class CZ_Gate
{
public:
void addVerticeAndEdge(QuantumProgMap &,qsize_t qubit1, qsize_t qubit2);
};

class CNOT_Gate
{
public:
void addVerticeAndEdge(QuantumProgMap &,qsize_t qubit1, qsize_t qubit2);
};

class PMeasure_Gate
{
public:
vector<pair<qsize_t, double>> getPMeasureResult(QuantumProgMap &,vector<qsize_t> &);
};
class X_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit);
};

class Y_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit);
};

class Z_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit);
};

class T_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit);
};

class S_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit);
};

class U4_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit,double,double,double,double);
};


class U1_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit,double angle);
};



class ISWAP_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit1, qsize_t qubit2);
};



*/
/*
class RX_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit,double angle);
};

class RY_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit,double angle);
};

class RZ_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit,double angle);
};

class CU_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit1,qsize_t qubit2,double angle);
};

class CNOT_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit1,qsize_t qubit2);
};



class ISWAPTheta_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit1,qsize_t qubit2, double angle);
};



class SQISWAP_Gate
{
public:
    void addVerticeAndEdge(qsize_t qubit1, qsize_t qubit2);
};
*/
#endif // !QUANTUMGATES_H
