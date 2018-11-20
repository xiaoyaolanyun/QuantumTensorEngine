#ifndef TENSOR_ENGINE_H
#define TENSOR_ENGINE_H
#include "TensorNode.h"
#include <thread>


typedef struct QubitVertice
{
    qsize_t m_qubit_id;
    qsize_t m_num;
    qsize_t m_max;
    std::thread m_t;
    QubitVertice(const QubitVertice & old)
    {
        m_num = old.m_num;
        m_qubit_id = old.m_qubit_id;
        m_max = old.m_max;
    }

    QubitVertice operator = (const QubitVertice & old)
    {
        m_num = old.m_num;
        m_qubit_id = old.m_qubit_id;
        m_max = old.m_max;
        return *this;
    }
    QubitVertice() :m_num(0), m_qubit_id(0), m_max(0)
    {}
} qubit_vertice_t;


class TensorEngine
{
public:
    static qubit_vertice_t getNoValueVertice(QuantumProgMap & prog_map);
    static qubit_vertice_t getNoValueAndContectEdgeMaxVertice
                                (QuantumProgMap & prog_map);
    static void Merge(QuantumProgMap & prog_map);
    /*
    qcomplex_data_t split(QuantumProgMap & prog_map,
                          qubit_vertice_t *);
                          */
    static qcomplex_data_t Merge(QuantumProgMap & prog_map,
                          qubit_vertice_t *);
    static qcomplex_data_t computing(QuantumProgMap & prog_map);
    /*
    qsize_t getSubQuantumProgMap(QuantumProgMap &,
                                 QuantumProgMap &, 
                                 qubit_vertice_t &,
                                 int);
    */
    static void MergeQuantumProgMap(QuantumProgMap &,
                                qubit_vertice_t &);
    
    static void dimDecrementbyValue(QuantumProgMap &,
        qubit_vertice_t &,int value);
};

#endif // !TENSOR_ENGINE_H


