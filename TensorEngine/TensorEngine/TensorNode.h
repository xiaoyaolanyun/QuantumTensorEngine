#ifndef TENSOR_NODE_H
#define TENSOR_NODE_H
#include <complex>
#include <vector>
#include <map>
#include <exception>
using std::map;
using std::complex;
using std::vector;
using std::exception;
using std::pair;

typedef double qdata_t;
typedef size_t qsize_t;
typedef size_t qsize_t;
typedef complex<qdata_t> qcomplex_data_t;
typedef vector<qcomplex_data_t> qstate_t;
typedef vector<pair<qsize_t, qsize_t>> qubit_vector_t;
class Edge;
typedef map<size_t, Edge> EdgeMap;



class ComplexTensor
{
public:
    ComplexTensor();
    ~ComplexTensor();

    friend class ComplexTensor;
    int getRank() const;
    qcomplex_data_t getElem(size_t num);
    friend ComplexTensor matrixMultiplication(const ComplexTensor &matrix_left,
        const ComplexTensor &matrix_right);
    void dimIncrement() ;
    void getSubTensor(size_t num,int value);
    void dimDecrement(size_t num) ;
    void swap(qsize_t, qsize_t);
	ComplexTensor & operator * (ComplexTensor &);
    ComplexTensor operator + (ComplexTensor &);
    ComplexTensor& operator = (ComplexTensor &);
    ComplexTensor(const ComplexTensor& old) 
        :m_rank(old.m_rank),
         m_tensor(old.m_tensor)
    {
    }

    ComplexTensor(int rank,qstate_t tensor)
        :m_rank(rank),
         m_tensor(tensor)
    {
    }

private:
    int m_rank;
    qstate_t m_tensor;
};

inline ComplexTensor::ComplexTensor():m_rank(0)
{
}

inline ComplexTensor::~ComplexTensor()
{
}

class VerticeMatrix;
class QuantumProgMap;
class Edge
{
public:
    friend class Edge;
    inline Edge(qsize_t qubit_count,
                ComplexTensor &tensor,
                vector<pair<qsize_t, qsize_t>> &contect_vertice) noexcept :
                m_tensor(tensor),
                m_contect_vertice(contect_vertice),
                m_qubit_count(qubit_count)
    {
    }
    ~Edge() {};

    inline Edge(const Edge & old) :
        m_tensor(old.getComplexTensor()),
        m_qubit_count(old.m_qubit_count),
        m_contect_vertice(old.m_contect_vertice)
    {

    }

    void earseContectVertice(qsize_t qubit, size_t num);
    qsize_t getQubitCount()const noexcept;
    void premultiplication(Edge &);
    bool mergeEdge(Edge & );
    void dimDecrementbyValue(qsize_t ,qsize_t,
                             int value);
    void dimDecrement(qsize_t qubit, qsize_t num);
    void dimIncrementByEdge(Edge &);
    void swapByEdge(Edge &);
    int getRank() const noexcept;
    ComplexTensor getComplexTensor() const noexcept;
    qcomplex_data_t getElem(VerticeMatrix &vertice);
    void setComplexTensor( ComplexTensor &) noexcept;
    vector<pair<qsize_t, qsize_t>> getContectVertice() const noexcept;
    void setContectVerticeVector(const qubit_vector_t &) noexcept;
    void setContectVertice(qsize_t qubit, qsize_t src_num, qsize_t des_num);
private:
    qsize_t m_qubit_count;
    ComplexTensor m_tensor;
    qubit_vector_t m_contect_vertice;
};



class Vertice
{
public:
    inline Vertice(int value, vector<qsize_t> &contect_edge) noexcept 
                   :m_contect_edge(contect_edge), m_value(value)
    {}

    inline Vertice() : m_value(-1)
    {}
    ~Vertice()
    {
    }
    
    inline Vertice operator = (const Vertice & old)
    {
        m_value = old.getValue();
        auto temp = old.getContectEdge();
        m_contect_edge.resize(0);
        for (auto aiter : temp)
        {
            m_contect_edge.push_back(aiter);
        }

        return *this;
    }
    
    inline Vertice(const Vertice & old) 
    {
        m_value = old.getValue();
        auto temp = old.getContectEdge();
        m_contect_edge.resize(0);
        for (auto aiter : temp)
        {
            m_contect_edge.push_back(aiter);
        }
    }

    vector<qsize_t> getContectEdge()const noexcept;
    void addContectEdge(qsize_t) noexcept;
    void setContectEdge(qsize_t, qsize_t);
    void setContectEdgebyID(qsize_t id, qsize_t value);
    void deleteContectEdge(qsize_t);
    int getValue() const noexcept;
    void setValue(int) noexcept;

private:
    vector<qsize_t> m_contect_edge;
    int m_value;
};


typedef map<qsize_t, Vertice> vertice_map_t;
typedef vector<vertice_map_t> vertice_matrix_t;
class VerticeMatrix
{
public:
    friend class VerticeMatrix;
    VerticeMatrix();
    VerticeMatrix(const VerticeMatrix &);
    VerticeMatrix operator = (const VerticeMatrix &);

    qsize_t getQubitCount() const;
    qsize_t getVerticeCount() const;
    void subVerticeCount();
    qsize_t addVertice(qsize_t);
    qsize_t addVertice(qsize_t, qsize_t);
    qsize_t addVertice(qsize_t,qsize_t, Vertice &);
    int getVerticeValue(qsize_t, qsize_t);
    qsize_t getEmptyVertice();
    void setVerticeValue(qsize_t, qsize_t, int);
    void initVerticeMatrix(qsize_t);
    map<qsize_t, Vertice>::iterator deleteVertice(qsize_t, qsize_t);
    qsize_t getQubitVerticeLastID(qsize_t);
    vector<qsize_t> getContectEdge(qsize_t, qsize_t) const;
    vertice_matrix_t::iterator begin()noexcept;
    vertice_matrix_t::iterator end()noexcept;
    vertice_matrix_t::iterator getQubitMapIter(qsize_t qubit)noexcept;
    vertice_map_t::iterator getQubitMapIterBegin(qsize_t qubit);
    vertice_map_t::iterator getQubitMapIterEnd(qsize_t qubit);
    vertice_map_t::iterator getVertice(qsize_t, qsize_t);
    void addContectEdge(qsize_t, qsize_t, qsize_t);
    void changeContectEdge(qsize_t, qsize_t, qsize_t,qsize_t);
    void deleteContectEdge(qsize_t, qsize_t, qsize_t);
    void clearVertice() noexcept;
    ~VerticeMatrix();
private:
    qsize_t m_qubit_count;
    qsize_t m_vertice_count;
    vector<map<qsize_t, Vertice>> m_vertice_matrix;
};

class QuantumProgMap
{
private:
    VerticeMatrix * m_vertice_matrix;
    EdgeMap * m_edge_map;
public:
    friend class QuantumProgMap;
    QuantumProgMap()
    {
        m_vertice_matrix = new VerticeMatrix();
        m_edge_map = new EdgeMap;
    }
    ~QuantumProgMap()
    {
        delete m_vertice_matrix;
        delete m_edge_map;
    }

    inline QuantumProgMap(const QuantumProgMap & old)
    {
        m_vertice_matrix = new VerticeMatrix(*(old.m_vertice_matrix));
        m_edge_map = new EdgeMap(*(old.m_edge_map));
    }

    inline VerticeMatrix * getVerticeMatrix()
    {
        return m_vertice_matrix;
    }

    inline EdgeMap * getEdgeMap()
    {
        return m_edge_map;
    }

    inline void clearVerticeValue() noexcept
    {
        m_vertice_matrix->clearVertice();
    }
    qsize_t getVerticeCount() const;

};



#endif // !TENSOR_NODE_H
