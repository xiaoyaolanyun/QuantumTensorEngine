#include "TensorNode.h"
#include <algorithm>
#include <omp.h>
int ComplexTensor::getRank() const
{
    return m_rank;
}

qcomplex_data_t ComplexTensor::getElem(size_t num)
{
    try
    {
        return m_tensor[num];
    }
    catch (const std::exception&e)
    {
        throw e;
    }
}


static bool isPerfectSquare(int number)
{
    for (int i = 1; number > 0; i += 2)
    {
        number -= i;
    }
    return  0 == number;
}


ComplexTensor matrixMultiplication(const ComplexTensor & tensor_left, 
                                const ComplexTensor & tensor_right)
{

    auto matrix_left = tensor_left.m_tensor;
    auto matrix_right = tensor_right.m_tensor;

    if ((matrix_left.size() != matrix_right.size())  // insure dimension of the two matrixes is same
        || (!isPerfectSquare((int)matrix_left.size())))
    {
        throw exception();
    }

    int size = (int)matrix_left.size();
    qstate_t matrix_result(size, 0);
    int dimension = (int)sqrt(size);

    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            complex<double> temp = 0;
            for (int k = 0; k < dimension; k++)
            {
                temp += matrix_left[i*dimension + k] * matrix_right[k*dimension + j];
            }
            matrix_result[i*dimension + j] = temp;
        }
    }

    ComplexTensor temp(dimension, matrix_result);
    return temp;
}


void ComplexTensor::dimIncrement()
{
	m_rank = m_rank + 1;
	auto size = m_tensor.size();

	for (long long i = 0; i < size; i++)
	{
		m_tensor.push_back(m_tensor[i]);
	}
}

void ComplexTensor::getSubTensor(size_t num, int value)
{
    if (num > m_rank)
    {
        throw exception();
    }

	auto size = m_tensor.size();
    auto sub = m_rank - num;
    qsize_t step = 1ull << sub;
	size_t j;
#pragma omp parallel for private(j)
    for (long long i = 0; i < size; i += step * 2)
    {
        for ( j = i; j < i + step; j++)
        {

            if (0 == value)
            {
                m_tensor[j] = m_tensor[j];
            }
            else if (1 == value)
            {
				m_tensor[j] = m_tensor[j+step];
            }
            else
            {
                throw exception();
            }
            
        }
    }
	m_rank = this->m_rank - 1;
	m_tensor.resize(1ull << m_rank);

}


void ComplexTensor::dimDecrement(size_t num)
{
    if ((num > m_rank)||(m_rank == 0))
    {
        throw exception();
    }
	qsize_t j;
	auto size = m_tensor.size();
    auto sub = m_rank - num;
    qsize_t step = 1ull << sub;
#pragma omp parallel for private(j)
    for (long long i = 0; i < size; i += step * 2)
    {
        for (j = i; j < i + step; j++)
        {
            m_tensor[j] = (m_tensor[j] + m_tensor[j+ step]);
        }
    }
	m_rank = this->m_rank - 1;
	m_tensor.resize(1ull << m_rank, 0);
}

void ComplexTensor::swap(qsize_t src, qsize_t des)
{
    if ((src > m_rank) || (des > m_rank ))
    {
        throw exception();
    }

    auto size = m_tensor.size();
    qsize_t step_src = 1ull << (m_rank - src);
    qsize_t step_des = 1ull << (m_rank - des);
    qsize_t step_temp = 0;

    if (step_src < step_des)
    {
        step_temp = step_src;
        step_src = step_des;
        step_des = step_temp;
    }

    step_temp = step_src - step_des;

    qcomplex_data_t a, b;
	size_t j, k;
#pragma omp parallel for private(j,k,a,b)
    for (long long i = 0; i < size; i = i + 2 * step_src)
    {
        for (j = i + step_des; j < i + step_src; j = j + 2 * step_des)
        {
            for (k = j; k < j + step_des; k++)
            {
                a = m_tensor[k];
                b = m_tensor[k + step_temp];

                m_tensor[k] = b;
                m_tensor[k + step_temp] = a;
            }
        }
    }
    
}

ComplexTensor & ComplexTensor::operator*(ComplexTensor & old)
{
    if (this->m_rank != old.m_rank)
    {
        throw exception();
    }
	qcomplex_data_t complex_temp;
#pragma omp parallel for private(complex_temp)
    for (long long i = 0; i < this->m_tensor.size(); i++)
    {
        complex_temp = this->m_tensor[i] * old.m_tensor[i];
        m_tensor[i] = complex_temp;
    }

	return *this;
}

ComplexTensor ComplexTensor::operator+(ComplexTensor & old)
{
    if (this->m_rank != old.m_rank)
    {
        throw exception();
    }

    ComplexTensor temp;
    for (size_t i = 0; i < this->m_tensor.size(); i++)
    {
        auto complex_temp = this->m_tensor[i] + old.m_tensor[i];
        temp.m_tensor.push_back(complex_temp);
    }
    temp.m_rank = old.m_rank;

    return temp;
}

ComplexTensor & ComplexTensor::operator=(ComplexTensor & old)
{
    m_rank = old.m_rank;
    m_tensor.resize(0);
    for (auto  iter : old.m_tensor)
    {
        m_tensor.push_back(iter);
    }
    return *this;
}

qsize_t Edge::getQubitCount() const noexcept
{
    return m_qubit_count;
}

void Edge::premultiplication(Edge & src_edge)
{
    ComplexTensor left = src_edge.getComplexTensor();
    ComplexTensor right = this->m_tensor;
    m_tensor = matrixMultiplication(left , right);

    auto size = m_contect_vertice.size();
    for (auto iter : src_edge.m_contect_vertice)
    {
        auto find_result = std::find(m_contect_vertice.begin(),
            m_contect_vertice.end(), iter);
        if (find_result == m_contect_vertice.end())
        {
            m_contect_vertice.push_back(iter);
        }
        else
        {
            m_contect_vertice.erase(find_result);
        }
    }
}

bool Edge::mergeEdge(Edge & edge)
{
    dimIncrementByEdge(edge);
    edge.dimIncrementByEdge(*this);
    edge.swapByEdge(*this);
	m_tensor = m_tensor * edge.m_tensor;
    return true;
}

void Edge::dimDecrementbyValue(qsize_t qubit,qsize_t num, int value)
{
    if ((0 != value) && (1 != value))
    {
        throw exception();
    }

    int i = 1;
    for (auto iter = m_contect_vertice.begin();
        iter != m_contect_vertice.end();
        ++iter)
    {
        if (((*iter).first == qubit) &&
            ((*iter).second == num))
        {
            m_tensor.getSubTensor(i,value);
            iter = m_contect_vertice.erase(iter);
            return;
        }
        i++;
    }
}


void Edge::earseContectVertice(qsize_t qubit, size_t num)
{
    int i = 1;
    for (auto iter = m_contect_vertice.begin();
        iter != m_contect_vertice.end();
        ++iter)
    {
        if (((*iter).first == qubit) && ((*iter).second == num))
        {
            iter = m_contect_vertice.erase(iter);
            return;
        }
        i++;
    }
}

void Edge::dimDecrement(qsize_t qubit, qsize_t num)
{
    int i = 1;
    for (auto iter = m_contect_vertice.begin();
         iter != m_contect_vertice.end();
         ++iter)
    {
        if (((*iter).first == qubit)&&((*iter).second == num))
        {
            m_tensor.dimDecrement(i);
            iter = m_contect_vertice.erase(iter);
            return;
        }
        i++;
    }
}

void Edge::dimIncrementByEdge(Edge & edge)
{
    for (auto iter : edge.m_contect_vertice)
    {
        bool isTrue = false;
        for (auto this_iter : m_contect_vertice)
        {
            if ((iter.first == this_iter.first) &&
                (iter.second == this_iter.second))
            {
                isTrue = true;
                break;
            }
        }

        if (!isTrue)
        {
            m_contect_vertice.insert(m_contect_vertice.begin(),iter);
            m_tensor.dimIncrement();
        }
    }
}

void Edge::swapByEdge(Edge & edge)
{
    int i = 0;
    auto size = m_contect_vertice.size();
    for (auto iter : edge.m_contect_vertice)
    {
        if ((iter.first != m_contect_vertice[i].first) ||
            (iter.second != m_contect_vertice[i].second))
        {
            for (size_t j = i+1; j < size; j++)

            {
                if ((iter.first == m_contect_vertice[j].first) &&
                    (iter.second == m_contect_vertice[j].second))
                {
                    m_tensor.swap(i+1,j+1);
                    std::swap(m_contect_vertice[i], m_contect_vertice[j]);
                    break;
                }
            }
        }
        i++;
    }
}

int Edge::getRank() const noexcept
{
    return m_tensor.getRank();
}

ComplexTensor Edge::getComplexTensor() const noexcept
{
    return m_tensor;
}

qcomplex_data_t Edge::getElem(VerticeMatrix &vertice)
{
    auto vertice_vector = getContectVertice();
    auto size = vertice_vector.size();
	if (m_tensor.getRank() == 0)
	{
		return m_tensor.getElem(0);
	}
    size_t M = 0;
    for (size_t i = 0; i < size; i++)
    {
        M += vertice.getVerticeValue(vertice_vector[i].first,
            vertice_vector[i].second) <<
            (size - i - 1);
    }
    auto result = m_tensor.getElem(M);
    return result;
}

void Edge::setComplexTensor( ComplexTensor &tensor) noexcept
{
    m_tensor = tensor;
}

vector<pair<qsize_t, qsize_t>> Edge::getContectVertice() const noexcept
{
    return m_contect_vertice;
}

void Edge::setContectVerticeVector
           (const qubit_vector_t & contect_vertice) noexcept
{
    m_contect_vertice.resize(0);
    for (auto aiter : contect_vertice)
    {
        m_contect_vertice.push_back(aiter);
    }
}

void Edge::setContectVertice(qsize_t qubit, qsize_t src_num, qsize_t des_num)
{
    try
    {
        for (size_t i = 0; i < m_contect_vertice.size(); i++)
        {
            if ((m_contect_vertice[i].first == qubit) && 
                 (m_contect_vertice[i].second == src_num))
            {
                m_contect_vertice[i].second = des_num;
            }
        }
    }
    catch (const std::exception&e)
    {
        throw e;
    }
}


vector<qsize_t> Vertice::getContectEdge() const noexcept
{
    return m_contect_edge;
}

void Vertice::addContectEdge(qsize_t edge_id) noexcept
{
    m_contect_edge.push_back(edge_id);
}

void Vertice::setContectEdge(qsize_t id, qsize_t value)
{
    try
    {
        for (auto i = m_contect_edge.begin();i != m_contect_edge.end(); ++i)
        {
            if (*i == id)
            {
                (*i) = value;
            }
        }
    }
    catch (const std::exception &e)
    {
        throw e;
    }

}

void Vertice::setContectEdgebyID(qsize_t id, qsize_t value)
{
    try
    {
        m_contect_edge[id] = value;
    }
    catch (const std::exception &e)
    {
        throw e;
    }
}

void Vertice::deleteContectEdge(qsize_t edge_num)
{
    for (auto iter = m_contect_edge.begin();
        iter != m_contect_edge.end();
        ++iter)
    {
        if (edge_num == *iter)
        {
            m_contect_edge.erase(iter);
            return;
        }
    }
}

int Vertice::getValue() const noexcept
{
    return m_value;
}

void Vertice::setValue(int value)noexcept
{
    m_value = value;
}

VerticeMatrix::VerticeMatrix():m_qubit_count(0), m_vertice_count(0)
{

}

VerticeMatrix::VerticeMatrix(const VerticeMatrix & old)
{
    for (auto vector_iter : old.m_vertice_matrix)
    {
        vertice_map_t temp;
        for (auto map_iter : vector_iter)
        {
            temp.insert(map_iter);
        }
        m_vertice_matrix.push_back(temp);
    }

    m_qubit_count = old.m_qubit_count;
    m_vertice_count = old.m_vertice_count;
}

VerticeMatrix VerticeMatrix::operator=(const VerticeMatrix &old)
{
    for (auto vector_iter : old.m_vertice_matrix)
    {
        vertice_map_t temp;
        for (auto map_iter : vector_iter)
        {
            temp.insert(map_iter);
        }
        m_vertice_matrix.push_back(temp);
    }
    m_qubit_count = old.m_qubit_count;
    m_vertice_count = old.m_vertice_count;
    return *this;
}

qsize_t VerticeMatrix::getQubitCount() const
{
    return m_qubit_count;
}

qsize_t VerticeMatrix::getVerticeCount() const
{
    return m_vertice_count;
}

void VerticeMatrix::subVerticeCount()
{
    if (m_vertice_count > 0)
    {
        m_vertice_count--;
    }
}

void VerticeMatrix::initVerticeMatrix(qsize_t qubit_num)
{
    m_qubit_count = qubit_num;
    m_vertice_count = qubit_num;
    for (qsize_t i = 0; i < qubit_num; i++)
    {
        Vertice temp;
        temp.setValue(0);
        map<qsize_t, Vertice> vertice_vector;
        vertice_vector.insert(pair<qsize_t, Vertice>(0, temp));
        m_vertice_matrix.push_back(vertice_vector);
    }
    //m_vertice_matrix[0][0].setValue(1);
}

map<qsize_t,Vertice>::iterator VerticeMatrix::deleteVertice(qsize_t qubit, 
                                                            qsize_t num)
{
    try
    {
        auto aiter = m_vertice_matrix[qubit].find(num);
        return m_vertice_matrix[qubit].erase(aiter);
    }
    catch (const std::exception& e)
    {
        throw e;
    }
    m_vertice_count--;
}


qsize_t VerticeMatrix::addVertice(qsize_t qubit)
{
    try
    {
        Vertice temp;
        qsize_t size = m_vertice_matrix[qubit].size();
        m_vertice_matrix[qubit].insert(pair<qsize_t,Vertice>(size,temp));
        m_vertice_count++;
        return size;
    }
    catch (const std::exception& e)
    {
        throw e;
    }
}

qsize_t VerticeMatrix::addVertice(qsize_t qubit, qsize_t num)
{
    Vertice temp;
    return addVertice(qubit,num,temp);
}


qsize_t VerticeMatrix::addVertice(qsize_t qubit, qsize_t num,Vertice & vertice)
{
    try
    {
        auto is_true = m_vertice_matrix[qubit].
                       insert(pair<qsize_t, Vertice>(num, vertice));

        if (!is_true.second)
        {
            m_vertice_matrix[qubit][num] = vertice;
        }
        else
        {
            m_vertice_count++;
        }

        return num;
    }
    catch (const std::exception& e)
    {
        throw e;
    }
}

int VerticeMatrix::getVerticeValue(qsize_t qubit, qsize_t num)
{
    try
    {
        auto iter = m_vertice_matrix[qubit].find(num);
        return (*iter).second.getValue();
    }
    catch (const std::exception& e)
    {
        throw e;
    }
}

qsize_t VerticeMatrix::getEmptyVertice()
{
    qsize_t result = 0;
    for (auto matrix_iter : m_vertice_matrix)
    {
        for (auto map_iter : matrix_iter)
        {
            auto vertice = map_iter.second.getContectEdge();
            if (0 == vertice.size())
            {
                result++;
            }
        }
    }
    return result;
}

void VerticeMatrix::setVerticeValue(qsize_t qubit, qsize_t num, int value)
{
    if (value > 1 && value < -1)
    {
        throw exception();
    }

    try
    {
        return m_vertice_matrix[qubit][num].setValue(value);
    }
    catch (const std::exception& e)
    {
        throw e;
    }
}

qsize_t VerticeMatrix::getQubitVerticeLastID(qsize_t qubit)
{
    try
    {
        auto iter = m_vertice_matrix[qubit].end();
        iter--;
        return (*iter).first;
    }
    catch (const std::exception&e)
    {
        throw e;
    }
}

vector<qsize_t> VerticeMatrix::getContectEdge(qsize_t qubit, 
                                              qsize_t num) const
{
    try
    {
        auto iter = m_vertice_matrix[qubit].find(num);
        return (*iter).second.getContectEdge();
    }
    catch (const std::exception&e)
    {
        throw e;
    }

}

vertice_matrix_t::iterator  VerticeMatrix::begin() noexcept
{
    return  m_vertice_matrix.begin();
}
vertice_matrix_t::iterator VerticeMatrix::end()noexcept
{
    return  m_vertice_matrix.end();
}

vertice_matrix_t::iterator VerticeMatrix::getQubitMapIter(qsize_t qubit) noexcept
{
    return m_vertice_matrix.begin()+qubit;
}

vertice_map_t::iterator VerticeMatrix::getQubitMapIterBegin(qsize_t qubit)
{
    try
    {
        return m_vertice_matrix[qubit].begin();
        
    }
    catch (const std::exception &e)
    {
        throw e;
    }
}

vertice_map_t::iterator VerticeMatrix::getQubitMapIterEnd(qsize_t qubit)
{
    try
    {
        return m_vertice_matrix[qubit].end();
    }
    catch (const std::exception &e)
    {
        throw e;
    }
}

vertice_map_t::iterator  VerticeMatrix::getVertice(qsize_t qubit, qsize_t num)
{
    try
    {
        return m_vertice_matrix[qubit].find(num);
    }
    catch (const std::exception&e)
    {
        throw e;
    }
}

void VerticeMatrix::addContectEdge(qsize_t qubit, qsize_t vertice_num, qsize_t edge_id)
{
    try
    {
        auto contect_edge = m_vertice_matrix[qubit][vertice_num].getContectEdge();
        for (auto edge : contect_edge)
        {
            if (edge_id == edge)
            {
                return;
            }
        }
        m_vertice_matrix[qubit][vertice_num].addContectEdge(edge_id);
    }
    catch (const std::exception&e)
    {
        throw e;
    }
}

void VerticeMatrix::changeContectEdge(qsize_t qubit,
                                      qsize_t vertice_num,
                                      qsize_t contect_edge_vector_num,
                                      qsize_t contect_edge)
{
    try
    {
        m_vertice_matrix[qubit][vertice_num].
            setContectEdge(contect_edge_vector_num,
                contect_edge);
    }
    catch (const std::exception&e)
    {
        throw e;
    }
}

void VerticeMatrix::deleteContectEdge(qsize_t qubit,
                                      qsize_t vertice_num,
                                      qsize_t contect_edge_vector_num )
{
    m_vertice_matrix[qubit][vertice_num].
        deleteContectEdge(contect_edge_vector_num);
}

void VerticeMatrix::clearVertice() noexcept
{
    for (auto iter = m_vertice_matrix.begin();
        iter != m_vertice_matrix.end();
        ++iter)
    {
        for (auto map_iter = (*iter).begin();
            map_iter != (*iter).end();
            ++map_iter)
        {
            (*map_iter).second.setValue(-1);
        }
    }
}

VerticeMatrix::~VerticeMatrix()
{
}

qsize_t QuantumProgMap::getVerticeCount() const
{
    return m_vertice_matrix->getVerticeCount();
}


