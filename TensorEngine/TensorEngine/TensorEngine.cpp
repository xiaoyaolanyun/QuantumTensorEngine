#include "TensorEngine.h"
#include <algorithm>
static
void mergeVerticeAndEdge(QuantumProgMap * prog_map,qsize_t qubit)
{
    auto vertice = prog_map->getVerticeMatrix();
    auto vertice_map_iter = vertice->getQubitMapIter(qubit);
    auto i = vertice->getQubitMapIterBegin(qubit);
    while (i != vertice->getQubitMapIterEnd(qubit))
    {
        auto contect_edge = vertice->getContectEdge(qubit, (*i).first);
        auto contect_count = contect_edge.size();
        auto edge_map = prog_map->getEdgeMap();

        if (2 == contect_count)
        {
            auto edge_iter_second = edge_map->find(contect_edge[1]);
            auto edge_second = (*edge_iter_second).second;
            auto edge_iter = edge_map->find(contect_edge[0]);
            auto edge_first = (*edge_iter).second;

            if ((edge_second.getRank() != 2) ||
                (edge_first.getRank() != 2))
            {
                ++i;
                continue;
            }
            if(((*edge_iter).second.getQubitCount() != 1) ||
                ((*edge_iter_second).second.getQubitCount() != 1))
            {
                ++i;
                continue;
            }
            (*edge_iter).second.premultiplication((*edge_iter_second).second);
            /*
            (*edge_iter).second.mergeEdge((*edge_iter_second).second);
            (*edge_iter).second.setContectVerticeVector(
                (*edge_iter_second).second.getContectVertice());
            (*edge_iter).second.dimDecrement(qubit, (*i).first);
            */
            edge_map->erase(edge_iter_second);
            i = (*vertice_map_iter).erase(i);
            vertice->subVerticeCount();
            (*i).second.setContectEdgebyID(0, (*edge_iter).first);
        }
        else
        {
            ++i;
        }
    }
}

void getNoValueMaxRankVertice(QuantumProgMap * prog_map,
                              qubit_vertice_t * qubit_vertice)
{
    if (nullptr == prog_map)
    {
        throw exception();
    }
    auto vertice = prog_map->getVerticeMatrix();
    auto vertice_matrix_iter = 
                        vertice ->getQubitMapIter(qubit_vertice->m_qubit_id);
    auto edge = prog_map->getEdgeMap();
    size_t min = 12345678;
    size_t num = 0;
    qsize_t target_id = 0;
    for (auto aiter : (*vertice_matrix_iter))
    {
        size_t temp = aiter.second.getContectEdge().size();
        bool is_true = true;
        if (min > temp)
            if (aiter.second.getValue() < 0)
            {
				min = temp;
				target_id = aiter.first;
            }
    }

    qubit_vertice->m_num = target_id;
    qubit_vertice->m_max = min;
}

qubit_vertice_t TensorEngine::getNoValueVertice(QuantumProgMap & prog_map)
{
    qubit_vertice_t qubit_vertice;
    qubit_vertice.m_max = 0;
    qubit_vertice.m_num = 0;
    qubit_vertice.m_qubit_id = 0;
    auto vertice = prog_map.getVerticeMatrix();
    int i = 0;
    for (auto iter = vertice->begin(); iter != vertice->end(); iter++)
    {
        for (auto map_iter = (*iter).begin(); map_iter != (*iter).end(); ++map_iter)
        {
            if (-1 == (*map_iter).second.getValue())
            {
                qubit_vertice.m_num = (*map_iter).first;
                qubit_vertice.m_qubit_id = i;
            }
        }
        i++;
    }

    return qubit_vertice;
}

qubit_vertice_t TensorEngine::getNoValueAndContectEdgeMaxVertice
                                (QuantumProgMap & prog_map)
{
    auto vertice = prog_map.getVerticeMatrix();

    vector<qubit_vertice_t> qubit_vertice;
    qubit_vertice.resize(vertice->getQubitCount());

    qsize_t i = 0;
    for (auto iter = vertice->begin(); iter != vertice->end(); iter++)
    {
        qubit_vertice[i].m_qubit_id = i;
        qubit_vertice[i].m_t = std::thread (getNoValueMaxRankVertice,
                                            &prog_map,&qubit_vertice[i]);
        i++;
    }

    qubit_vertice_t temp;
    temp.m_qubit_id = 0;
    temp.m_max = 1234567;
    temp.m_num = 0;
    for (auto aiter = qubit_vertice.begin();
         aiter != qubit_vertice.end();
         ++aiter)
    {
        (*aiter).m_t.join();

        if ((*aiter).m_max < temp.m_max)
        {
            temp.m_num = (*aiter).m_num;
            temp.m_qubit_id = (*aiter).m_qubit_id;
            temp.m_max = (*aiter).m_max;
        }
    }
    return temp;
} 

void TensorEngine::Merge(QuantumProgMap & prog_map)
{
    auto vertice = prog_map.getVerticeMatrix();
    auto count = vertice->getQubitCount();
    vector<std::thread> thread_vector;
    thread_vector.resize(count);

    for (auto i = 0; i < count; i++)
    {
        thread_vector[i] = std::thread(mergeVerticeAndEdge,&prog_map,i);
    }

    for (auto i = 0; i < count; i++)
    {
        thread_vector[i].join();
    }
}

static void split(QuantumProgMap * prog_map,
                             qubit_vertice_t * qubit_vertice,
                             qcomplex_data_t * result)
{
    qubit_vertice_t  temp;
    if ((nullptr == prog_map) || (nullptr == result))
    {
        throw exception();
    }
    if (nullptr == qubit_vertice)
    {
        temp = TensorEngine::getNoValueVertice(*prog_map);
        split(prog_map, &temp,result);
    }
    else
    {
        if ((0 == qubit_vertice->m_qubit_id)&&(0 == qubit_vertice->m_num))
        {
            (*result) = TensorEngine::computing(*prog_map);
        }
        else
        {

            QuantumProgMap *new_map = new QuantumProgMap(*prog_map);
 
            TensorEngine::dimDecrementbyValue(*prog_map, *qubit_vertice, 0);
            temp = TensorEngine::getNoValueVertice(*prog_map);
            qcomplex_data_t result_zero(0);
            std::thread thread = std::thread(split, prog_map, &temp, &result_zero);

            TensorEngine::dimDecrementbyValue(*new_map, *qubit_vertice, 1);
            qcomplex_data_t result_one(0);
            split(new_map, &temp,&result_one);
            thread.join();
            delete new_map;
            *result = result_one + result_zero;
        }
    }
}

qcomplex_data_t TensorEngine::Merge(QuantumProgMap & prog_map,
                                    qubit_vertice_t *qubit_vertice)
{
    if(nullptr == qubit_vertice)
    {
        qubit_vertice_t temp = getNoValueAndContectEdgeMaxVertice(prog_map);
        return Merge(prog_map, &temp);
    }
    else
    {
        if ((qubit_vertice->m_qubit_id != 0) || (qubit_vertice->m_num != 0))
        {
			MergeQuantumProgMap(prog_map,
				*qubit_vertice);
			qubit_vertice_t temp = getNoValueAndContectEdgeMaxVertice(prog_map);
			auto result = Merge(prog_map, &temp);
			return result;
        }
        else
        {
            qcomplex_data_t result(0);
			result = TensorEngine::computing(prog_map);
            return result;
        }
       

    }
}

qcomplex_data_t TensorEngine::computing(QuantumProgMap & prog_map)
{
    auto edge_map = prog_map.getEdgeMap();
    qcomplex_data_t result = 1;

    for (auto iter = edge_map->begin(); iter != edge_map->end(); ++iter)
    {
        result *= (*iter).second.getElem(*prog_map.getVerticeMatrix());
    }
    return result;
}

#include <iostream>
void TensorEngine::MergeQuantumProgMap(QuantumProgMap & prog_map,
                                          qubit_vertice_t & qubit_vertice)
{
    auto vertice = prog_map.getVerticeMatrix();
    auto edge_map = prog_map.getEdgeMap();

    auto contect_edge = vertice->getContectEdge(qubit_vertice.m_qubit_id,
                                                qubit_vertice.m_num);
    try
    {
        qsize_t i = 0;
        auto edge_iter = edge_map->find((*contect_edge.begin()));
        vector<pair<qsize_t, qsize_t>> vertice_vector;
        
        auto first_edge = edge_map->find(contect_edge[0]);
        for (size_t i = 1; i < contect_edge.size(); i++)
        {
            auto edge = edge_map->find(contect_edge[i]);
            (*first_edge).second.mergeEdge((*edge).second);
        }
        (*first_edge).second.dimDecrement(qubit_vertice.m_qubit_id,
                                          qubit_vertice.m_num);
        for (auto contect_edge_iter : contect_edge)
        {
            auto iter = edge_map->find(contect_edge_iter);
            auto contect_vertice = (*iter).second.getContectVertice();

            for (auto contect_vertice_iter : contect_vertice)
            {
                if ((contect_vertice_iter.first != qubit_vertice.m_qubit_id) ||
                    (contect_vertice_iter.second != qubit_vertice.m_num))
                {
                    vertice->deleteContectEdge(contect_vertice_iter.first,
                        contect_vertice_iter.second,
                        contect_edge_iter);
                    vertice->addContectEdge(contect_vertice_iter.first,
                        contect_vertice_iter.second,
                        (*edge_iter).first);
                }
            }
            if (0 != i)
            {
                edge_map->erase(iter);
            }
            i++;
        }
        vertice->deleteVertice(qubit_vertice.m_qubit_id, qubit_vertice.m_num);
    }
    catch (const std::exception&e)
    {
        throw e;
    }
}

void TensorEngine::dimDecrementbyValue(QuantumProgMap & prog_map,
    qubit_vertice_t & qubit_vertice,int value)
{
    auto vertice = prog_map.getVerticeMatrix();
    auto edge_map = prog_map.getEdgeMap();
    auto contect_edge = vertice->getContectEdge(qubit_vertice.m_qubit_id,
                                                qubit_vertice.m_num);

    for (auto iter : contect_edge)
    {
        auto find_iter =edge_map->find(iter);
        if (find_iter != edge_map->end())
        {
            (*find_iter).second.dimDecrementbyValue(qubit_vertice.m_qubit_id,
                                                    qubit_vertice.m_num,
                                                    value);
        }
    }

    vertice->deleteVertice(qubit_vertice.m_qubit_id,
                           qubit_vertice.m_num);
}




