#include "QuantumGates.h"
#include "TensorEngine.h"

#define SQRT2 1.4142135623731
static qsize_t edge_count = 0;



static void addSingleGateNonDiagonalVerticeAndEdge(QuantumProgMap & prog_map,
                                                   qstate_t &gate_tensor,
                                                   qsize_t qubit)
{
    EdgeMap * edge_map = prog_map.getEdgeMap();
    ComplexTensor temp(2, gate_tensor);

    VerticeMatrix * vertice_matrix = prog_map.getVerticeMatrix();
    auto vertice_id = vertice_matrix->getQubitVerticeLastID(qubit);
    auto vertice_id2 = vertice_matrix->addVertice(qubit);

    vector<pair<qsize_t, qsize_t>> contect_vertice = 
                       { { qubit,vertice_id },
                         { qubit,vertice_id2 } };
    edge_count++;
    Edge edge(1, temp, contect_vertice);
    edge_map->insert(pair<qsize_t,Edge>(edge_count, edge));
    vertice_matrix->addContectEdge(qubit, vertice_id, edge_count);
    vertice_matrix->addContectEdge(qubit, vertice_id2, edge_count);
}

static void addSingleGateDiagonalVerticeAndEdge(QuantumProgMap & prog_map,
    qstate_t &gate_tensor,
    qsize_t qubit)
{
    EdgeMap * edge_map = prog_map.getEdgeMap();
    ComplexTensor temp(1, gate_tensor);

    VerticeMatrix * vertice_matrix = prog_map.getVerticeMatrix();
    auto vertice_id = vertice_matrix->getQubitVerticeLastID(qubit);

    vector<pair<qsize_t, qsize_t>> contect_vertice =
    { { qubit,vertice_id } };
    edge_count++;
    Edge edge(1, temp, contect_vertice);
    edge_map->insert(pair<qsize_t, Edge>(edge_count, edge));
    vertice_matrix->addContectEdge(qubit, vertice_id, edge_count);
}




static void addDoubleDiagonalGateVerticeAndEdge(QuantumProgMap & prog_map,
                                                qstate_t &gate_tensor,
                                                qsize_t qubit1,
                                                qsize_t qubit2)
{
    EdgeMap * edge_map = prog_map.getEdgeMap();
    ComplexTensor temp(2, gate_tensor);
    VerticeMatrix * vertice_matrix = prog_map.getVerticeMatrix();
    auto vertice_qubit1_id = vertice_matrix->getQubitVerticeLastID(qubit1);

    auto vertice_qubit2_id = vertice_matrix->getQubitVerticeLastID(qubit2);

    vector<pair<qsize_t, qsize_t>> contect_vertice =
                { { qubit1,vertice_qubit1_id },
                  { qubit2,vertice_qubit2_id } };
    edge_count++;
    Edge edge(2, temp, contect_vertice);
    edge_map->insert(pair<qsize_t, Edge>(edge_count, edge));
    vertice_matrix->addContectEdge(qubit1, vertice_qubit1_id, edge_count);
    vertice_matrix->addContectEdge(qubit2, vertice_qubit2_id, edge_count);
}



static void addDoubleNonDiagonalGateVerticeAndEdge(QuantumProgMap & prog_map,
                                                   qstate_t &gate_tensor,
                                                   qsize_t qubit1,
                                                   qsize_t qubit2)
{
    EdgeMap * edge_map = prog_map.getEdgeMap();
    ComplexTensor temp(4, gate_tensor);
    VerticeMatrix * vertice_matrix = prog_map.getVerticeMatrix();
    auto vertice_qubit1_id = vertice_matrix->getQubitVerticeLastID(qubit1);
    auto vertice_qubit1_id2 = vertice_matrix->addVertice(qubit1);

    auto vertice_qubit2_id = vertice_matrix->getQubitVerticeLastID(qubit2);
    auto vertice_qubit2_id2 = vertice_matrix->addVertice(qubit2);

    vector<pair<qsize_t, qsize_t>> contect_vertice 
        = { { qubit1,vertice_qubit1_id },
            { qubit2,vertice_qubit2_id },
            { qubit1,vertice_qubit1_id2 },
            { qubit2,vertice_qubit2_id2 } };
    edge_count++;
    Edge edge(2, temp, contect_vertice);
    edge_map->insert(pair<qsize_t, Edge>(edge_count, edge));
    vertice_matrix->addContectEdge(qubit1, vertice_qubit1_id, edge_count);
    vertice_matrix->addContectEdge(qubit1, vertice_qubit1_id2, edge_count);

    vertice_matrix->addContectEdge(qubit2, vertice_qubit2_id, edge_count);
    vertice_matrix->addContectEdge(qubit2, vertice_qubit2_id2, edge_count);
}


void H(QuantumProgMap & prog_map,qsize_t qubit)
{
    qstate_t gate_tensor(4, 0);
    gate_tensor[0] = 1 / SQRT2;
    gate_tensor[1] = 1 / SQRT2;
    gate_tensor[2] = 1 / SQRT2;
    gate_tensor[3] = -1 / SQRT2;
    addSingleGateNonDiagonalVerticeAndEdge(prog_map, gate_tensor, qubit);
}

void X(QuantumProgMap & prog_map, qsize_t qubit)
{
    qstate_t gate_tensor(4, 0);
    gate_tensor[1] = 1;
    gate_tensor[2] = 1;
    addSingleGateNonDiagonalVerticeAndEdge(prog_map, gate_tensor, qubit);
}

void RX(QuantumProgMap & prog_map, qsize_t qubit,double angle)
{
	qstate_t gate_tensor(4, 0);
	gate_tensor[0] = cos(angle / 2);
	gate_tensor[1].imag(-1 * sin(angle / 2));
	gate_tensor[2].imag(-1 * sin(angle / 2));
	gate_tensor[3] = cos(angle / 2);
	addSingleGateNonDiagonalVerticeAndEdge(prog_map, gate_tensor, qubit);
}

void Y(QuantumProgMap & prog_map, qsize_t qubit)
{
    qstate_t gate_tensor(4, 0);
    gate_tensor[1].imag(-1) ;
    gate_tensor[2].imag(1);
    addSingleGateNonDiagonalVerticeAndEdge(prog_map, gate_tensor, qubit);
}

void RY(QuantumProgMap & prog_map, qsize_t qubit, double angle)
{
	qstate_t gate_tensor(4, 0);
	gate_tensor[0] = cos(angle / 2);
	gate_tensor[1] = -sin(angle / 2);
	gate_tensor[2] = sin(angle / 2);
	gate_tensor[3] = cos(angle / 2);
	addSingleGateNonDiagonalVerticeAndEdge(prog_map, gate_tensor, qubit);
}

void X1(QuantumProgMap & prog_map, qsize_t qubit)
{
    qstate_t gate_tensor(4, 0);
    gate_tensor[0] = 1 / SQRT2;
    gate_tensor[1] = qcomplex_data_t(0, -1 / SQRT2);
    gate_tensor[2] = qcomplex_data_t(0, -1 / SQRT2);
    gate_tensor[3] = 1 / SQRT2;
    addSingleGateNonDiagonalVerticeAndEdge(prog_map, gate_tensor, qubit);
}

void Y1(QuantumProgMap & prog_map, qsize_t qubit)
{
    qstate_t gate_tensor(4, 0);
    gate_tensor[0] = 1 / SQRT2;
    gate_tensor[1] = -1 / SQRT2;
    gate_tensor[2] = 1 / SQRT2;
    gate_tensor[3] = 1 / SQRT2;
    addSingleGateNonDiagonalVerticeAndEdge(prog_map, gate_tensor, qubit);
}



void Z(QuantumProgMap & prog_map, qsize_t qubit)
{
    qstate_t gate_tensor(2, 0);
    gate_tensor[0] = 1;
    gate_tensor[1] = -1;
    addSingleGateDiagonalVerticeAndEdge(prog_map, gate_tensor, qubit);
}

void RZ(QuantumProgMap & prog_map, qsize_t qubit, double angle)
{
	qstate_t gate_tensor(2, 0);
	gate_tensor[0].real(cos(angle / 2));
	gate_tensor[0].imag(-1 * sin(angle / 2));
	gate_tensor[1].real(cos(angle / 2));
	gate_tensor[1].imag(1 * sin(angle / 2));
	addSingleGateDiagonalVerticeAndEdge(prog_map, gate_tensor, qubit);
}

void Z1(QuantumProgMap & prog_map, qsize_t qubit)
{
    qstate_t gate_tensor(2, 0);
    gate_tensor[0] = qcomplex_data_t(1 / SQRT2, -1 / SQRT2);;
    gate_tensor[1] = qcomplex_data_t(1 / SQRT2, 1 / SQRT2);;
    addSingleGateDiagonalVerticeAndEdge(prog_map, gate_tensor, qubit);
}

void CZ(QuantumProgMap & prog_map,
                                qsize_t qubit1,
                                qsize_t qubit2)
{
    qstate_t gate_tensor(4, 1);
    gate_tensor[0] = 1 ;
    gate_tensor[1] = 1;
    gate_tensor[2] = 1;
    gate_tensor[3] = -1;
    addDoubleDiagonalGateVerticeAndEdge(prog_map,
                                        gate_tensor,
                                        qubit1,
                                        qubit2);
}

void CNOT(QuantumProgMap &prog_map,
                                  qsize_t qubit1,
                                  qsize_t qubit2)
{
    qstate_t gate_tensor(16, 0);
    gate_tensor[0] = 1;
    gate_tensor[5] = 1;
    gate_tensor[11] = 1;
    gate_tensor[14] = 1;
    addDoubleNonDiagonalGateVerticeAndEdge(prog_map,
                                           gate_tensor,
                                           qubit1,
                                           qubit2);
}

#include <iostream>
vector<pair<qsize_t, double>> PMeasure
                    (QuantumProgMap &prog_map,
                     vector<qsize_t> & qubit_vector)
{
    vector<pair<qsize_t, double>> temp;
    auto vertice = prog_map.getVerticeMatrix();


	qubit_vertice_t qubit_vertice_end;
	qubit_vertice_t qubit_vertice_begen;
	for (size_t i = 0; i < qubit_vector.size(); i++)
	{
		auto iter = vertice->getQubitMapIter(qubit_vector[i]);
		auto vertice_map_iter = (*iter).end();
		auto vertice_map_iter_b = (*iter).begin();
		vertice_map_iter--;
		(*vertice_map_iter).second.setValue(1);

		qubit_vertice_end.m_qubit_id = i;
		qubit_vertice_end.m_num = (*vertice_map_iter).first;
		TensorEngine::dimDecrementbyValue(prog_map, qubit_vertice_end, 0);
		qubit_vertice_begen.m_qubit_id = i;
		qubit_vertice_begen.m_num = (*vertice_map_iter_b).first;
		TensorEngine::dimDecrementbyValue(prog_map, qubit_vertice_begen, 0);
	}

    

    /*

		auto iter = vertice->getQubitMapIter(qubit_vector[0]);
	auto vertice_map_iter = (*iter).end();
	vertice_map_iter--;
	(*vertice_map_iter).second.setValue(1);

	iter = vertice->getQubitMapIter(qubit_vector[1]);
	vertice_map_iter = (*iter).end();
	vertice_map_iter--;
	(*vertice_map_iter).second.setValue(0);


	iter = vertice->getQubitMapIter(qubit_vector[2]);
	vertice_map_iter = (*iter).end();
	vertice_map_iter--;
	(*vertice_map_iter).second.setValue(1);

	iter = vertice->getQubitMapIter(qubit_vector[3]);
	vertice_map_iter = (*iter).end();
	vertice_map_iter--;
	(*vertice_map_iter).second.setValue(0);

	iter = vertice->getQubitMapIter(qubit_vector[4]);
	vertice_map_iter = (*iter).end();
	vertice_map_iter--;
	(*vertice_map_iter).second.setValue(1);

	iter = vertice->getQubitMapIter(qubit_vector[5]);
	vertice_map_iter = (*iter).end();
	vertice_map_iter--;
	(*vertice_map_iter).second.setValue(0);

	iter = vertice->getQubitMapIter(qubit_vector[6]);
	vertice_map_iter = (*iter).end();
	vertice_map_iter--;
	(*vertice_map_iter).second.setValue(1);

	iter = vertice->getQubitMapIter(qubit_vector[7]);
	vertice_map_iter = (*iter).end();
	vertice_map_iter--;
	(*vertice_map_iter).second.setValue(0);

	iter = vertice->getQubitMapIter(qubit_vector[8]);
	vertice_map_iter = (*iter).end();
	vertice_map_iter--;
	(*vertice_map_iter).second.setValue(1);

	iter = vertice->getQubitMapIter(qubit_vector[9]);
	vertice_map_iter = (*iter).end();
	vertice_map_iter--;
	(*vertice_map_iter).second.setValue(0);





    */

    

    
    TensorEngine t;
    auto a = t.Merge(prog_map,nullptr);
    auto result = a.real() * a.real() + a.imag() * a.imag();
    std::cout << result << std::endl;
    return temp;
}

/*



void Y_Gate::addVerticeAndEdge(qsize_t qubit)
{
    MatrixXcd gate_matrix(2, 2);
    gate_matrix(0, 0) = 0;
    gate_matrix(0, 1).real(0);
    gate_matrix(0, 1).imag(-1);
    gate_matrix(1, 0).real(0);
    gate_matrix(1, 0).imag(1);
    gate_matrix(1, 1) = 0;
    addSingleGateNonDiagonalVerticeAndEdge(gate_matrix, qubit);
}

void Z_Gate::addVerticeAndEdge(qsize_t qubit)
{
    MatrixXcd gate_matrix(2, 2);
    gate_matrix(0, 0) = 1;
    gate_matrix(0, 1) = 0;
    gate_matrix(1, 0) = 0;
    gate_matrix(1, 1) = -1;
    addSingleGateNonDiagonalVerticeAndEdge(gate_matrix, qubit);
}

void T_Gate::addVerticeAndEdge(qsize_t qubit)
{
    MatrixXcd gate_matrix(2, 2);
    gate_matrix(0, 0) = 0;
    gate_matrix(0, 1) = 0;
    gate_matrix(1, 0) = 0;
    gate_matrix(1, 1).real(1 / SQRT2);
    gate_matrix(1, 1).imag(1 / SQRT2);
    addSingleGateNonDiagonalVerticeAndEdge(gate_matrix, qubit);
}

void S_Gate::addVerticeAndEdge(qsize_t qubit)
{
    MatrixXcd gate_matrix(2, 2);
    gate_matrix(0, 0) = 1;
    gate_matrix(0, 1) = 0;
    gate_matrix(1, 0) = 0;
    gate_matrix(1, 1).real(0);
    gate_matrix(1, 1).imag(1);
    addSingleGateNonDiagonalVerticeAndEdge(gate_matrix, qubit);
}

void U4_Gate::addVerticeAndEdge(qsize_t qubit, double alpha, double beta, double gamma, double delta)
{
    MatrixXcd gate_matrix(2, 2);
    gate_matrix(0,0) = (qcomplex_data_t(cos(alpha - beta / 2 - delta / 2)*cos(gamma / 2),
        sin(alpha - beta / 2 - delta / 2)*cos(gamma / 2)));
    gate_matrix(0, 1) = (qcomplex_data_t(-cos(alpha - beta / 2 + delta / 2)*sin(gamma / 2),
        -sin(alpha - beta / 2 + delta / 2)*sin(gamma / 2)));
    gate_matrix(1, 0) = (qcomplex_data_t(cos(alpha + beta / 2 - delta / 2)*sin(gamma / 2),
        sin(alpha + beta / 2 - delta / 2)*sin(gamma / 2)));
    gate_matrix(1, 1) = (qcomplex_data_t(cos(alpha + beta / 2 + delta / 2)*cos(gamma / 2),
        sin(alpha + beta / 2 + delta / 2)*cos(gamma / 2)));
    addSingleGateNonDiagonalVerticeAndEdge(gate_matrix, qubit);
}

void U1_Gate::addVerticeAndEdge(qsize_t qubit, double angle)
{
    MatrixXcd gate_matrix(2, 2);
    gate_matrix(0, 0) = 1;
    gate_matrix(0, 1) = 0;
    gate_matrix(1, 0) = 0;
    gate_matrix(1, 1).real(cos(angle));
    gate_matrix(1, 1).imag(1 * sin(angle));
    addSingleGateNonDiagonalVerticeAndEdge(gate_matrix, qubit);
}

void ISWAP_Gate::addVerticeAndEdge(qsize_t qubit1, qsize_t qubit2)
{
    MatrixXcd gate_matrix = MatrixXcd::Zero(4, 4);
    gate_matrix(0,0) = 1;
    gate_matrix(1,1) = 0;
    gate_matrix(1,2).imag(-1);
    gate_matrix(2,1).imag(-1);
    gate_matrix(2,2) = 0;
    gate_matrix(3, 3) = 1;
    addDoubleNonDiagonalGateVerticeAndEdge(gate_matrix, qubit1, qubit2);
}







*/
