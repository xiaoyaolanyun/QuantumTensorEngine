#include <iostream>
#include "QuantumGates.h"
#include "TensorEngine.h"
#include "TensorNode.h"
#include "QRunesParserAPI.h"
using namespace std;

void init(QuantumProgMap & prog_map, qsize_t qubit)
{
    VerticeMatrix  *vertice_matrix = prog_map.getVerticeMatrix();
    vertice_matrix->initVerticeMatrix(qubit);
}

int main()
{
    QuantumProgMap prog_map;
	string filename = "QRunes.txt";
	QList mQList;
	int iResult = QRunesParser(filename, mQList);
    init(prog_map, iResult);

	QNode *ptQNode = mQList.qList.head->next;

	while (nullptr != ptQNode)
	{
	
		if (ptQNode->gateSpecifier == RX_GATE)
		{
			SingleAngleGateNode * pSingleAngleGateNode = (SingleAngleGateNode *)ptQNode;
			RX(prog_map, pSingleAngleGateNode->opQubit, pSingleAngleGateNode->angle);
		}
		else if (ptQNode->gateSpecifier == RY_GATE)
		{
			SingleAngleGateNode * pSingleAngleGateNode = (SingleAngleGateNode *)ptQNode;
			RY(prog_map, pSingleAngleGateNode->opQubit, pSingleAngleGateNode->angle);
		}
		else if (ptQNode->gateSpecifier == RZ_GATE)
		{
			SingleAngleGateNode * pSingleAngleGateNode = (SingleAngleGateNode *)ptQNode;
			RZ(prog_map, pSingleAngleGateNode->opQubit, pSingleAngleGateNode->angle);
		}
		else if (ptQNode->gateSpecifier == CR_GATE)
		{
			DoubleAngleGateNode * pDoubleAngle = (DoubleAngleGateNode *)ptQNode;
			CZ(prog_map, pDoubleAngle->opQubit1, pDoubleAngle->opQubit2);
		}

		ptQNode = ptQNode->next;
	}
	


	auto start = clock();
    TensorEngine::Merge(prog_map);
    vector<qsize_t> vector;
    for (size_t i = 0; i < iResult; i++)
    {
        vector.push_back(i);
    }
    PMeasure(prog_map,vector);
	auto end = clock();
	cout << "end - start = " << end - start << endl;
    system("pause");
    return 0;
}
