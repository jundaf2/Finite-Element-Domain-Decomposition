/*
 * main.cpp
 *
 *  Created on: Apr 15, 2021
 *      Author: ubuntu
 */
#include<iostream>
#include<fstream>
#include<cmath>
#include<assert.h>
#include<stdio.h>
#include<cmath>
#include<sys/time.h>
#include <string>
#include <sstream>

#include "Serial_CGM.h"
using namespace std;

const string file_name = "naive_fe_ddm_";
float_t Hmax;

const int max_asmiter_num = 100;
const int Num_Subdomain = 4;
const float_t electric_charge = 1;
const float_t epsilon = 1;

// determine which neighboring rank to receive from
const int Nodes_Interface_InWhichSubDomain[4][3] = {{4,3,2},{4,3,1},{4,2,1},{3,2,1}};

//
int WriteOutput(vector<float_t>& phi){
	stringstream buf;
	buf<<Hmax;
	string str;
	str = buf.str();

	ofstream outfile(file_name+str+"out.txt",fstream::out);

	if (outfile.is_open())
   {
		for(auto a:phi) outfile << " " << a;
		outfile.close();
		return 0;
	}
	else{
		cout << "Write Output txt file failed" << endl;
		return 1;
	}
}

int ReadSubdomainInfo(int_t* Num_Nodes,int_t* Num_Elements,vector<float_t>* coords,vector<int_t>* Element2NodeList,vector<int_t>* Node_Source,vector<int_t>* Elements_Node_Source,vector<vector<int_t>>* Nodes_InteriorSubDomain,vector<vector<vector<int_t>>>* Nodes_InterfaceSubDomain,vector<vector<int_t>>* Nodes_Dirichlet,vector<vector<int_t>>* Elements_Subdomain,vector<int_t>* Nodes_onBoundary){
	// subdomain geometry information generated by matlab
	stringstream buf;
	buf<<Hmax;
	string str;
	str = buf.str();

	ifstream infile(file_name+str+".txt",fstream::in); // load file
	if (!infile.is_open())
	{
		cout << "read failed" << endl;
		return 0;
	}
	string one_line; // one line in the .txt file
	int_t node_id, element_id;  // temp variable
	float_t coord; // temp variable
	vector<int_t> lines; // temp variable
	vector<vector<int_t>> lines_feat; // temp variable


	getline(infile, one_line);
	stringstream stringin(one_line);
	stringin >> Num_Nodes[0];
	stringin >> Num_Elements[0];

	getline(infile, one_line);
	stringin.clear();
	stringin.str(one_line);
	while (stringin >> coord) {
		coords->push_back(coord);
	}
	getline(infile, one_line);
	stringin.clear();
	stringin.str(one_line);
	while (stringin >> node_id) {
		Element2NodeList->push_back(node_id);
	}

	getline(infile, one_line);
	stringin.clear();
	stringin.str(one_line);
	while (stringin >> node_id) {
		Node_Source->push_back(node_id);
	}

	getline(infile, one_line);
	stringin.clear();
	stringin.str(one_line);
	while (stringin >> element_id) {
		Elements_Node_Source->push_back(element_id);
	}

	for(int i=0;i<Num_Subdomain;i++){
		getline(infile, one_line);
		stringin.clear();
		stringin.str(one_line);
		lines.clear();
		while (stringin >> node_id) {
			lines.push_back(node_id);
		}
		Nodes_InteriorSubDomain->push_back(lines);
	}

	for(int i=0;i<Num_Subdomain;i++){
		lines_feat.clear();
		for(int j=0;j<Num_Subdomain-1;j++){
			getline(infile, one_line);
			stringin.clear();
			stringin.str(one_line);
			lines.clear();
			while (stringin >> node_id) {
				lines.push_back(node_id);
			}
			lines_feat.push_back(lines);
		}
		Nodes_InterfaceSubDomain->push_back(lines_feat);
	}

	// on each subdomain
	for(int i=0;i<Num_Subdomain;i++){
		getline(infile, one_line);
		stringin.clear();
		stringin.str(one_line);
		lines.clear();
		while (stringin >> node_id) {
			lines.push_back(node_id);
		}
		Nodes_Dirichlet->push_back(lines);
	}

	int Num_Subdomain_Elements;
	getline(infile, one_line); // no need to store the nums for now

	for(int i=0;i<Num_Subdomain;i++){
		getline(infile, one_line);
		stringin.clear();
		stringin.str(one_line);
		lines.clear();
		while (stringin >> element_id) {
			lines.push_back(element_id);
		}
		Elements_Subdomain->push_back(lines);
	}

	getline(infile, one_line);
	stringin.clear();
	stringin.str(one_line);
	while (stringin >> node_id) {
		Nodes_onBoundary->push_back(node_id);
	}

	if(!infile.eof()) cout << "... read sequence has some problem ..." << endl;
	else cout << "... read sequence has no problem ..." << endl;

	infile.close();
	return 0;
}

// select sub matrix of original matrix
vector<vector<float_t>> SubMatrix(vector<vector<float_t>>&A, vector<int_t>&rows,vector<int_t>&columns){
	vector<vector<float_t>> Asub;
	vector<float_t> Arow;
	for(auto row:rows){
		Arow.clear();
		for(auto column:columns){
			Arow.push_back(A[row][column]);
		}
		Asub.push_back(Arow);
	}

	return Asub;
}

// select sub array of original array
vector<float_t> SubArray(vector<float_t>&x, vector<int_t>&indices){
	vector<float_t> xsub;
	for(auto index:indices){
		xsub.push_back(x[index]);
	}

	return xsub;
}

// assign values of sub array to original array
void AssignSubArray(vector<float_t>*x, vector<float_t>&xsub, vector<int_t>&indices){
	int_t i= 0;
	for(auto index:indices){
		x->at(index) = xsub[i++];
	}
}

int main(int argv,char **args){
	string s  = args[1];
	Hmax = stof(s);

	// geometry
	int_t Num_Nodes,Num_Elements;
	vector<float_t> coords;
	vector<int_t> Element2NodeList;
	vector<int_t> Node_Source;
	vector<int_t> Elements_Node_Source;
	vector<vector<int_t>> Nodes_InteriorSubDomain;
	vector<vector<vector<int_t>>> Nodes_InterfaceSubDomain;
	vector<vector<int_t>> Nodes_Dirichlet; // on each subdomain
	vector<vector<int_t>> Elements_Subdomain;
	vector<int_t> Nodes_onBoundary;

	try {
		ReadSubdomainInfo(&Num_Nodes,&Num_Elements,&coords, &Element2NodeList, &Node_Source,&Elements_Node_Source,&Nodes_InteriorSubDomain,&Nodes_InterfaceSubDomain,&Nodes_Dirichlet,&Elements_Subdomain,&Nodes_onBoundary);
		cout << "Num_Nodes " << Num_Nodes << endl;
		cout << "Num_Elements " << Num_Elements << endl;
		cout << "coords size " << coords.size() << endl;
		cout << "Element2NodeList size " << Element2NodeList.size() << endl;
	}
	catch (std::exception& e) {
	  cout << e.what() << endl;
	}

	// timing
	timeval  start,end;

	// the large total system, solve this for testing the FEM part
	vector<vector<float_t>> K(Num_Nodes,vector<float_t>(Num_Nodes,0));
	vector<vector<float_t>> Ke(3,vector<float_t>(3,0));
	vector<float_t> b(Num_Nodes,0);
	vector<float_t> eps(Num_Elements,epsilon);
	vector<float_t> phi(Num_Nodes,0); // solution for the entire domain
	vector<float_t> phi_asm(Num_Nodes,0);
	// start timing
	gettimeofday(&start,NULL);

	float_t ie[3], xe[3], ye[3], be[3], ce[3], Delta_e;
	for(int_t e=0;e<Num_Elements;e++){
		for(int ne = 0; ne < 3; ++ne) {
			ie[ne] = Element2NodeList [3 * e + ne];
			xe[ne] = coords[2 * Element2NodeList [3 * e + ne]];
			ye[ne] = coords[2 * Element2NodeList [3 * e + ne] + 1];
		}
		be[0]=ye[1]-ye[2];
		be[1]=ye[2]-ye[0];
		be[2]=ye[0]-ye[1];

		ce[0]=xe[2]-xe[1];
		ce[1]=xe[0]-xe[2];
		ce[2]=xe[1]-xe[0];
		Delta_e=(be[0]*ce[1]-be[1]*ce[0])/2;
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				Ke[i][j]=eps[e]/(4*Delta_e)*(be[i]*be[j]+ce[i]*ce[j]);
			}
		}
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				K[ie[i]][ie[j]]+=Ke[i][j];
			}
		}
	}

	// Impose the Dirichlet boundary condition
	int_t node_id, element_id;
	for(int i =0;i<Nodes_onBoundary.size();i++){
		node_id = Nodes_onBoundary[i];
		b[Nodes_onBoundary[i]]=0;
		K[node_id][node_id]=1;
		for(int n=0;n<Num_Nodes;n++){
			if(n==node_id) continue;
			b[n]-=K[n][node_id]*0;
			K[n][node_id]=0;
			K[node_id][n]=0;
		}
	}

	// Impose Source
	for(int e =0;e<Elements_Node_Source.size();e++){
		element_id = Elements_Node_Source[e];
		for(int ne = 0; ne < 3; ++ne) {
			ie[ne] = Element2NodeList [3 * element_id + ne];
			xe[ne] = coords[2 * Element2NodeList [3 * element_id + ne]];
			ye[ne] = coords[2 * Element2NodeList [3 * element_id + ne] + 1];
		}

		be[0]=ye[1]-ye[2];
		be[1]=ye[2]-ye[0];
		be[2]=ye[0]-ye[1];

		ce[0]=xe[2]-xe[1];
		ce[1]=xe[0]-xe[2];
		ce[2]=xe[1]-xe[0];
		Delta_e=(be[0]*ce[1]-be[1]*ce[0])/2;

		for(int n=0;n<Node_Source.size();n++){
			node_id=Node_Source[n];
			if(ie[0]==node_id||ie[1]==node_id||ie[2]==node_id)
				b[node_id]+=electric_charge*1/3*Delta_e;
		}
	}
/*
	// traditional FEM solving
	//displayMat(K);
	//displayVec(b);
	phi = Serial_CG(K,b);
	//displayVec(phi);

	// end timing
	gettimeofday(&end,NULL);
	float_t run_time = (end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec);
	run_time /= 1000000;
	cout<<" ... traditional serial fem costs time "<<run_time<<" s ... "<<endl;
// start timing
	gettimeofday(&start,NULL);
	*/

	// dividing the global matrix to sub domain matrices

	// four linear systems for sub domains
	// the linear system has a length of the interior node numbers, do not include the artificial boundary
	// there will be a separate array for artificial boundary nodes of each sub domain. They will be used as the MPI receiving buffer as well as the array to enforce artificial boundary.
	// after receive the other ranks, the
	// each rank also acts as the sender, so the rank need to know which nodal coefficients to send


	// although the linear system has a length of interior node numbers, the assembly procedure should be done for the discretized region including the interface boundary.

	// Nodes_InteriorSubDomain,Nodes_InterfaceSubDomain

	vector<float_t> b1, b2, b3, b4;
	vector<float_t> phi1, phi2, phi3, phi4;
	vector<vector<float_t>> K1, K2, K3, K4, K14, K13, K12, K24, K23, K21, K34, K32, K31, K43, K42, K41;
	b1 = SubArray(b,Nodes_InteriorSubDomain[0]);
	b2 = SubArray(b,Nodes_InteriorSubDomain[1]);
	b3 = SubArray(b,Nodes_InteriorSubDomain[2]);
	b4 = SubArray(b,Nodes_InteriorSubDomain[3]);

	phi1 = SubArray(phi_asm,Nodes_InteriorSubDomain[0]);
	phi2 = SubArray(phi_asm,Nodes_InteriorSubDomain[1]);
	phi3 = SubArray(phi_asm,Nodes_InteriorSubDomain[2]);
	phi4 = SubArray(phi_asm,Nodes_InteriorSubDomain[3]);
	cout << " ... SubArray success ..." << endl;

	K1 = SubMatrix(K,Nodes_InteriorSubDomain[0],Nodes_InteriorSubDomain[0]);
	K2 = SubMatrix(K,Nodes_InteriorSubDomain[1],Nodes_InteriorSubDomain[1]);
	K3 = SubMatrix(K,Nodes_InteriorSubDomain[2],Nodes_InteriorSubDomain[2]);
	K4 = SubMatrix(K,Nodes_InteriorSubDomain[3],Nodes_InteriorSubDomain[3]);

	K14 = SubMatrix(K,Nodes_InteriorSubDomain[0],Nodes_InterfaceSubDomain[0][0]);
	K13 = SubMatrix(K,Nodes_InteriorSubDomain[0],Nodes_InterfaceSubDomain[0][1]);
	K12 = SubMatrix(K,Nodes_InteriorSubDomain[0],Nodes_InterfaceSubDomain[0][2]);

	K24 = SubMatrix(K,Nodes_InteriorSubDomain[1],Nodes_InterfaceSubDomain[1][0]);
	K23 = SubMatrix(K,Nodes_InteriorSubDomain[1],Nodes_InterfaceSubDomain[1][1]);
	K21 = SubMatrix(K,Nodes_InteriorSubDomain[1],Nodes_InterfaceSubDomain[1][2]);

	K34 = SubMatrix(K,Nodes_InteriorSubDomain[2],Nodes_InterfaceSubDomain[2][0]);
	K32 = SubMatrix(K,Nodes_InteriorSubDomain[2],Nodes_InterfaceSubDomain[2][1]);
	K31 = SubMatrix(K,Nodes_InteriorSubDomain[2],Nodes_InterfaceSubDomain[2][2]);

	K43 = SubMatrix(K,Nodes_InteriorSubDomain[3],Nodes_InterfaceSubDomain[3][0]);
	K42 = SubMatrix(K,Nodes_InteriorSubDomain[3],Nodes_InterfaceSubDomain[3][1]);
	K41 = SubMatrix(K,Nodes_InteriorSubDomain[3],Nodes_InterfaceSubDomain[3][2]);
	cout << " ... SubMatrix success ..." << endl;

	vector<float_t>* bdm[4] = { &b1, &b2, &b3, &b4};
	vector<float_t>* phidm[4] = { &phi1, &phi2, &phi3, &phi4};
	vector<vector<float_t>>* Kdm[4] = { &K1, &K2, &K3, &K4};
	vector<vector<float_t>>* Kdm_slf_adj[4][3] = {{&K14, &K13, &K12},{&K24, &K23, &K21},{&K34, &K32, &K31},{&K43, &K42, &K41}};


	cout << " ... serial variables define success ..." << endl;
	
	vector<float_t> phi_temp;
	vector<float_t> enforceABC;
	vector<float_t> r;
	float_t rho, error;
	float_t norm_b = sqrt(INNER_PRODUCT(b,b));
	for(int asm_iter=0; asm_iter<max_asmiter_num; asm_iter++){
		// serial version, for each sub domain
		for(int slf = 0; slf<Num_Subdomain;slf++){
			vector<float_t> rhs(Nodes_InteriorSubDomain[slf].size(),0);
			// calculate RHS enforced by artificial boundary condition
			for(int i=0;i<rhs.size();i++){
				rhs[i]+=(*bdm[slf])[i];
			}
			for(int adj = 0; adj<Num_Subdomain-1;adj++){
				phi_temp = SubArray(phi_asm,Nodes_InterfaceSubDomain[slf][adj]);
				enforceABC = MATRIX_VECTOR_PRODUCT(*Kdm_slf_adj[slf][adj], phi_temp);
				for(int i=0;i<rhs.size();i++){
					rhs[i] -= enforceABC[i];
				}
			}

			// solve the sub domain linear system using the interior K and enforced RHS
			phi_temp = Serial_CG(*Kdm[slf],rhs);

			// serial: update the global phi_asm
			AssignSubArray(&phi_asm,phi_temp,Nodes_InteriorSubDomain[slf]);
			// parallel: send the inside nodal coefficient to the neighboring ranks to act as artificial boundary.
			//			 recv from the neighboring ranks and use that for enforcing RHS

		}
		// check convergence on the total domain
		//r = Residual(K,phi_asm,b);
		//rho = INNER_PRODUCT(r,r);
		//error = sqrt(rho)/norm_b;
		// cout << error << endl;
		//if(error < TOL*10){
		//	break;
		//}
	}
	// end timing
	gettimeofday(&end,NULL);
	float_t run_time = (end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec);
	run_time /= 1000000;
	cout<<" ... ASM serial version costs time "<<run_time<<" s ... "<<endl;


	//for(int i=0;i<Num_Nodes;i++){
	//	cout << abs(phi_asm[i]-phi[i]) << " " << endl;
	//}
	return WriteOutput(phi_asm);
}

