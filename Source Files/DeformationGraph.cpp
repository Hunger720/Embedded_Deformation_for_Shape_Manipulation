#include <cmath>
#include "DeformationGraph.h"


//void DeformationGraph::neighbours(const int j, int &n_neighbours, int *output){
//	n_neighbours = 0;
//}


void DeformationGraph::computeWeights(const float *vertex, float *weights, int *idx){
	float *d=new float[n_nodes], temp[3], sum = 0;
	//find k+1-nearest nodes
	for(int j=0; j<n_nodes; j++){
		for(int i=0; i<3; i++)
			temp[i] = vertex[i]-nodes[j][i];
		d[j] = sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]);
	}
	for(int i=0; i<k+1; i++){   //selection sort
		idx[i] = i;
		for(int j=n_nodes-1; j>i; j--)
			if(d[j]<d[idx[i]]) idx[i] = j;
		temp[0] = d[i];
		d[i] = d[idx[i]];
		d[idx[i]] = temp[0];
	}
	//compuate wj(v)
	for(int i=0; i<3; i++)
		temp[i] = vertex[i]-nodes[idx[k]][i];
	float dmax =  sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]);
	for(int j=0; j<k; j++){
		for(int i=0; i<3; i++)
			temp[i] = vertex[i]-nodes[idx[j]][i];
		weights[j] = pow(1-sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2])/dmax,2);
		sum += weights[j];
	}
	//normalize to sum to 1
	for(int j=0; j<k; j++)
		weights[j] /= sum;
	delete[] d;
}


void DeformationGraph::predict(const float *input, float *output){
	output[0] = 0;output[1] = 0;output[2] = 0;
	int *idx = new int[k+1];
	float *weights=new float[k], temp[3];
	computeWeights(input,weights,idx);
	for(int j=0; j<k; j++){
		for(int i=0; i<3; i++)
			temp[i] = input[i] - nodes[idx[j]][i];
		for(int i=0; i<3; i++)
			output[i] += weights[j]*(rot[idx[j]][i]*temp[0]+rot[idx[j]][i+3]*temp[1]+rot[idx[j]][i+6]*temp[2]+nodes[idx[j]][i]+trans[idx[j]][i]);
	}
	delete[] idx;
	delete[] weights;
}


float DeformationGraph::Erot(){
	e_rot = 0;
	for(int j=0; j<n_nodes; j++){
		e_rot += pow(rot[j][0]*rot[j][3]+rot[j][1]*rot[j][4]+rot[j][2]*rot[j][5],2);
		e_rot += pow(rot[j][0]*rot[j][6]+rot[j][1]*rot[j][7]+rot[j][2]*rot[j][8],2);
		e_rot += pow(rot[j][3]*rot[j][6]+rot[j][4]*rot[j][7]+rot[j][5]*rot[j][8],2);
		e_rot += pow(rot[j][0]*rot[j][0]+rot[j][1]*rot[j][1]+rot[j][2]*rot[j][2]-1,2);
		e_rot += pow(rot[j][3]*rot[j][3]+rot[j][4]*rot[j][4]+rot[j][5]*rot[j][5]-1,2);
		e_rot += pow(rot[j][6]*rot[j][6]+rot[j][7]*rot[j][7]+rot[j][8]*rot[j][8]-1,2);
	}
	return e_rot;
}


float DeformationGraph::Ereg(){
	float temp[3];
	e_reg = 0;
	for(int j=0; j<n_nodes; j++){
		for(int n=0; n<n_nodes; n++){
			if(edges[j][n]){
				for(int i=0; i<3; i++)
					temp[i] = nodes[n][i]-nodes[j][i];
				for(int i=0; i<3; i++)
					e_reg += pow(rot[j][i]*temp[0]+rot[j][i+3]*temp[1]+rot[j][i+6]*temp[2]+nodes[j][i]+trans[j][i]-nodes[n][i]-trans[n][i],2);
			}
		}
	}
	return e_reg;
}


float DeformationGraph::Econ(const int p, const float **v, const float **q){
	float _v[3];
	e_con = 0;
	for(int l=0; l<p; l++){
		predict(v[l],_v);
		for(int i=0; i<3; i++){
			e_con += pow(_v[i]-q[l][i],2);
		}
	}
	return e_con;
}