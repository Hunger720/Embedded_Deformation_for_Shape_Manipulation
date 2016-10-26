#include <cmath>
#include "DeformationGraph.h"



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
}


void DeformationGraph::predict(const float *input, float *&output){
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
}