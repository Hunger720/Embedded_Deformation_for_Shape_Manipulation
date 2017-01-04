#include <cmath>
#include "DeformationGraph.h"


DeformationGraph::DeformationGraph()
	:n_nodes(0),
	 n_edges(0),
	 k_nearest(1),
	 w_rot(1),
	 w_reg(1),
	 w_con(1){
	nodes = NULL;
	edges = NULL;
	rot = NULL;
	trans = NULL;
}


DeformationGraph::DeformationGraph(int n_nodes, Node *nodes, int n_vertices, double **vertices, int k)
	:n_nodes(n_nodes),
	 n_edges(0),
	 k_nearest(k),
	 w_rot(1),
	 w_reg(1),
	 w_con(1){
	this->nodes = new Node[n_nodes];
	this->edges = new bool*[n_nodes];
	rot = new Rotation[n_nodes];
	trans = new Translation[n_nodes];
	for(int j=0; j<n_nodes; j++){
		//initialize nodes
		for(int i=0; i<3; i++) this->nodes[j][i] = nodes[j][i];
		//initialize edges
		this->edges[j] = new bool[n_nodes];
		for(int i=0; i<n_nodes; ++i) edges[j][i] = false;
		//initialize rotation matrixes
		for(int i=0; i<9; i++) rot[j][i] = 0.1;
		rot[j][0] = 0.9; rot[j][4] = 0.9; rot[j][8] = 0.9;
		//initialize translation vectors
		for(int i=0; i<3; i++) trans[j][i] = 0;
	}
	//compute graph connectivity
	int *idx = new int[k_nearest];
	for(int i=0; i<n_vertices; ++i){
		findNearestNodes(vertices[i],k_nearest,idx);
		for(int j=0; j<k_nearest; ++j)
			for(int jj=j+1; jj<k_nearest; ++jj){
				edges[idx[j]][idx[jj]] = true;
				edges[idx[jj]][idx[j]] = true;
			}
	}
	for(int i=0; i<n_nodes; ++i)
		for(int j=0; j<n_nodes; ++j)
			if(edges[i][j]) ++n_edges;
	n_edges /= 2;
}


DeformationGraph::DeformationGraph(int n_nodes, Node *nodes, int n_edges, bool **edges)
	:n_nodes(n_nodes),
	 n_edges(n_edges),
	 k_nearest(1),
	 w_rot(1),
	 w_reg(1),
	 w_con(1){
	this->nodes = new Node[n_nodes];
	this->edges = new bool*[n_nodes];
	rot = new Rotation[n_nodes];
	trans = new Translation[n_nodes];
	for(int j=0; j<n_nodes; j++){
		//initialize nodes
		for(int i=0; i<3; i++) this->nodes[j][i] = nodes[j][i];
		//initialize edges
		this->edges[j] = new bool[n_nodes];
		for(int i=0; i<n_nodes; i++) this->edges[j][i] = edges[j][i];
		//initialize rotation matrixes
		for(int i=0; i<9; i++) rot[j][i] = 0.1;
		rot[j][0] = 0.9; rot[j][4] = 0.9; rot[j][8] = 0.9;
		//initialize translation vectors
		for(int i=0; i<3; i++) trans[j][i] = 0;
	}
}


DeformationGraph::~DeformationGraph(){
	delete[] nodes; nodes = NULL;
	delete[] rot; rot = NULL;
	delete[] trans; trans = NULL;
	for(int j=0; j<n_nodes; j++){
		delete[] edges[j];
		edges[j] = NULL;
	}
	delete[] edges; edges = NULL;
}


void DeformationGraph::findNearestNodes(const double *vertex, int k, int *idx)const{
	double *d=new double[n_nodes], temp[3];
	int *index = new int[n_nodes], idx_min, t;
	//find k nearest nodes
	for(int j=0; j<n_nodes; j++){
		index[j] = j;
		for(int i=0; i<3; i++) temp[i] = vertex[i]-nodes[j][i];
		d[j] = sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]);
	}
	for(int i=0; i<k; ++i){   //selection sort
		idx_min = i;
		for(int j=n_nodes-1; j>i; j--) if(d[j]<d[idx_min]) idx_min = j;
		idx[i] = index[idx_min];
		temp[0] = d[i]; t=index[i];
		d[i] = d[idx_min]; index[i]=index[idx_min];
		d[idx_min] = temp[0]; index[idx_min]=t;
	}
}


void DeformationGraph::computeWeights(const double *vertex, double *weights, int *idx) const{
	double *d=new double[n_nodes], dist, temp[3], sum = 0;
	int *index = new int[n_nodes], idx_min, t;
	//find k+1-nearest nodes
	for(int j=0; j<n_nodes; j++){
		index[j] = j;
		for(int i=0; i<3; i++) temp[i] = vertex[i]-nodes[j][i];
		d[j] = sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]);
	}
	for(int i=0; i<k_nearest+1; i++){   //selection sort
		idx_min = i;
		for(int j=n_nodes-1; j>i; j--) if(d[j]<d[idx_min]) idx_min = j;
		idx[i] = index[idx_min];
		temp[0] = d[i]; t=index[i];
		d[i] = d[idx_min]; index[i]=index[idx_min];
		d[idx_min] = temp[0]; index[idx_min]=t;
	}
	//compuate wj(v)
	for(int i=0; i<3; i++) temp[i] = vertex[i]-nodes[idx[k_nearest]][i];
	double dmax =  sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]);
	for(int j=0; j<k_nearest; j++){
		for(int i=0; i<3; i++) temp[i] = vertex[i]-nodes[idx[j]][i];
		dist = sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2]);
		weights[j] = pow(1-dist/dmax,2);
		sum += weights[j];
	}
	//normalize to sum to 1
	if(k_nearest==1) weights[0]=1;
	else for(int j=0; j<k_nearest; j++) weights[j] /= sum;
	delete[] d; d = NULL;
	delete[] index; index = NULL;
}


void DeformationGraph::predict(const double *input, double *output) const{
	for(int i=0; i<3; ++i) output[i]=0;
	int    *idx     = new int   [k_nearest+1];
	double *weights = new double[k_nearest], temp[3];
	computeWeights(input,weights,idx);
	for(int j=0; j<k_nearest; j++){
		for(int i=0; i<3; i++)
			temp[i] = input[i]-nodes[idx[j]][i];
		for(int i=0; i<3; i++)
			output[i] += weights[j]*(rot[idx[j]][i]*temp[0]+rot[idx[j]][i+3]*temp[1]+rot[idx[j]][i+6]*temp[2]+nodes[idx[j]][i]+trans[idx[j]][i]);
	}
	delete[] idx;     idx = NULL;
	delete[] weights; weights = NULL;
}


double DeformationGraph::Erot(){
	double e_rot = 0;
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


double DeformationGraph::Ereg(){
	double temp[3], e_reg = 0;
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


double DeformationGraph::Econ(const int p, double **v, double **q){
	double _v[3], e_con = 0;
	for(int l=0; l<p; l++){
		predict(v[l],_v);
		for(int i=0; i<3; i++) e_con += pow(_v[i]-q[l][i],2);
	}
	return e_con;
}