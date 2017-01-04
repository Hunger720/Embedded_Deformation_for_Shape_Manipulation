#ifndef __DEFORMATION_GRAPH_H__
#define __DEFORMATION_GRAPH_H__

#ifndef __OPENCV_HPP__
#define __OPENCV_HPP__
#include <opencv2/opencv.hpp>
using namespace cv;
#endif

typedef double Node[3];
typedef double Rotation[9];
typedef double Translation[3];


class DeformationGraph{
public:
	DeformationGraph();
	//create a deformation graph according to nodes data, mesh vertices data and k
	DeformationGraph(int n_nodes, Node *nodes, int n_vertices, double **vertices, int k);
	//input nodes and edges data to create a deformation graph
	DeformationGraph(int n_nodes, Node *nodes, int n_edges, bool **edges);
	
	~DeformationGraph();

	void findNearestNodes(const double *vertex, int k, int *idx)const;

	//calculate the weight of a vertex
	//weights is a double array, contains k-nearest nodes' weights
	//idx is an int array, contain k+1-nearest nodes' indexes
	void computeWeights(const double *vertex, double *weights, int *idx) const;

	// set k of k-nearest node(s)
	void setK(int k){ k_nearest = k; }

	//set weights of Erot, Ereg and Econ
	inline void setWeights(double w_rot, double w_reg, double w_con){
		this->w_rot = w_rot;
		this->w_reg = w_reg;
		this->w_con = w_con;
	}

	//predict the position of a vertex or a node
	//if the number of nodes of deformation graph is 0, the output will be (0,0,0)
	void predict(const double *v, double *_v) const;

	double Erot();  //compute the rotation error
	double Ereg();  //compute the regularization error
	double Econ(const int p, double **v, double **q);  //compute the constraints error

	friend void computef(const DeformationGraph &dg, int &p, double **v, double **q, OutputArray f);
	friend void computeJ(const DeformationGraph &dg, const int &p, double **v, double **q, OutputArray J);
	friend double F(DeformationGraph &dg, Mat &x, int p, double **v, double **q);   //x is used to update dg
	friend void gaussNewton(DeformationGraph &dg, int p, double **v, double **q);

private:
	int n_nodes;         //number of nodes
	int n_edges;         //number of edges
	int k_nearest;       //number of k (k-nearest nodes)
	Node *nodes;         //nodes array
	bool **edges;        //edges, a adjacency matrix, n_nodes * n_nodes dimensions
	Rotation *rot;       //rotation matrixes
	Translation *trans;  //translation vectors
	double w_rot;        //rotation error weight
	double w_reg;        //regularization error weight
	double w_con;        //constraints error weig
};

#endif