#ifndef __DEFORMATION_GRAPH_H__
#define __DEFORMATION_GRAPH_H__

typedef float Node[3];
typedef float Rotation[9];
typedef float Translation[3];


class DeformationGraph{
public:
	//findout the neighbours of a node
	//don't forget to delete output when you finish with it
	//void neighbours(const int j, int &n_neighbouts, int *output);

	//calculate the weight of a vertex
	//weights is a float array, contains k-nearest nodes' weights
	//idx is an int array, contain k+1-nearest nodes' indexes
	void computeWeights(const float *vertex, float *weights, int *idx);

	//predict the position of a vertex or a node
	//if the number of nodes of deformation graph is 0, the output will be (0,0,0)
	void predict(const float *input, float *output);

	float Erot();  //compute the rotation error
	float Ereg();  //compute the regularization error
	float Econ(const int p, const float **v, const float **q);  //compute the constraints error

	//compute the value of the objective function
	inline float F(const float w_rot, const float w_reg, const float w_con){
		return w_rot*e_rot+w_reg*e_reg+w_con*e_con;
	}
private:
	int n_nodes;   //number of nodes
	int n_edges;   //number of edges
	int k;         //number of k( k-nearest nodes)
	Node *nodes;   //nodes array
	bool **edges;  //edges, a adjacency matrix, n_nodes * n_nodes dimensions
	Rotation *rot; //rotation matrixes
	Translation *trans;   //translation vectors
	float e_rot;   //rotation error
	float e_reg;   //regularization error
	float e_con;   //constraints error
};

#endif