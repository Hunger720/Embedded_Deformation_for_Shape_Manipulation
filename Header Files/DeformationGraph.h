typedef float Node[3];
typedef float Rotation[9];
typedef float Translation[3];


class DeformationGraph{
public:
	//findout the neighbours of a node

	//calculate the weight of a vertex
	//weights is a float array, contains k-nearest nodes' weights
	//idx is an int array, contain k+1-nearest nodes' indexes
	void computeWeights(const float *vertex, float *weights, int *idx);

	//predict the position of a vertex or a node
	//if the number of nodes of deformation graph is 0, the output will be (0,0,0)
	void predict(const float *input, float *&output);
private:
	int n_nodes;   //number of nodes
	int n_edges;   //number of edges
	int k;         //number of k( k-nearest nodes)
	Node *nodes;   //nodes array
	float **edges; //edges, a adjacency matrix
	Rotation *rot;   //rotation matrixes
	Translation *trans;   //translation vectors
};