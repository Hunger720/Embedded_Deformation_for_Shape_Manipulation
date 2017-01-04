#include <cmath>
#include "DeformationGraph.h"


void computef(const DeformationGraph &dg, int &p, double **v, double **q, OutputArray f){
	double _v[3];
	int idx = 0;
	Mat fx = Mat::zeros(6*dg.n_nodes+6*dg.n_edges+3*p,1,CV_64F);
	//Erot
	for(int j=0; j<dg.n_nodes; j++){
		fx.at<double>(idx++) = (dg.rot[j][0]*dg.rot[j][3]+dg.rot[j][1]*dg.rot[j][4]+dg.rot[j][2]*dg.rot[j][5])*sqrt(dg.w_rot);
		fx.at<double>(idx++) = (dg.rot[j][0]*dg.rot[j][6]+dg.rot[j][1]*dg.rot[j][7]+dg.rot[j][2]*dg.rot[j][8])*sqrt(dg.w_rot);
		fx.at<double>(idx++) = (dg.rot[j][3]*dg.rot[j][6]+dg.rot[j][4]*dg.rot[j][7]+dg.rot[j][5]*dg.rot[j][8])*sqrt(dg.w_rot);
		fx.at<double>(idx++) = (dg.rot[j][0]*dg.rot[j][0]+dg.rot[j][1]*dg.rot[j][1]+dg.rot[j][2]*dg.rot[j][2]-1)*sqrt(dg.w_rot);
		fx.at<double>(idx++) = (dg.rot[j][3]*dg.rot[j][3]+dg.rot[j][4]*dg.rot[j][4]+dg.rot[j][5]*dg.rot[j][5]-1)*sqrt(dg.w_rot);
		fx.at<double>(idx++) = (dg.rot[j][6]*dg.rot[j][6]+dg.rot[j][7]*dg.rot[j][7]+dg.rot[j][8]*dg.rot[j][8]-1)*sqrt(dg.w_rot);
	}
	//Ereg
	for(int j=0; j<dg.n_nodes; j++){
		for(int k=0; k<dg.n_nodes; k++){
			if(dg.edges[j][k]){
				for(int i=0; i<3; i++)
					fx.at<double>(idx++) = (dg.rot[j][i]*(dg.nodes[k][0]-dg.nodes[j][0])
											+dg.rot[j][i+3]*(dg.nodes[k][1]-dg.nodes[j][1])
											+dg.rot[j][i+6]*(dg.nodes[k][2]-dg.nodes[j][2])
											+dg.nodes[j][i]+dg.trans[j][i]-dg.nodes[k][i]-dg.trans[k][i])*sqrt(dg.w_reg);
			}
		}
	}
	//Econ
	for(int l=0; l<p; l++){
		dg.predict(v[l],_v);
		for(int i=0; i<3; i++) fx.at<double>(idx++) = (_v[i]-q[l][i])*sqrt(dg.w_con);
	}
	//f(x)
	fx.copyTo(f);
}


void computeJ(const DeformationGraph &dg, const int &p, double **v, double **q, OutputArray J){
	double *w = new double[dg.k_nearest];
	int idx = 0, *index = new int[dg.k_nearest+1];
	Mat Jacobi= Mat::zeros(6*dg.n_nodes+6*dg.n_edges+3*p, 12*dg.n_nodes, CV_64F);
	//Erot
	for(int j=0; j<dg.n_nodes; j++){
		Jacobi.at<double>(idx,0+12*j) = dg.rot[j][3]*sqrt(dg.w_rot);Jacobi.at<double>(idx,1+12*j) = dg.rot[j][4]*sqrt(dg.w_rot);Jacobi.at<double>(idx,2+12*j) = dg.rot[j][5]*sqrt(dg.w_rot);
		Jacobi.at<double>(idx,3+12*j) = dg.rot[j][0]*sqrt(dg.w_rot);Jacobi.at<double>(idx,4+12*j) = dg.rot[j][1]*sqrt(dg.w_rot);Jacobi.at<double>(idx++,5+12*j) = dg.rot[j][2]*sqrt(dg.w_rot);

		Jacobi.at<double>(idx,0+12*j) = dg.rot[j][6]*sqrt(dg.w_rot);Jacobi.at<double>(idx,1+12*j) = dg.rot[j][7]*sqrt(dg.w_rot);Jacobi.at<double>(idx,2+12*j) = dg.rot[j][8]*sqrt(dg.w_rot);
		Jacobi.at<double>(idx,6+12*j) = dg.rot[j][0]*sqrt(dg.w_rot);Jacobi.at<double>(idx,7+12*j) = dg.rot[j][1]*sqrt(dg.w_rot);Jacobi.at<double>(idx++,8+12*j) = dg.rot[j][2]*sqrt(dg.w_rot);

		Jacobi.at<double>(idx,3+12*j) = dg.rot[j][6]*sqrt(dg.w_rot);Jacobi.at<double>(idx,4+12*j) = dg.rot[j][7]*sqrt(dg.w_rot);Jacobi.at<double>(idx,5+12*j) = dg.rot[j][8]*sqrt(dg.w_rot);
		Jacobi.at<double>(idx,6+12*j) = dg.rot[j][3]*sqrt(dg.w_rot);Jacobi.at<double>(idx,7+12*j) = dg.rot[j][4]*sqrt(dg.w_rot);Jacobi.at<double>(idx++,8+12*j) = dg.rot[j][5]*sqrt(dg.w_rot);

		for(int i=0; i<9; i++){
			if(i == 3 || i == 6) ++idx;
			Jacobi.at<double>(idx,i+12*j) = 2*dg.rot[j][i]*sqrt(dg.w_rot);
		}
		++idx;
	}
	//Ereg
	for(int j=0; j<dg.n_nodes; j++){
		for(int k=0; k<dg.n_nodes; ++k){
			if(dg.edges[j][k]){  //k-th node is j-th node's neighbour
				for(int i = 0; i<3; i++){
					for(int ii=0; ii<3; ++ii)
						Jacobi.at<double>(idx,12*j+3*ii+i) = (dg.nodes[k][ii]-dg.nodes[j][ii])*sqrt(dg.w_reg);
					Jacobi.at<double>(idx,12*j+i+9) = sqrt(dg.w_reg);
					Jacobi.at<double>(idx++,12*k+i+9) = -sqrt(dg.w_reg);
				}
			}
		}
	}
	//Econ
	for(int l=0; l<p; l++){
		dg.computeWeights(v[l],w,index);
		for(int i=0; i<3; i++){
			for(int j=0; j<dg.k_nearest; j++){
				for(int ii=0; ii<3; ++ii) Jacobi.at<double>(idx,ii*3+i+12*index[j]) = w[j]*(v[l][ii]-dg.nodes[index[j]][ii])*sqrt(dg.w_con);
				Jacobi.at<double>(idx,9+i+12*index[j]) = w[j]*sqrt(dg.w_con);
			}
			++idx;
		}
	}
	//J
	Jacobi.copyTo(J);
	delete[] w; w = NULL;
	delete[] index; index = NULL;
}


double F(DeformationGraph &dg, Mat &x, int p, double **v, double **q){
	for(int j=0; j<dg.n_nodes; j++){
		for(int i=0; i<9; ++i) dg.rot[j][i] = x.at<double>(12*j+i);
		for(int i=0; i<3; ++i) dg.trans[j][i] = x.at<double>(12*j+i+9);
	}
	double e_rot = dg.Erot();
	double e_reg = dg.Ereg();
	double e_con = dg.Econ(p,v,q);
	std::cout<<"Erot = "<<e_rot<<'\n'
			 <<"Ereg = "<<e_reg<<'\n'
			 <<"Econ = "<<e_con<<'\n';
	return dg.w_rot*e_rot+dg.w_reg*e_reg+dg.w_con*e_con;
}


void gaussNewton(DeformationGraph &dg, int p, double **v, double **q){
	const double epsilon = 0.0000001;
	double Fk=0, Fk1;
	int iter_times = 10;
	Mat f,J,h,x(12*dg.n_nodes,1,CV_64F);
	for(int j=0; j<dg.n_nodes; j++){
		for(int i=0; i<9; i++) x.at<double>(i+j*12) = dg.rot[j][i];
		for(int i=0; i<3; i++) x.at<double>(i+9+j*12) = dg.trans[j][i];
	}
	Fk1 = F(dg,x,p,v,q);
	for(int iter_count=0; (iter_count<iter_times)&&(abs(Fk1-Fk) >= epsilon*(1+Fk1)); ++iter_count){
		std::cout<<"iteration "<<iter_count<<": E = "<<Fk1<<'\n';
		Fk = Fk1;
		computef(dg,p,v,q,f);
		computeJ(dg,p,v,q,J);
		Mat JtJ = J.t()*J;
		Mat I = Mat::eye(JtJ.size(),CV_64F);
		JtJ += 0.00000001*I;
		bool check = solve(JtJ,-J.t()*f,h,DECOMP_CHOLESKY);
		if(!check){ std::cerr<<"solve linear system of equations failed!\n"; break; }
		x += h;
		Fk1 = F(dg,x,p,v,q);
	}
	std::cout<<"Iteration finished.\n";
}