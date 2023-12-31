#include "LMPC.h"
#include <qpOASES.hpp>


using namespace Eigen;

//get return like this
//auto [Ad, Bd] = discretizeState(i, Ac, Bc);
//Ref : Bilinear transformation of continusous time state space system
//https://dsp.stackexchange.com/questions/45042/bilinear-transformation-of-continuous-time-state-space-system
auto LMPC::discretizeState(const int& state_num, const MatrixXd& Ac, const  MatrixXd& Bc)->std::tuple<MatrixXd, MatrixXd>
{
	int Acrow = Ac.rows();	int Accol = Ac.cols();
	int Bcrow = Bc.rows();	int BCcol = Bc.cols();
	
	MatrixXd I = MatrixXd::Identity(state_num, state_num);
	MatrixXd M(Acrow, Accol);
	MatrixXd Ad(Acrow, Accol);
	MatrixXd Bd(Bcrow, BCcol);
	
	M.setZero(); Ad.setZero(); Bd.setZero();

	M = (I - 0.5 * _Ts * Ac).inverse();
	Ad= (I + 0.5 * _Ts * Ac) * M;
	Bd = 0.5 * _Ts * (Ad + I) * Bc;

	return std::make_tuple(Ad, Bd);
}


auto LMPC::stackingMatrix(const int &state_num, const MatrixXd &Ad , const MatrixXd &Bd) -> std::tuple<std::vector<Eigen::MatrixXd>, std::vector<Eigen::MatrixXd>>
{
	//vector type Eigen::Matrix return
	std::vector<Eigen::MatrixXd> Ps_Matrices; // projection state (state묶음 ex- alpha velocitu alphadot)
	std::vector<Eigen::MatrixXd> Pu_Matrices; // projection input

	//Return Matrices initalize
	for (int i = 0; i < state_num; i++)
	{
		Ps_Matrices.push_back(Eigen::MatrixXd::Zero(_Dim_Preview, state_num));
		Pu_Matrices.push_back(Eigen::MatrixXd::Zero(_Dim_Preview, _Dim_Preview));
	}
	
	//unnecessary for acceleration input, necessary for torque input
	Eigen::MatrixXd Ds(state_num, state_num);
	Eigen::MatrixXd Du(Bd.rows(), _NumTarget_Input);
	Ds = Ad;
	Du = Bd;

	//stacking projection matrix
	for (int i = 0; i < _Dim_Preview; i++)
	{
		
		for (int n = 0; n < state_num; n++)
		{	
			//stacking Ps matrix
			Ps_Matrices[n].row(i) = Ds.row(n);

			//stacking Pu matrix
			for (int j = i; j < _Dim_Preview; j++)
			{
				int p = _NumTarget_Input * (j - i);

				for (int k = 0; k < _NumTarget_Input; k++)
				{
					Pu_Matrices[n](j,p+k) = Du.row(n)[k];
				}
			}
		}
		Ds = Ad * Ds;
		Du = Ad * Du;
	}
	return std::make_tuple(Ps_Matrices, Pu_Matrices);
}


auto LMPC::assemble_HMatrix(const int& state_num, std::vector<Eigen::MatrixXd>& Pus, const Eigen::VectorXd& Weight_vec) -> Eigen::MatrixXd
{
	//Np x Np
	Eigen::MatrixXd H(_Dim_Preview, _Dim_Preview); 
	Eigen::MatrixXd E_MAT = Eigen::MatrixXd::Identity(_Dim_Preview, _Dim_Preview);
	H.setZero();
	
	for (int i = 0; i < (state_num + 1) ; i++)
	{
		//input weight
		if (i == state_num)
		{
			H += Weight_vec[i] * E_MAT;
			break;
		}
		
		H  += Weight_vec[i] * (Pus[i].transpose() * Pus[i]);

	}
	return std::move(H);
}


auto LMPC::assemble_gMatrix(const int& state_num, const std::vector<Eigen::MatrixXd> &Pss, std::vector<Eigen::MatrixXd>& Pus,
	Eigen::VectorXd& state_vec, const Eigen::VectorXd& Weight_vec, const std::vector<Eigen::VectorXd*> &refs) -> Eigen::VectorXd
{
	Eigen::VectorXd g(_Dim_Preview);
	g.setZero();
	int count{ 0 };
		
	for(int i = 0 ; i < state_num ; i++)
	{
		g += Weight_vec(i) * (((Pss[i] * state_vec - *refs[i]).transpose()) * Pus[i]);
	}
	return g;
}


auto LMPC::constraint_AMatrix() ->Eigen::MatrixXd
{

	//기존 : lower 과 upper를 하나로 묶음
	//Eigen::MatrixXd C_input(2 * Dim, Dim);
	//C_input.block(0, 0, Dim, Dim) = Eigen::MatrixXd::Identity(Dim, Dim);
	//C_input.block(Dim, 0, Dim, Dim) = -Eigen::MatrixXd::Identity(Dim, Dim);

	//수정 :lower과 upper를 따로 : (2 * Dim 할 이유 없음)
	Eigen::MatrixXd C_input = Eigen::MatrixXd::Identity(_Dim_Preview, _Dim_Preview);

	return C_input;
}


auto LMPC::Bound_AMatrix(const int& upper, const int& lower) ->std::tuple<Eigen::VectorXd, Eigen::VectorXd>
{
	//기존
	//Eigen::VectorXd Abound(2*_Dim_Preview);
	//Abound.segment(0, _Dim_Preview).setConstant(upper);
	//Abound.segment(_Dim_Preview, _Dim_Preview).setConstant(lower);

	//수정
	Eigen::VectorXd A_upbound(_Dim_Preview);
	Eigen::VectorXd A_lwbound(_Dim_Preview);

	A_upbound.setConstant(upper);
	A_lwbound.setConstant(lower);

	return std::make_tuple(A_upbound, A_lwbound);
}

void LMPC::Get_Reference(const unsigned int& cur_time_stp, const vector<double>& from, Eigen::VectorXd& dest)
{
	// Assuming dest has been resized to _Dim_Preview before this function call.
	if (cur_time_stp + _Dim_Preview < _Dim_Total) {
		
		std::copy(from.begin() + cur_time_stp + 1, from.begin() + cur_time_stp + 1 + _Dim_Preview, dest.data());
	}
}

//col major Eigen to row major real_t of qpOASES
//Eigen: colmajor , qpOASES : rowmajor
qpOASES::real_t* LMPC::Convert2RealT(const Eigen::MatrixXd& mat)
{
	qpOASES::real_t* qp = new qpOASES::real_t[_Dim_Preview * _Dim_Preview];
	
	auto matT = mat.transpose();

	for (int i = 0; i < mat.cols() * mat.rows(); ++i)
	{
		qp[i] = matT(i);
	}

	return qp;
}

qpOASES::real_t* LMPC::ConvertEigenToRealT(const Eigen::MatrixXd& mat)
{
	if (mat.size() == 0) return nullptr;

	// Allocate memory for the new array
	auto realTArray = std::make_unique<qpOASES::real_t[]>(mat.size());

	// Copy data from the Eigen matrix to the real_t array in row-major order
	for (int i = 0; i < mat.rows(); ++i) {
		for (int j = 0; j < mat.cols(); ++j) {
			realTArray[i * mat.cols() + j] = static_cast<qpOASES::real_t>(mat(i, j));
		}
	}
	// Return the pointer to the new array
	return realTArray.release();
}