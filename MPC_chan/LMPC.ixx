module;
#include <iostream>
#include <thread>
#include <future>
#include <qpOASES.hpp>
#include <tuple>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <initializer_list>
#include <assert.h>
#include <thread>
#include <future>
#include "CommUtil.h"
#include "Robotics_func.h"

export module LMPC;
import std.core;
	/////////////////////////////MPC관련 상수///////////////////////////////////////////////////
	constexpr double _dT = 0.030;			// Sampling Time [s]
	constexpr double _Time = 10.0;			// Total Time (총 작동시간)
	constexpr double _Time_Preview = 1.5;	// Horizon [s]

	int _NumTarget_Input = 1;						//system의 입력변수 u
	const double _Ts = _dT;							// sampling time ex) 0.01 = 100Hz

	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrixXd;

export namespace LMPC
{
	///////////////////////////////////////////////////////////////////////////////////////////////
	constexpr int _Dim_Time = int(_Time / _dT);					// 전체 Time Dimension
	constexpr int _Dim_Preview = int(_Time_Preview / _dT);		// Horizon Dimension(Np)
	constexpr int _Dim_Total = _Dim_Time + _Dim_Preview;		// 전체 Time Dimension (마지막 Horizon Dimension의 MPC 계산을 위해)
	//////////////////////////////////////////////////////////////////////////////////////
	
	// get return like this
	// auto [Ad, Bd] = discretizeState(Ac, Bc, i);
	auto discretizeState(const int& state_num, const MatrixXd& Ac, const  MatrixXd& Bc)->std::tuple<MatrixXd, MatrixXd>
	{
		int Acrow = Ac.rows();	int Accol = Ac.cols();
		int Bcrow = Bc.rows();	int BCcol = Bc.cols();

		MatrixXd I = MatrixXd::Identity(state_num, state_num);
		MatrixXd M(Acrow, Accol);
		MatrixXd Ad(Acrow, Accol);
		MatrixXd Bd(Bcrow, BCcol);

		M.setZero(); Ad.setZero(); Bd.setZero();

		M = (I - 0.5 * _Ts * Ac).inverse();
		Ad = (I + 0.5 * _Ts * Ac) * M;
		Bd = 0.5 * _Ts * (Ad + I) * Bc;

		return std::make_tuple(Ad, Bd);
	}

	// return multiple Ps & Pu at once
	// auto [Pss, Pus] = stackingMatrix(const int &stateA_num, const MatrixXd &Ad , const MatrixXd &Bd); 
	// auto Ps_alpha = Pss[0], auto Ps_beta = Pss[1]....
	auto stackingMatrix(const int& state_num, const Eigen::MatrixXd& Ad, const Eigen::MatrixXd& Bd) -> std::tuple<std::vector<Eigen::MatrixXd>, std::vector<Eigen::MatrixXd>>
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
						Pu_Matrices[n](j, p + k) = Du.row(n)[k];
					}
				}
			}
			Ds = Ad * Ds;
			Du = Ad * Du;
		}
		return std::make_tuple(Ps_Matrices, Pu_Matrices);
	}
	
	//use output of Pss & Pus from auto stackingMatrix.
	//Weighting_vec must included Weighting of INPUT at LAST
	auto assemble_HMatrix(const int& state_num, std::vector<Eigen::MatrixXd>& Pus, const Eigen::VectorXd& Weight_vec) -> Eigen::MatrixXd
	{
		//Np x Np
		Eigen::MatrixXd H(_Dim_Preview, _Dim_Preview);
		Eigen::MatrixXd E_MAT = Eigen::MatrixXd::Identity(_Dim_Preview, _Dim_Preview);
		H.setZero();

		for (int i = 0; i < state_num + 1; i++)
		{
			//input weight
			if (i == state_num)
			{
				H += Weight_vec.coeff(Weight_vec.size() - 1) * E_MAT;
				break;
			}

			H += Weight_vec[0] * (Pus[0].transpose()) * Pus[0];

		}
		return std::move(H);
	}

	//grdient matirx 설정
	//reference matrix를 받아와서 g 계산
	//state_vec은 뒤에서도 씀
	//일단 다 넣어놨음 -> 뺄 수 있는 것 빼기
	auto assemble_gMatrix(const int& state_num, const std::vector<Eigen::MatrixXd>& Pss, std::vector<Eigen::MatrixXd>& Pus,
		Eigen::VectorXd& state_vec, const Eigen::VectorXd& Weight_vec, const std::vector<Eigen::VectorXd*>& refs) -> Eigen::VectorXd
	{
		Eigen::VectorXd g(_Dim_Preview);
		g.setZero();
		int count = 0;

		//if (count == 0)
		//{
		//	for (int i = 0; i < 1; i++)
		//	{
		//		cout << "!!" << endl;
		//		//문제없음
		//		//cout << Weight_vec[i] << endl; 
		//		//cout << Weight_vec.size() << endl;
		//		cout << "Pss.size() : " << Pss.size() << endl;
		//		cout <<"Pss[i].size() : " << Pss[i].size() << endl;
		//		/*cout << Pss[0] << endl;*/
		//		cout << "Pus : " << Pus.size() << endl;
		//		cout << "Pus[0] : " << Pus[0].size() << endl;
		//		
		//		cout << " refs[i]->size() :  " << refs[i]->size() << endl;
		//		cout << " refs[i] :  " << refs[i] << endl;
		//		//cout << " *refs[i] :  " << *refs[i] << endl;
		//		//cout << "A : " << A.size() << endl;
		//		//cout << "A[0] : " << A[0].size() << endl;
		//		//cout << A[0] << endl;
		//		//cout << Pss[state_num] << endl << endl;
		//	}
		//	count = 1;
		//}

		for (int i = 0; i < state_num; i++)
		{
			g += Weight_vec[i] * ((Pss[i] * state_vec - *refs[i]).transpose()) * Pus[i];
		}
		return std::move(g.transpose());
	}

	//A : constraint matrix (for TWIP)
	auto constraint_AMatrix() ->Eigen::MatrixXd
	{

		//기존 : lower 과 upper를 하나로 묶음
		//Eigen::MatrixXd C_input(2 * Dim, Dim);
		//C_input.block(0, 0, Dim, Dim) = Eigen::MatrixXd::Identity(Dim, Dim);
		//C_input.block(Dim, 0, Dim, Dim) = -Eigen::MatrixXd::Identity(Dim, Dim);

		//수정 :lower과 upper를 따로 : (2 * Dim 할 이유 없음)
		Eigen::MatrixXd C_input = Eigen::MatrixXd::Identity(_Dim_Preview, _Dim_Preview);


		return C_input;
	}

	//A matrix의 constraint 설정
	auto Bound_AMatrix(const int& upper, const int& lower) ->std::tuple<Eigen::VectorXd, Eigen::VectorXd>
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

		return std::tuple(A_upbound, A_lwbound);
	}

	//col major Eigen to row major real_t of qpOASES
	qpOASES::real_t* Convert2RealT(const Eigen::MatrixXd& mat)
	{
		// Create a row-major map of the matrix
		Eigen::Map<const RowMajorMatrixXd> rowMajorMap(mat.data(), mat.rows(), mat.cols());

		// Return the internal data pointer of the map
		return const_cast<qpOASES::real_t*>(rowMajorMap.data());
	}

	//////////////////////////////////////////////getter & setter///////////////////////////////////////////////////	
	int getNumTarget_input() { return _NumTarget_Input; }
	void setNumTarget_input(int& num) { _NumTarget_Input = num; }

	double getSamplingTime() { return _dT; }
	double getOperationTime() { return _Time; }
	double getPreviewTime() { return _Time_Preview; }
	int getHorizonDim()  { return _Dim_Preview; }
	const int getTotalHorizon()  { return _Dim_Total; }
	int getDimTIME()  { return _Dim_Time; }
}

