#pragma once
#include <iostream>
#include <chrono>
#include <ctime>
#include <math.h>
#include <time.h>
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

class LMPC
{
protected:
	
	/////////////////////////////MPC���� ���///////////////////////////////////////////////////
	static constexpr double _dT = 0.030;			// Sampling Time [s]
	static constexpr double _Time = 10.0;			// Total Time (�� �۵��ð�)
	static constexpr double _Time_Preview = 1.5;	// Horizon [s]

	int _NumTarget_Input = 1;						//system�� �Էº��� u
	const double _Ts = _dT;							// sampling time ex) 0.01 = 100Hz
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	static constexpr int _Dim_Time = int(_Time / _dT);					// ��ü Time Dimension
	static constexpr int _Dim_Preview = int(_Time_Preview / _dT);		// Horizon Dimension(Np)
	static constexpr int _Dim_Total = _Dim_Time + _Dim_Preview ;		// ��ü Time Dimension (������ Horizon Dimension�� MPC ����� ����)
	//////////////////////////////////////////////////////////////////////////////////////

public:
	LMPC()
	{}
	~LMPC()
	{}

	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrixXd;

	// get return like this
	// auto [Ad, Bd] = discretizeState(Ac, Bc, i);
	auto discretizeState(const int& state_num, const Eigen::MatrixXd& Ac, const Eigen::MatrixXd& Bc) -> std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>;

	// return multiple Ps & Pu at once
	// auto [Pss, Pus] = stackingMatrix(const int &stateA_num, const MatrixXd &Ad , const MatrixXd &Bd); 
	// auto Ps_alpha = Pss[0], auto Ps_beta = Pss[1]....
	auto stackingMatrix(const int &state_num, const Eigen::MatrixXd &Ad, const Eigen::MatrixXd &Bd) -> std::tuple<std::vector<Eigen::MatrixXd>, std::vector<Eigen::MatrixXd>>;

	//use output of Pss & Pus from auto stackingMatrix.
	//Weighting_vec must included Weighting of INPUT at LAST
	auto assemble_HMatrix(const int &state_num , std::vector<Eigen::MatrixXd> &Pus , const Eigen::VectorXd &Weight_vec) -> Eigen::MatrixXd;
	
	//grdient matirx ����
	//reference matrix�� �޾ƿͼ� g ���
	//state_vec�� �ڿ����� ��
	//�ϴ� �� �־���� -> �� �� �ִ� �� ����
	auto assemble_gMatrix(const int& state_num, const std::vector<Eigen::MatrixXd> &Pss, std::vector<Eigen::MatrixXd>& Pus,
		Eigen::VectorXd& state_vec, const Eigen::VectorXd& Weight_vec, const std::vector<Eigen::VectorXd*> &refs) -> Eigen::VectorXd;

	//A : constraint matrix (for TWIP)
	auto constraint_AMatrix() -> Eigen::MatrixXd;

	//A matrix�� constraint ����
	auto Bound_AMatrix(const int& upper, const int& lower) -> std::tuple<Eigen::VectorXd, Eigen::VectorXd>;

	//Change from Eigen to QP real_t
	qpOASES::real_t* Convert2RealT(const Eigen::MatrixXd& mat);
	qpOASES::real_t* Convert2RealT2(const Eigen::MatrixXd& mat);

	//////////////////////////////////////////////getter & setter///////////////////////////////////////////////////	
	int getNumTarget_input() { return _NumTarget_Input; }
	void setNumTarget_input(int& num) { _NumTarget_Input = num; }

	double getSamplingTime()		  { return _dT; }
	double getOperationTime()		  { return _Time; }
	double getPreviewTime()			  { return _Time_Preview; }
	int getHorizonDim() const		  { return _Dim_Preview; }
	const int getTotalHorizon() const { return _Dim_Total; }
	int getDimTIME() const			  { return _Dim_Time; }
	
};


