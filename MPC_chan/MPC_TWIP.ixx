module;
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
export module MPC_TWIP;

import LMPC;

////////////////////////////////MPC 가중치 in TWIP//////////////////////////////////////////////
const double W_phi{ 1000.0 };
const double W_phidot{ 100.0 };
const double W_yaw{ 0.01 };

const double W_alpha{ 400.0 };
const double W_velocity{ 100.0 };
const double W_alphadot{ 1.0 };
const double W_pitch{ 0.0001 };  //pdf 상 smooth input 

const double W_beta{ 1000.0 };
const double W_betadot{ 1.0 };
const double W_roll{ 0.01 };

////////////////////// From mother class for convenience /////////////////////////////////
int Dim;		//Np
int DimTotal;
int Dim_TIME;	//작동시간 (Sampling time)
double dT;



Eigen::Vector<double, 3> W_YAW;
Eigen::Vector<double, 4> W_PITCH;
Eigen::Vector<double, 3> W_ROLL;

///////////////////////////////TWIP constant/////////////////////////////////////////////////////////
const int DOF_P{ 3 };									// x_P = ( alpha, v_M, alphadot )
const int DOF_Y{ 2 };									// x_Y = ( phi_M, phidot_M )
const int DOF_R{ 2 };

const double l_G{ 0.69 };								//	CoM height
const double r_W{ 0.075 };								//  Wheel radius
const double W_b{ 0.215 };								//  Wheel Track / 2
const double m_w{ 1.5 };								//  Wheel mass
const double J_wa{ 0.004219 };
const double J_wd{ 0.002222 };

const double m_p{ 40.0 };								//	Pendulum mass
const double I_Px{ 6.48 };
const double I_Py{ 7.01 };
const double I_Pz{ 1.81 };

//////////////////////////////////constants/////////////////////////////////////////////
const double g_const = 9.81;
const int upperbound = 100;
const int lowerbound = -100;
int global_indx = 0;
double main_time = 0.0;

double Trgt_a_yaw{ 0.0 };
double Trgt_a_pitch{ 0.0 };
double Trgt_a_roll{ 0.0 };

//////////////////////////matirx of continous state/////////////////////////////////
Eigen::MatrixXd Ac_pitch = Eigen::MatrixXd::Zero(DOF_P, DOF_P);
Eigen::Vector3d Bc_pitch = Eigen::MatrixXd::Zero(DOF_P, 1);
Eigen::MatrixXd Ac_yaw = Eigen::MatrixXd::Zero(DOF_Y, DOF_Y);
Eigen::VectorXd Bc_yaw = Eigen::MatrixXd::Zero(DOF_Y, 1);
Eigen::MatrixXd Ac_roll = Eigen::MatrixXd::Zero(DOF_R, DOF_R);
Eigen::VectorXd Bc_roll = Eigen::MatrixXd::Zero(DOF_R, 1);

////////////////////state for first order term of QP//////////////////////////////////

Eigen::VectorXd x_Yaw = Eigen::VectorXd::Zero(DOF_Y);
Eigen::VectorXd x_Pitch = Eigen::VectorXd::Zero(DOF_P);
Eigen::VectorXd x_Roll = Eigen::VectorXd::Zero(DOF_R);

////////////////////////reference for first order term of QP////////////////////////////
Eigen::VectorXd Preview_PHI_ref, Preview_PHIDOT_ref, Preview_ALPHA_ref, Preview_VELOCITY_ref, Preview_ALPHADOT_ref, Preview_BETA_ref, Preview_BETADOT_ref;
qpOASES::real_t Yaw_Opt[LMPC::_Dim_Preview]{ 0.0 }, Pitch_Opt[LMPC::_Dim_Preview]{ 0.0 }, Roll_Opt[LMPC::_Dim_Preview]{ 0.0 };

export namespace MPC_TWIP
{
	//생성자가 없네...
	void someFunction() {
		Dim = LMPC::getHorizonDim(); 
		DimTotal = LMPC::getTotalHorizon();
		dT = LMPC::getSamplingTime();
		Dim_TIME = LMPC::getDimTIME()
	}
	

	//continuous ver.
	void setStateMatrix()
	{
		std::cout << "set state matrix" << endl;

		double p1 = m_w + J_wa / SQR(r_W);
		double mu2 = I_Py + m_p * l_G * (l_G + r_W);
		double lambda2 = g_const * m_p * l_G / mu2;
		double gamma2 = -(m_p * l_G + r_W * (m_p + 2 * p1)) / mu2;

		///Pitch
		Ac_pitch(0, 2) = 1;
		Ac_pitch(2, 0) = lambda2;

		Bc_pitch(1, 0) = 1;
		Bc_pitch(2, 0) = gamma2;

		///Yaw
		Ac_yaw(0, 0) = 0.0;  Ac_yaw(0, 1) = 1.0;
		Ac_yaw(1, 0) = 0.0;  Ac_yaw(1, 1) = 0.0;

		Bc_yaw[0] = 0.0;	 Bc_yaw[1] = 1.0;

		///Roll
		Ac_roll(0, 0) = 0.0;  Ac_roll(0, 1) = 1.0;
		Ac_roll(1, 0) = 0.0;  Ac_roll(0, 1) = 0.0;

		Bc_roll[0] = 0.0;	  Bc_roll[1] = 1.0;
	}



}