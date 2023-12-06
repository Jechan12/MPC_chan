#pragma once
#include "LMPC.h"
#include "Traj_MPC.h"



//MPC in TWIP with PFL
class MPC_TWIP :public LMPC, public Traj_MPC
{

	////////////////////////////////MPC 가중치 in TWIP//////////////////////////////////////////////
	double W_phi   { 100.0 };
	double W_phidot{ 20.0 };
	double W_yaw   { 0.001 };
	
	double W_alpha   { 1000.0 };
	double W_velocity{ 800.0 };
	double W_alphadot{ 200.0 };
	double W_pitch   { 0.0001 };  //pdf 상 smooth input 

	double W_beta    { 1000.0 };
	double W_betadot { 1.0 };
	double W_roll    { 0.01 };

	///////////////////////////////TWIP constant/////////////////////////////////////////////////////////
	const double l_G{ 0.69 };								//	CoM height
	const double r_W{ 0.075 };							//  Wheel radius
	const double W_b{ 0.215 };							//  Wheel Track / 2
	const double m_w{ 1.5 };								//  Wheel mass
	const double J_wa{ 0.004219 };
	const double J_wd{ 0.002222 };

	const double m_p{ 40.0 };								//	Pendulum mass
	const double I_Px{ 6.48 };
	const double I_Py{ 7.01 };
	const double I_Pz{ 1.81 };

	///////////////////////////////FOR MPC THREAD/////////////////////////////////////////////////////////
	Uint m_TimerID;
	Uint mmTimePeriod; //[ms]


public:
	///////////////////////////////TWIP system DOF/////////////////////////////////////////////////////////
	const int DOF_P{ 3 };									// x_P = ( alpha, v_M, alphadot )
	const int DOF_Y{ 2 };									// x_Y = ( phi_M, phidot_M )
	const int DOF_R{ 2 };									// x_R = ( beta_M, betadot_M )

	//////////////////////////////////constants/////////////////////////////////////////////
	const double g_const = 9.81;
	const int upperbound = 100;
	const int lowerbound = -100;

	double Trgt_a_yaw = 0.0;
	double Trgt_a_pitch = 0.0;
	double Trgt_a_roll = 0.0;

	//////////////////////////matirx of continous state/////////////////////////////////
	Eigen::MatrixXd Ac_pitch = Eigen::MatrixXd::Zero(DOF_P, DOF_P);
	Eigen::Vector3d Bc_pitch = Eigen::MatrixXd::Zero(DOF_P, 1);
	Eigen::MatrixXd Ac_yaw   = Eigen::MatrixXd::Zero(DOF_Y, DOF_Y);
	Eigen::VectorXd Bc_yaw   = Eigen::MatrixXd::Zero(DOF_Y, 1);
	Eigen::MatrixXd Ac_roll  = Eigen::MatrixXd::Zero(DOF_R, DOF_R);
	Eigen::VectorXd Bc_roll  = Eigen::MatrixXd::Zero(DOF_R, 1);

	//......................................................................................................//
	std::vector<Eigen::MatrixXd> Pss_yaw;   std::vector<Eigen::MatrixXd> Pus_yaw;
	std::vector<Eigen::MatrixXd> Pss_pitch; std::vector<Eigen::MatrixXd> Pus_pitch;
	std::vector<Eigen::MatrixXd> Pss_roll;  std::vector<Eigen::MatrixXd> Pus_roll;

	Eigen::MatrixXd H_yaw;   Eigen::MatrixXd cA_yaw;   Eigen::VectorXd Aub_yaw;   Eigen::VectorXd Alb_yaw;
	Eigen::MatrixXd H_pitch; Eigen::MatrixXd cA_pitch; Eigen::VectorXd Aub_pitch; Eigen::VectorXd Alb_pitch;
	Eigen::MatrixXd H_roll;  Eigen::MatrixXd cA_roll;  Eigen::VectorXd Aub_roll;  Eigen::VectorXd Alb_roll;

	qpOASES::real_t* H_yaw_qp;	 qpOASES::real_t* cA_yaw_qp;	qpOASES::real_t* Aub_yaw_qp;	qpOASES::real_t* Alb_yaw_qp;
	qpOASES::real_t* H_pitch_qp; qpOASES::real_t* cA_pitch_qp;	qpOASES::real_t* Aub_pitch_qp;	qpOASES::real_t* Alb_pitch_qp;
	qpOASES::real_t* H_roll_qp;	 qpOASES::real_t* cA_roll_qp;	qpOASES::real_t* Aub_roll_qp;	qpOASES::real_t* Alb_roll_qp;

	qpOASES::real_t* g_yaw_qp;  qpOASES::real_t* g_pitch_qp;
	Eigen::VectorXd g_pitch;

	


	///////////// Vector of Weighting //////////////////////////////
	Eigen::VectorXd W_YAW = Eigen::VectorXd::Zero(3);
	Eigen::VectorXd W_PITCH = Eigen::VectorXd::Zero(4);
	Eigen::VectorXd W_ROLL = Eigen::VectorXd::Zero(3);

	
	////////////////////state for first order term of QP//////////////////////////////////
	Eigen::VectorXd x_Yaw = Eigen::VectorXd::Zero(DOF_Y);   //2 (phi, phidot)
	Eigen::VectorXd x_Pitch = Eigen::VectorXd::Zero(DOF_P); //3 (alpha, velocity ,alphadot)
	Eigen::VectorXd x_Roll = Eigen::VectorXd::Zero(DOF_R);  //2 (beta, betadot)

	////////////////////////reference for first order term of QP(얘도 Traj_MPC로?)////////////////////////////
	Eigen::VectorXd Preview_PHI_ref, Preview_PHIDOT_ref, Preview_ALPHA_ref, Preview_VELOCITY_ref, Preview_ALPHADOT_ref, Preview_BETA_ref, Preview_BETADOT_ref;
	
	//Optimal input 
	qpOASES::real_t Yaw_Opt[_Dim_Preview]{0.0}, Pitch_Opt[_Dim_Preview]{ 0.0 }, Roll_Opt[_Dim_Preview]{0.0};

	std::chrono::steady_clock::time_point AfterInput;
	std::chrono::steady_clock::time_point lastUpdate;
	std::chrono::milliseconds minInterval;

	////////////////////// From mother class for convenience /////////////////////////////////
	int Dim;		//Np
	int DimTotal;
	int Dim_TIME;	//작동시간 (Sampling time)
	double dT;
	double OperationTime;
	double StartTime;
	int global_indx = 0;
	double main_time = 0.0;

	MPC_TWIP()
	{
		Dim       = LMPC::getHorizonDim();		//_Dim_Preview (N_p)
		DimTotal  = LMPC::getTotalHorizon();	// 전체 Time Dimension (마지막 Horizon Dimension의 MPC 계산을 위해)
		dT        = LMPC::getSamplingTime() ;	//_dT
		Dim_TIME  = LMPC::getDimTIME();			//_Dim_Time // 전체 Time Dimension
		StartTime = LMPC::getStratTime();
		OperationTime = LMPC::getOperationTime();
		mmTimePeriod = LMPC::getSamplingTime() * 1000; //[s -> ms]
	}
	
	~MPC_TWIP()
	{
		/////	Stop and destroy the multimedia timer
		//if (m_TimerID != 0) {
		//	timeKillEvent(m_TimerID);
		//}

		delete[] H_yaw_qp, cA_yaw_qp , Aub_yaw_qp , Alb_yaw_qp;
		delete[] H_pitch_qp , cA_pitch_qp , Aub_pitch_qp , Alb_pitch_qp;
		delete[] H_roll_qp , cA_roll_qp , Aub_roll_qp , Alb_roll_qp;
		cout << "MPC_TWIP destrucor" << endl;
		cout << "global index : " << global_indx << endl;
	}

	Traj_MPC traj;

	//g , A , lower and upper bound depends on system.  -> set after inheritance
	void setStateMatrix();

	void Weight2vec();

	void initializeRefVec(int& size);
	
	//Get_Reference 모음
	void solveRef();
	
	//루프돌기 전
	void OutLoopSetYaw();
	void OutLoopSetPitch();
	void OutLoopSetRoll();

	//루프 내부 QP계산	
	void InLoopSolvYaw();
	void InLoopSolvPitch();
	void InLoopSolvRoll();

	//QP계산용 루프열기
	void OpenMPCThread();


	void setTIME(double &steptime);

};
