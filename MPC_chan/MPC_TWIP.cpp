#pragma once
#include "MPC_TWIP.h"


//continuous ver.
void MPC_TWIP::setStateMatrix()
{
	std::cout << "Set State Matrix" << endl;
	
	double p1 = m_w + J_wa / SQR<const double>(r_W);
	double mu2 = I_Py + m_p * l_G * (l_G + r_W);  //28.124
	double lambda2 = g_const * m_p * l_G / mu2;   //9.627222301
	double gamma2 = -(m_p * l_G + r_W * (m_p + 2 * p1)) / mu2; 

	////////////////  Pitch  ///////////////////	
	Ac_pitch(0, 2) = 1;
	Ac_pitch(2, 0) = lambda2;

	Bc_pitch[1] = 1.0;
	Bc_pitch[2] = gamma2;

	////////////////  Yaw  ///////////////////	
	Ac_yaw(0, 0) = 0.0;  Ac_yaw(0, 1) = 1.0;
	Ac_yaw(1, 0) = 0.0;  Ac_yaw(1, 1) = 0.0;

	Bc_yaw[0] = 0.0;	 Bc_yaw[1] = 1.0;
	
	////////////////  Roll  ///////////////////		  
	Ac_roll(0, 0) = 0.0;  Ac_roll(0, 1) = 1.0;
	Ac_roll(1, 0) = 0.0;  Ac_roll(0, 1) = 0.0;

	Bc_roll[0] = 0.0;	  Bc_roll[1] = 1.0;

}

void MPC_TWIP::Weight2vec()
{
	cout << "W2V" << endl;

	W_YAW << W_phi, W_phidot, W_yaw;

	W_PITCH << W_alpha, W_velocity, W_alphadot, W_pitch;

	W_ROLL << W_beta, W_betadot, W_roll;
}

void MPC_TWIP::initializeRefVec(int& size)
{
	cout << "Ref init" << endl;

	Preview_PHI_ref.resize(size); Preview_PHIDOT_ref.resize(size); Preview_ALPHA_ref.resize(size);
	Preview_VELOCITY_ref.resize(size); Preview_ALPHADOT_ref.resize(size); 
	Preview_BETA_ref.resize(size); Preview_BETADOT_ref.resize(size);
	
	Preview_PHI_ref.setZero(); Preview_PHIDOT_ref.setZero();
	Preview_ALPHA_ref.setZero(); Preview_VELOCITY_ref.setZero(); Preview_ALPHADOT_ref.setZero();
	Preview_BETA_ref.setZero(); Preview_BETADOT_ref.setZero();
}




//루프 돌기 전 세팅
void MPC_TWIP::OutLoopSetYaw()
{
	//discretize
	auto [Ad_yaw, Bd_yaw] = discretizeState(DOF_Y, Ac_yaw, Bc_yaw);
	//stacking prediction matrix 
	std::tie(Pss_yaw, Pus_yaw) = stackingMatrix(DOF_Y, Ad_yaw, Bd_yaw);

	//get H matirx for QP
	H_yaw = assemble_HMatrix(DOF_Y, Pus_yaw, W_YAW);

	//constraint A matrix for QP
	cA_yaw = constraint_AMatrix();

	//constraint A bound for QP
	std::tie(Aub_yaw, Alb_yaw) = Bound_AMatrix(upperbound, lowerbound);

	//convert Eigen to qpOASES
	H_yaw_qp   = Convert2RealT(H_yaw);
	cA_yaw_qp  = Convert2RealT(cA_yaw);
	Aub_yaw_qp = Convert2RealT(Aub_yaw);
	Alb_yaw_qp = Convert2RealT(Alb_yaw);
}


void MPC_TWIP::OutLoopSetPitch()
{
	//discretize
	auto [Ad_pitch, Bd_pitch] = discretizeState(DOF_P, Ac_pitch, Bc_pitch);
	//stacking prediction matrix 
	std::tie(Pss_pitch, Pus_pitch) = stackingMatrix(DOF_P, Ad_pitch, Bd_pitch);
	
	//get H matirx for QP
	H_pitch = assemble_HMatrix(DOF_P, Pus_pitch, W_PITCH);
	//constraint A matrix for QP
	cA_pitch = constraint_AMatrix();
	//constraint A bound for QP
	std::tie(Aub_pitch, Alb_pitch) = Bound_AMatrix(upperbound, lowerbound);

	//convert Eigen to qpOASES
	H_pitch_qp = Convert2RealT(H_pitch);
	cA_pitch_qp = Convert2RealT(cA_pitch);
	Aub_pitch_qp = Convert2RealT(Aub_pitch);
	Alb_pitch_qp = Convert2RealT(Alb_pitch);
}

void MPC_TWIP::OutLoopSetRoll()
{
	//discretize
	auto [Ad_roll, Bd_roll] = discretizeState(DOF_R, Ac_roll,Bc_roll);
	//stacking prediction matrix 
	std::tie(Pss_roll, Pus_roll) = stackingMatrix(DOF_R, Ad_roll, Bd_roll);

	//get H matirx for QP
	H_roll = assemble_HMatrix(DOF_R, Pus_roll, W_ROLL);
	//constraint A matrix for QP
	cA_roll = constraint_AMatrix();
	//constraint A bound for QP
	std::tie(Aub_roll, Alb_roll) = Bound_AMatrix(upperbound, lowerbound);

	//convert Eigen to qpOASES
	H_roll_qp   = Convert2RealT(H_roll);
	cA_roll_qp  = Convert2RealT(cA_roll);
	Aub_roll_qp = Convert2RealT(Aub_roll);
	Alb_roll_qp = Convert2RealT(Alb_roll);
}

void MPC_TWIP::InLoopSolvYaw()
{
	Eigen::VectorXd g_yaw = assemble_gMatrix(DOF_Y, Pss_yaw, Pus_yaw, x_Yaw, W_YAW, { &(Preview_PHI_ref),  &Preview_PHIDOT_ref });
	g_yaw_qp = ConvertEigenToRealT(g_yaw);

	//.................................................................................//
	qpOASES::int_t nWSR = 1000;
	qpOASES::SQProblem yawExample(MPC_TWIP::Dim, 1);
	yawExample.init(H_yaw_qp, g_yaw_qp, cA_yaw_qp, 0, 0, Alb_yaw_qp, Aub_yaw_qp, nWSR, 0);
	yawExample.getPrimalSolution(Yaw_Opt);
	//.................................................................................//

	//select first acc
	Trgt_a_yaw = Yaw_Opt[0];
	

	delete[] g_yaw_qp;
}

void MPC_TWIP::InLoopSolvPitch()
{
	Eigen::VectorXd g_pitch = assemble_gMatrix(MPC_TWIP::DOF_P, Pss_pitch, Pus_pitch, MPC_TWIP::x_Pitch, W_PITCH,{ &(Preview_ALPHA_ref),&(Preview_VELOCITY_ref),&(Preview_ALPHADOT_ref) });
	g_pitch_qp = ConvertEigenToRealT(g_pitch.transpose());

	//.................................................................................//
	qpOASES::int_t nWSR = 1000;
	qpOASES::SQProblem pitchExample(MPC_TWIP::Dim, MPC_TWIP::Dim);
	pitchExample.init(H_pitch_qp, g_pitch_qp, cA_pitch_qp, Alb_pitch_qp, Aub_pitch_qp, 0, 0, nWSR, 0);
	pitchExample.getPrimalSolution(Pitch_Opt);
	//.................................................................................//

	//select first acc
	Trgt_a_pitch = Pitch_Opt[0]; 
	
	delete[] g_pitch_qp;
}

void MPC_TWIP::InLoopSolvRoll()
{
	
	Eigen::VectorXd g_roll = assemble_gMatrix(MPC_TWIP::DOF_R, Pss_roll, Pus_roll, x_Roll, W_ROLL, { &(Preview_BETA_ref), &(Preview_BETADOT_ref) });
	//qpOASES::real_t* g_roll_qp = Convert2RealT(g_roll);
	qpOASES::real_t* g_roll_qp = ConvertEigenToRealT(g_roll);
	qpOASES::int_t nWSR = 200;
	qpOASES::SQProblem rollExample(MPC_TWIP::Dim, MPC_TWIP::Dim);
	//qp푸는 함수 추가

	delete[] g_roll_qp;
}


