#pragma once
#include "MPC_TWIP.h"

//continuous ver.
void MPC_TWIP::setStateMatrix()
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

void MPC_TWIP::Weight2vec()
{

	W_YAW.setZero();  W_PITCH.setZero();   W_ROLL.setZero();
	
	W_YAW << W_phi, W_phidot, W_yaw;

	W_PITCH << W_alpha, W_velocity, W_alphadot, W_pitch;

	W_ROLL << W_beta, W_betadot, W_roll;
}

void MPC_TWIP::initializeRefVec(int& size)
{
	Preview_PHI_ref.resize(size), Preview_PHIDOT_ref.resize(size), Preview_ALPHA_ref.resize(size), 
	Preview_VELOCITY_ref.resize(size), Preview_ALPHADOT_ref.resize(size), Preview_BETA_ref.resize(size), Preview_BETADOT_ref.resize(size);
}

//전체 traj에서 horizon 만큼 잘라서 보낸다.
void MPC_TWIP::solveRef()
{
	
	Get_Reference(global_indx, Traj_MPC::traj.phi, Preview_PHI_ref);
	Get_Reference(global_indx, Traj_MPC::traj.phidot, Preview_PHIDOT_ref);

	Get_Reference(global_indx, Traj_MPC::traj.alpha, Preview_ALPHA_ref);
	Get_Reference(global_indx, Traj_MPC::traj.velocity, Preview_VELOCITY_ref);
	Get_Reference(global_indx, Traj_MPC::traj.alphadot, Preview_ALPHADOT_ref);

	Get_Reference(global_indx, Traj_MPC::traj.beta, Preview_BETA_ref);
	Get_Reference(global_indx, Traj_MPC::traj.betadot, Preview_BETADOT_ref);
}


void MPC_TWIP::solveOptimization(qpOASES::SQProblem& example,
	qpOASES::real_t* H,
	qpOASES::real_t* g,
	qpOASES::real_t* cA,
	qpOASES::real_t* lbA,
	qpOASES::real_t* ubA,
	qpOASES::real_t* lb,
	qpOASES::real_t* ub,
	int nWSR,
	qpOASES::real_t* Opt) {
	// Initialize the problem
	example.init(H, g, cA, lb, ub, lbA, ubA, nWSR, 0);

	// Solve the problem and store the solution in Opt
	example.getPrimalSolution(Opt);
}

//루프돌기전 세팅 마무리?

void MPC_TWIP::OutLoopSetYaw()
{
	//discretize
	auto [Ad_yaw, Bd_yaw] = discretizeState(DOF_Y, Ac_yaw, Bc_yaw);
	//stacking prediction matrix 
	auto [Pss_yaw, Pus_yaw] = stackingMatrix(DOF_Y, Ad_yaw, Bd_yaw);
	
	//get H matirx for QP
	auto H_yaw = assemble_HMatrix(DOF_Y, Pus_yaw, W_YAW);
	//constraint A matrix for QP
	auto cA_yaw = constraint_AMatrix();
	//constraint A bound for QP
	auto [Aub_yaw, Alb_yaw] = Bound_AMatrix(upperbound, lowerbound);

	qpOASES::real_t* H_yaw_qp   = Convert2RealT(H_yaw);
	qpOASES::real_t* cA_yaw_qp  = Convert2RealT(cA_yaw);
	qpOASES::real_t* Aub_yaw_qp = Convert2RealT(Aub_yaw);
	qpOASES::real_t* Alb_yaw_qp = Convert2RealT(Alb_yaw);
}


void MPC_TWIP::OutLoopSetPitch()
{
	auto [Ad_pitch, Bd_pitch] = discretizeState(DOF_P, Ac_pitch, Bc_pitch);

	auto [Pss_pitch, Pus_pitch] = stackingMatrix(DOF_P, Ad_pitch, Bd_pitch);

	auto H_pitch = assemble_HMatrix(DOF_P, Pus_pitch, W_PITCH);

	auto cA_pitch = constraint_AMatrix();

	auto [Aub_pitch, Alb_pitch] = Bound_AMatrix(upperbound, lowerbound);

	qpOASES::real_t* H_pitch_qp = Convert2RealT(H_pitch);
	qpOASES::real_t* cA_pitch_qp = Convert2RealT(cA_pitch);
	qpOASES::real_t* Aub_pitch_qp = Convert2RealT(Aub_pitch);
	qpOASES::real_t* Alb_pitch_qp = Convert2RealT(Alb_pitch);
}

void MPC_TWIP::OutLoopSetRoll()
{
}
