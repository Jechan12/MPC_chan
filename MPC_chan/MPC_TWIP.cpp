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


/////////////명준씨 ver////////////////
//void MPC_TWIP::TrajectoryPlanning(Traj &traj)
//{
//
//	//std::cout << "planning" << endl;
//	volatile double local_t = 0.0;
//	volatile int local_indx = 0;
//	//std::cout << "local_indx : " << local_indx << std::endl;
//
//	while (local_indx < DimTotal) //전체 스텝동안 가야할 trajectory
//	{
//
//		//std::cout << "local_indx : " << local_indx << std::endl;
//		local_t = local_indx * MPC_TWIP::dT;
//
//		traj.phi[local_indx] = 1.0 * tanh(local_t);
//		traj.phidot[local_indx] = 0.0;
//
//		traj.alpha[local_indx] = 0.0;
//		traj.velocity[local_indx] = 18.0 * tanh(local_t);
//		traj.alphadot[local_indx] = 0.0;
//
//		traj.beta[local_indx] = sin(local_t);
//		traj.betadot[local_indx] = cos(local_t);
//
//		local_indx++;
//	}
//	//cout << "traj loop out" << endl;
//
//}

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

