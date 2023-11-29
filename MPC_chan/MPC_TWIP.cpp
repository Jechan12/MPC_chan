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
	std::tie(Pss_yaw, Pus_yaw) = stackingMatrix(DOF_Y, Ad_yaw, Bd_yaw);
	//cout << Pss_yaw[0] << endl;
	//cout << Pss_yaw[1] << endl;
	//get H matirx for QP
	H_yaw = assemble_HMatrix(DOF_Y, Pus_yaw, W_YAW);
	//constraint A matrix for QP
	cA_yaw = constraint_AMatrix();
	//constraint A bound for QP
	std::tie(Aub_yaw, Alb_yaw) = Bound_AMatrix(upperbound, lowerbound);

	//cout << "Ad_yaw.size() : " << Ad_yaw.size() << endl;
	//cout << "H_yaw.size() : " << H_yaw.size() << endl;

	H_yaw_qp   = Convert2RealT(H_yaw);
	cA_yaw_qp  = Convert2RealT(cA_yaw);
	Aub_yaw_qp = Convert2RealT(Aub_yaw);
	Alb_yaw_qp = Convert2RealT(Alb_yaw);
}


void MPC_TWIP::OutLoopSetPitch()
{
	auto [Ad_pitch, Bd_pitch] = discretizeState(DOF_P, Ac_pitch, Bc_pitch);

	std::tie(Pss_pitch, Pus_pitch) = stackingMatrix(DOF_P, Ad_pitch, Bd_pitch);

	H_pitch = assemble_HMatrix(DOF_P, Pus_pitch, W_PITCH);
	//cout << "H_pitch" << endl;
	//cout << H_pitch<< endl;
	cA_pitch = constraint_AMatrix();

	std::tie(Aub_pitch, Alb_pitch) = Bound_AMatrix(upperbound, lowerbound);

	//cout << "Pss_pitch.size() : " << Pss_pitch.size() << endl;

	H_pitch_qp   = Convert2RealT(H_pitch);
	cA_pitch_qp  = Convert2RealT(cA_pitch);
	Aub_pitch_qp = Convert2RealT(Aub_pitch);
	Alb_pitch_qp = Convert2RealT(Alb_pitch);
}

void MPC_TWIP::OutLoopSetRoll()
{
	auto [Ad_roll, Bd_roll] = discretizeState(DOF_R, Ac_roll,Bc_roll);

	std::tie(Pss_roll, Pus_roll) = stackingMatrix(DOF_R, Ad_roll, Bd_roll);
	//cout << "ROLL" << endl;
	//cout << Pss_roll[0] << endl;
	H_roll = assemble_HMatrix(DOF_R, Pus_roll, W_ROLL);
	//cout << "H_roll" << endl;
	//cout << H_roll << endl;
	cA_roll = constraint_AMatrix();

	std::tie(Aub_roll, Alb_roll) = Bound_AMatrix(upperbound, lowerbound);

	H_roll_qp   = Convert2RealT(H_roll);
	cA_roll_qp  = Convert2RealT(cA_roll);
	Aub_roll_qp = Convert2RealT(Aub_roll);
	Alb_roll_qp = Convert2RealT(Alb_roll);
}

void MPC_TWIP::InLoopSolvYaw()
{
	qpOASES::int_t nWSR = 100;
	qpOASES::SQProblem yawExample(MPC_TWIP::Dim, 1);
	Eigen::VectorXd g_yaw = assemble_gMatrix(MPC_TWIP::DOF_Y, Pss_yaw, Pus_yaw, MPC_TWIP::x_Yaw, W_YAW, { &(Preview_PHI_ref),  &Preview_PHIDOT_ref });
	qpOASES::real_t* g_yaw_qp = Convert2RealT(g_yaw);
	solveOptimization(yawExample, H_yaw_qp, g_yaw_qp, cA_yaw_qp, Alb_yaw_qp, Aub_yaw_qp, nullptr, nullptr, nWSR, Yaw_Opt);
	Trgt_a_yaw = Yaw_Opt[0];
}

void MPC_TWIP::InLoopSolvPitch()
{
	qpOASES::int_t nWSR = 100;
	qpOASES::SQProblem pitchExample(MPC_TWIP::Dim, 1);
	Eigen::VectorXd g_pitch = assemble_gMatrix(MPC_TWIP::DOF_P, Pss_pitch, Pus_pitch, MPC_TWIP::x_Pitch, W_PITCH, { &(Preview_ALPHA_ref),&(Preview_VELOCITY_ref),&(Preview_ALPHADOT_ref) });
	qpOASES::real_t* g_pitch_qp = Convert2RealT(g_pitch);

	//cout << "sizeof(H_pitch_qp) : " << sizeof(&H_pitch_qp) << endl;

	//cout << "sizeof(pitch_g_qp) : " << sizeof(&pitch_g_qp) << endl;

	//cout << "sizeof(cA_pitch_qp) : " << sizeof(&cA_pitch_qp) << endl;

	//cout << "sizeof(Alb_pitch_qp) : " << sizeof(&Alb_pitch_qp) << endl;

	//cout << "sizeof(Aub_pitch_qp) : " << sizeof(&Aub_pitch_qp) << endl;

	solveOptimization(pitchExample, H_pitch_qp, g_pitch_qp, cA_pitch_qp, Alb_pitch_qp, Aub_pitch_qp, nullptr, nullptr, nWSR, Pitch_Opt);
	Trgt_a_pitch = Pitch_Opt[0]; //첫번째 가속도
}

void MPC_TWIP::InLoopSolvRoll()
{
	qpOASES::int_t nWSR = 100;
	
	qpOASES::SQProblem rollExample(MPC_TWIP::Dim, 1);

	Eigen::VectorXd g_roll = assemble_gMatrix(MPC_TWIP::DOF_R, Pss_roll, Pus_roll, x_Roll, W_ROLL, { &(Preview_BETA_ref), &(Preview_BETADOT_ref) });
	
	qpOASES::real_t* g_roll_qp = Convert2RealT(g_roll);

	solveOptimization(rollExample, H_roll_qp, g_roll_qp, cA_roll_qp, Alb_roll_qp, Aub_roll_qp, nullptr, nullptr, nWSR, Roll_Opt);
	
	Trgt_a_roll = Roll_Opt[0];
}
