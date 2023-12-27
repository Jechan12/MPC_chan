#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Dense>

struct Traj
{
	std::vector<double> phi, phidot;
	std::vector<double> alpha, alphadot, velocity;
	std::vector<double> beta, betadot;
	std::vector<double> x_d, y_d;
};

class Traj_MPC
{
	////////////// Desired velocities of each trajectories///////////
	// DIRECT
	double v_r_Direct = 6.0;    double phidot_r_Direct = 0.1;    double phi_r_Direct = 0.0;
	// CIRCLE
	double magR = 5.0;  int NoCycle = 1;
	// SLARUM

	// RAMP_velovity
	double rampUPDuration = 2.0;		double rampDOWNDuration = v_r_Direct;  // [s]

public:
	Traj traj;

	double TIME_TASK_EXECUTION{0.0};

	int TrjType;
	int	KinCtrlType;
	

	double	   x_r{ 0.0 }, y_r{ 0.0 }, alpha_r{ 0.0 };
	double	xdot_r{ 0.0 },  ydot_r{ 0.0 },  alphadot_r{ 0.0 };
	double xddot_r{ 0.0 }, yddot_r{ 0.0 };
	//phi_r{ 0.0 },phidot_r{ 0.0 },
	


	//Total �۵� �ð���ŭ�� step �����Ҵ�
	void setsize_Traj(const int size);

	// traj �� �⺻ �ʱ�ȭ + planning
	void TrajectoryPlanning(const int& trajtype , const int& total_step_size, const double& steptime, const double& execution_time);

	//Kinematic Controller

	
	//traj�ڸ��� (horizon ������ŭ �ڸ���) -> LMPC�� ����?
	
	enum Trajectory_Types
	{
		STAY = 0,
		DIRECT,
		CIRCLE,
		TEST
	};


	enum KinematicCtrl_Types
	{
		NO_CTRL = 0,
		LQR_CTRL,
		LINEAR_CTRL,
		NONLINEAR_CTRL
	};


};

