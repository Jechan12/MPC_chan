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
public:
	Traj traj;


	double TIME_TASK_EXECUTION = 0.0;

	int TrjType;
	int	KinCtrlType;
	
	double	   x_r{ 0.0 },	   y_r{ 0.0 },    phi_r{ 0.0 },    alpha_r{ 0.0 },  v_r{ 0.0 };
	double	xdot_r{ 0.0 },  ydot_r{ 0.0 }, phidot_r{ 0.0 }, alphadot_r{ 0.0 };
	double xddot_r{ 0.0 }, yddot_r{ 0.0 };

	//Total 작동 시간만큼의 step 동적할당
	void setsize_Traj(const int size);

	// traj 별 기본 초기화 + planning
	void TrajectoryPlanning(const int& trajtype , const int& total_step_size, const double& steptime, const double& execution_time);

	//Kinematic Controller

	
	//traj자르기 (horizon 갯수만큼 자르기) -> LMPC에 있음?
	
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

