#pragma once
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

struct Traj
{
	std::vector<double> phi, phidot;
	std::vector<double> alpha, alphadot, velocity;
	std::vector<double> beta, betadot;
};

enum Trajectory_Types
{
	STAY = 0,
	DIRECT,
	CIRCLE,
};

enum KinematicCtrl_Types
{
	NO_CTRL = 0,
	LQR_CTRL,
	LINEAR_CTRL,
	NONLINEAR_CTRL
};

class Traj_MPC
{
private:
	Traj traj;

public:
	int TrjType;
	int	KinCtrlType;


	//Total 작동 시간만큼의 step 동적할당
	void setsize_Traj(int size);

	// traj 모양 별 기본 초기화
	void initializeTrj(const int& trajtype, int total_step_size, int steptime);

	//Kinematic Controller

	
	//traj자르기 (horizon 갯수만큼 자르기)
	//void getTraj();






};

