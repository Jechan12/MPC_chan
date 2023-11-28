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


	//Total �۵� �ð���ŭ�� step �����Ҵ�
	void setsize_Traj(int size);

	// traj ��� �� �⺻ �ʱ�ȭ
	void initializeTrj(const int& trajtype, int total_step_size, int steptime);

	//Kinematic Controller

	
	//traj�ڸ��� (horizon ������ŭ �ڸ���)
	//void getTraj();






};

