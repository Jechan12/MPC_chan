#include "MPC_TWIP.h"

int main()
{

	std::cout << "main in MPC_TWIP.cpp" << endl;

	MPC_TWIP twip;
	int trajtype = twip.DIRECT;
	
	//FILE* fout;
	std::chrono::system_clock::time_point start_time3 = std::chrono::system_clock::now();
	std::cout << "Planning START" << endl;
	twip.TrajectoryPlanning(trajtype , twip.DimTotal ,twip.dT, twip.OperationTime);
	std::cout << "Planning END" << endl;

	//errno_t err = fopen_s(&fout, "MPC.txt", "wt");

	
	twip.setStateMatrix();
	twip.Weight2vec();

	auto [Ad_yaw, Bd_yaw] = twip.discretizeState(twip.DOF_Y, twip.Ac_yaw, twip.Bc_yaw);
	auto [Ad_pitch, Bd_pitch] = twip.discretizeState(twip.DOF_P, twip.Ac_pitch, twip.Bc_pitch);
	auto [Ad_roll, Bd_roll] = twip.discretizeState(twip.DOF_R, twip.Ac_roll, twip.Bc_roll);

	//stacking prediction matrix
	auto [Pss_yaw, Pus_yaw] = twip.stackingMatrix(twip.DOF_Y, Ad_yaw, Bd_yaw);
	auto [Pss_pitch, Pus_pitch] = twip.stackingMatrix(twip.DOF_P, Ad_pitch, Bd_pitch);
	auto [Pss_roll, Pus_roll] = twip.stackingMatrix(twip.DOF_R, Ad_roll, Bd_roll);

	//Pus 바로 H matirx로 대입
	auto H_yaw = twip.assemble_HMatrix(twip.DOF_Y, Pus_yaw, twip.W_YAW);
	auto H_pitch = twip.assemble_HMatrix(twip.DOF_P, Pus_pitch, twip.W_PITCH);
	auto H_roll = twip.assemble_HMatrix(twip.DOF_R, Pus_roll, twip.W_ROLL);

	//constraint A Matrix
	auto cA_yaw = twip.constraint_AMatrix();
	auto cA_pitch = twip.constraint_AMatrix();
	auto cA_roll = twip.constraint_AMatrix();

	//constraint A bound
	auto [Aub_yaw, Alb_yaw] = twip.Bound_AMatrix(twip.upperbound, twip.lowerbound);
	auto [Aub_pitch, Alb_pitch] = twip.Bound_AMatrix(twip.upperbound, twip.lowerbound);
	auto [Aub_roll, Alb_roll] = twip.Bound_AMatrix(twip.upperbound, twip.lowerbound);
	cout << "Ad_yaw.size() : " << Ad_yaw.size() << endl;
	cout << "Pss_pitch.size() : " << Pss_pitch.size() << endl;
	cout << "H_yaw.size() : " << H_yaw.size() << endl;

	////real_t로 바꿀 수 있는 건 바꿔보자
	qpOASES::real_t* H_yaw_qp = twip.Convert2RealT(H_yaw);
	qpOASES::real_t* H_pitch_qp = twip.Convert2RealT(H_pitch);
	qpOASES::real_t* H_roll_qp = twip.Convert2RealT(H_roll);

	qpOASES::real_t* cA_yaw_qp = twip.Convert2RealT(cA_yaw);
	qpOASES::real_t* cA_pitch_qp = twip.Convert2RealT(cA_pitch);
	qpOASES::real_t* cA_roll_qp = twip.Convert2RealT(cA_roll);

	qpOASES::real_t* Aub_yaw_qp = twip.Convert2RealT(Aub_yaw);
	qpOASES::real_t* Alb_yaw_qp = twip.Convert2RealT(Alb_yaw);
	qpOASES::real_t* Aub_pitch_qp = twip.Convert2RealT(Aub_pitch);
	qpOASES::real_t* Alb_pitch_qp = twip.Convert2RealT(Alb_pitch);
	qpOASES::real_t* Aub_roll_qp = twip.Convert2RealT(Aub_roll);
	qpOASES::real_t* Alb_roll_qp = twip.Convert2RealT(Alb_roll);

	//ref initialize.
	twip.initializeRefVec(twip.Dim);

	std::cout << "loopstart" << std::endl;
	std::chrono::system_clock::time_point end_time3 = std::chrono::system_clock::now();

	while (twip.global_indx < twip.Dim_TIME)
	{
		USING_NAMESPACE_QPOASES

			twip.main_time = double(twip.global_indx) * twip.dT;

		std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
		
		twip.solveRef();

		//std::chrono::system_clock::time_point start_time2 = std::chrono::system_clock::now();
		//Eigen::VectorXd yaw_g = twip.assemble_gMatrix(twip.DOF_Y, Pss_yaw, Pus_yaw, twip.x_Yaw, twip.W_YAW, { &(twip.Preview_PHI_ref),  &twip.Preview_PHIDOT_ref });
		//Eigen::VectorXd roll_g = twip.assemble_gMatrix(twip.DOF_R, Pss_roll, Pus_roll, twip.x_Roll, twip.W_ROLL, { &(twip.Preview_BETA_ref), &(twip.Preview_BETADOT_ref) });
		//Eigen::VectorXd pitch_g = twip.assemble_gMatrix(twip.DOF_P, Pss_pitch, Pus_pitch, twip.x_Pitch, twip.W_PITCH, { &(twip.Preview_ALPHA_ref),&(twip.Preview_VELOCITY_ref),&(twip.Preview_ALPHADOT_ref) });

		//qpOASES::real_t* yaw_g_qp   = twip.Convert2RealT(yaw_g);
		//qpOASES::real_t* pitch_g_qp = twip.Convert2RealT(pitch_g);
		//qpOASES::real_t* roll_g_qp  = twip.Convert2RealT(roll_g);
		//std::chrono::system_clock::time_point end_time2 = std::chrono::system_clock::now();

		int_t nWSR = 100;
		//..............................< ORIGINAL >.................................//
		//////////////////example.init(H, g, A, lb, ub, lbA, ubA, nWSR);///////////////////////
		//qpOASES::SQProblem example(twip.Dim , 1);
		//
		//example.init(H_yaw_qp, yaw_g_qp, cA_yaw_qp, 0, 0, Alb_yaw_qp, Aub_yaw_qp, nWSR, 0);
		//example.getPrimalSolution(twip.Yaw_Opt);
		//
		//example.init(H_pitch_qp, pitch_g_qp, cA_pitch_qp, 0, 0, Alb_pitch_qp, Aub_pitch_qp, nWSR, 0);
		//example.getPrimalSolution(twip.Pitch_Opt);
		//
		//example.init(H_roll_qp, roll_g_qp, cA_roll_qp, 0, 0, Alb_roll_qp, Aub_roll_qp, nWSR, 0);
		//example.getPrimalSolution(twip.Roll_Opt);
		//
		//twip.Trgt_a_yaw   = twip.Yaw_Opt[0];
		//twip.Trgt_a_pitch = twip.Pitch_Opt[0];
		//twip.Trgt_a_roll  = twip.Roll_Opt[0];
		//................................................................................//

		//...........................<Task Parallelism>.....................................//
		qpOASES::SQProblem yawExample(twip.Dim, 1);
		qpOASES::SQProblem pitchExample(twip.Dim, 1);
		qpOASES::SQProblem rollExample(twip.Dim, 1);

		auto futurePitch = std::async(std::launch::async, [&]() {
			Eigen::VectorXd pitch_g = twip.assemble_gMatrix(twip.DOF_P, Pss_pitch, Pus_pitch, twip.x_Pitch, twip.W_PITCH, { &(twip.Preview_ALPHA_ref),&(twip.Preview_VELOCITY_ref),&(twip.Preview_ALPHADOT_ref) });
			qpOASES::real_t* pitch_g_qp = twip.Convert2RealT(pitch_g);
			//qpOASES::real_t* pitch_g_qp = twip.Convert2RealT2(pitch_g);
			twip.solveOptimization(pitchExample, H_pitch_qp, pitch_g_qp, cA_pitch_qp, Alb_pitch_qp, Aub_pitch_qp, nullptr, nullptr, nWSR, twip.Pitch_Opt);
			twip.Trgt_a_pitch = twip.Pitch_Opt[0];
			//delete[] pitch_g_qp;
			//pitch_g_qp = nullptr;
			});


		auto futureYaw = std::async(std::launch::async, [&]() {
			Eigen::VectorXd yaw_g = twip.assemble_gMatrix(twip.DOF_Y, Pss_yaw, Pus_yaw, twip.x_Yaw, twip.W_YAW, { &(twip.Preview_PHI_ref),  &twip.Preview_PHIDOT_ref });
			qpOASES::real_t* yaw_g_qp = twip.Convert2RealT(yaw_g);
			//qpOASES::real_t* yaw_g_qp = twip.Convert2RealT2(yaw_g);
			twip.solveOptimization(yawExample, H_yaw_qp, yaw_g_qp, cA_yaw_qp, Alb_yaw_qp, Aub_yaw_qp, nullptr, nullptr, nWSR, twip.Yaw_Opt);
			twip.Trgt_a_yaw = twip.Yaw_Opt[0];
			//delete[] yaw_g_qp;
			//yaw_g_qp = nullptr;
			});


		auto futureRoll = std::async(std::launch::async, [&]() {
			Eigen::VectorXd roll_g = twip.assemble_gMatrix(twip.DOF_R, Pss_roll, Pus_roll, twip.x_Roll, twip.W_ROLL, { &(twip.Preview_BETA_ref), &(twip.Preview_BETADOT_ref) });
			qpOASES::real_t* roll_g_qp = twip.Convert2RealT(roll_g);
			//qpOASES::real_t* roll_g_qp = twip.Convert2RealT2(roll_g);
			twip.solveOptimization(rollExample, H_roll_qp, roll_g_qp, cA_roll_qp, Alb_roll_qp, Aub_roll_qp, nullptr, nullptr, nWSR, twip.Roll_Opt);
			twip.Trgt_a_roll = twip.Roll_Opt[0];
			//delete[] roll_g_qp;
			//roll_g_qp = nullptr;
			});

		futureYaw.get();
		futureRoll.get();
		futurePitch.get();
		//.........................................................................................//


		std::chrono::system_clock::time_point end_time = std::chrono::system_clock::now();
		std::chrono::microseconds micro = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
		cout << "micro : " << micro << endl;
		//std::chrono::microseconds micro2 = std::chrono::duration_cast<std::chrono::microseconds>(end_time2 - start_time2);
		//cout << "micro2 : " << micro2 << endl;
		twip.global_indx++;

	};

	std::chrono::microseconds micro3 = std::chrono::duration_cast<std::chrono::microseconds>(end_time3 - start_time3);
	cout << "micro3 : " << micro3 << endl;

	std::cout << "loop END " << endl;
	std::cout << "global_index : " << twip.global_indx << endl;

	return 0;
}