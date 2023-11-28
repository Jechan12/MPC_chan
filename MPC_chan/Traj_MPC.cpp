#include "Traj_MPC.h"

void Traj_MPC::setsize_Traj(int total_step_size)
{
    traj.phi.resize(total_step_size, 0.0);
    traj.phidot.resize(total_step_size, 0.0);

    traj.alpha.resize(total_step_size , 0.0);
    traj.velocity.resize(total_step_size, 0.0);
    traj.alphadot.resize(total_step_size, 0.0);

    traj.beta.resize(total_step_size, 0.0);
    traj.betadot.resize(total_step_size, 0.0);
}

//traj 유형 , total stepsize , step별 크기 , ...
void Traj_MPC::initializeTrj(const int& trajtype , int total_step_size , int steptime)
{
    setsize_Traj(total_step_size);

    TrjType = trajtype;
    double local_t = 0.0;
    int local_index = 0;

    switch (TrjType) {
    case STAY:
        while (local_index < total_step_size)
        {







            local_index++;
        }

    case DIRECT:
        while (local_index < total_step_size)
        {







            local_index++;
        }

    case CIRCLE:
        while (local_index < total_step_size)
        {







            local_index++;
        }
    }
}
