#include "Traj_MPC.h"

void Traj_MPC::setsize_Traj(const int total_step_size)
{
    traj.phi.resize(total_step_size, 0.0);
    traj.phidot.resize(total_step_size, 0.0);

    traj.alpha.resize(total_step_size , 0.0);
    traj.velocity.resize(total_step_size, 0.0);
    traj.alphadot.resize(total_step_size, 0.0);

    traj.beta.resize(total_step_size, 0.0);
    traj.betadot.resize(total_step_size, 0.0);

    traj.x_d.resize(total_step_size, 0.0);
    traj.y_d.resize(total_step_size, 0.0);
}

//traj 유형 , total stepsize , step별 크기 , 총 작동 시간 ...
void Traj_MPC::TrajectoryPlanning(const int& trajtype , const int& total_step_size , const double& steptime , const double& execution_time)
{
    TIME_TASK_EXECUTION = execution_time;
    setsize_Traj(total_step_size);
    
    TrjType = trajtype;
    double local_t = 0.0;
    int local_index = 0;

    if (TrjType == STAY)
    {
        while (local_index < total_step_size)
        {
            // all desired trajctories are already ZERO
            break;
        }
    }
    else if (TrjType == DIRECT) //Given v_r & phidot_r
    {
        
        
        double phidot_r = phidot_r_Direct;
        double phi_r = phi_r_Direct;
        double v_r = v_r_Direct;

        while (local_index < total_step_size)
        {
            local_t = local_index * steptime;
            
            if (local_t < rampUPDuration)
            { //Ramp-up phase
                v_r = (local_t / rampUPDuration) * v_r_Direct;
            }
            else if (local_t > (TIME_TASK_EXECUTION - rampDOWNDuration ))
            { //Ramp-down phase
                v_r = ((TIME_TASK_EXECUTION - local_t) / rampDOWNDuration) * v_r_Direct;
                phidot_r = ((TIME_TASK_EXECUTION - local_t) / rampDOWNDuration) * phidot_r_Direct;
            }
            else
            { //constant velocity
                v_r = v_r_Direct;
            }

            xdot_r = v_r * cos(phi_r);
            ydot_r = v_r * sin(phi_r);

            //kinematic controller에서 활용?
            x_r += xdot_r * steptime;
            y_r += ydot_r * steptime;
            phi_r += phidot_r * steptime;

            //position (수정예정)
            traj.x_d[local_index] = traj.x_d[local_index] + x_r;
            traj.y_d[local_index] = traj.y_d[local_index] + y_r;

            //Yaw
            traj.phi[local_index] = phi_r;
            traj.phidot[local_index] = phidot_r;

            //Pitch
            traj.alpha[local_index] = 0.0;
            traj.velocity[local_index] = v_r;
            traj.alphadot[local_index] = 0.0;

            //Roll
            traj.beta[local_index] = 0.0;
            traj.betadot[local_index] = 0.0;

            local_index++;
        }
    }
    else if (TrjType == CIRCLE) //	Given linear velocity ref. v_r ( Thus phidot_r = v_r / radius )
    {
        //유지보수를 위해 헤더로 보내는게 좋을 듯

        double v_r = 2 * M_PI * NoCycle * magR / TIME_TASK_EXECUTION;
        double phidot_r = v_r / magR;

        while (local_index < total_step_size)
        {
            local_t = local_index * steptime;






            local_index++;
        }
    }
    else if (TrjType == TEST)
    {
        //명준씨 ver
        while (local_index < total_step_size)
        {
            local_t = local_index * steptime;

            traj.phi[local_index] = 1.0 * tanh(local_t);
            traj.phidot[local_index] = 0.0;

            traj.alpha[local_index] = 0.0;
            traj.velocity[local_index] = 18.0 * tanh(local_t);
            traj.alphadot[local_index] = 0.0;

            traj.beta[local_index] = sin(local_t);
            traj.betadot[local_index] = cos(local_t);

            local_index++;
        }
    }
}


