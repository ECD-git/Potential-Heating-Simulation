#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>

// Real world consts
float M = 1.00784*1.66*pow(10,-27); //[kg] //using atomic mass of H
float K_B = 1.380649*pow(10, -23); //boltzmann const
// USEFUL FUNCTIONS
float KE_Conv(float v)
{
    // converts a velocity in [m/s] to a kinetic energy [J]
    return 0.5*M*pow(v,2);
};
float Vel_Conv(float KE)
{
    // converts a KE [J] to a velocity [m/s]
    return pow((2*KE)/M, 0.5);
};
float eV_Conv(float KE)
{
    // Converts an energy [J] to an energy [eV]
    return KE/(1.6*pow(10,-19));
}
float Kel_E_Conv(float K)
{
    // converts a temperature [K] to an energy [J]
    return K_B*K;
}
float E_Kel_Conv(float E)
{
    return E/(K_B);
}

std::vector<float> LogSpace(int start, int stop, int step)
{
    // generates a log space of values between index start and stop 
    // and seperated by intex 10^step, in integer steps
    std::vector<float> result;
    for(int i = 0; start + i <= stop; i+=step)
    {
        result.push_back(pow(10, start + i));
    }
    return result;
}

// Particle definitions and initials
float K_i = 0; // kelvin
float v_0 = Vel_Conv(Kel_E_Conv(K_i)); //[m/s] // 4ms corresponds to about 1 milikelvin
float x_0 = -0.25; //[m] //particle starts at far left side 

// Spring Potential Consts
float springPEMax = 0.0005; // Kelvin, this is the max value of the potential in the well
float springK = 2*Kel_E_Conv(springPEMax)/pow(x_0/2,2); // [N/m]  

// vectors for storing positions and velocity
/* storing in the form (a,v,x):
calc each step in order v -> x -> a
take 1/2 of each velocity index to get acc value
*/
std::vector<float> X;
std::vector<float> V;
std::vector<float> A;

// vector of total energies against time step
std::vector<float> E;

float Pendulum_Force(float x)
{
    /*
    Models the force experienced in a PE = 1/2kx^2 potential
    F = -kx
    Offset by half the length of the well so its centered
    this potential is also symmetric -> periodic conditions
    */
   return -1*springK*(x);
}
float Pendulum_Potential(float x)
{
    /*
    Potential ver of Pendulum_Force
    In case I need it, will be removed if obselete in final ver
    */
   return 0.5*springK*pow(x,2);
}

// time step and sim length
std::vector<float> timeSteps = LogSpace(-3, 3, 1); //[s]
float N = 1000; // number of oscillations. (roughly, time is rounded to nearest num of steps)

float stopTime = N * pow(M/springK, 0.5) * 2 * M_PI; //[s] 


// MAIN
int main()
{
    // open results file
    std::ofstream result;
    result.open("results.dat");

    float E_i;

    for(const float& step : timeSteps)
    {
        int N_t = std::round(stopTime/step); // number of time steps
        // define i=0 positions from initials
        //X.push_back(x_0);
        //V.push_back(v_0);
        //A.push_back(Pendulum_Force(X[0])/M);

        float x_i = x_0;
        float v_i = v_0;
        float a_i = Pendulum_Force(x_i)/M;

        for (int i=1; i<N_t; i++)
        {
            // time = i * timestep
            float v_i = v_i + a_i*step;
            float x_i = x_i + v_i*step;
            float a_i = Pendulum_Force(x_i)/M;
        }
        // get total energy at the last step
        E_i = KE_Conv(v_i) + Pendulum_Potential(x_i); // [j]
        // push to results with the step
        result<<step<<','<<E_Kel_Conv(E_i)<<'\n';
        X.clear(); V.clear(); A.clear(); //clear all 3 vectors for reuse
    }

    result.close();
    std::cout<<"EXECUTE SUCESSFUL";
    return 0;
}