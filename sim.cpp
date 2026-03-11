#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>

float M = 1.00784*1.66*pow(10,-27); //[kg] //using atomic mass of H

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
    return 1.380649*pow(10, -23)*K;
}
float E_Kel_Conv(float E)
{
    return E/(1.380649*pow(10, -23));
}

// well definitions (rough area the particle is confined to)
float wellLength = -0.25; //[m]
float wellDepth = 0.5; // [kelvin]

// Particle definitions and initials
float K_i = 0.000; // kelvin
float v_0 = Vel_Conv(Kel_E_Conv(K_i)); //[m/s] // 4ms corresponds to about 1 milikelvin
float x_0 = wellLength; //[m] //particle starts at far left side 

// bump consts
float bumpAmp = 0.0005; // [kelvin]
float bumpLength = 0.05; //[m]
float bumpPos = (wellLength/2) - (bumpLength/2); //[m] //places bump in middle
float bumpOmega = 0.5*(M_PI*v_0/wellLength); // [rad/s] //driving frequency of potential 'bump'
float bumpPhase = M_PI; //[rads]
// Spring Potential Consts
float springPEMax = 0.0005; // Kelvin, this is the max value of the potential in the well
float springK = 2*Kel_E_Conv(springPEMax)/pow(wellLength/2,2); // [N/m]  
// time step and sim length
float timeStep = 0.0001; //[s]
float stopTime = 1; //[s] //Overall Length of sim (will be rounded if not divisible)
int N_t = std::round(stopTime/timeStep); // number of time steps
// vectors for storing positions and velocity
/* storing in the form (a,v,x):
calc each step in order v -> x -> a
take 1/2 of each velocity index to get acc value
*/
std::vector<float> X;
std::vector<float> V;
std::vector<float> A;

// Potential/force functions
float Bump_Potential(float x, float time)
{
    /* Potential function for a confined bump over place in the middle of a well
    modelled as a abs(cos(X) - 1) function for a smooth transition to the flat well
    then off set by a cos(omega t) to oscillate the bump up and down
    PE = A/2(cos(pi/l * x)-1)*cos(omega*t) */
    // THIS POTENTIAL IS SYMMETRIC, CAN USE PERIODIC BOUNDARY CONDS
    // check if over bump
    if (x > bumpPos && x < bumpPos + bumpLength)
    {
        std::cout<<"Over bump ";
        return Kel_E_Conv(bumpAmp)*0.5*(std::cos((M_PI/bumpLength)*(x-bumpPos))-1)*std::cos(bumpOmega*time + bumpPhase);
    }
    else
    {
        return 0;
    }
}
float Bump_Force(float x, float time)
{
    /*
    Force version of the bump potential, just the derivative F = -dU/dx
    */
    if (x > bumpPos && x < bumpPos + bumpLength)
    {
        std::cout<<"Over bump ";
        return Kel_E_Conv(bumpAmp)*(M_PI/bumpLength)*0.5*std::sin((M_PI/bumpLength)*(x-bumpPos))*std::cos(bumpOmega*time + bumpPhase);
    }
    else
    {
        return 0;
    }
}
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

// MAIN
int main()
{
    // open results file
    std::ofstream result;
    result.open("results.dat");

    // define i=0 positions from initials
    X.push_back(x_0);
    V.push_back(v_0);
    A.push_back(Pendulum_Force(X[0])/M);

    // push initial time and KE
    result<<0<<','<<E_Kel_Conv(KE_Conv(V[0]))<<'\n';

    for (int i=1; i<N_t; i++)
    {
    // time = i * timestep
    float v_i = V[i-1] + A[i-1]*timeStep;
    V.push_back(v_i);
    float x_i = X[i-1] + V[i]*timeStep;
    X.push_back(x_i);
    float a_i = Pendulum_Force(X[i])/M;
    A.push_back(a_i);

    // push back resulting ke and time from this step
    result<<i*timeStep<<','<<E_Kel_Conv(KE_Conv(V[i]))<<'\n';
    
    
    }

    result.close();
    std::cout<<"EXECUTE SUCESSFUL";
    return 0;
}