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

// time step and sim length
float timeStep = 0.0001; //[s]
float stopTime = 0.1; //[s] //Overall Length of sim (will be rounded if not divisible)
int N_t = std::round(stopTime/timeStep); // number of time steps
// Particle definitions and initials
float K_i = 0.001; // kelvin
float v_0 = Vel_Conv(Kel_E_Conv(K_i)); //[m/s] // 4ms corresponds to about 1 milikelvin
float x_0 = 0; //[m] //particle starts at far left side 

// Spring Potential Consts
float springPEMax = 0.0005; // Kelvin, this is the max value of the potential in the well
float springK = 2*Kel_E_Conv(springPEMax)/pow(x_0/2,2); // [N/m]  

// Standing wave values
float standingPEMax = 0.0001; // Kelvin
float lambda = 500 * pow(10,-9); // [m] typical light wavelength
float tau_On = 0.001; // [s] tuned so we get about 10 time steps with the force on
float k = M_PI/(tau_On*v_0);
// make tau_Off some random variable in the scale tau_On;

// vectors for storing positions and velocity
/* storing in the form (a,v,x):
calc each step in order v -> x -> a
take 1/2 of each velocity index to get acc value
*/
std::vector<float> X;
std::vector<float> V;
std::vector<float> A;

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
float Standing_Force(float x, bool on)
{
    // force of standing wave implementing force on particle, tuned to A_0
    if (on == 1)
    {
        return Kel_E_Conv(standingPEMax) *k* std::sin(2*k*x);
    }
    else
    {
        return 0;
    }
}
float Standing_Potential(float x, bool on)
{
    // potential of a standing wave implementing a force on our particle
    if (on == 1)
    {
        return Kel_E_Conv(standingPEMax) * std::cos(k*x);
    }
    else
    {
        return 0;
    }
}

// function to grab the potential at any time t or position x
// TODO design this for the turning off / on of a standing wave

// MAIN
int main()
{
    // open results file
    std::ofstream result;
    result.open("results.dat");
    std::ofstream walk;
    walk.open("walk.dat");

    // define i=0 positions from initials
    bool PE_on = true; // potential is initially on
    float cooldown; // a cooldown timer for switching the potential on/off
    if (PE_on==1)
    {
        cooldown = tau_On;
    }
    else
    {
        cooldown = tau_On; // TAU OFF
    }
    X.push_back(x_0);
    V.push_back(v_0);
    A.push_back(Standing_Force(X[0], PE_on)/M);
    cooldown -= timeStep;

    // push initial time and KE
    result<<0<<','<<E_Kel_Conv(KE_Conv(V[0]))<<'\n';
    walk<<0<<','<<V[0]*M<<'\n';

    for (int i=1; i<N_t; i++)
    {
    // time = i * timestep
    float v_i = V[i-1] + A[i-1]*timeStep;
    V.push_back(v_i);
    float x_i = X[i-1] + V[i]*timeStep;
    X.push_back(x_i);
    if(cooldown <= 0)
    {
        if(PE_on==1)
        {
            PE_on = false;
            cooldown = tau_On/2; // GET RANDOM TAU OFF HERE
            walk<<i*timeStep<<','<<V[i]*M<<'\n';
        }
        else 
        {
            PE_on = true;
            cooldown = tau_On/2; // TAU OFF
        }
    }
    float a_i = Standing_Force(X[i], PE_on)/M;
    std::cout<<a_i<<','<<v_i<<','<<PE_on<<std::endl; // Debug purposes
    A.push_back(a_i);
    cooldown -= timeStep;

    // push back resulting ke and time from this step
    result<<i*timeStep<<','<<E_Kel_Conv(KE_Conv(V[i]))<<'\n';
    }

    result.close();
    walk.close();
    std::cout<<"EXECUTE SUCESSFUL";
    return 0;
}