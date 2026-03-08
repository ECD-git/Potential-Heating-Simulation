#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>

// Particle definitions and initials
float v_0 = 0.1; //[m/s] 
float M = 1.00784*1.66*pow(10,-27); //[kg] //using atomic mass of H
float x_0 = 0; //[m] //particle starts at far left side 
// well definitions 
float wellLength = 0.25; //[m]
float wellDepth = 0.5; // [kelvin]
// bump constraints
float bumpAmp = 0.0000001; // [kelvin]
float bumpLength = 0.05; //[m]
float bumpPos = (wellLength/2) - (bumpLength/2); //[m] //places bump in middle
float bumpOmega = 0.5*(M_PI*v_0/wellLength); // [rad/s] //driving frequency of potential 'bump'
float bumpPhase = M_PI; //[rads]
// time step and sim length
float timeStep = 0.01; //[s]
float stopTime = 20; //[s] //Overall Length of sim

// USEFUL FUNCTIONS
float KEConv(float v)
{
    // converts a velocity in [m/s] to a kinetic energy [J]
    return 0.5*M*pow(v,2);
};
float velConv(float KE)
{
    // converts a KE [J] to a velocity [m/s]
    return pow((2*KE)/M, 0.5);
};
float eVConv(float KE)
{
    // Converts an energy [J] to an energy [eV]
    return KE/(1.6*pow(10,-19));
}
float KelEConv(float K)
{
    // converts a temperature [K] to an energy [J]
    return 1.380649*pow(10, -23)*K;
}
float BumpPotential(float x, float time)
{
    /* Potential function for a confined bump over place in the middle of a well
    modelled as a cos^2 - 1 function for a smooth transition to the flat well
    then off set by a cos(omega t) to oscillate the bump up and down
    PE = A(cos^2(pi/l * x)-1)*cos(omega*t) */

    // check if over bump
    if (x > bumpPos && x < bumpPos + bumpLength)
    {
        std::cout<<"Over bump ";
        return KelEConv(bumpAmp)*(pow(std::cos((M_PI/bumpLength)*(x-bumpPos)),2)-1)*std::cos(bumpOmega*time);
    }
    else
    {
        return 0;
    }
}

// MAIN
int main()
{
    float time = 0;
    float x = x_0;
    float v = v_0;
    int i = 1; // inital velocity is to the right
    float KE = KEConv(v_0); // init KE [j]
    float ET = KE + BumpPotential(x,time); //[J] Get E_total at x_0 and time=0

    //open file write our initials 
    std::ofstream result;
    result.open("results.dat");
    result<<time<<','<<KE<<'\n';

    while(time < stopTime)
    {   
        x += v*timeStep; // new x pos
        // check if out of bounds
        if (x < 0)
        {
            x = x*-1;
            i = 1; // vel to the right
        }
        if (x > wellLength)
        {
            x = 2*wellLength - x;
            i = -1; // vel to the left
        } 
        KE = ET - BumpPotential(x,time);
        v = velConv(KE)*i;
        // work out new KE and export result
        
        result<<time<<','<<KE<<'\n';

        // update velocity and time
        time += timeStep;
        std::cout<<x<<','<<v<<','<<time<<'\n';
    };

    result.close();
    std::cout<<"EXECUTE SUCESSFUL";
    return 0;
}