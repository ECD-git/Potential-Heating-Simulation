#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>

// Define our constants
float wellLength = 0.1; //[m]
float wellDepth = 0.5; // [kelvin]

float bumpAmp = 0.01; // [kelvin]
float bumpLength = 0.01; //[m]
float bumpPos = (wellLength/2) - (bumpLength/2); //[m] //places bump in middle
// this will need to be varied
float bumpOmega = 0; // [rad/s] //driving frequency of potential 'bump'
float bumpPhase = 0; //[rads]

float v_0 = 0.02; //[m/s] //a few cm a second seems about right
float M = 1.00784*1.66*pow(10,-27); //[kg] //using atomic mass of H
float x_0 = 0; //[m] //particle starts at far left side 

float timeStep = 0.1; //[s]
float stopTime = 10; //[s] //Overall Length of sim

// useful functions
float KEConv(float v, float M)
{
    return 0.5*M*pow(v,2);
};
float velConv(float KE, float M)
{
    return pow((2*KE)/M, 0.5);
};
float eVConv(float KE)
{
    return KE/(1.6*pow(10,-19));
}
float KelEConv(float K)
{
    return 1.380649*pow(10, -23)*K;
}
float BumpPotential(float x, float time)
{
    // modeled as a sin wave starting from the left side of the bump pos
    // with a half wavelength equal to the bump length
    // designed to take raw x value, auto offsets by the bump position
    return KelEConv(bumpAmp) * std::sin((M_PI/bumpLength)*(x-bumpPos) + (bumpOmega*time) + bumpPhase);
}

float GetPotential(float x, float time)
{
    // check if over bump
    if (x > bumpPos && x < bumpPos + bumpLength)
    {
        std::cout<<"Over bump ";
        return BumpPotential(x,time);
    }
    else
    {
        return 0;
    }
}

int main()
{
    float time = 0;
    float x = x_0;
    float v = v_0;
    float KE = KEConv(v_0,M); // init KE [j]
    float ET = KE + GetPotential(x,time); //[J] Get E_total at x_0 and time=0

    //open file write our initials 
    std::ofstream result;
    result.open("results.dat");
    result<<time<<','<<eVConv(KE)<<'\n';

    while(time < stopTime)
    {   
        x += v*timeStep; // new x pos
        // check if out of bounds
        if (x < 0)
        {
            x = x*-1;
            v = abs(v);
        }
        if (x > wellLength)
        {
            x = 2*wellLength - x;
            v = abs(v)*-1;
        }
        // work out new KE and export result
        KE = ET - GetPotential(x,time);
        result<<time<<','<<eVConv(KE)<<'\n';

        time += timeStep;
        std::cout<<x<<','<<time<<'\n';
    };

    result.close();
    std::cout<<"EXECUTE SUCESSFUL";
    return 0;
}