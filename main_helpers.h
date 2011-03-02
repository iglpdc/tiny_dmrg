/** 
 * @file main_helpers.cpp
 * @brief A couple of small functions needed in the main DMRG file
 * 
 * @author Roger Melko 
 * @author Ivan Gonzalez
 * @date February 9th, 2011
 * 
 */
#ifndef MAIN_HELPERS_H
#define MAIN_HELPERS_H  

#include<iostream>

/**
 * @brief A function to calculate the minimum size of the enviroment
 *
 * @param m the number of states you are keeping
 * @param numberOfSites the number of sites in the whole chain
 *
 * When you are sweeping there is no need to go to very small environment
 * size because you can solve exactly when the environment is small
 * enough.
 */
inline double calculateMinEnviromentSize(int m, int numberOfSites)
{
    int result;
    for (result=3; result<numberOfSites; result++)
	if (powf(2.0,result) >= 2.0*m) break;
    return result;
}

/**
 * @brief A function to print the energy in a fancy way
 *
 * Just prints the thing: sites in the left block, sites in the right
 * block, and energy per site
 */
inline void printGroundStateEnergy(int sitesInLeft, int sitesInRight, 
	double groundStateEnergy)
{
    std::cout<<std::setprecision(16);
    std::cout<<sitesInLeft<<" "<<sitesInRight\
	<<" "<<groundStateEnergy/(sitesInLeft+sitesInRight)<<std::endl;
}
#endif //MAIN_HELPERS_H
