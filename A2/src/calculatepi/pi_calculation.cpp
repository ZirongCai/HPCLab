#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <chrono>
#include <omp.h>

#include "Stopwatch.h"
#define PRECESION 1000
#define PI 3.14159265358979323846  
#define LOOP 1000

//#define NUM_THREAD 8
using namespace std::chrono;

double pi_integrate_serial(double start, double end, int precision);
double pi_integrate_critical(double start, double end, int precision);
double pi_integrate_reduction(double start, double end, int precision);

bool checkPi(double result)
{
    if(abs(result - PI) < 0.001)
    {
        return true;
    }
    return false;
}

int main(int argc, char **argv)
{
    int precesion = PRECESION;
    if (argc > 1)
	{
		precesion = atoi(argv[1]);
	}

    double aver_time = 0;
    auto start = steady_clock::now();
    for (size_t i = 0; i < LOOP; i++)
    {
        if(!checkPi(pi_integrate_serial(0, 1.0, precesion)))
        {
            printf("Calculation wrong");
        }
    }
    auto end = steady_clock::now();
    duration<double> elapse_time = end - start;
    aver_time = elapse_time.count()/LOOP;
    std::cout <<"time spent in serial : " << aver_time << " sec" << std::endl;

    start = steady_clock::now();
    for (size_t i = 0; i < LOOP; i++)
    {
        if(!checkPi(pi_integrate_critical(0, 1.0, precesion)))
        {
            printf("Calculation wrong");
        }
    }
    end = steady_clock::now();
    elapse_time = end - start;
    aver_time = elapse_time.count()/LOOP;
    std::cout << "time spent in critical : " << aver_time << " sec" << std::endl;

    
    start = steady_clock::now();
    for (size_t i = 0; i < LOOP; i++)
    {
        if(!checkPi(pi_integrate_reduction(0, 1.0, precesion)))
        {
            printf("Calculation wrong");
        }       
    }
    end = steady_clock::now();
    elapse_time = end - start;
    aver_time = elapse_time.count()/LOOP;
    std::cout << "time spent in reduction : " << aver_time << " sec" << std::endl;



}

double pi_func(double x)
{
    return 1 / (1 + x * x );
}

double pi_integrate_serial(double start, double end, int precision)
{
    double length = (end - start ) / precision;
    double x = length/2.0;
    double sum = 0;
    for (size_t i = 0; i < precision; i++)
    {
        double tmp = pi_func(x + i*length);
        sum += tmp;
    }
    return 4 * sum * length;
    
}

double pi_integrate_critical(double start, double end, int precision)
{
    double length = (end - start ) / precision;
    double x = length/2.0;
    double sum = 0;
    //omp_set_num_threads(NUM_THREAD);
    #pragma omp parallel for shared(x, length, sum)
    for (size_t i = 0; i < precision; i++)
    {
        double tmp = pi_func(x + i*length);
        #pragma omp critical
        sum += tmp;
    }
    return 4 * sum * length;
    
}

double pi_integrate_reduction(double start, double end, int precision)
{
    double length = (end - start ) / precision;
    double x = length/2.0;
    double sum = 0;
    //omp_set_num_threads(NUM_THREAD);
    #pragma omp parallel for shared(x, length), reduction(+: sum)
    for (size_t i = 0; i < precision; i++)
    {
        double tmp = pi_func(x + i*length);
        sum += tmp;
    }
    return 4 * sum * length;
    
}