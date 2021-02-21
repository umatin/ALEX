#pragma once

#include <vector>


/*
    Fits a linear interpolation function to 2 points
*/
void fit_line(std::vector<double>& data, long fr, long to, double & a, double & b){
    a = (double)(to-fr) / (double)(data[to] - data[fr]);
    b = (double) fr - ((double)data[fr]) * a;
}

/*
    Fits a linear interpolation function to 2 points
*/
void fit_line(std::vector<double>& datax, std::vector<double>& datay, long fr, long to, double & a, double & b){
    a = (double)(datay[to] - datay[fr]) / (double)(datax[to] - datax[fr]);
    b = (double) datay[fr] - ((double)datax[fr]) * a;
}
