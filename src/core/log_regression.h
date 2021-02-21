#pragma once
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <random>
#include <cassert>
#include <limits>

#include "error_function.h"
#include "fast_random.h"
#include "interpolation.h"


/*
    A brute force algorithm to create the optimal regression function for a dataset
    Only works for error functions in O(error)
*/
template<ERROR_TYPE E, bool ROUND = true, bool BOUNDED = true>
double create_regression_optimal(std::vector<double>& data, double & a, double & b){
    double e = INFINITY;
    for(long fr = 0; fr < data.size(); fr++){
        for(long to = fr+1; to < data.size(); to++){
            double a_tmp;
            double b_tmp;
            fit_line(data, fr, to, a_tmp, b_tmp);
            double error_tmp = calculate_error<E,ROUND,BOUNDED>(data, a_tmp, b_tmp);
            if(error_tmp < e){
                a = a_tmp;
                b = b_tmp;
                e = error_tmp;
            }
        }
    }
    return e;
}

/*
    A recursive helper function for create_regression_tournament_selection
*/
template<ERROR_TYPE E, bool ROUND = true, bool BOUNDED = true>
void create_regression_tournament_selection_recursive(std::vector<double>& data, std::vector<int>& tournament, int fr, int to, double & a, double & b){
    assert(to - fr >= 2);
    //std::cout << "fr = " << fr << "to = " << to << std::endl;
    if (to - fr == 2){
        fit_line(data, tournament[fr], tournament[to-1], a, b);
        return;
    }
    double al, bl, ar, br;
    create_regression_tournament_selection_recursive<E,ROUND,BOUNDED>(data, tournament, fr, (fr+to)/2, al, bl);
    create_regression_tournament_selection_recursive<E,ROUND,BOUNDED>(data, tournament, (fr+to)/2, to, ar, br);
    double errorl = 0;
    for (int i = fr; i < to; i++){
        errorl += calculate_error_single_element<E,ROUND,BOUNDED> (data, al, bl, tournament[i]);
    }
    double errorr = 0;
    for (int i = fr; i < to; i++){
        errorr += calculate_error_single_element<E,ROUND,BOUNDED> (data, ar, br, tournament[i]);
    }
    if (errorl < errorr){
        a = al;
        b = bl;
    }
    else{
        a = ar;
        b = br;
    }
}

/*
    A randomized approximation algorihm for the logarithm error regression.
    Only works properly for Logarithmic errors.
    Works by facing candidtate regression functions off in a tounament to obtain a good candidate.
    Runtime is in O(nlogn)
    The function rounds the size or the data up to the next highest power of 2.
*/
template<ERROR_TYPE E, bool ROUND = true, bool BOUNDED = true>
void create_regression_tournament_selection(std::vector<double>& data, double & a, double & b){
    long s = std::pow(2, (int)std::ceil(log2(data.size())));
    s *= 2;
    std::vector<int> tournament;
    tournament.resize(s);
    for (int i = 0; i < s; i++){
        tournament[i] = rand()%data.size();
    }
    create_regression_tournament_selection_recursive<E,ROUND,BOUNDED>(data, tournament, 0, s, a, b);
}

/*
    A repeated execution of create_regression_tournament_selection
*/
template<ERROR_TYPE E, bool ROUND = true, bool BOUNDED = true>
void create_regression_tournament_selection_repeated(std::vector<double>& data, double & a, double & b, int n){
    double min_error = std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; i++){
        double candidate_a, candidate_b;
        create_regression_tournament_selection<E,ROUND,BOUNDED>(data, candidate_a, candidate_b);
        double candidate_error = calculate_error<E,ROUND,BOUNDED>(data, candidate_a, candidate_b);
        if(candidate_error < min_error){
            min_error = candidate_error;
            a = candidate_a;
            b = candidate_b;
        }
    }
}

template<ERROR_TYPE E, bool ROUND = true, bool BOUNDED = true>
void create_regression_tournament_selection_fast(std::vector<double> & data, int & fr, int & to, int height){
    if(height == 0){
        fr = xorshf96()%data.size();
        to = xorshf96()%data.size();
        while(to == fr){
            to = xorshf96()%data.size();
        }
        return;
    }
    int fr_l, to_l;
    int fr_r, to_r;
    create_regression_tournament_selection_fast<E,ROUND,BOUNDED>(data, fr_l, to_l, height-1);
    create_regression_tournament_selection_fast<E,ROUND,BOUNDED>(data, fr_r, to_r, height-1);
    double a_l,b_l;
    double a_r,b_r;
    fit_line(data,fr_l,to_l,a_l,b_l);
    fit_line(data,fr_r,to_r,a_r,b_r);
    double err_l = 0;
    double err_r = 0;
    err_l += calculate_error_single_element<E,ROUND,BOUNDED>(data, a_l, b_l, fr_r);
    err_l += calculate_error_single_element<E,ROUND,BOUNDED>(data, a_l, b_l, to_r);
    err_r += calculate_error_single_element<E,ROUND,BOUNDED>(data, a_r, b_r, fr_l);
    err_r += calculate_error_single_element<E,ROUND,BOUNDED>(data, a_r, b_r, to_l);
    
    for(int i = 0; i < (2<<(height-1)); i++){
        int ele = xorshf96()%data.size();
        err_l += calculate_error_single_element<E,ROUND,BOUNDED>(data, a_l, b_l, ele);
        err_r += calculate_error_single_element<E,ROUND,BOUNDED>(data, a_r, b_r, ele);
    }

    if(err_l < err_r){
        fr = fr_l;
        to = to_l;
    }
    else{
        fr = fr_r;
        to = to_r;
    }
}

template<ERROR_TYPE E, bool ROUND = true, bool BOUNDED = true>
void create_regression_tournament_selection_fast(std::vector<double> & data, double & a, double & b){
    int fr,to;
    int height = std::ceil(log2(data.size()));
    create_regression_tournament_selection_fast<E,ROUND,BOUNDED>(data, fr, to, height);
    fit_line(data,fr,to,a,b);
}

template<ERROR_TYPE E, bool ROUND = true, bool BOUNDED = true>
void create_regression_tournament_selection(std::vector<double> & datax, std::vector<double> & datay, int & fr, int & to, int height, int min_, int max_){
    if(height == 0){
        fr = xorshf96()%datax.size();
        to = xorshf96()%datax.size();
        while(to == fr){
            to = xorshf96()%datax.size();
        }
        return;
    }
    int fr_l, to_l;
    int fr_r, to_r;
    create_regression_tournament_selection<E,ROUND,BOUNDED>(datax, datay, fr_l, to_l, height-1, min_, max_);
    create_regression_tournament_selection<E,ROUND,BOUNDED>(datax, datay, fr_r, to_r, height-1, min_, max_);
    double a_l,b_l;
    double a_r,b_r;
    fit_line(datax,datay,fr_l,to_l,a_l,b_l);
    fit_line(datax,datay,fr_r,to_r,a_r,b_r);
    double err_l = 0;
    double err_r = 0;
    err_l += calculate_error_single_element<E,ROUND,BOUNDED>(datax, datay, a_l, b_l, fr_r, min_, max_);
    err_l += calculate_error_single_element<E,ROUND,BOUNDED>(datax, datay, a_l, b_l, to_r, min_, max_);
    err_r += calculate_error_single_element<E,ROUND,BOUNDED>(datax, datay, a_r, b_r, fr_l, min_, max_);
    err_r += calculate_error_single_element<E,ROUND,BOUNDED>(datax, datay, a_r, b_r, to_l, min_, max_);
    
    for(unsigned int i = 0; i < (unsigned int)(2<<(height-1)) && i < datax.size(); i++){
        int ele = xorshf96()%datax.size();
        err_l += calculate_error_single_element<E,ROUND,BOUNDED>(datax, datay, a_l, b_l, ele, min_, max_);
        err_r += calculate_error_single_element<E,ROUND,BOUNDED>(datax, datay, a_r, b_r, ele, min_, max_);
    }

    if(err_l < err_r){
        fr = fr_l;
        to = to_l;
    }
    else{
        fr = fr_r;
        to = to_r;
    }
}
