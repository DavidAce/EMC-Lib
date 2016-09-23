#ifndef RANDOMFUNCTIONS_HPP   // if x.h hasn't been included yet
#define RANDOMFUNCTIONS_HPP   //  #define this so the compiler knows it has been included
#include <random>
#include <iostream>
#include "constants.hpp"

namespace EMC_rnd{

    //Random functions
    extern std::mt19937 rng;
    extern std::uniform_int_distribution<>  rand_int_1;
    extern std::uniform_real_distribution<> rand_real_1;


    inline int uniform_integer_1(){
        return rand_int_1(rng);
    }


    inline double uniform_double_1(){
        return rand_real_1(rng);
    }

    template <typename inType>
    inline inType uniform_integer(const inType lowerLimit, const inType upperLimit ){
        std::uniform_int_distribution<>  rand_int(lowerLimit,upperLimit);
        return rand_int(rng);
    }

    template <typename inType>
    inline inType uniform_double(const inType lowerLimit, const inType upperLimit ){
        std::uniform_real_distribution<>  rand_real(lowerLimit,upperLimit);
        return rand_real(rng);
    }

    template <typename inType1, typename inType2>
    inline inType1 gaussian_truncated(const inType1 lowerLimit, const inType1 upperLimit, const inType2 mean, const inType2 std){
        std::normal_distribution<double> distribution(mean,std);
        inType1 ul = fmaxl(lowerLimit, upperLimit);
        inType1 ll = fminl(lowerLimit, upperLimit);
        inType1 number;
        while (true) {
//            std::cout << "low: "<< lowerLimit << " high: " << upperLimit << " mean: " << mean << " std: "<< std << endl;

            number = distribution(rng);
            if (number >= ll && number <= ul) {
                return number;
            }
        }
    }


    template <typename inType1, typename inType2>
    inline double gaussian(const inType1 mean, const inType2 std){
        std::normal_distribution<double> distribution(mean,std);
        return distribution(rng);
    }

    inline int triangular_increasing(const int lowerLimit, const int upperLimit ){
        return (int) ( lowerLimit + (upperLimit - lowerLimit)*sqrt(uniform_double_1()));
    }

    inline int triangular_decreasing(const int lowerLimit, const int upperLimit ){
        return (int) ( upperLimit - (upperLimit-lowerLimit) * sqrt(1-uniform_double_1()));
    }



};

#endif