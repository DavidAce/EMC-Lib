//
// Created by david on 9/2/16.
//

#include "source/EMC.h"

ArrayXd myfunction(ArrayXXd &dos, ArrayXd &E, ArrayXd &M, ArrayXd &T){

    return T.cwiseAbs2();
}


int main(){

    ArrayXd min_bound(1);
    ArrayXd max_bound(1);

    min_bound <<  -5;
    max_bound << 5;
    int parameters = 1;
    minimization(myfunction,parameters, min_bound, max_bound);
    return 0;
}