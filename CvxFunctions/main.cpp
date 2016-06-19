//
//  main.cpp
//  CvxFunctions
//
//  Created by Sabarna Choudhuri on 6/5/16.
//  Copyright © 2016 Sabarna Choudhuri. All rights reserved.
//

#include <iostream>
#include <vector>
#include <math.h>
#include "calc.h"

double check_test (std::vector<double> vec )
{
    return  vec[0] * vec[1] + 2 * vec[2] * vec[0];
}






int main()
{
    std::vector<double> vector;
    vector.push_back(1);
    vector.push_back(2);
    vector.push_back(3);
    
    std::vector <std::vector<double> > hess = getHessian <std::vector <std::vector<double> >> (check_test,vector);
    
    for(unsigned int i = 0 ; i < hess.size(); ++i)
    {
        for(unsigned int j = 0 ; j < hess[i].size(); ++j)
        {
            std::cout<< "i " << i <<"j " << j <<std::endl;
            std::cout << hess[i][j] << std::endl;
        }
        
    }
    
    
    return 0;
}


