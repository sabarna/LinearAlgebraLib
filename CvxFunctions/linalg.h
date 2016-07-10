//
//  linalg.h
//  CvxFunctions
//
//  Created by Sabarna Choudhuri on 7/9/16.
//  Copyright © 2016 Sabarna Choudhuri. All rights reserved.
//

#include "linalg.cpp"


#ifndef linalg_h
#define linalg_h



/**
 Usage :
 Define a double differentiable single variable function as :
 -----------------------------------
 double check_test (double x)
 {
    return  x*x + 3 *x*x*x + 2 ;
 }
 
 
 
 Get numerically solved double differential,evaluated at x, of the function as :
 double x = 6;
 double z = doubleDifferential<double>(check_test,x);
 **/

template<typename F>
inline double doubleDifferential(F fxn, const double x);


/**
 Usage :
 Define a differentiable single variable function as :
 -----------------------------------
 double check_test (double x)
 {
 return  x*x + 3 *x*x*x + 2 ;
 }
 
 
 
 Get numerically solved differential,evaluated at x, of the function as :
 double x = 6;
 double z = doubleDifferential<double>(check_test,x);
 **/

template<typename F>
inline double differential(F fxn , const double x);

/**
 Usage :
 Define a differentiable single variable function of a vector as :
 -----------------------------------
 double check_test (vector<double> x)
 {
    return  x[0]*x[1] + 3 *x[2] *x[1]*x[1] + 2*x[0] ;
 }
 
 Get numerically solved differential,evaluated at x, of the function as :
 double x = 6;
 double z = gradient<double>(check_test,x);
 **/

template<typename F>
inline double gradient(F fxn , std::vector<double> x);


/**
 Usage :
 double x;
 bool sign = sgn(x);
 **/
template<typename T>
int sgn( T val);


/**
 Usage :
 double check_test (double x)
 {
    return 2*x*x + 3x + 2 ;
 }
 double from = -5;
 double to = 5;
 double z = bisectionMethod<double>(check_test, from, to);
 **/

template<typename F>
inline double bisectionMethod(F fxn, double from , double to );

/**
 Usage :
 double check_test (double x)
 {
    return 2*x*x + 3x + 2 ;
 }
 double intiVal = -5
 double z = bisectionMethod<double>(check_test, intiVal);
 **/


template<typename F>
inline double newtonRootMethod(F fxn, double initVal);



/**
 Usage :
 Define a multivariable function as :
 -----------------------------------
 double check_test (std::vector<double> vec )
 {
 return  vec[0] * vec[1] + 2 * vec[2] * vec[0];// xy + 2xz
 }
 
 Create a vector variable a = [x,y,z]
 -------------------------------------
 std::vector<double> vector;
 vector.push_back(1);
 vector.push_back(2);
 vector.push_back(3);
 
 Get Hessian matrix hess of the function evaluated at [1,2,3
 
std::vector <std::vector<double> > hess = getHessian <std::vector <std::vector<double> >> (check_test,vector);
 
 **/


template<typename F>
void hessian1(F fxn,std::vector<double> x, std::vector <std::vector<double> >& hess );

template<typename F,typename F1>
F getHessian(F1 fxn,std::vector<double> x);

/**
 Usage:
std::vector <std::vector<double> > checkvec;
double matDim = 4  ;
checkvec.resize(matDim);
for(unsigned int i = 0; i < checkvec.size(); ++i)
checkvec[i].resize(4,0.0);

checkvec[0][0] = 1;
checkvec[0][1] = 2;
checkvec[0][2] = 3;
checkvec[0][3] = 4;


checkvec[1][0] = 5;
checkvec[1][1] = 6;
checkvec[1][2] = 7;
checkvec[1][3] = 8;

checkvec[2][0] = 9;
checkvec[2][1] = 10;
checkvec[2][2] = 11;
checkvec[2][3] = 12;

checkvec[3][0] = 13;
checkvec[3][1] = 14;
checkvec[3][2] = 15;
checkvec[3][3] = 16;

getTransposeSqMat<std::vector <std::vector<double>>>(checkvec);**/

template<typename F>
F getTransposeSqMat(std::vector <std::vector<double> > ipVec);


/**
 Usage :
 std::vector <std::vector<double> > checkvec;
 double matDim = 4  ;
 checkvec.resize(matDim);
 for(unsigned int i = 0; i < checkvec.size(); ++i)
 checkvec[i].resize(4,0.0);
 
 checkvec[0][0] = 18;
 checkvec[0][1] = 22;
 checkvec[0][2] = 54;
 checkvec[0][3] = 42;
 
 
 checkvec[1][0] = 22;
 checkvec[1][1] = 70;
 checkvec[1][2] = 86;
 checkvec[1][3] = 62;
 
 checkvec[2][0] = 54;
 checkvec[2][1] = 86;
 checkvec[2][2] = 174;
 checkvec[2][3] = 134;
 
 checkvec[3][0] = 42;
 checkvec[3][1] = 62;
 checkvec[3][2] = 134;
 checkvec[3][3] = 106;
 
 
 getCholesky<std::vector <std::vector<double> >>(checkvec);**/
template<typename F>
F getCholesky(std::vector <std::vector<double> > ipMat);



/*   Usage:
 std::vector <std::vector<double> > checkvec;
 double matDim = 4  ;
 checkvec.resize(matDim);
 for(unsigned int i = 0; i < checkvec.size(); ++i)
 checkvec[i].resize(4,0.0);
 
 checkvec[0][0] = 3;
 checkvec[0][1] = 0;
 checkvec[0][2] = 2;
 checkvec[0][3] = -1;
 
 
 checkvec[1][0] = 1;
 checkvec[1][1] = 2;
 checkvec[1][2] = 0;
 checkvec[1][3] = -2;
 
 checkvec[2][0] = 4;
 checkvec[2][1] = 0;
 checkvec[2][2] = 6;
 checkvec[2][3] = -3;
 
 checkvec[3][0] = 5;
 checkvec[3][1] = 0;
 checkvec[3][2] = 2;
 checkvec[3][3] = 0;
 
 std::cout << detFour(checkvec) <<std::endl;*/

double detTwo(std::vector <std::vector<double> > ipVec);

double detThree(std::vector <std::vector<double> > ipVec);


double detFour(std::vector <std::vector<double> > ipVec);




#endif /* linalg_h */
