//
//  Header.h
//  CvxFunctions
//
//  Created by Sabarna Choudhuri on 6/5/16.
//  Copyright © 2016 Sabarna Choudhuri. All rights reserved.
//

#ifndef Header_h
#define Header_h

template<typename F>
inline double doubleDifferential(F fxn, const double x)
{
    double eps = 0.001;
    
    auto differential = [&](F fxn , const double x)
    {
        
        
        return (4/3) * (fxn(x + eps) - fxn(x - eps))/(2*eps) - (1/3) * (fxn(x + 2 * eps) - fxn(x - 2 * eps))/(4*eps);
        
    };
   
     return (4/3) * (differential(fxn, x + eps) - differential(fxn, x - eps))/(2*eps) - (1/3) * (differential(fxn, x + 2 * eps) - differential(fxn,x - 2 * eps))/(4*eps) ;
    
   
}

template<typename F>
inline double differential(F fxn , const double x)
{
    double eps = 0.001;
    return (4/3) * (fxn(x + eps) - fxn(x - eps))/(2*eps) - (1/3) * (fxn(x + 2 * eps) - fxn(x - 2 * eps))/(4*eps);
}

template<typename T>
int sgn( T val)
{
    return (T(0) < val) - (val < T(0));
}

template<typename F>
inline double bisectionMethod(F fxn, double from = -1e7, double to = 1e7)
{
    double eps = 1e-7;
    double tol = 1e-7;
    if(from >= to)
    {
        std::cout<<"lower boundary of search space >= upper boundary of search space"<<std::endl;
        return 99999999999;
    }
    if((fxn(from) >0 && fxn(to) > 0) || (fxn(from) <0 && fxn(to) < 0))
    {
        std::cout<<"Invalid search space"<<std::endl;
        return 99999999999;
    }
    
    int maxIter = 1000000;
    double a  = from;
    double b = to;
    double c = 0;
    double checkTol = 0;
    double iter = 1;
    while(iter <= maxIter)
    {
        c = (a+b)/2;
        checkTol = (b - a)/2;
        double fc = fxn(c);
        double fa = fxn(a);
        if(fc == 0 || checkTol < tol)
            break ;// solution found
        iter ++;
        
        if (sgn(fc) == sgn(fa))
            a = c;
        else
            b = c;
    }
    return c;
}

template<typename F>
inline double newtonRootMethod(F fxn, double initVal)
{
    
    double tol = 1e-7;
    double x{0};
    double diffF = differential(fxn , initVal);
    x = initVal - fxn(initVal)/diffF;
    while (fxn(x) - 0 > tol)
    {
        diffF = differential(fxn , x);
        x = x - fxn(x)/diffF;
    }
    return x;
}



template<typename F>

void doubleDifferential1(F fxn,std::vector<double> x, std::vector <std::vector<double> >& hess )
{
    
    double eps = 0.001;
    
    size_t length = x.size();
    
    for (unsigned int i = 0; i < length ; i++ )
    {
        auto differential1 = [&](F fxn , const std::vector<double> x)
        {
            auto a1Plus = x;
            auto a1Minus = x;
            auto b1Plus = x;
            auto b1Minus = x;
            
            a1Plus[i] = a1Plus[i] + eps;
            a1Minus[i] = a1Minus[i] - eps;
            b1Plus[i] = b1Plus[i] + 2 * eps;
            b1Minus[i] = b1Minus[i] - 2 * eps;
            
            return (4/3) * (fxn(a1Plus) - fxn(a1Minus))/(2*eps) - (1/3) * (fxn(b1Plus) - fxn(b1Minus))/(4*eps);
        };
        
        for (unsigned int j = i ; j < length ; j++ )
        {
            auto a2Plus = x;
            auto a2Minus = x;
            auto b2Plus = x;
            auto b2Minus = x;
            
            a2Plus[j] = a2Plus[j] + eps;
            a2Minus[j] = a2Minus[j] - eps;
            b2Plus[j] = b2Plus[j] + 2 * eps;
            b2Minus[j] = b2Minus[j] - 2 * eps;
            
            hess[i][j] = (4/3) * (differential1(fxn, a2Plus) - differential1(fxn, a2Minus))/(2*eps) - (1/3) * (differential1(fxn, b2Plus) - differential1(fxn,b2Minus))/(4*eps);
            hess[j][i] = hess[i][j];
        }
        
    }
    return ;
}


template<typename F,typename F1>
F getHessian(F1 fxn,std::vector<double> x)
{
    std::vector <std::vector<double> > hess;
    double matDim = x.size()  ;
    hess.resize(matDim);
    for(unsigned int i = 0 ; i < hess.size(); ++i)
        hess[i].resize(matDim);
    doubleDifferential1(fxn,x, hess);
    return hess;
}









#endif /* Header_h */
