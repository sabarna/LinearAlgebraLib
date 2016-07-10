//
//  Header.h
//  CvxFunctions
//
//  Created by Sabarna Choudhuri on 6/5/16.
//  Copyright © 2016 Sabarna Choudhuri. All rights reserved.
//

#ifndef Header_h
#define Header_h

    /**
     * returns the double differential of a function fxn w.r.t. variable x
     * input parameters are the function fxn and variable x
     **/
    
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
    
    
    /**
     * returns the differential of a function fxn w.r.t. variable x
     * input parameters are the function fxn and variable x
     **/
    
    template<typename F>
    inline double differential(F fxn , const double x)
    {
        double eps = 0.001;
        return (4/3) * (fxn(x + eps) - fxn(x - eps))/(2*eps) - (1/3) * (fxn(x + 2 * eps) - fxn(x - 2 * eps))/(4*eps);
    }
    
    
    /**
     * returns the gradient of a function fxn w.r.t. variable x; where x is a vector
     * input parameters are the function fxn and variable x
     **/
    
    template<typename F>
    inline double gradient(F fxn , std::vector<double> x)
    {
        double eps = 0.001;
        std::vector<double> incr_vec;
        incr_vec.resize(x.size(),0.0);
        std::vector<double> decr_vec;
        decr_vec.resize(x.size(),0.0);
        std::vector<double> incr_2vec;
        incr_2vec.resize(x.size(),0.0);
        std::vector<double> decr_2vec;
        decr_2vec.resize(x.size(),0.0);
        
        for (unsigned int i = 0; i < x.size(); i++)
        {
            double incr = 0 ;
            double decr = 0 ;
            double db_incr = 0;
            double db_decr = 0;
            incr = x[i] + eps;
            decr = x[i] - eps;
            db_incr = x[i] + 2 *eps;
            db_decr = x[i] - 2 * eps;
            
            incr_vec.push_back(incr);
            decr_vec.push_back(decr);
            
            incr_2vec.push_back(db_incr);
            decr_2vec.push_back(db_decr);
            
        }
        
        return (4/3) * (fxn(incr_vec) - fxn(decr_vec))/(2*eps) - (1/3) * (fxn(incr_2vec) - fxn(decr_2vec))/(4*eps);
    }
    
    
    /**
     * returns 1 if input value is positive else 0
     * input parameter is the input value
     **/
    
    template<typename T>
    int sgn( T val)
    {
        return (T(0) < val) - (val < T(0));
    }
    
    /**
     * returns the root of a function
     * input parameters are the function fxn and the search space limits
     **/
    
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
    
    /**
     * returns the root of a function fxn using Newton Raphson's method
     * input parameters are the function fxn and the initial search value
     * output parameter is the root value
     **/
    
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
    
    /**
     * populates the Hessian matrix of a function fxn
     * input parameters are the function fxn, variable vector x and the Hessian matrix of the function
     **/
    
    template<typename F>
    void hessian1(F fxn,std::vector<double> x, std::vector <std::vector<double> >& hess )
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
    
    /**
     * returns the Hessian matrix of a function fxn
     * input parameters are the function fxn and variable vector x
     **/
    
    template<typename F,typename F1>
    F getHessian(F1 fxn,std::vector<double> x)
    {
        std::vector <std::vector<double> > hess;
        double matDim = x.size()  ;
        hess.resize(matDim);
        for(unsigned int i = 0 ; i < hess.size(); ++i)
            hess[i].resize(matDim);
        hessian1(fxn,x, hess);
        return hess;
    }
    
    
    /**
     * returns the Cholesky decomposition of a matrix
     * input parameter is the matrix whose Cholesky decomposition needs to be performed
     **/
    
    template<typename F>
    F getCholesky(std::vector <std::vector<double> > ipMat)
    {
        std::vector <std::vector<double> > decomp;
        double matDim = ipMat.size()  ;
        decomp.resize(matDim);
        for(unsigned int i = 0; i < decomp.size(); ++i)
            decomp[i].resize(matDim,0.0);
        
        for(unsigned int j = 0; j < ipMat.size(); ++j)
        {
            for (unsigned int i = 0; i < ipMat.size(); ++i)
            {
                std::cout << "i = " << i << "j = " << j << std::endl;
                if(i == j)
                {
                    double lSum1 = 0.0;
                    for (unsigned int k = 0 ; k < j ;k++ )
                        lSum1 += decomp[i][k] * decomp[i][k] ;
                    decomp[i][j] = sqrt(ipMat[i][j] - lSum1);
                    std::cout << decomp[i][j] << std::endl;
                }
                else if (j > i)
                {
                    decomp[i][j] = 0.0;
                    std::cout << decomp[i][j] << std::endl;
                }
                else
                {
                    double lSum2 = 0.0;
                    for(unsigned int k = 0; k < j ; k++ )
                        lSum2 += decomp[i][k] * decomp[j][k];
                    decomp[i][j] = (1 / decomp[j][j]) * (ipMat[i][j] - lSum2);
                    std::cout << decomp[i][j] << std::endl;
                }
                
            }
        }
        
        for(unsigned int i = 0; i < decomp.size(); ++i)
        {
            for (unsigned int j = 0; j < decomp.size(); ++j)
            {
                std::cout << decomp[i][j] << std::endl;
            }
        }
        
        return decomp;
        
    }
    
    /**
     * returns the determinant of a 2*2 matrix
     * input parameter is the matrix
     **/
    
    double detTwo(std::vector <std::vector<double> > ipVec)
    {
        return ipVec[0][0] * ipVec[1][1] - ipVec[0][1] * ipVec[1][0] ;
    }
    
    
    /**
     * returns the determinant of a 3*3 matrix
     * input parameter is the matrix
     **/
    double detThree(std::vector <std::vector<double> > ipVec)
    {
        
        std::vector <std::vector<double> > ipVecTemp1;
        ipVecTemp1.resize(2);
        for(unsigned int iterTemp = 0 ; iterTemp < 2; iterTemp++)
            ipVecTemp1[iterTemp].resize(2);
        
        std::vector <std::vector<double> > ipVecTemp2;
        ipVecTemp2.resize(2);
        for(unsigned int iterTemp = 0 ; iterTemp < 2; iterTemp++)
            ipVecTemp2[iterTemp].resize(2);
        
        std::vector <std::vector<double> > ipVecTemp3;
        ipVecTemp3.resize(2);
        for(unsigned int iterTemp = 0 ; iterTemp < 2; iterTemp++)
            ipVecTemp3[iterTemp].resize(2);
        
        
        ipVecTemp1[0][0] = ipVec[1][1];
        ipVecTemp1[0][1] = ipVec[1][2];
        ipVecTemp1[1][0] = ipVec[2][1];
        ipVecTemp1[1][1] = ipVec[2][2];
        
        ipVecTemp2[0][0] = ipVec[1][0];
        ipVecTemp2[0][1] = ipVec[1][2];
        ipVecTemp2[1][0] = ipVec[2][0];
        ipVecTemp2[1][1] = ipVec[2][2];
        
        ipVecTemp3[0][0] = ipVec[1][0];
        ipVecTemp3[0][1] = ipVec[1][1];
        ipVecTemp3[1][0] = ipVec[2][0];
        ipVecTemp3[1][1] = ipVec[2][1];
        
        
        return ipVec[0][0] * detTwo(ipVecTemp1) - ipVec[0][1] * detTwo(ipVecTemp2) + ipVec[0][2] * detTwo(ipVecTemp3);
        
    }
    
    /**
     * returns the determinant of a 4*4 matrix
     * input parameter is the matrix
     **/
    
    double detFour(std::vector <std::vector<double> > ipVec)
    {
        std::vector <std::vector<double> > ipVecTemp1;
        ipVecTemp1.resize(3);
        for(unsigned int iterTemp = 0 ; iterTemp < 3; iterTemp++)
            ipVecTemp1[iterTemp].resize(3);
        
        std::vector <std::vector<double> > ipVecTemp2;
        ipVecTemp2.resize(3);
        for(unsigned int iterTemp = 0 ; iterTemp < 3; iterTemp++)
            ipVecTemp2[iterTemp].resize(3);
        
        std::vector <std::vector<double> > ipVecTemp3;
        ipVecTemp3.resize(3);
        for(unsigned int iterTemp = 0 ; iterTemp < 3; iterTemp++)
            ipVecTemp3[iterTemp].resize(3);
        
        std::vector <std::vector<double> > ipVecTemp4;
        ipVecTemp4.resize(3);
        for(unsigned int iterTemp = 0 ; iterTemp < 3; iterTemp++)
            ipVecTemp4[iterTemp].resize(3);
        
        ipVecTemp1[0][0] = ipVec[1][1];
        ipVecTemp1[0][1] = ipVec[1][2];
        ipVecTemp1[0][2] = ipVec[1][3];
        ipVecTemp1[1][0] = ipVec[2][1];
        ipVecTemp1[1][1] = ipVec[2][2];
        ipVecTemp1[1][2] = ipVec[2][3];
        ipVecTemp1[2][0] = ipVec[3][1];
        ipVecTemp1[2][1] = ipVec[3][2];
        ipVecTemp1[2][2] = ipVec[3][3];
        
        
        ipVecTemp2[0][0] = ipVec[1][0];
        ipVecTemp2[0][1] = ipVec[1][2];
        ipVecTemp2[0][2] = ipVec[1][3];
        ipVecTemp2[1][0] = ipVec[2][0];
        ipVecTemp2[1][1] = ipVec[2][2];
        ipVecTemp2[1][2] = ipVec[2][3];
        ipVecTemp2[2][0] = ipVec[3][0];
        ipVecTemp2[2][1] = ipVec[3][2];
        ipVecTemp2[2][2] = ipVec[3][3];
        
        
        ipVecTemp3[0][0] = ipVec[1][0];
        ipVecTemp3[0][1] = ipVec[1][1];
        ipVecTemp3[0][2] = ipVec[1][3];
        ipVecTemp3[1][0] = ipVec[2][0];
        ipVecTemp3[1][1] = ipVec[2][1];
        ipVecTemp3[1][2] = ipVec[2][3];
        ipVecTemp3[2][0] = ipVec[3][0];
        ipVecTemp3[2][1] = ipVec[3][1];
        ipVecTemp3[2][2] = ipVec[3][3];
        
        
        ipVecTemp4[0][0] = ipVec[1][0];
        ipVecTemp4[0][1] = ipVec[1][1];
        ipVecTemp4[0][2] = ipVec[1][2];
        ipVecTemp4[1][0] = ipVec[2][0];
        ipVecTemp4[1][1] = ipVec[2][1];
        ipVecTemp4[1][2] = ipVec[2][2];
        ipVecTemp4[2][0] = ipVec[3][0];
        ipVecTemp4[2][1] = ipVec[3][1];
        ipVecTemp4[2][2] = ipVec[3][2];
        
        return ipVec[0][0] * detThree(ipVecTemp1) - ipVec[0][1] * detThree(ipVecTemp2) + ipVec[0][2] * detThree(ipVecTemp3) - ipVec[0][3] * detThree(ipVecTemp4);
        
    }
    
    
    
    /**
     * returns the transpose of a square matrix
     * input parameter is the matrix
     **/
    
    template<typename F>
    F getTransposeSqMat(std::vector <std::vector<double> > ipVec)
    {
        std::vector <std::vector<double> > opVec;
        opVec.resize(ipVec.size());
        for (unsigned int iterTemp = 0 ; iterTemp < ipVec.size(); iterTemp++)
            opVec[iterTemp].resize(ipVec.size());
        
        for(unsigned int i = 0 ; i < opVec.size(); i++)
        {
            for(unsigned int j = 0 ; j < opVec.size(); j++)
                opVec[i][j] = ipVec[j][i];
        }
        
        for(unsigned int i = 0 ; i < opVec.size(); i++)
        {
            for(unsigned int j = 0 ; j < opVec.size(); j++)
                std::cout <<opVec[i][j] <<std::endl;
        }
        
        
        return opVec;
        
    }
    
    /*
     template<typename F>
     F getMatProd(std::vector <std::vector<double> > ipVec1, std::vector <std::vector<double> > ipVec2)
     {
     size_t nRows1 = ipVec1.size();
     size_t nCols1 = ipVec1[0].size();
     
     size_t nRows2 = ipVec2.size();
     size_t nCols2 = ipVec2[0].size();
     
     if(nCols1 != nRows2)
     {
     std::cout << "Dimension Mismatch error" << std::endl;
     exit(1);
     }
     
     
     
     
     }
     
     std::vector<double> update_var(std::vector<double> x)
     {
     std::vector<double> step_size;
     step_size.resize(x.size(),0);
     
     }
     
     
     template<typename F>
     inline double doubleDifferential(F fxn, std::vector<double> init_vec)
     {
     int iter = 0;
     double eps  = 0.000001;
     double grd = gradient(fxn , init_vec);
     std::vector<double> search_var;
     search_var = init_vec;
     std::vector <std::vector<double> > hess ;
     double desc = 0.0001;
     
     while(desc != 0)
     {
     
     grd = gradient(fxn , search_var);
     hess = getHessian <std::vector <std::vector<double> >> (fxn,search_var);
     //desc = ;
     
     
     }
     
     
     
     
     
     }*/


#endif /* Header_h */
