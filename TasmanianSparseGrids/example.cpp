#include <iostream>
#include "math.h"

#include "TasmanianSparseGrid.hpp"

using namespace std;

int main( int argc, const char ** argv ){
        
        TasGrid::TasmanianSparseGrid grid;
        
        { // EXAMPLE 1:
                // make a classical Smolyak grid using Clenshaw-Curtis quadrature
                // integrate the function f(x,y) = exp( -x^2 ) * cos( y )
                // the exact answer is: 2.513723354063905e+00
                cout << endl << "EXAMPLE 1" << endl;
                const int dimension = 2;
                const int outputs = 0;
                const int level = 7;
                grid.makeGlobalGrid( dimension, outputs, level, TasGrid::type_level, TasGrid::rule_clenshawcurtis );
                
                double *points = 0;
                double *weights = 0;
                grid.getPoints( points );
                grid.getWeights( weights );
                
                double sum = 0.0;
                for( int i=0; i<grid.getNumPoints(); i++ ){
                        sum += weights[i] * exp( - points[2*i] * points[2*i] ) * cos( points[2*i+1] );
                }
                cout.precision(17);
                cout << "grid has: " << grid.getNumPoints() << " points." << endl;
                cout << "The integral is: " << std::scientific << sum << endl;
                cout.precision(4);
                cout << "          error: " << std::scientific << fabs( sum - 2.513723354063905e+00 ) << endl;
                delete[] points;
                delete[] weights;
        }
        
        { // EXAMPLE 2:
                // make a Clenshaw-Curtis rule that interpolates exactly polynomials of order up to 10
                // integrate the function f(x,y) = exp( -x^2 ) * cos( y )
                // exact values is: 6.990131267703512e-01
                cout << endl << "EXAMPLE 2" << endl;
                const int dimension = 2;
                const int outputs = 1;
                const int precision = 10;
                grid.makeGlobalGrid( dimension, outputs, precision, TasGrid::type_basis, TasGrid::rule_clenshawcurtis );
                
                int N = grid.getNumNeededPoints();
                double *points = 0;
                grid.getNeededPoints( points );
                
                double *values = new double[N];
                
                for( int i=0; i<N; i++ ){
                        values[i] = exp( - points[2*i] * points[2*i] ) * cos( points[2*i+1] );
                }
                grid.loadNeededPoints( values );
                
                double I;
                
                double x[2] = { 0.3, 0.7 }; // a point of interest
                grid.evaluate( x, &I );
                cout.precision(17);
                cout << "grid has: " << grid.getNumPoints() << " points." << endl;
                cout << "The value of the interpolant at (0.3,0.7) is: " << std::scientific << I << endl;
                cout.precision(4);
                cout << "                                       error: " << std::scientific << fabs( I - 6.990131267703512e-01 ) << endl;
                
                grid.integrate( &I );
                cout.precision(17);
                cout << "The integral is: " << std::scientific << I << endl;
                cout.precision(4);
                cout << "          error: " << std::scientific << fabs( I - 2.513723354063905e+00 ) << endl;
                delete[] values;
                delete[] points;
        }
        
        { // EXAMPLE 3:
                // make a quadrature that used Gauss-Legendre points and integrates exactly polynomials up to order 10
                // integrate the function f(x,y) = exp( -x^2 ) * cos( y )
                cout << endl << "EXAMPLE 3" << endl;
                const int dimension = 2;
                const int outputs = 0;
                const int precision = 10;
                grid.makeGlobalGrid( dimension, outputs, precision, TasGrid::type_basis, TasGrid::rule_gausslegendre );
                
                double *points = 0;
                double *weights = 0;
                grid.getPoints( points );
                grid.getWeights( weights );
                
                double sum = 0.0;
                for( int i=0; i<grid.getNumPoints(); i++ ){
                        sum += weights[i] * exp( - points[2*i] * points[2*i] ) * cos( points[2*i+1] );
                }
                cout.precision(17);
                cout << "grid has: " << grid.getNumPoints() << " points." << endl;
                cout << "The integral is: " << std::scientific << sum << endl;
                cout.precision(4);
                cout << "          error: " << std::scientific << fabs( sum - 2.513723354063905e+00 ) << endl;
                delete[] points;
                delete[] weights;
        }
        
        { // EXAMPLE 4:
                // make a quadrature that uses Gauss-Gegenbauer points and uses 8 times more points in y direction
                // integrate the function f(x,y) = ( x - 2 )^3 * exp( -y^2 ) * ( 1 - x^2 )^0.4 * ( 1 - y^2 )^0.4
                // exact integral is: -2.029979511486524e+01
                cout << endl << "EXAMPLE 4" << endl;
                const int dimension = 2;
                const int outputs = 0;
                const int depth = 16;
                const int anisotropic_weights[2] = { 8, 1 };
                const double alpha = 0.4;
                grid.makeGlobalGrid( dimension, outputs, depth, TasGrid::type_level, TasGrid::rule_gaussgegenbauer, anisotropic_weights, &alpha );
                
                double *points = 0;
                double *weights = 0;
                grid.getPoints( points );
                grid.getWeights( weights );
                
                double sum = 0.0;
                for( int i=0; i<grid.getNumPoints(); i++ ){
                        sum += weights[i] * ( points[2*i] - 2.0 )*( points[2*i] - 2.0 )*( points[2*i] - 2.0 )* exp( -points[2*i+1]*points[2*i+1] );
                }
                cout.precision(17);
                cout << "grid has: " << grid.getNumPoints() << " points." << endl;
                cout << "The integral is: " << std::scientific << sum << endl;
                cout.precision(4);
                cout << "          error: " << std::scientific << fabs( sum + 2.029979511486524e+01 ) << endl;
                delete[] points;
                delete[] weights;
        }
        
        { // EXAMPLE 5:
                // interpolates the function f(x,y) = exp( -x^2 ) * cos( y ) 
                // using adaptive piece-wise local quadratic polynomials over [0,1] x [0,1]
                // integrate the interpolant
                // the exact integral is: 6.990131267703512e-01
                cout << endl << "EXAMPLE 5" << endl;
                const int dimension = 2;
                const int outputs = 1;
                const int initial_level = 4;
                const int order = 2;
                grid.makeLocalPolynomialGrid( dimension, outputs, initial_level, order, TasGrid::rule_pwpolynomial );
                
                const double a[2] = { 0.0, 0.0 };
                const double b[2] = { 1.0, 1.0 };
                grid.setTransformAB( a, b ); // the transformation is simply the new boundaries of the domain
                
                const double tolerance = 1.E-6;
                
                int iteration = 0;
                
                while ( grid.getNumNeededPoints() > 0 ){
                        int N = grid.getNumNeededPoints();
                        
                        double *points = 0;
                        grid.getNeededPoints( points );
                        
                        double *values = new double[N];
                        
                        for( int i=0; i<N; i++ ){
                                values[i] = exp( - points[2*i] * points[2*i] ) * cos( points[2*i+1] );
                        }
                        
                        grid.loadNeededPoints( values );
                        
                        double I;
                        double x[2] = { 0.3, 0.7 };
                        grid.evaluate( x, &I );
                        if ( grid.getNumPoints() > 99 ){ // lazy formatting, I should probably use setw() here
                                cout << "Iteration: " << iteration << "  number of samples: " << grid.getNumPoints() << "  value: " << I << "  error: " << fabs( I - 6.990131267703512e-01 ) << endl;
                        }else{
                                cout << "Iteration: " << iteration << "  number of samples:  " << grid.getNumPoints() << "  value: " << I << "  error: " << fabs( I - 6.990131267703512e-01 ) << endl;
                        }
                        
                        grid.setRefinement( tolerance, TasGrid::refine_classic );
                        
                        iteration++;
                        delete[] values;
                        delete[] points;
                }        
        }
        
        cout << endl;
        
        return 0;
}
