/*
 * Code Author: Miroslav Stoyanov, Mar 2013
 *
 * Copyright (C) 2013  Miroslav Stoyanov
 *
 * This file is part of
 * Toolkit for Adaprive Stochastic Modeling And Non-Intrusive Approximation
 *              a.k.a. TASMANIAN
 *
 * TASMANIAN is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TASMANIAN is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TASMANIAN.  If not, see <http://www.gnu.org/licenses/>
 *
 */

#ifndef __TASMANIAN_SPARSE_GRID_RULE_GJ_CPP
#define __TASMANIAN_SPARSE_GRID_RULE_GJ_CPP

#include "tsgRuleGaussJacobi.hpp"

namespace TasGrid{

RuleGaussJacobi::RuleGaussJacobi( const int level, double walpha, double wbeta ) : max_level(level), tol(NUM_TOL), alpha(walpha), beta(wbeta), OneDRule(){
        int total_points = 0;
        levels = new int[max_level+1]; levels[0] = 0;
        for( int l=0; l<max_level; l++ ){
                levels[l+1] = levels[l] + getNumPoints(l);
                total_points += getNumPoints(l);
        }

        level_points = new int[total_points];
        weights = new double[total_points];

        int num_known_points = 0;
        double *known_x = new double[total_points];

        double *x = 0, *w = 0;

        for( int l=0; l<max_level; l++ ){
                int num_points = getNumPoints( l );
                buildOneLevel( l, w, x );
                for( int i=0; i<num_points; i++ ){
                        weights[ levels[l] + i ] = w[i];

                        int point = -1;
                        for( int j=0; j<num_known_points; j++ ){
                                if ( fabs( x[i] - known_x[j] ) < tol ){
                                        point = j;
                                        break;
                                }
                        }
                        if ( point == - 1){ // new point found
                                known_x[num_known_points] = x[i];
                                point = num_known_points;
                                num_known_points++;
                        }
                        level_points[levels[l] + i] = point;
                }
        }

        nodes = new double[num_known_points];
        tcopy( num_known_points, known_x, nodes );

        delete[] known_x;
        delete[] x;
        delete[] w;
}

RuleGaussJacobi::~RuleGaussJacobi(){
        if ( nodes != 0 ){ delete[] nodes; }
        if ( levels != 0 ){ delete[] levels; }
        if ( weights != 0 ){ delete[] weights; }
        if ( level_points != 0 ){ delete[] level_points; }
}

int RuleGaussJacobi::getMaxLevel() const{ return max_level; }

TypeOneDRule RuleGaussJacobi::getType() const{
        return rule_gaussgegenbauer;
}

int RuleGaussJacobi::getNumPoints( int level ) const{
        return level+1;
}

int RuleGaussJacobi::getBasisLevel( int level ) const{
        return 2*level;
}

void RuleGaussJacobi::getPoints( int level, int* &pnts ) const{
        if ( pnts != 0 ){ delete[] pnts; }
        pnts = new int[getNumPoints( level )];
        tcopy( getNumPoints( level ), &( level_points[levels[level]] ), pnts );
}

const char* RuleGaussJacobi::getDescription() const{
        return "Gauss-Jacobi points and weights and Lagrange Polynomials";
}

double RuleGaussJacobi::getX( int point ) const{ return nodes[point]; }

double RuleGaussJacobi::getWeight( int level, int point ) const{
        for( int i=levels[level]; i<levels[level+1]; i++ ){
                if ( level_points[i] == point ){ return weights[i]; }
        }
        return 0.0; // this should never happen
}

double RuleGaussJacobi::eval( int level, int point, double x ) const{
        double value = 1.0, d = nodes[point];
        for( int i=levels[level]; i<levels[level+1]; i++ ){
                value *= ( level_points[i] != point ) ? ( x - nodes[level_points[i]] ) / ( d - nodes[level_points[i]] ) : 1.0;
        }
        return value;
}

void RuleGaussJacobi::buildOneLevel( int level, double* &w, double* &x ){
        // get Gauss-Jacobi quadrature points
        int m = getNumPoints( level );
        if ( w != 0 ){ delete[] w; }
        if ( x != 0 ){ delete[] x; }
        w = new double[m];
        x = new double[m];

        double *s = new double[m];

        for( int i=0; i<m; i++ ){ x[i] = w[i] = s[i] = 0.0; }

        double ab = alpha + beta;

        w[0] = sqrt( pow( 2.0, 1.0 + ab ) * tgamma( alpha + 1.0 ) * tgamma( beta + 1.0 ) / tgamma( 2.0 + ab ) );

        x[0] = ( beta - alpha ) / ( 2.0 + ab );
        s[0] = sqrt( 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) / (( 3.0 + ab) * ( 2.0 + ab ) * ( 2.0 + ab )) );
        for( int i=1; i<m; i++ ){
                double di = (double) (i+1);
                x[i] = ( beta*beta - alpha*alpha ) / ( (2.0*di + ab -2.0)*(2.0*di + ab) );
                s[i] = sqrt( 4.0 * di * (di + alpha ) * (di + beta) * ( di + ab )/ ( ( (2.0*di + ab)*(2.0*di + ab) - 1.0 ) * ( 2.0*di + ab ) * ( 2.0*di + ab ) ) );
        }
        s[m-1] = 0.0;

        decompose( m, x, s, w );

        delete[] s;
}

}

#endif

