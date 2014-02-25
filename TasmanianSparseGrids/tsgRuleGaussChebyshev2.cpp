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

#ifndef __TASMANIAN_SPARSE_GRID_RULE_GC2_CPP
#define __TASMANIAN_SPARSE_GRID_RULE_GC2_CPP

#include "tsgRuleGaussChebyshev2.hpp"

namespace TasGrid{

RuleGaussChebyshev2::RuleGaussChebyshev2( const int level ) : max_level(level), tol(NUM_TOL), OneDRule(){
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
                buildOneLevel( l, w, x);
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

RuleGaussChebyshev2::~RuleGaussChebyshev2(){
        if ( nodes != 0 ){ delete[] nodes; }
        if ( levels != 0 ){ delete[] levels; }
        if ( weights != 0 ){ delete[] weights; }
        if ( level_points != 0 ){ delete[] level_points; }
}

int RuleGaussChebyshev2::getMaxLevel() const{ return max_level; }

TypeOneDRule RuleGaussChebyshev2::getType() const{
        return rule_gausschebyshev2;
}

int RuleGaussChebyshev2::getNumPoints( int level ) const{
        return level+1;
}

int RuleGaussChebyshev2::getBasisLevel( int level ) const{
        return 2*level;
}

void RuleGaussChebyshev2::getPoints( int level, int* &pnts ) const{
        if ( pnts != 0 ){ delete[] pnts; }
        pnts = new int[getNumPoints( level )];
        tcopy( getNumPoints( level ), &( level_points[levels[level]] ), pnts );
}

const char* RuleGaussChebyshev2::getDescription() const{
        return "Gauss-Chebychev points and weights of type 2, and Lagrange Polynomials";
}

double RuleGaussChebyshev2::getX( int point ) const{ return nodes[point]; }

double RuleGaussChebyshev2::getWeight( int level, int point ) const{
        for( int i=levels[level]; i<levels[level+1]; i++ ){
                if ( level_points[i] == point ){ return weights[i]; }
        }
        return 0.0; // this should never happen
}

double RuleGaussChebyshev2::eval( int level, int point, double x ) const{
        double value = 1.0, d = nodes[point];
        for( int i=levels[level]; i<levels[level+1]; i++ ){
                value *= ( level_points[i] != point ) ? ( x - nodes[level_points[i]] ) / ( d - nodes[level_points[i]] ) : 1.0;
        }
        return value;
}

void RuleGaussChebyshev2::buildOneLevel( int level, double* &w, double* &x ){
        // get Gauss-Chebyshev-type2 quadrature points
        int m = getNumPoints( level );
        if ( w != 0 ){ delete[] w; }
        if ( x != 0 ){ delete[] x; }
        w = new double[m];
        x = new double[m];

        for( int i=0; i<m; i++ ){
                double theta = M_PI*( (double)(i+1) )/( (double)(m+1) );
                x[m-i-1] = cos( theta );
                w[i] = (M_PI / ( (double)(m+1) )) * sin( theta )* sin( theta );
        }
}

}

#endif

