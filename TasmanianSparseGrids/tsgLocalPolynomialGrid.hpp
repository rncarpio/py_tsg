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

#ifndef __TASMANIAN_SPARSE_GRID_LOCAL_POLYNOMIAL_GRID_HPP
#define __TASMANIAN_SPARSE_GRID_LOCAL_POLYNOMIAL_GRID_HPP

#include "tsgEnumerate.hpp"
#include "tsgHelperFunctions.hpp"

#include "tsgBase1DRule.hpp"
#include "tsgRulePieceWiseLocal.hpp"
#include "tsgRulePieceWiseLocalZero.hpp"

#include "tsgBaseGrid.hpp"


namespace TasGrid{

class LocalPolynomialGrid : public Grid{
public:
        LocalPolynomialGrid();
        LocalPolynomialGrid( int dimensions, int outputs, int depth, int order, TypeOneDRule boundary );
        ~LocalPolynomialGrid();

        void reset( int dimensions, int outputs, int depth, int order, TypeOneDRule boundary );
        //virtual void reset( int outputs );
        void updateOrder( int new_order ); // resets the order of the underlying grid (recomputes surpluses)
        int getOrder() const;

        int getNumDimensions() const;
        int getNumOutputs() const;
        TypeOneDRule getOneDRule() const;
        const char *getOneDRuleDescription() const;

        void write( std::ofstream &ofs ) const; // write the grid to a file
        bool read( std::ifstream &ifs ); // reads the grid

        int getNumPoints() const;

        void getPoints( double* &pnts ) const;
        void getWeights( double* &weights ) const;
        void getInterpolantWeights( const double x[], double* &weights ) const;

        int getNumNeededPoints() const;
        void getNeededPoints( double* &pnts ) const;
        void loadNeededPoints( const double vals[] );
        void loadNeededPoints( const IndexSet *data );

        void evaluate( const double x[], double y[] ) const;
        void integrate( double y[] ) const;

        // refinement functions
        const IndexSet* getState() const; // returns enough information to rebuild the grid (i.e. the tensor list or the points list)
        void setState( const IndexSet* state ); // copy the tensors and/or points
        void getData( IndexSet* &data ); // returns a list of the set points and their values
        void getUpdateState( IndexSet* &update, double tol, TypeRefinement criteria ) const; // give the new set of points or tensors
        void setUpdate( const IndexSet *update ); // creates a grid with the data updates

protected:
        void clear();

        void addChild( const int point[], int direction, IndexSet *destination, IndexSet *exclude = 0 ) const;
        // adds children in direction to destination only if they are not part of exclude
        bool addParent( const int point[], int direction, IndexSet *destination, IndexSet *exclude = 0 ) const;
        // adds the parent if it has not been excluded and returns true if anything has been added

        void recomputeSurpluses();
        void applySurplusMapTransposed( double w[] ) const; // needed to compute the weights

        int makeLevelMap( int* &map ) const; // returns the max level
        double evalBasis( const int p[], const double x[] ) const;
        double evalIntegral( const int p[] ) const;

        // the map has dimensions num_points x num_dimensions, for each point and each direction, it flags wheather it should be refined or not
        void buildUpdateMap( int* &map, double tol, TypeRefinement criteria ) const; // use int for the map so we can flag more than true/false (-1 do not refine, 0 not set, 1 refine)

        bool sameLine( const int a[], const int b[], int direction ) const;

        // returns the L-\infty norm in each direction, if the norm is less than a tolerance, it is replaced by 1 indicating that we will work with absolute error in that direction
        void computeOutputNormalization( double* &norm ) const;

        //bool isChildless( const int point[], int direction ) const; //returnstrue if the children of this point are missing from points

private:
        OneDHierarchicalRule *rule1D;

        RulePieceWiseLocal pwp;
        RulePieceWiseLocalZero pwp0;

        TypeOneDRule rule; // keeps track whether we use zero boundary or regular rule

        int num_dimensions, num_outputs;

        double *surplus;

        IndexSet *points;
        IndexSet *needed_points;
};


}
#endif

