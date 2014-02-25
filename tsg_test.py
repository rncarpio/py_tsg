#!/usr/bin/env python

import scipy, scipy.integrate, itertools

import pyublas
import _py_tsg as tsg

# examples copied from example.cpp in TasmianianSparseGrid distribution
def test1():
	grid = tsg.TSG()
	
	# EXAMPLE 1:
	# make a classical Smolyak grid using Clenshaw-Curtis quadrature
	# integrate the function f(x,y) = exp( -x^2 ) * cos( y )
	# the exact answer is: 2.513723354063905e+00
	def fn1(x):		return scipy.exp(-x[0]*x[0])*scipy.cos(x[1])
	dimension = 2
	outputs = 0
	level = 7
	grid.make_global_grid( dimension, outputs, level, tsg.TypeDepth.type_level, tsg.TypeOneDRule.rule_clenshawcurtis, scipy.array([], dtype=int), 0, 0 )
	points = grid.get_points()
	weights = grid.get_weights()
	sum = scipy.sum( [w*fn1(x) for (x,w) in zip(points, weights)] )
	print("\nExample 1")
	print("grid has: %d points" % grid.get_num_points())
	print("integral is: %.17f" % sum)
	print("error: %.17f" % scipy.fabs( sum - 2.513723354063905e+00 ))

	# EXAMPLE 2:
	# make a Clenshaw-Curtis rule that interpolates exactly polynomials of order up to 10
	# integrate the function f(x,y) = exp( -x^2 ) * cos( y )
	# exact values is: 6.990131267703512e-01	
	dimension = 2
	outputs = 1
	precision = 10
	grid.make_global_grid( dimension, outputs, precision, tsg.TypeDepth.type_basis, tsg.TypeOneDRule.rule_clenshawcurtis, scipy.array([], dtype=int), 0, 0 )
	points = grid.get_needed_points()
	values = fn1(points.T)
	grid.load_needed_points(values)
	x = scipy.array([0.3, 0.7])
	f = grid.evaluate(x)
	print("\nExample 2")
	print("grid has: %d points" % grid.get_num_points())
	print("value at %s is %.17f" % (x, f))
	print("error: %.17f" % scipy.fabs( f - 6.990131267703512e-01 ))
	sum = grid.integrate()
	print("integral is: %.17f" % sum)
	print("error: %.17f" % scipy.fabs( sum - 2.513723354063905e+00 ))

	# EXAMPLE 3:
	# make a quadrature that used Gauss-Legendre points and integrates exactly polynomials up to order 10
	# integrate the function f(x,y) = exp( -x^2 ) * cos( y )	
	dimension = 2
	outputs = 0
	precision = 10
	grid.make_global_grid( dimension, outputs, precision, tsg.TypeDepth.type_basis, tsg.TypeOneDRule.rule_gausslegendre, scipy.array([], dtype=int), 0, 0 )
	points = grid.get_points()
	weights = grid.get_weights()
	sum = scipy.sum( [w*fn1(x) for (x,w) in zip(points, weights)] )
	print("\nExample 3")
	print("grid has: %d points" % grid.get_num_points())
	print("integral is: %.17f" % sum)
	print("error: %.17f" % scipy.fabs( sum - 2.513723354063905e+00 ))
	
	# EXAMPLE 4:
	# make a quadrature that uses Gauss-Gegenbauer points and uses 8 times more points in y direction
	# integrate the function f(x,y) = ( x - 2 )^3 * exp( -y^2 ) * ( 1 - x^2 )^0.4 * ( 1 - y^2 )^0.4
	# exact integral is: -2.029979511486524e+01	
	def fn2(x):	return (x[0]-2.0)**3 * scipy.exp(-x[1]**2)
	dimension = 2
	outputs = 0
	depth = 16
	anisotropic_weights = scipy.array([8, 1])
	alpha = 0.4
	grid.make_global_grid( dimension, outputs, depth, tsg.TypeDepth.type_level, tsg.TypeOneDRule.rule_gaussgegenbauer, anisotropic_weights, alpha, 0 )
	points = grid.get_points()
	weights = grid.get_weights()
	sum = scipy.sum( [w*fn2(x) for (x,w) in zip(points, weights)] )
	print("\nExample 4")
	print("grid has: %d points" % grid.get_num_points())
	print("integral is: %.17f" % sum)
	print("error: %.17f" % scipy.fabs( sum - -2.029979511486524e+01 ))
	
	# EXAMPLE 5:
	# interpolates the function f(x,y) = exp( -x^2 ) * cos( y ) 
	# using adaptive piece-wise local quadratic polynomials over [0,1] x [0,1]
	# integrate the interpolant
	# the exact integral is: 6.990131267703512e-01	
	grid = tsg.TSG()
	dimension = 2
	outputs = 1
	initial_level = 4
	order = 2
	grid.make_local_polynomial_grid( dimension, outputs, initial_level, order, tsg.TypeOneDRule.rule_pwpolynomial )
	a = scipy.array([0.0, 0.0])
	b = scipy.array([1.0, 1.0])
	grid.set_transform_AB(a, b)

	tolerance = 1e-6	
	iteration = 0
	print("\nExample 5")
	while (grid.get_num_needed_points() > 0):
		N = grid.get_num_needed_points()
		points = grid.get_needed_points()
		values = fn1(points.T)
		#BREAK()
		# print("loading %d values" % len(values))
		# for i in range(len(values)):
			# print("%.5f, %.5f -> %.10f" % (points[i,0], points[i,1], values[i]))
		grid.load_needed_points(values)
		x = scipy.array([0.3, 0.7])
		f = grid.evaluate(x)		
		print("iteration: %d number of samples: %d value: %.17f error: %.17f" % (iteration, grid.get_num_points(), f, scipy.fabs( f - 6.990131267703512e-01 )))
		grid.set_refinement(tolerance, tsg.TypeRefinement.refine_classic)
		iteration += 1

if __name__ == "__main__":
	test1()

	
	
	
	
	
	
	
	
	
	
	
