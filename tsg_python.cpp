//
// Copyright (c) 2014 Ronaldo Carpio
//                                     

/*
	This file is part of py_tsg.

    py_tsg is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/                                                  

#include <TasmanianSparseGrids/TasmanianSparseGrid.hpp>
#include <Python.h>
#include <pyublas/numpy.hpp>
#include <tuple>
#include <string>
#include <vector>

using namespace TasGrid;
namespace bpl = boost::python;
using std::vector;

#define bpl_assert(exp, s)	\
	if (!(exp)) {				\
      PyErr_SetString(PyExc_ValueError, (s));	\
      boost::python::throw_error_already_set();	\
	}	

typedef pyublas::numpy_vector<double> dPyArr;
typedef pyublas::numpy_vector<int> iPyArr;
typedef vector<double> dVec;
typedef vector<int> iVec;
typedef vector<dPyArr> dPyArrVector;

// convert any type with begin() and end() -> a Python list
template <class IterableT>
struct forward_iterable_to_list { 
  static PyObject* convert(IterableT const &cont) {
    bpl::list result;
    for (auto iter=cont.begin(); iter != cont.end(); iter++) {
      result.append(*iter);
    }                      
    return boost::python::incref(result.ptr());
  }                                            
}; 

dPyArr double_carray_to_dPyArr(int n, const double* p) {
  dPyArr result(n);
  std::copy(p, p+n, result.begin());
  return result;
}

class TSG_Wrap : public TasmanianSparseGrid {
public:
  void make_global_grid(int dimensions, int outputs, int depth, TypeDepth type, TypeOneDRule oned, iPyArr const &a_weights, double alpha, double beta) {
    double alpha_beta[2] = {alpha, beta};
    iVec a_weights2;	
	const int* a_weights3 = NULL;
	if (a_weights.size() > 0) {
	  bpl_assert(a_weights.size() == dimensions, "anisotropic weights must match dimensions");
	  a_weights2.resize(dimensions);
	  std::copy(a_weights.begin(), a_weights.end(), a_weights2.begin());
	  a_weights3 = &a_weights2[0];
	}
    this->makeGlobalGrid(dimensions, outputs, depth, type, oned, a_weights3, alpha_beta);
  }
  
  void make_full_tensor_grid(int dimensions, int outputs, iPyArr order, TypeOneDRule oned, double alpha, double beta) {
    bpl_assert(dimensions == order.size(), "dimensions should equal len(order)");
    double alpha_beta[2] = {alpha, beta};
	vector<int> order2(order.begin(), order.end());
    this->makeFullTensorGrid(dimensions, outputs, &order2[0], oned, alpha_beta);
  }
  
  void recycle_global_grid(int depth, TypeDepth type) {
    this->recycleGlobalGrid(depth, type);
  }

  void recycle_full_tensor_grid(iPyArr order) {
    vector<int> order2(order.begin(), order.end());
    this->recycleFullTensorGrid(&order2[0]);
  }
  
/*  
  std::string write_string() const {
    std::ostringstream result;
	this->write(result);
	return result;
  }
  
  bool read_string(std::string const &s) {
    std::istringstream source(s);
	return this->read(source);
  }
*/  
  void write_file(std::string const &filename) const {
    this->write(filename.c_str());
  }
  
  bool read_file(std::string const &filename) {
    return this->read(filename.c_str());
  }
  
  void set_transform_AB(dPyArr const &a_array, dPyArr const &b_array) {
    int dims = this->getNumDimensions();
	bpl_assert(dims == a_array.size() && dims == b_array.size(), "dimensions must match array length");
	dVec a(a_array.begin(), a_array.end());
	dVec b(b_array.begin(), b_array.end());
	this->setTransformAB(&a[0], &b[0]);
  }
  
  vector<dPyArr> get_transform_AB() const {
    vector<dPyArr> result(2);
	int dims = this->getNumDimensions();
    double *a=NULL, *b=NULL;
	this->getTransformAB(a, b);
	result[0] = double_carray_to_dPyArr(dims, a);
	result[1] = double_carray_to_dPyArr(dims, b);
	delete a;
	delete b;
	return result;
  }
  
  dPyArr get_points() const {
    int dims = this->getNumDimensions();
	int n_points = this->getNumPoints();
	int array_dims[] = {n_points, dims};	
	double* output = NULL;
	this->getPoints(output);
	dPyArr result = double_carray_to_dPyArr(dims*n_points, output);
	delete output;
	result.reshape(2, array_dims);
	return result;
  }
  
  dPyArr get_weights() const {
    int n_points = this->getNumPoints();
	double* output = NULL;
	this->getWeights(output);
	dPyArr result = double_carray_to_dPyArr(n_points, output);
	delete output;
	return result;
  }
  
  dPyArr get_interpolant_weights(dPyArr const &x) const {
    int n_points = this->getNumPoints();
    vector<double> x2(x.begin(), x.end());
	double* output = NULL;
	this->getInterpolantWeights(&x2[0], output);
	dPyArr result = double_carray_to_dPyArr(n_points, output);
	delete output;
	return result;
  }
  
  dPyArr get_needed_points() const {
    int dims = this->getNumDimensions();
	int n_points = this->getNumNeededPoints();
	int array_dims[] = {n_points, dims};	
	double* output = NULL;
	this->getNeededPoints(output);
	dPyArr result;
	if (output != NULL) {
	  result = double_carray_to_dPyArr(dims*n_points, output);
	  delete output;
	  result.reshape(2, array_dims);
	}
	return result;
  }
  
  void load_needed_points(dPyArr const &vals) {
    int outputs = this->getNumOutputs();
	int n_points = this->getNumNeededPoints();
    bpl_assert(vals.size() == outputs*n_points, "vals has wrong size");
    vector<double> vals2(vals.begin(), vals.end());
    this->loadNeededPoints(&vals2[0]);
  }

  dPyArr evaluate_wrap(dPyArr const &x) const {
    vector<double> x2(x.begin(), x.end());
    dPyArr result(this->getNumOutputs());
    this->evaluate(&x2[0], &result[0]);	
	return result;
  }
  
  dPyArr integrate_wrap() const {
    dPyArr result(this->getNumOutputs());
    this->integrate(&result[0]);	
	return result;
  }
    
};

BOOST_PYTHON_MODULE(_py_tsg)
{
  bpl::to_python_converter<dPyArrVector, forward_iterable_to_list<dPyArrVector>>();
  
  // functions that don't need a wrapper
// getVersion()
// getLicense()
// makeLocalPolynomialGrid
// makeWaveletGrid
// recycleLocalPolynomialGrid
// recycleWaveletGrid
// clearTransformAB
// getNumDimensions
// getNumOutputs
// getOneDRule
// getOneDRuleDescription
// getNumPoints
// getNumNeededPoints
// printStats
// setRefinement

  bpl::class_<TSG_Wrap, boost::noncopyable>("TSG", bpl::init<>())
		.def("get_version", &TSG_Wrap::getVersion)
		.def("get_license", &TSG_Wrap::getLicense)
		.def("make_global_grid", &TSG_Wrap::make_global_grid)
		.def("make_local_polynomial_grid", &TSG_Wrap::makeLocalPolynomialGrid)
		.def("make_wavelet_grid", &TSG_Wrap::makeWaveletGrid)
		.def("make_full_tensor_grid", &TSG_Wrap::make_full_tensor_grid)
		.def("recycle_global_grid", &TSG_Wrap::recycle_global_grid)
		.def("recycle_local_polynomial_grid", &TSG_Wrap::recycleLocalPolynomialGrid)
		.def("recycle_wavelet_grid", &TSG_Wrap::recycleWaveletGrid)
		.def("recycle_full_tensor_grid", &TSG_Wrap::recycle_full_tensor_grid)
		//.def("write_string", &TSG_Wrap::write_string)
		//.def("read_string", &TSG_Wrap::read_string)
		.def("write_file", &TSG_Wrap::write_file)
		.def("read_file", &TSG_Wrap::read_file)
		.def("set_transform_AB", &TSG_Wrap::set_transform_AB)
		.def("clear_transform_AB", &TSG_Wrap::clearTransformAB)		
		.def("get_transform_AB", &TSG_Wrap::get_transform_AB)		
		.def("get_num_dimensions", &TSG_Wrap::getNumDimensions)
		.def("get_num_outputs", &TSG_Wrap::getNumOutputs)		
		.def("get_oned_outputs", &TSG_Wrap::getOneDRule)				
		.def("get_oned_rule_description", &TSG_Wrap::getOneDRuleDescription)				
		.def("get_num_points", &TSG_Wrap::getNumPoints)				
		.def("get_points", &TSG_Wrap::get_points)				
		.def("get_weights", &TSG_Wrap::get_weights)				
		.def("get_interpolant_weights", &TSG_Wrap::get_interpolant_weights)				
		.def("get_num_needed_points", &TSG_Wrap::getNumNeededPoints)				
		.def("get_needed_points", &TSG_Wrap::get_needed_points)				
		.def("load_needed_points", &TSG_Wrap::load_needed_points)				
		.def("evaluate", &TSG_Wrap::evaluate_wrap)				
		.def("integrate", &TSG_Wrap::integrate_wrap)				
		.def("print_stats", &TSG_Wrap::printStats)				
		.def("set_refinement", &TSG_Wrap::setRefinement)	
  ;
  
  bpl::enum_<TypeOneDRule>("TypeOneDRule")
	.value("rule_base", rule_base)
	.value("rule_clenshawcurtis", rule_clenshawcurtis)
	.value("rule_chebyshev", rule_chebyshev)
	.value("rule_gausslegendre", rule_gausslegendre)
	.value("rule_gausschebyshev1", rule_gausschebyshev1)
	.value("rule_gausschebyshev2", rule_gausschebyshev2)
	.value("rule_chebyshevN2P", rule_chebyshevN2P)
	.value("rule_fejer2", rule_fejer2)
	.value("rule_gaussgegenbauer", rule_gaussgegenbauer)
	.value("rule_gaussjacobi", rule_gaussjacobi)
	.value("rule_gausslaguerre", rule_gausslaguerre)
	.value("rule_gausshermite", rule_gausshermite)
	.value("rule_pwpolynomial", rule_pwpolynomial)
	.value("rule_pwpolynomial0", rule_pwpolynomial0)
	.value("rule_wavelet", rule_wavelet)
	.value("rule_fulltensor", rule_fulltensor)       
  ;  

  bpl::enum_<TypeDepth>("TypeDepth")  		
	.value("type_level", type_level) 
	.value("type_basis", type_basis)
	.value("type_hyperbolic", type_hyperbolic)
  ;

  bpl::enum_<TypeRefinement>("TypeRefinement")
    .value("refine_classic", refine_classic)
	.value("refine_parents_first", refine_parents_first)
	.value("refine_direction_selective", refine_direction_selective)
	.value("refine_fds", refine_fds)
  ;

}


