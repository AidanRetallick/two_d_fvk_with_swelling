//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC//
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC//
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC//
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC//
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC//
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//====================================================================
#include <fenv.h>
//Generic routines
#include "generic.h"

// The equations
#include "C1_foeppl_von_karman.h"

// The mesh
#include "meshes/triangle_mesh.h"
#include "C1_basis/my_geom_object.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

/*                      OUTLINE OF PROBLEM CONSTRUCTION                       */
// The basic constuction is much the same as the usual order of things in a
// problem. Underneath is the order of actions (with stars next to non actions
// that are unique to these types of problems).
// 1.  Setup mesh parameters
// 2.  Build the mesh
// 3.* Upgrade Elements
//     We upgrade edge elements on relevant boundaries to be curved C1 elements.
//     This involves working out which edge is to be upgraded and then passing
//     information about the global curve and start and end points of the
//     element edge on that curve to the element.
// 4.* Rotate edge degrees of freedom.
//     We rotate the Hermite dofs that lie on the edge into the normal -
//     tangential basis so that we can set physical boundary conditions like
//     clamping or resting conditions.
// 5.  Complete problem Setup and set Boundary conditions.

/*                            REQUIRED DEFINITIONS                            */
// Per Curve Section we will need:
// 1.  A parametric function defining the curve section.
// 2.  The tangential derivative of the parametric function defining
//     the curve section.
// 3.* (For 5 order boundary representation) The second tangential derivative
//     of the parametric function defining the curve section.
// 4.  A unit normal and tangent to each curve section and corresponding
//     derivatives, to allow the rotation of boundary coordinates.
// It also convenient to define:
// 1.  An inverse function (x,y) -> s (the arc coordinate) to help in setting
//     up the nodal positions in terms of this parametric coordinate.

namespace Params
{
 //Shape of the domain
 double A = 1.0;
 double B = 1.0;
 // The coupling of the stretching energy
 double p_mag = 0;
 double nu = 0.0;
 double mu = 0.0;
 double eta = 0.1;//12*(1-nu*nu);
 
 // Dimensions
 double thickness = 0.0054;
 double L_dim = 5.0e-3;
 double E_dim = 1.2139e6;
 double P_dim = E_dim * thickness / eta;

 // Mesh parameters
 double element_area=0.1;

 // Output directory
 string output_dir="RESLT";


 /*                     PARAMETRIC BOUNDARY DEFINITIONS                        */
 // Here we create the geom objects for the Parametric Boundary Definition
 CurvilineCircleTop parametric_curve_top;
 CurvilineCircleBottom parametric_curve_bottom;

 // The normal and tangential directions. We need the derivatives so we can form
 // The Hessian and the Jacobian of the rotation
 void get_normal_and_tangent(const Vector<double>& x, Vector<double>& n,
			     Vector<double>& t, DenseMatrix<double>& Dn, DenseMatrix<double>& Dt)
 {
  // Fill in the normal and derivatives of the normal
  n[0] = x[0]/sqrt(x[0]*x[0]+x[1]*x[1]);
  n[1] = x[1]/sqrt(x[0]*x[0]+x[1]*x[1]);

  // The (x,y) derivatives of the (x,y) components
  Dn(0,0) = x[1]*x[1] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
  Dn(1,0) =-x[1]*x[0] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
  Dn(0,1) =-x[0]*x[1] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
  Dn(1,1) = x[0]*x[0] * pow(x[0]*x[0]+x[1]*x[1],-1.5);

  // Fill in the tangent and derivatives of the tangent
  t[0] =-x[1]/sqrt(x[0]*x[0]+x[1]*x[1]);
  t[1] = x[0]/sqrt(x[0]*x[0]+x[1]*x[1]);

  Dt(0,0) =-Dn(1,0);
  Dt(1,0) = Dn(0,0);
  Dt(0,1) =-Dn(1,1);
  Dt(1,1) = Dn(0,1);
 }

 /*                           PROBLEM DEFINITIONS                              */
 // Assigns the value of pressure depending on the position (x,y)
 void get_pressure(const Vector<double>& x, double& pressure)
 {
  // Constant pressure
  pressure = p_mag;
 }

 // Pressure wrapper so we can output the pressure function
 void get_pressure(const Vector<double>& X, Vector<double>& pressure)
 {
  pressure.resize(1);
  get_pressure(X,pressure[0]);
 }

 // Temperature wrapper so we can output the temperature function
 void get_swelling(const Vector<double>& X, Vector<double>& swelling)
 {
  swelling.resize(1);
  swelling[0]=0;
 }

 // Assigns the value of in plane forcing depending on the position (x,y)
 void get_in_plane_force(const Vector<double>& x, Vector<double>& grad)
 {
  // No in plane force
  grad[0]=0.0;
  grad[1]=0.0;
 }

 // This metric will flag up any non--axisymmetric parts
 void axiasymmetry_metric(const Vector<double>& x, const
			  Vector<double>& u, const Vector<double>& u_exact, double& error, double& norm)
 {
  // We use the theta derivative of the out of plane deflection
  error = pow((-x[1]*u[1] + x[0]*u[2])/sqrt(x[0]*x[0]+x[1]*x[1]),2);
  norm  = pow(( x[0]*u[1] + x[1]*u[2])/sqrt(x[0]*x[0]+x[1]*x[1]),2);
 }

 // Get the exact solution
 void get_null_fct(const Vector<double>& X, double& exact_w)
 {
  exact_w = 0.0;
 }

 // Get the exact solution(s)
 void dummy_exact_w(const Vector<double>& x, Vector<double>& exact_w)
 {
  /* Do nothing -> no exact solution in this case */
 }

}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredFvKProblem : public virtual Problem
{

public:

 /// Constructor
 UnstructuredFvKProblem();

 /// Destructor
 ~UnstructuredFvKProblem()
 {
  Trace_file_dim.close();
  Trace_file_nondim.close();
  delete (Surface_mesh_pt);
  delete (Bulk_mesh_pt);
  // Clean up memory
  delete Outer_boundary_pt;
  delete Outer_boundary_ellipse_pt;
  delete Outer_curvilinear_boundary_pt[0];
  delete Outer_curvilinear_boundary_pt[1];
  delete Inner_open_boundaries_pt[0];
  delete Inner_open_boundaries_pt[1];
  delete Boundary2_pt;
  delete Boundary3_pt;
 };

 /// Setup and build the mesh
 void build_mesh();

 /// Return the number of newton iterations in the last solve attempt
 unsigned nnewton_iter_taken()
 {
  return Nnewton_iter_taken;
 }
 
 /// Update after solve (empty)
 void actions_after_newton_solve()
 {
  /* No actions before newton solve */
 }

 /// Information about the previous newton step
 void actions_after_newton_step()
 {
  DoubleVector dx;
  get_residuals(dx);
  double maxres = dx.max();
  
  ofstream newton_stream;
  char filename[100];
  string tempstring =
   Params::output_dir + "/newton_info" + std::to_string(Doc_info.number()) + ".dat";
  strcpy(filename, tempstring.c_str());
  newton_stream.open(filename);
  newton_stream << Doc_info.number() << " "
		<< maxres            << std::endl;
  newton_stream.close();
 }
 
 /// Pin the in-plane displacements and set to zero at centre
 void pin_in_plane_displacements_at_centre_node();

 /// Update the problem specs before solve: Re-apply boundary conditions
 /// Empty as the boundary conditions stay fixed
 void actions_before_newton_solve()
 {
  /* Reapply boundary conditions */
  apply_boundary_conditions();
 }

 /// Doc the solution
 void doc_solution(const std::string& comment="");

 /// \short Overloaded version of the problem's access function to
 /// the mesh. Recasts the pointer to the base Mesh object to
 /// the actual mesh type.
 TriangleMesh<ELEMENT>* mesh_pt()
 {
  return Bulk_mesh_pt;
 }

 /// Doc info object for labeling output
 DocInfo Doc_info;

private:

 // Triangle Mesh Parameter Data
 // This is the data used to set-up the mesh, we need to store the pointers
 // HERE otherwise we will not be able to clean up the memory once we have
 // finished the problem.
 Ellipse* Outer_boundary_ellipse_pt;
 // The outer curves
 Vector<TriangleMeshCurveSection*> Outer_curvilinear_boundary_pt;
 // The Internal curves
 Vector<TriangleMeshOpenCurve *> Inner_open_boundaries_pt;
 // The close outer boundary
 TriangleMeshClosedCurve* Outer_boundary_pt;
 // The first of the internal boundaries
 TriangleMeshPolyLine* Boundary2_pt;
 // The second of the internal boundaries
 TriangleMeshPolyLine* Boundary3_pt;

 /// Actions to be performed after read in of meshes
 void actions_after_read_unstructured_meshes()
 {
  // Curved Edges need to be upgraded after the rebuild
  upgrade_edge_elements_to_curve(0,Bulk_mesh_pt);
  upgrade_edge_elements_to_curve(1,Bulk_mesh_pt);
  // Rotate degrees of freedom
  rotate_edge_degrees_of_freedom(Bulk_mesh_pt);
  // Make the problem fully functional
  complete_problem_setup();
  // Apply any boundary conditions
  apply_boundary_conditions();
 }

 /// Helper function to apply boundary conditions
 void apply_boundary_conditions();

 /// \short Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup();

 /// Temporal error norm function.
 double global_temporal_error_norm()
 {
  double global_error = 0.0;
  
  //Find out how many nodes there are in the problem
  unsigned n_node = mesh_pt()->nnode();
 
  //Loop over the nodes and calculate the estimated error in the values
  for(unsigned i=0;i<n_node;i++)
   {
    double error = 0;
    Node* node_pt = mesh_pt()->node_pt(i);
    // Get error in solution: Difference between predicted and actual
    // value for nodal value 2 (only if we are at a vertex node)
    if(node_pt->nvalue()>2)
     {
      error = node_pt->time_stepper_pt()->
       temporal_error_in_value(mesh_pt()->node_pt(i),2);
     }
    //Add the square of the individual error to the global error
    global_error += error*error;
   }
    
  // Divide by the number of nodes
  global_error /= double(n_node);
 
  // Return square root...
  return sqrt(global_error);
 }
 
  /// Trace file to document norm of solution
  ofstream Trace_file_dim, Trace_file_nondim;

  // Keep track of boundary ids
  enum
  {
   Outer_boundary0 = 0,
   Outer_boundary1 = 1,
   Inner_boundary0 = 2,
   Inner_boundary1 = 3
  };

 double Element_area;

 /// \short Loop over all curved edges, then loop over elements and upgrade
 /// them to be curved elements
 void upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const &
				     bulk_mesh_pt);

 /// \short Loop over all edge elements and rotate the Hermite degrees of freedom
 /// to be in the directions of the two in-plane vectors specified in Params
 void rotate_edge_degrees_of_freedom(Mesh* const &bulk_mesh_pt);

 /// \short Delete traction elements and wipe the surface mesh
 void delete_traction_elements(Mesh* const &surface_mesh_pt);

 /// Pointer to "bulk" mesh
 TriangleMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to "surface" mesh
 Mesh* Surface_mesh_pt;

}; // end_of_problem_class

/// Constructor definition
template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem()
 :
 Element_area(Params::element_area)
{
 add_time_stepper_pt(new BDF<2>(true));
 Problem::Always_take_one_newton_step = true;
 
 Problem::Always_take_one_newton_step = true;
 // Build the mesh
 build_mesh();

 // Curved Edge upgrade
 upgrade_edge_elements_to_curve(0,Bulk_mesh_pt);
 upgrade_edge_elements_to_curve(1,Bulk_mesh_pt);

 // Rotate degrees of freedom
 rotate_edge_degrees_of_freedom(Bulk_mesh_pt);

 // Store number of bulk elements
 complete_problem_setup();

 char filename[100];
 ofstream Param_file;
 strcpy(filename, (Params::output_dir + "/parameters.dat").c_str());
 Param_file.open(filename);

 // Output plate parameters
 Param_file << "thickness    " << Params::thickness << std::endl
	    << "nu           " << Params::nu        << std::endl
	    << "eta          " << Params::eta       << std::endl
	    << std::endl
	    << "L_dim        " << Params::L_dim     << std::endl
	    << "E_dim        " << Params::E_dim     << std::endl
	    << "P_dim        " << Params::P_dim     << std::endl
	    << std::endl
	    << "Element area " << Params::element_area << std::endl;
 
 strcpy(filename, (Params::output_dir + "/trace_nondim.dat").c_str());
 Trace_file_nondim.open(filename);
 strcpy(filename, (Params::output_dir + "/trace_dim.dat").c_str());
 Trace_file_dim.open(filename);

 oomph_info << "Number of equations: "
	    << assign_eqn_numbers() << '\n';
} // end Constructor

/// Set up and build the mesh
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::build_mesh()
{
 Vector<double> zeta(1);
 Vector<double> posn(2);
 
 //Outer boundary
 //--------------
 
 double A = 1.0;
 double B = 1.0;
 Outer_boundary_ellipse_pt = new Ellipse(A, B);
 
 //First bit
 double zeta_start = 0.0;
 double zeta_end = MathematicalConstants::Pi;
 unsigned nsegment = (int)(MathematicalConstants::Pi/sqrt(Element_area));
 
 Outer_curvilinear_boundary_pt.resize(2);
 Outer_curvilinear_boundary_pt[0] =
  new TriangleMeshCurviLine(Outer_boundary_ellipse_pt, zeta_start,
			    zeta_end, nsegment, Outer_boundary0);
 
 //Second bit
 zeta_start = MathematicalConstants::Pi;
 zeta_end = 2.0*MathematicalConstants::Pi;
 nsegment = (int)(MathematicalConstants::Pi/sqrt(Element_area));
 Outer_curvilinear_boundary_pt[1] =
  new TriangleMeshCurviLine(Outer_boundary_ellipse_pt, zeta_start,
			    zeta_end, nsegment, Outer_boundary1);
 
 Outer_boundary_pt =
  new TriangleMeshClosedCurve(Outer_curvilinear_boundary_pt);
 
 // Internal open boundaries
 // Total number of open curves in the domain
 unsigned n_open_curves = 2;
 // We want internal open curves
 Inner_open_boundaries_pt.resize(n_open_curves);
 
 // Internal bit - this means we can have a boundary which is just the centre
 // We start by creating the internal boundaries
 // The boundary 2 is defined by its two vertices
 // Open curve 1
 Vector<Vector<double> > vertices(2,Vector<double>(2,0.0));
 vertices[0][0] =-0.5;
 vertices[0][1] = 0.0;

 vertices[1][0] = 0.5;
 vertices[1][1] = 0.0;
 unsigned boundary_id = Inner_boundary0;

 Boundary2_pt =
  new TriangleMeshPolyLine(vertices, boundary_id);
 // Open Curve 2
 vertices[0][0] = 0.0;
 vertices[0][1] =-0.5;

 vertices[1][0] = 0.0;
 vertices[1][1] = 0.5;
 boundary_id = Inner_boundary1;

 Boundary3_pt =
  new TriangleMeshPolyLine(vertices, boundary_id);

 // Each internal open curve is defined by a vector of
 // TriangleMeshCurveSection,
 // on this example we only need one curve section for each internal boundary
 Vector<TriangleMeshCurveSection *> internal_curve_section1_pt(1);
 internal_curve_section1_pt[0] = Boundary2_pt;

 Vector<TriangleMeshCurveSection *> internal_curve_section2_pt(1);
 internal_curve_section2_pt[0] = Boundary3_pt;

 // The open curve that define this boundary is composed of just one
 // curve section
 Inner_open_boundaries_pt[0] =
  new TriangleMeshOpenCurve(internal_curve_section1_pt);

 Inner_open_boundaries_pt[1] =
  new TriangleMeshOpenCurve(internal_curve_section2_pt);

 //Create the mesh
 //---------------
 //Create mesh parameters object
 TriangleMeshParameters mesh_parameters(Outer_boundary_pt);
 
 mesh_parameters.element_area() = Element_area;
 
 // Specify the internal open boundaries
 mesh_parameters.internal_open_curves_pt() = Inner_open_boundaries_pt;
 
 // Build an assign bulk mesh
 Bulk_mesh_pt=new TriangleMesh<ELEMENT>(mesh_parameters, time_stepper_pt());
 
 // Create "surface mesh" that will contain only the prescribed-traction
 // elements. The constructor creates the mesh without adding any nodes
 // elements etc.
 Surface_mesh_pt =  new Mesh;
 
 //Add two submeshes to problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);
 
 // Combine submeshes into a single Mesh
 build_global_mesh();
}// end build_mesh

//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of
/// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::pin_in_plane_displacements_at_centre_node()
{
 // Pin the node that is at the centre in the domain!
 // Get the num of nods on internal_boundary 2
 unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(2);
 for (unsigned inod=0;inod<num_int_nod;inod++)
  {
   // Get node point
   Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(2,inod);
   // If the node is on the other internal boundary too
   if( nod_pt->is_on_boundary(3))
    {
     // Pin it! It's the centre of the domain!
     // In-plane dofs are always 0 and 1 - on vertex nodes we also have
     // out-of-plane dofs from 2-7.
     nod_pt->pin(0);
     nod_pt->set_value(0,0.0);
     nod_pt->pin(1);
     nod_pt->set_value(1,0.0);
    }
  }

}

//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of
/// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{
 // Set the boundary conditions for problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.

 // The node at the very centre should be pinned when we have "do-nothing"
 // conditions on the in-plane displacements (i.e stress free b/c)
 pin_in_plane_displacements_at_centre_node();

 // Complete the build of all elements so they are fully functional
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   
   // Resize the internal data to make room for past timesteps
   //   unsigned n_bdof = el_pt->ninternal_dofs();
   unsigned index_of_internal_data = el_pt->index_of_internal_data();
   // for(unsigned i_bdof=0; i_bdof<n_bdof; i_bdof++)
   //  {
   oomph_info << "About to resize for el " << e << std::endl; //" and dof " << i_bdof << std::endl; 
     el_pt->internal_data_pt(index_of_internal_data)
      ->set_time_stepper(time_stepper_pt(), false);
   //  }
    
   //Set the pressure function pointers and the physical constants
   el_pt->pressure_fct_pt() = &Params::get_pressure;
   el_pt->in_plane_forcing_fct_pt() = &Params::get_in_plane_force;
   // There is no error metric in this case
   el_pt->error_metric_fct_pt() = &Params::axiasymmetry_metric;
   el_pt->nu_pt() = &Params::nu;
   el_pt->mu_pt() = &Params::mu;
   el_pt->eta_pt() = &Params::eta;
  }
 // Set the boundary conditions
 apply_boundary_conditions();
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{
 // Set the boundary conditions
 //Just loop over outer boundary since inner boundary doesn't have boundary
 //conditions
 unsigned nbound = Outer_boundary1 + 1;
 for(unsigned b=0;b<nbound;b++)
  {
   const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
   for(unsigned e=0;e<nb_element;e++)
    {
     // Get pointer to bulk element adjacent to b
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
					     Bulk_mesh_pt->boundary_element_pt(b,e));
     // A true clamp, so we set everything except the second normal to zero
     for(unsigned idof=0; idof<6; ++idof)
      {
       // Cannot set second normal derivative
       if(idof!=3)
	{
	 el_pt->fix_out_of_plane_displacement_dof(idof,b,Params::get_null_fct);
	}
      }
    }
  }
} // end set bc

/// A function that upgrades straight sided elements to be curved. This involves
// Setting up the parametric boundary, F(s) and the first derivative F'(s)
// We also need to set the edge number of the upgraded element and the positions
// of the nodes j and k (defined below) and set which edge (k) is to be exterior
/*            @ k                                                             */
/*           /(                                                               */
/*          /. \                                                              */
/*         /._._)                                                             */
/*      i @     @ j                                                           */
// For RESTING or FREE boundaries we need to have a C2 CONTINUOUS boundary
// representation. That is we need to have a continuous 2nd derivative defined
// too. This is well discussed in by [Zenisek 1981] (Aplikace matematiky ,
// Vol. 26 (1981), No. 2, 121--141). This results in the necessity for F''(s)
// as well.
//==start_of_upgrade_edge_elements============================================
template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT >::
upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const &bulk_mesh_pt) 
{
 // How many bulk elements adjacent to boundary b
 unsigned n_element = bulk_mesh_pt-> nboundary_element(b);
 // These depend on the boundary we are on
 CurvilineGeomObject* parametric_curve_pt; 

 // Define the functions for each part of the boundary
 switch (b)
  {
  // Upper boundary
  case 0:
   parametric_curve_pt = &Params::parametric_curve_top;
  break;

  // Lower boundary
  case 1:
   parametric_curve_pt = &Params::parametric_curve_bottom;
  break;

  default:
   throw OomphLibError(
    "I have encountered a boundary number that I wasn't expecting. Please fill \
me in if you want additional curved boundaries..",
    "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
    OOMPH_EXCEPTION_LOCATION);
  break;
 }
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
   
   // Loop over (vertex) nodes
   unsigned nnode=3; //This should always be = 3 for triangles
   unsigned index_of_interior_node=3;

   // Enum for the curved edge
   MyC1CurvedElements::Edge edge(MyC1CurvedElements::none);

   // Vertices positions
   Vector<Vector<double> > xn(3,Vector<double>(2,0.0));
 
   // Get vertices for debugging
   Vector<Vector<double> > verts(3,Vector<double>(2,0.0));
   // Loop nodes
   for(unsigned n=0;n<nnode;++n)
    {
     // If it is on boundary
     Node* nod_pt = bulk_el_pt->node_pt(n);
     verts[n][0]=nod_pt->x(0);
     verts[n][1]=nod_pt->x(1);

     // Check if it is on the outer boundaries
     if(nod_pt->is_on_boundary(Outer_boundary0) 
         || nod_pt->is_on_boundary(Outer_boundary1))
      {
       xn[n][0]=nod_pt->x(0);
       xn[n][1]=nod_pt->x(1);
      }
     // The edge is denoted by the index of the  opposite (interior) node
     else {index_of_interior_node = n;}
    }
   // Initialise s_ubar s_obar (start and end respectively)
   double s_ubar, s_obar;

   // s at the next (cyclic) node after interior
   s_ubar = parametric_curve_pt->get_zeta(xn[(index_of_interior_node+1) % 3]);
   // s at the previous (cyclic) node before interior
   s_obar = parametric_curve_pt->get_zeta(xn[(index_of_interior_node+2) % 3]);

   // Assign edge case
   switch(index_of_interior_node)
    {
     case 0: edge= MyC1CurvedElements::zero; 
      break;
     case 1: edge= MyC1CurvedElements::one; 
      break;
     case 2: edge= MyC1CurvedElements::two; 
      break;
     // Should break it here HERE
     default: edge= MyC1CurvedElements::none; 
      throw OomphLibError(
       "The edge number has been set to a value greater than two: either we have\
 quadrilateral elements or more likely the index_of_interior_node was never set\
 and remains at its default value.",
       "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
       OOMPH_EXCEPTION_LOCATION);
      break;
     }
   // Check for inverted elements HERE
   if (s_ubar>s_obar)
    {
     oomph_info <<"Apparent clockwise direction of parametric coordinate."
                <<"This will probably result in an inverted element."
                <<"s_start="<<s_ubar<<"; s_end ="<<s_obar<<std::endl;
     throw OomphLibError(
       "The Edge coordinate appears to be decreasing from s_start to s_end. \
Either the parametric boundary is defined to be clockwise (a no-no) or \
the mesh has returned an inverted element (less likely)",
       "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
       OOMPH_EXCEPTION_LOCATION);
    }

   // Upgrade it
   bulk_el_pt->upgrade_element_to_curved(edge,s_ubar,s_obar,
    parametric_curve_pt,3);     
  }
}// end upgrade elements
// Function to set up rotated nodes on the boundary: necessary if we want to set
// up physical boundary conditions on a curved boundary with Hermite type dofs.
// For example if we know w(n,t) = f(t) (where n and t are the
// normal and tangent to a boundary) we ALSO know dw/dt and d2w/dt2.
// NB no rotation is needed if the edges are completely free!
// begin rotate_edge_degrees_of_freedom
template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
rotate_edge_degrees_of_freedom( Mesh* const &bulk_mesh_pt)
{
 // How many bulk elements
 unsigned n_element = bulk_mesh_pt-> nelement();
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
 
   // Loop nodes
   unsigned nnode =3;
   unsigned nbnode=0;
   // Count the number of boundary nodes
   for (unsigned n=0; n<nnode;++n)
     {
      // Check it isn't on an internal boundary
      bool on_boundary_2=el_pt->node_pt(n)->is_on_boundary(2);
      bool on_boundary_3=el_pt->node_pt(n)->is_on_boundary(3);
      // If it isn't on an internal boundary but it is on an external boundary
      if(!(on_boundary_2 || on_boundary_3))
       {nbnode+=unsigned(el_pt->node_pt(n)->is_on_boundary());}
     }

   // Now if we have nodes on boundary in this element 
   if(nbnode>0)
    {
     // Set up vector
     Vector<unsigned> bnode (nbnode,0);
     unsigned inode(0);

     // Fill in the bnode Vector
     for (unsigned n=0; n<nnode;++n)
      {
       // Check it isn't on an internal boundary
       bool on_boundary_2=el_pt->node_pt(n)->is_on_boundary(2);
       bool on_boundary_3=el_pt->node_pt(n)->is_on_boundary(3);
       if(!(on_boundary_2 || on_boundary_3))
       {
       // If it is on the boundary
       if(el_pt->node_pt(n)->is_on_boundary())
        {
         // Set up the Vector holding the boundary nodes
         bnode[inode]=n;
         ++inode;
        }
       }
      }
    // Now rotate the nodes by passing the index of the nodes and the 
    // normal / tangent vectors to the element
    el_pt->set_up_rotated_dofs(nbnode,bnode,&Params::get_normal_and_tangent);
   }
 }
}// end rotate_edge_degrees_of_freedom

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(const
						   std::string& comment)
{
 ofstream some_file;
 char filename[100];
 string tempstring;

 // Number of plot points
 unsigned npts = 2;

 tempstring = Params::output_dir + "/coarse_soln" + std::to_string(Doc_info.number()) + ".dat";
 strcpy(filename, tempstring.c_str());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
	   << comment << "\"\n";
 some_file.close();

 // Number of plot points
 npts = 10;

 tempstring = Params::output_dir + "/soln" + std::to_string(Doc_info.number()) + ".dat";
 strcpy(filename, tempstring.c_str());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
	   << comment << "\"\n";
 some_file.close();


 // Write the pressure, degree of swelling,
 // deflection and time to the trace file
 //-----------------------------------------------------
 // Get the centre deflection first
 Vector<double> w_centre(1,0.0);

 // Get the num of nods on internal_boundary 2
 unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(2);
 for (unsigned inod=0;inod<num_int_nod;inod++)
  {
   // Get node point
   Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(2,inod);
   // If the node is on the other internal boundary too
   if( nod_pt->is_on_boundary(3))
    {
     // Then the node is in the centre!
     w_centre[0]=nod_pt->value(2);
     

 Trace_file_nondim << Doc_info.number()              << " "
		   << Params::p_mag                  << " "
		   << w_centre[0]                    << " "
		   << time()                         << " "
		   << nnewton_iter_taken()           << endl;
 Trace_file_dim    << Doc_info.number()              << " "
		   << Params::p_mag*Params::P_dim    << " "
		   << w_centre[0]*Params::L_dim      << " "
		   << time()                         << " "
		   << nnewton_iter_taken()           << endl;
    }
  }
 // Increment the doc_info number
 Doc_info.number()++;


 //  // Store a profile deflection slice through x=0
 // //-----------------------------------------------------
 // unsigned n_ppoint=101;
 // Vector<double> ppoint(2,0.0), ppoint_s(2,0.0), z(8,0.0);
 // GeomObject* ppoint_element_pt;
 // double h = 2.0/(n_ppoint-1.0);
  
 // tempstring = Params::output_dir + "/profile" + std::to_string(Doc_info.number()) + ".dat";
 // strcpy(filename, tempstring.c_str());
 // some_file.open(filename);

 // ppoint[0]=-1.0;
 // ppoint[1]=-1.0;
 // for(unsigned j=0; j<n_ppoint; j++)
 //  {
 //   ppoint[0]+=h;
 //   ppoint[1]+=h;
 //   for(unsigned i=0; i<n_element; i++)
 //    {
 //     try
 //      {
 //       dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i))
 // 	->locate_zeta(ppoint, ppoint_element_pt, ppoint_s);
 //      }
 //     catch(...)
 //      {
 //       oomph_info << "HEY, I CAUGHT AN ERROR:" << std::endl;
 //       oomph_info << "NEVERMIND, CARRY ON." << std::endl;
 //      }
 //     if(ppoint_element_pt!=NULL)
 //      {
 //       z = dynamic_cast<ELEMENT*>(ppoint_element_pt)
 // 	->interpolated_u_foeppl_von_karman(ppoint_s);
 //       some_file << ppoint[0]+z[6] << " " << ppoint[1]+z[7] << " " << z[0] << std::endl;
 //       break;
 //      }
 //    }  
 //  }
 // some_file.close();


} // end of doc

//============start_of_delete_flux_elements==============================
/// Delete Poisson Flux Elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>
::delete_traction_elements(Mesh* const &surface_mesh_pt)
{
 // How many surface elements are in the surface mesh
 unsigned n_element = surface_mesh_pt->nelement();

 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete surface_mesh_pt->element_pt(e);
  }

 // Wipe the mesh
 surface_mesh_pt->flush_element_and_node_storage();

} // end of delete_flux_elements

//=======start_of_main========================================
///Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{
 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 // Directory for solution
 CommandLineArgs::specify_command_line_flag("--dir", &Params::output_dir);

 // Poisson Ratio
 CommandLineArgs::specify_command_line_flag("--nu", &Params::nu);

 // Dampening coefficient
 CommandLineArgs::specify_command_line_flag("--mu", &Params::mu);

 double p_inc = 10.0;
 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--p", &p_inc);

 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--eta", &Params::eta);

 // Pin u_alpha everywhere
 CommandLineArgs::specify_command_line_flag("--pininplane");

 // Element Area (no larger element than 0.09)
 Params::element_area=0.1;
 CommandLineArgs::specify_command_line_flag("--element_area", &Params::element_area);

 // Parse command line
 CommandLineArgs::parse_and_assign();

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
 UnstructuredFvKProblem<ThermoFoepplVonKarmanC1CurvedBellElement<4>>
  problem;
 
 if(CommandLineArgs::command_line_flag_has_been_set("--pininplane"))
  {
   cout << "gonna pin em" << endl;
   // Pin the in-plane displacements
   unsigned nnode = problem.mesh_pt()->nnode();
   for(unsigned inode=0; inode<nnode; inode++)
    {
     problem.mesh_pt()->node_pt(inode)->pin(0);
     problem.mesh_pt()->node_pt(inode)->pin(1);
    }
   cout << "successfully pinned" << endl;
   problem.assign_eqn_numbers();
  }

 // Set ICs of each type to be zero everywhere for all time.
 // First each timestep
 unsigned ntime = problem.time_stepper_pt()->ntstorage();
 for(unsigned itime=0; itime<ntime; itime++)
  {
   // Next at each node
   unsigned nnode = problem.mesh_pt()->nnode();
   for(unsigned inode=0; inode<nnode; inode++)
    {
     // And finally for every dof
     Node* node_pt = problem.mesh_pt()->node_pt(inode);
     unsigned ndof = node_pt->nvalue();
     for(unsigned idof=0; idof<ndof; idof++)
      {
       node_pt->set_value(itime, idof, 0.0);
      }
    }
  }
 
 // // Check ICs
 // ofstream iccheck;
 // iccheck.open("RESLT/IC_file.dat");
 // for(unsigned in=0; in<problem.mesh_pt()->nnode(); in++)
 //  {
 //   iccheck << "Node: " << in << std::endl;
 //   Node* nod_pt = problem.mesh_pt()->node_pt(in);
 //   unsigned nval = nod_pt->nvalue();
 //   unsigned ntime = nod_pt->ntstorage();
 //   for
 //   iccheck << " " << std::endl;
 //  }
  
 // Set up some problem paramters
 problem.newton_solver_tolerance()=1e-10;
 problem.max_residuals()=1e4;
 problem.max_newton_iterations()=30;

 double dt = 1.0;
 double t_max = 100*dt;
 double epsilon = 1.0e-4;
 problem.initialise_dt(dt);

 //double c_inc = 0.001;

 cout << "about to steady solve" << endl;
 // INITIAL SOLVE
 problem.steady_newton_solve(); // SOLVE
 problem.doc_solution(); // AND DOCUMENT

 // INFLATION
 Params::p_mag+=p_inc;   // 10
 
 // unsigned nstep = unsigned(t_max/dt);
 while(problem.time() < t_max)
  {
   oomph_info << "STARTING A NEW TIMESTEP :)" << std::endl;
   oomph_info << endl << Params::eta << " " << Params::mu << endl;
   double dt_next=problem.adaptive_unsteady_newton_solve(dt,epsilon);
   oomph_info << dt_next << ", " << dt << std::endl;
   dt=dt_next;
   problem.doc_solution();
  }

} //End of main
