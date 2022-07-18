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
#include "my_exact_solutions.h"

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

namespace TestSoln
{
//Shape of the domain
double A = 1.0;
double B = 1.0;
// The coupling of the stretching energy
double p_mag = 1; 
double nu = 0.5;
double eta = 12*(1-nu*nu);
double eta_david =1.0;

/*                         VALIDATION DEFINITIONS                             */
// Validation cases as enum
enum Validation_case {no_validate=-1, one=1, two=2};
// Classes that contain the two validation solutions - these are in the 'DGR 
// nondimensionalisation' 
Validation_case validation_case = no_validate;
// Classes that contain the two validation solutions
ManufacturedSolutionWithLinearDisplacements 
  manufactured_solution_linear_u(&nu,&eta_david);
ManufacturedSolutionDecoupledExtension 
  manufactured_solution_decoupled_extension(&nu,&p_mag);

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
 //pressure = p_mag; //constant pressure - can make a function of x
 // NB if you are going to use 1/r make sure to take out internal boundaries
 // Or there will be a div by zero error here
 
 // Dependent on case call the forcing
 if(validation_case== one)
  manufactured_solution_linear_u.get_pressure(x,pressure);
 else if(validation_case== two)
  manufactured_solution_decoupled_extension.get_pressure(x,pressure);
 // Default just use a constant force
 else
  pressure = p_mag;

 //Convert from `David Nondim' to `Airy Nondim'
 pressure *= 12*(1-nu*nu);
}

// Pressure wrapper so we can output the pressure function
void get_pressure(const Vector<double>& X, Vector<double>& pressure)
 {
  pressure.resize(1);
  get_pressure(X,pressure[0]);
 }

// Assigns the value of in plane forcing depending on the position (x,y)
void get_in_plane_force(const Vector<double>& x, Vector<double>& grad)
 {
 // Dependent on case call the forcing
 if(validation_case== one)
  manufactured_solution_linear_u.get_in_plane_force(x,grad);
 else if(validation_case== two)
  manufactured_solution_decoupled_extension.get_in_plane_force(x,grad);
 else
  {
  grad[0]=0.0;
  grad[1]=0.0;
  }
 }

// This metric will flag up any non--axisymmetric parts
void error_metric(const Vector<double>& x, const 
  Vector<double>& u, const Vector<double>& u_exact, double& error, double& norm)
{
 // We use the theta derivative of the out of plane deflection
 error = pow((-x[1]*u[1] + x[0]*u[2])/sqrt(x[0]*x[0]+x[1]*x[1]),2);
 norm  = pow(( x[0]*u[1] + x[1]*u[2])/sqrt(x[0]*x[0]+x[1]*x[1]),2);
}

// Get the exact solution, complete with radial and azimuthal derivatives
void get_exact_radial_w(const Vector<double>& x, Vector<double>& exact_w)
{
 // u_z and derivatives
 if(validation_case == one)
  manufactured_solution_linear_u.exact_radial_w_solution(x,exact_w);

 else if(validation_case == two)
  manufactured_solution_decoupled_extension.exact_radial_w_solution(x,exact_w);

 else
  { /* Do Nothing - no exact solution for this case*/}
}

// Get the exact solution(s)
void get_exact_w(const Vector<double>& x, Vector<double>& exact_w)
{
 // u_z and derivatives
 if(validation_case == one)
  manufactured_solution_linear_u.exact_solution(x,exact_w);

 else if(validation_case == two)
  manufactured_solution_decoupled_extension.exact_solution(x,exact_w);
 else 
 {
  /* Do nothing -> no exact solution in this case */
 }
}

// Get the exact u solution (wrapper)
void get_exact_ux(const Vector<double>& X, double& exact_w)
{
 Vector<double> w(8);
 get_exact_w(X,w);
 exact_w= w[6];
}

// Get the exact u solution (wrapper)
void get_exact_uy(const Vector<double>& X, double& exact_w)
{
 Vector<double> w(8);
 get_exact_w(X,w);
 exact_w= w[7];
}

// Get the exact solution
void get_null_fct(const Vector<double>& X, double& exact_w)
{
 exact_w = 0.0;
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
UnstructuredFvKProblem(double element_area = 0.09);

/// Destructor
~UnstructuredFvKProblem()
{
 Trace_file.close();
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

/// Update after solve (empty)
void actions_after_newton_solve()
{
/* No actions before newton solve */
}

/// Pin the in-plane displacements and set to zero at centre
void pin_in_plane_displacements_at_centre_node();

/// Set the initial values to the exact solution (useful for debugging)
void set_initial_values_to_exact_solution()
{
// Loop over all elements and reset the values
unsigned n_element = Bulk_mesh_pt->nelement();
for(unsigned e=0;e<n_element;e++)
{
 // Upcast from GeneralisedElement to the present element
 ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
 unsigned nnode = el_pt->nnode();
 for(unsigned i=0;i<3;++i)
  {
   // Containers
   Vector<double> x(2),w(8,0.0);
   // Get node pointer
   Node* nod_pt = el_pt->node_pt(i);
   // Get x 
   x[0]=nod_pt->x(0);
   x[1]=nod_pt->x(1);
   // Get test
   TestSoln::get_exact_w(x,w);
   // Set value
  for(unsigned l=0 ; l<6;++l)
   {
    nod_pt->set_value(2+l,w[l]);
   }
  }
 // Pin in Plane dofs
 for(unsigned i=0;i<nnode;++i)
  {
   // Containers
   Vector<double> x(2),w(8,0.0);
   // Get node pointer
   Node* nod_pt = el_pt->node_pt(i);
   // Get x 
   x[0]=nod_pt->x(0);
   x[1]=nod_pt->x(1);
   // Get test
   TestSoln::get_exact_w(x,w);
   // Set value
  for(unsigned l=0 ; l<6;++l)
   {
    nod_pt->set_value(0,w[6]);
    nod_pt->set_value(1,w[7]);
   }
  }
  // Pin the internal dofs
  const unsigned n_internal_dofs = el_pt->number_of_internal_dofs();
  for(unsigned i=0;i<n_internal_dofs;++i)
   {
   // Containers
   Vector<double> x(2,0.0), s(2,0.0), w(6,0.0);
   // Get internal data
   Data* internal_data_pt = el_pt->internal_data_pt(1);
   el_pt->get_internal_dofs_location(i,s);
   // Get test
   el_pt->get_coordinate_x(s,x); 
   TestSoln::get_exact_radial_w(x,w);
   // HERE need a function that can give us this lookup
   internal_data_pt->set_value(i,w[0]);
   }  

//Just loop over the boundary elements
unsigned nbound = Outer_boundary1 + 1;
for(unsigned b=0;b<nbound;b++)
 {
 const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<nb_element;e++)
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));
   // Loop over vertices of the element (i={0,1,2} always!)
   for(unsigned i=0;i<3;++i)
    {
     Node* nod_pt = el_pt->node_pt(i);
     // If it is on the bth boundary
     if(nod_pt -> is_on_boundary(b))
      {
       // Set the values of the unknowns that aren't pinned
       Vector<double> x(2,0.0),w(6,0.0);
       x[0] = nod_pt->x(0);
       x[1] = nod_pt->x(1);
       TestSoln::get_exact_radial_w(x,w);
       // Just set the nonzero values -> the others are pinned
       nod_pt->set_value(2+0,w[0]);
       nod_pt->set_value(2+1,w[1]);
       nod_pt->set_value(2+2,0.0);
       nod_pt->set_value(2+3,w[3]);
       nod_pt->set_value(2+4,0.0);
       nod_pt->set_value(2+5,0.0);
      }
    }
  }
 }
}// End loop over elements
}

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
return dynamic_cast<TriangleMesh<ELEMENT>*> (Problem::mesh_pt()); 
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

/// Trace file to document norm of solution
ofstream Trace_file;

// Keep track of boundary ids
enum
{
 Outer_boundary0 = 0,
 Outer_boundary1 = 1,
 Inner_boundary0 = 2,
 Inner_boundary1 = 3
};

double Element_area;

// The extra members required for flux type boundary conditions.
/// \short Number of "bulk" elements (We're attaching flux elements
/// to bulk mesh --> only the first Nkirchhoff_elements elements in
/// the mesh are bulk elements!)
// unsigned Nkirchhoff_elements;

/// \short Create bending moment elements on the b-th boundary of the
/// problems mesh 
void create_traction_elements(const unsigned &b, Mesh* const & bulk_mesh_py,
                            Mesh* const &surface_mesh_pt);

/// \short Loop over all curved edges, then loop over elements and upgrade
/// them to be curved elements
void upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const & 
 bulk_mesh_pt);

/// \short Loop over all edge elements and rotate the Hermite degrees of freedom
/// to be in the directions of the two in-plane vectors specified in TestSoln
void rotate_edge_degrees_of_freedom(Mesh* const &bulk_mesh_pt);

/// \short Delete traction elements and wipe the surface mesh
void delete_traction_elements(Mesh* const &surface_mesh_pt);

/// \short Set pointer to prescribed-flux function for all elements
/// in the surface mesh
void set_prescribed_traction_pt();

/// Pointer to "bulk" mesh
TriangleMesh<ELEMENT>* Bulk_mesh_pt;

/// Pointer to "surface" mesh
Mesh* Surface_mesh_pt;

}; // end_of_problem_class

/// Constructor
template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
:
Element_area(element_area)
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
unsigned nsegment = (int)(MathematicalConstants::Pi/sqrt(element_area));

Outer_curvilinear_boundary_pt.resize(2);
Outer_curvilinear_boundary_pt[0] =
new TriangleMeshCurviLine(Outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary0);

//Second bit
zeta_start = MathematicalConstants::Pi;
zeta_end = 2.0*MathematicalConstants::Pi;
nsegment = (int)(MathematicalConstants::Pi/sqrt(element_area));
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

mesh_parameters.element_area() = element_area;

// Specify the internal open boundaries
mesh_parameters.internal_open_curves_pt() = Inner_open_boundaries_pt;

// Build an assign bulk mesh
Bulk_mesh_pt=new TriangleMesh<ELEMENT>(mesh_parameters);

// Create "surface mesh" that will contain only the prescribed-traction
// elements. The constructor creates the mesh without adding any nodes
// elements etc.
Surface_mesh_pt =  new Mesh;

//Add two submeshes to problem
add_sub_mesh(Bulk_mesh_pt);
add_sub_mesh(Surface_mesh_pt);

// Combine submeshes into a single Mesh
build_global_mesh();

// Curved Edge upgrade
upgrade_edge_elements_to_curve(0,Bulk_mesh_pt);
upgrade_edge_elements_to_curve(1,Bulk_mesh_pt);
 
// Rotate degrees of freedom
rotate_edge_degrees_of_freedom(Bulk_mesh_pt);

// Store number of bulk elements
complete_problem_setup();

char filename[100];
sprintf(filename, "RESLT/trace.dat");
Trace_file.open(filename);

oomph_info << "Number of equations: "
        << assign_eqn_numbers() << '\n';

}

//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of 
/// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::pin_in_plane_displacements_at_centre_node()
{
// Pin the node that is at the centre in the domain!
// Get the num of nods on internal_boundary 0
unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(2);
for (unsigned inod=0;inod<num_int_nod;inod++)
{
 // Get node point
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(2,inod);
 // If the node is on the other internal boundary too
 if( nod_pt->is_on_boundary(3))
 {
  // Be verbose and let everyone know we have found the centre
  oomph_info<<"Found centre point\n";
  oomph_info<<"x = ("<<nod_pt->x(0)<<" "<<nod_pt->x(1)<<") \n";
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
if(TestSoln::validation_case != TestSoln::one)
 {
  pin_in_plane_displacements_at_centre_node();
 }

// Complete the build of all elements so they are fully functional
unsigned n_element = Bulk_mesh_pt->nelement();
for(unsigned e=0;e<n_element;e++)
{
// Upcast from GeneralisedElement to the present element
ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

//Set the pressure function pointers and the physical constants
el_pt->pressure_fct_pt() = &TestSoln::get_pressure;
el_pt->in_plane_forcing_fct_pt() = &TestSoln::get_in_plane_force;  
// There is no error metric in this case
if(TestSoln::validation_case == TestSoln::no_validate) 
 el_pt->error_metric_fct_pt() = &TestSoln::error_metric;
el_pt->nu_pt() = &TestSoln::nu;
el_pt->eta_pt() = &TestSoln::eta;
}

// Set the boundary conditions
apply_boundary_conditions();
/// Set pointer to prescribed traction function for traction elements
// /* No Traction */

// Re-apply Dirichlet boundary conditions (projection ignores
// boundary conditions!)
// /* No Projection */
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
   el_pt->fix_out_of_plane_displacement_dof(0,b,TestSoln::get_null_fct);
   el_pt->fix_out_of_plane_displacement_dof(1,b,TestSoln::get_null_fct);
   el_pt->fix_out_of_plane_displacement_dof(2,b,TestSoln::get_null_fct);
   el_pt->fix_out_of_plane_displacement_dof(4,b,TestSoln::get_null_fct);
   el_pt->fix_out_of_plane_displacement_dof(5,b,TestSoln::get_null_fct);
  if(TestSoln::validation_case == TestSoln::one)
   {
   el_pt->fix_in_plane_displacement_dof(0,b,TestSoln::get_exact_ux);
   el_pt->fix_in_plane_displacement_dof(1,b,TestSoln::get_exact_uy);
   }
  }
 }
} // end set bc

template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
create_traction_elements(const unsigned &b, Mesh* const &bulk_mesh_pt, 
                             Mesh* const &surface_mesh_pt)
{
/* NO TRACTION ELEMENTS YET DEFINED*/
}// end create traction elements


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
   parametric_curve_pt = &TestSoln::parametric_curve_top;
  break;

  // Lower boundary
  case 1:
   parametric_curve_pt = &TestSoln::parametric_curve_bottom;
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
    el_pt->set_up_rotated_dofs(nbnode,bnode,&TestSoln::get_normal_and_tangent);
   }
 }
}// end create traction elements

//==start_of_set_prescribed_traction_pt===================================
/// Set pointer to prescribed traction function for all elements in the 
/// surface mesh
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::set_prescribed_traction_pt()
{
}// end of set prescribed flux pt

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(const 
                                                    std::string& comment)
{ 
ofstream some_file;
char filename[100];

// Number of plot points
unsigned npts = 2;

sprintf(filename,"RESLT/coarse_soln%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output(some_file,npts); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

// Number of plot points 
//(Can output at extremely high resolution BUT to save diskspace limit to 5)
//(Try it with 50 on very low resolution to see HOW accurate we can be!)
npts = 5;

sprintf(filename,"RESLT/soln%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output(some_file,npts); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

//  Output exact solution
sprintf(filename,"%s/exact_soln%i-%f.dat","RESLT",Doc_info.number()
 ,Element_area);
some_file.open(filename);
Bulk_mesh_pt->output_fct(some_file,npts,TestSoln::get_exact_w); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

//  Output pressure function
sprintf(filename,"%s/pressure%i-%f.dat","RESLT",Doc_info.number()
 ,Element_area);
some_file.open(filename);
Bulk_mesh_pt->output_fct(some_file,npts,TestSoln::get_pressure); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

// Output boundaries
//------------------
sprintf(filename,"RESLT/boundaries%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output_boundaries(some_file);
some_file.close();

// Output regions
unsigned n_region = Bulk_mesh_pt->nregion();
if (n_region > 1)
{
for (unsigned r = 0; r < n_region; r++)
{
 //Attempt to output elements in different regions
 sprintf(filename,"RESLT/region%i%i-%f.dat",r,Doc_info.number(),
   Element_area);
 some_file.open(filename);
 unsigned nel = Bulk_mesh_pt->nregion_element(r);
 for (unsigned e = 0; e < nel; e++)
  {
   Bulk_mesh_pt->region_element_pt(r,e)->output(some_file,npts);
  }
 some_file.close();
}
}

// // Doc error and return of the square of the L2 error
// //---------------------------------------------------
// //double error,norm,dummy_error,zero_norm;
  double dummy_error,zero_norm;
 sprintf(filename,"RESLT/error%i-%f.dat",Doc_info.number(),Element_area);
 some_file.open(filename);
 
 Bulk_mesh_pt->compute_error(some_file,TestSoln::get_exact_w,
                        dummy_error,zero_norm);
 some_file.close();
 
 // Doc L2 error and norm of solution
 oomph_info << "Absolute norm of computed solution: " << sqrt(dummy_error) 
            << std::endl;
 
 oomph_info << "Norm of computed solution: " << sqrt(zero_norm)
            << std::endl;
 
 //Trace_file << TestSoln::p_mag << " " << "\n ";
 //
 // Find the solution at r=0
 //   // ----------------------
 MeshAsGeomObject* Mesh_as_geom_obj_pt=
  new MeshAsGeomObject(Bulk_mesh_pt);
 Vector<double> s(2);
 GeomObject* geom_obj_pt=0;
 Vector<double> r(2,0.0);
 Mesh_as_geom_obj_pt->locate_zeta(r,geom_obj_pt,s);
 oomph_info << "pointer: " << geom_obj_pt << std::endl;
 // The member function does not exist in this element
 // it is instead called interpolated_u_foeppl_von_karman and returns a vector of 
 // length 12 or 18 - the interface is pretty horrible so it may be something
 // we want to tidy up
 Vector<double> u_0(12,0.0);
 u_0=dynamic_cast<ELEMENT*>(geom_obj_pt)->interpolated_u_foeppl_von_karman(s);

 oomph_info << "w in the middle: " <<std::setprecision(15) << u_0[0] << std::endl;

 Trace_file << TestSoln::p_mag
            << " " << u_0[0] << '\n';

// Doc error and return of the square of the L2 error
//---------------------------------------------------
sprintf(filename,"RESLT/L2-norm%i-%f.dat",
        Doc_info.number(),
        Element_area);
some_file.open(filename);

some_file<<"### L2 Norm\n";
some_file<<"##  Format: err^2 norm^2 \n";
// Print error in prescribed format
some_file<< dummy_error <<" "<< zero_norm <<"\n";
some_file.close();

// Increment the doc_info number
Doc_info.number()++;

// Clean up
delete Mesh_as_geom_obj_pt;
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

 // The `restart' flag
 CommandLineArgs::specify_command_line_flag("--restart");

 // The `validate' flag
 CommandLineArgs::specify_command_line_flag("--validate");


 // The `problem dump' flag
 CommandLineArgs::specify_command_line_flag("--dump_at_every_step");

 // Directory for solution
 string output_dir="RESLT";
 CommandLineArgs::specify_command_line_flag("--dir", &output_dir);

 // Poisson Ratio
 CommandLineArgs::specify_command_line_flag("--nu", &TestSoln::nu);
 
 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--p", &TestSoln::p_mag);

 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--eta", &TestSoln::eta);

 // P_step
 double p_step=3;
 CommandLineArgs::specify_command_line_flag("--dp", &p_step);

 // Applied Pressure
 double p_max = 4;
 CommandLineArgs::specify_command_line_flag("--p_max", &p_max);

 // Element Area (no larger element than 0.09)
 double element_area=0.09;
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);

 // Which case are we using
 int iusage_case = -1;
 CommandLineArgs::specify_command_line_flag("--validation_case", &iusage_case);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // If restart flag has been set
 const bool do_restart=CommandLineArgs::
   command_line_flag_has_been_set("--restart");

 // If restart flag has been set
 const bool dump_at_every_step=CommandLineArgs::
   command_line_flag_has_been_set("--dump_at_every_step");

 // If validate flag provided
 const bool validate=CommandLineArgs::
   command_line_flag_has_been_set("--validate");

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Assign the boundary case
 if (iusage_case==-1 ||  (iusage_case>=1 && iusage_case <=2))
  {
  // Cast int to enum
  TestSoln::validation_case=(TestSoln::Validation_case)(iusage_case);
  }
 else // Default to what is set in TestSoln
  { 
   oomph_info<<"Boundary case \""<<iusage_case<<"\" not recognised.\n";
   return(-1);
  }

 // Validation loop
 if(validate)
  {
  //Validate 1
  {
  // Parameters for first test
  TestSoln::validation_case = TestSoln::one; 
  TestSoln::p_mag = 1;
  // Problem instance
  UnstructuredFvKProblem<FoepplVonKarmanC1CurvedBellElement<4> >
    problem(element_area);
  // Set max residuals
  problem.max_residuals()=1e3;
  // Set max residuals
  // Loop until target pressure
  // Newton solve
  problem.disable_info_in_newton_solve();
  problem.newton_solve();
  problem.doc_solution();
  }

  //Validate 2
  {
  // Parameters for first test
  TestSoln::validation_case = TestSoln::two; 
  TestSoln::p_mag = 0;
  p_max = 3;
  p_step = 1;
  // Problem instance
  UnstructuredFvKProblem<FoepplVonKarmanC1CurvedBellElement<4> >
    problem(element_area);
  // Set max residuals
  problem.max_residuals()=1e2;
  // Calculate the number of steps
  unsigned number_p_steps = ceil((p_max - TestSoln::p_mag)/p_step);
  // Loop until target pressure
  for(unsigned ip = 0; ip<number_p_steps;++ip)
   {
    // Newton solve
    TestSoln::p_mag += p_step; 
    problem.disable_info_in_newton_solve();
    problem.newton_solve();
    problem.doc_solution();
   }
   }
   // EXIT
   return 0;
  }

 // Problem instance
 UnstructuredFvKProblem<FoepplVonKarmanC1CurvedBellElement<4> >
   problem(element_area);

 // Set up some problem paramters
 problem.max_residuals()=1e3;
 problem.max_newton_iterations()=20;

 // If we are restarting
 if(do_restart)
 {
  oomph_info<<"Restarting\n";
  // Read the restart file
  std::ifstream filestream;
  filestream.open("fvk_problem_data.dump");
  problem.read(filestream);    
  filestream.close(); 

  // Check with a newton solve 
  problem.newton_solve(); 
  
  //Output solution
  problem.doc_solution();
 }

 // Now do simple loop until target pressure and solve
 while(TestSoln::p_mag<p_max)
  {
  // Do the newton solve
  oomph_info<<"Solving for p=" << TestSoln::p_mag<<"\n";
  problem.newton_solve();

  // Document
  problem.doc_solution();
  oomph_info << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << " Pcos (" << TestSoln::p_mag << ")" << std::endl;
  oomph_info << "Current dp  (" << p_step << ")" << std::endl;
  oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
  oomph_info << "Solution number (" <<problem.Doc_info.number()-1 << ")"
             << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << std::endl;

  // Dump the data 
  if(dump_at_every_step)
  {
   char filename[100];
   std::ofstream filestream;
   filestream.precision(15);
   sprintf(filename,"%s/fvk_circle_data%i-%f.dump",
           problem.Doc_info.directory().c_str(),
           problem.Doc_info.number(),
           TestSoln::p_mag
          );
   oomph_info<<"\nDumping to: \""<<filename<<"\"\n";
   filestream.open(filename);
   problem.dump(filestream);
   filestream.close();
  }
  // Increase the pressure 
  TestSoln::p_mag+=p_step;
  }
 // Print success
 oomph_info<<"Exiting Normally\n";

} //End of main

