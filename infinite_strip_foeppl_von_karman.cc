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

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;


namespace Params
{
 // Plate parameters
 double L1 = 1.0;
 double L2 = 1.0;
 double thickness = 0.0054;
 double nu = 0.495;
 double eta = 12*(1-nu*nu) / (thickness*thickness);
 // Control parameters
 double p_mag = 0.0;
 double c_swell = 0.0;

 // Dimensions
 double L_dim = 1.0;
 double E_dim = 1.21e6;
 double P_dim = E_dim * thickness / eta;

 // Output directory
 string output_dir="RESLT";
 unsigned n_edge_nodes=0;

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

 // Assigns the value of swelling depending on the position (x,y)
 void get_swelling(const Vector<double>& x, double& swelling)
 {
  // Almost uniform swelling with clamped boundaries
  swelling = c_swell;//*(1 - pow(x[0]/L1,10) - pow(x[1]/L2,10)
		     // + pow( (x[0]*x[1])/(L1*L2) , 10));
 }

 // Temperature wrapper so we can output the temperature function
 void get_swelling(const Vector<double>& X, Vector<double>& swelling)
 {
  swelling.resize(1);
  get_swelling(X,swelling[0]);
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
 UnstructuredFvKProblem(double element_area = 0.09);

 /// Destructor
 ~UnstructuredFvKProblem()
 {
  Trace_file_nondim.close();
  Trace_file_dim.close();
  delete (Surface_mesh_pt);
  delete (Bulk_mesh_pt);
  // Clean up memory
  delete Boundary_pt;
  //delete Inner_open_boundaries_pt[0];
  //delete Inner_open_boundaries_pt[1];
  delete Boundary0_pt;
  delete Boundary1_pt;
  delete Boundary2_pt;
  delete Boundary3_pt;
  //delete InnerBoudary0_pt;
  //delete InnerBoudary1_pt;
 };

 /// Setup and build the mesh
 void build_mesh();

 /// Update after solve (empty)
 void actions_after_newton_solve()
 {
  /* No actions before newton solve */
 }

 /// Pin the in-plane displacements and set to zero at centre
 void pin_in_plane_displacements_at_centre_node();

 /// Update the problem specs before solve: Re-apply boundary conditions
 /// Empty as the boundary conditions stay fixed
 void actions_before_newton_solve()
 {
  /* Reapply boundary conditions */
  apply_boundary_conditions();

  // Print some gubbins about parameters
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << "Solving for p=" << Params::p_mag
	     << "      c_swell=" << Params::c_swell << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
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
 // The initial (and maximum) element area
 double Element_area;
 // The mesh parameters
 TriangleMeshParameters* Triangle_mesh_parameters_pt;
 TriangleMeshClosedCurve* Boundary_pt;
 Vector<TriangleMeshOpenCurve *> Inner_open_boundaries_pt;
 TriangleMeshPolyLine* Boundary0_pt;
 TriangleMeshPolyLine* Boundary1_pt;
 TriangleMeshPolyLine* Boundary2_pt;
 TriangleMeshPolyLine* Boundary3_pt;
 TriangleMeshPolyLine* InnerBoudary0_pt;
 TriangleMeshPolyLine* InnerBoudary1_pt;
 
 /// Actions to be performed after read in of meshes
 void actions_after_read_unstructured_meshes()
 {
  // Curved Edges need to be upgraded after the rebuild
  //upgrade_edge_elements_to_curve(0,Bulk_mesh_pt);
  //upgrade_edge_elements_to_curve(1,Bulk_mesh_pt);
  // Rotate degrees of freedom
  // rotate_edge_degrees_of_freedom(Bulk_mesh_pt);
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
 ofstream Trace_file_dim, Trace_file_nondim;

 // Keep track of boundary ids
 enum
  {
   Outer_boundary0 = 0,
   Outer_boundary1 = 1,
   Outer_boundary2 = 2,
   Outer_boundary3 = 3,
   Inner_boundary0 = 4,
   Inner_boundary1 = 5
  }; 
 
 /// \short Loop over all curved edges, then loop over elements and upgrade
 /// them to be curved elements
 void upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const &
				     bulk_mesh_pt);

 /// \short Loop over all edge elements and rotate the Hermite degrees of freedom
 /// to be in the directions of the two in-plane vectors specified in TestSoln
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
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
 :
 Element_area(element_area)
{
 Problem::Always_take_one_newton_step = true;
 // Build the mesh
 build_mesh();

 // Store number of bulk elements
 complete_problem_setup();

 char filename[100];
 ofstream Param_file;
 Param_file.open("Parameters.dat");

 // Output plate parameters
 Param_file << "L1       " << Params::L1        << std::endl
	    << "L2       " << Params::L2        << std::endl
	    << "thickness" << Params::thickness << std::endl
	    << "nu       " << Params::nu        << std::endl
	    << "eta      " << Params::eta       << std::endl

	    << "L_dim    " << Params::L_dim     << std::endl
	    << "E_dim    " << Params::E_dim     << std::endl
	    << "P_dim    " << Params::P_dim     << std::endl;
  
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
 // Rectangular mesh boundary
 //
 //   V3        E2       V2
 //    O--------o--------O     ^
 //    |        V6       |     |
 //    |                 |     |
 // E3 oV7      oV8    V5o E1  L2
 //    |                 |     |
 //    |        V4       |     |
 //    O--------o--------O     v
 //   V0        E0       V1
 //    <--------L1------->
 
 double L1 = Params::L1;
 double L2 = Params::L2;

 
 unsigned N = Params::n_edge_nodes;
 double h = L1 / (N-1);
 Vector<Vector<double>> E0(N,Vector<double>(2,0.0)),
  E1(2,Vector<double>(2,0.0)),
  E2(N,Vector<double>(2,0.0)),
  E3(2,Vector<double>(2,0.0));
 for(unsigned i=0; i<N; i++)
  {
   E0[i][0] = -0.5*L1 + i*h;
   E0[i][1] = -0.5*L2;
   E2[i][0] =  0.5*L1 - i*h;
   E2[i][1] =  0.5*L2;
  }
 E1[0][0] =  0.5*L1;
 E1[0][1] = -0.5*L2;
 E1[1][0] =  0.5*L1;
 E1[1][1] =  0.5*L2;
 E3[0][0] = -0.5*L1;
 E3[0][1] =  0.5*L2;
 E3[1][0] = -0.5*L1;
 E3[1][1] = -0.5*L2;
  
 // Define boundaries from edges.
 Boundary0_pt = new TriangleMeshPolyLine(E0, 0);
 Boundary1_pt = new TriangleMeshPolyLine(E1, 1);
 Boundary2_pt = new TriangleMeshPolyLine(E2, 2);
 Boundary3_pt = new TriangleMeshPolyLine(E3, 3);
 
 Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);
 boundary_polyline_pt[0] = Boundary0_pt;
 boundary_polyline_pt[1] = Boundary1_pt;
 boundary_polyline_pt[2] = Boundary2_pt;
 boundary_polyline_pt[3] = Boundary3_pt;

 Boundary_pt = new TriangleMeshClosedCurve(boundary_polyline_pt);
 

 // Internal open boundaries
 // Total number of open curves in the domain
 unsigned n_open_curves = 2;
 // We want internal open curves
 Inner_open_boundaries_pt.resize(n_open_curves);
 
 // // Internal bit - this means we can have a boundary which is just the centre
 //InnerBoudary0_pt = new TriangleMeshPolyLine(E4, 4);
 //InnerBoudary1_pt = new TriangleMeshPolyLine(E5, 5);
 
 // // Connect the inner boundaries to the outer boundaries
 // InnerBoudary0_pt->connect_initial_vertex_to_polyline(Boundary1_pt,1);
 // InnerBoudary0_pt->connect_final_vertex_to_polyline(Boundary3_pt,1);
 // InnerBoudary1_pt->connect_initial_vertex_to_polyline(Boundary0_pt,1);
 // InnerBoudary1_pt->connect_final_vertex_to_polyline(Boundary2_pt,1);
 
 // // Each internal open curve is defined by a vector of
 // // TriangleMeshCurveSection,
 // on this example we only need one curve section for each internal boundary
 // Vector<TriangleMeshCurveSection *> internal_curve_section1_pt(1);
 // internal_curve_section1_pt[0] = Boundary2_pt;
 
 // Vector<TriangleMeshCurveSection *> internal_curve_section2_pt(1);
 // internal_curve_section2_pt[0] = Boundary3_pt;
 
 // // The open curve that define this boundary is composed of just one
 // // curve section
 //Inner_open_boundaries_pt[0] =
 // new TriangleMeshOpenCurve(internal_curve_section1_pt);
 
 //Inner_open_boundaries_pt[1] =
 // new TriangleMeshOpenCurve(internal_curve_section2_pt);

 
 TriangleMeshParameters Triangle_mesh_parameters(Boundary_pt);
 
 // // Specify the internal open boundaries
 //Triangle_mesh_parameters.internal_open_curves_pt() = Inner_open_boundaries_pt;

 // Set the maximum element area
 Triangle_mesh_parameters.element_area()=Element_area;
 
 // Build an assign bulk mesh
 Bulk_mesh_pt=new TriangleMesh<ELEMENT>(Triangle_mesh_parameters);
 
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
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{
 // Set the boundary conditions for problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.

 // The node at the very centre should be pinned when we have "do-nothing"
 // conditions on the in-plane displacements (i.e stress free b/c)
 // pin_in_plane_displacements_at_centre_node();

 // Complete the build of all elements so they are fully functional
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the pressure & temperature function pointers and the physical constants
   el_pt->pressure_fct_pt() = &Params::get_pressure;
   el_pt->swelling_fct_pt() = &Params::get_swelling;
   el_pt->in_plane_forcing_fct_pt() = &Params::get_in_plane_force;
   // There is no error metric in this case
   el_pt->error_metric_fct_pt() = &Params::axiasymmetry_metric;
   el_pt->nu_pt() = &Params::nu;
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
 unsigned n_bound = Bulk_mesh_pt->nboundary();
 for(unsigned b=0;b<n_bound;b++)
  {
   const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
   for(unsigned e=0; e<nb_element; e++)
    {
     // Get pointer to bulk element adjacent to b
     ELEMENT* el_pt =
      dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b,e));

     // Pin in-plane dofs
     el_pt->fix_in_plane_displacement_dof(0, b, Params::get_null_fct);
     if(b%2==0) // Only pin y for infinite edges
      {
       el_pt->fix_in_plane_displacement_dof(1, b, Params::get_null_fct);
      }
     
  
     for(unsigned idof=0; idof<6; ++idof)
      {
       // Symmetry conditions at x-normal edges (left and right)
       if(b%2==1)
	{
	 // Pin anything with an x derivative.
	 if(idof==1||idof==3||idof==4)
	  {
	   el_pt->fix_out_of_plane_displacement_dof(idof, b,
						    Params::get_null_fct);
	  }
	}
       // Clamp y-normal edges (top and bottom)
       if(b%2==0)
	{
	 // Pin all but the second y derivative.
	 if(idof!=5)
	  {
	   el_pt->fix_out_of_plane_displacement_dof(idof, b,
						    Params::get_null_fct);
	  }
	}
      }
    }
  }
} // end set bc


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
 npts = 20;

 tempstring = Params::output_dir + "/soln" + std::to_string(Doc_info.number()) + ".dat";
 strcpy(filename, tempstring.c_str());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
	   << comment << "\"\n";
 some_file.close();


 
 // Output regions
 unsigned n_region = Bulk_mesh_pt->nregion();
 if (n_region > 1)
  {
   for (unsigned r = 0; r < n_region; r++)
    {
     //Attempt to output elements in different regions
     tempstring = Params::output_dir + "/regions" + std::to_string(r) +  std::to_string(Doc_info.number()) + ".dat";
     strcpy(filename, tempstring.c_str());
     some_file.open(filename);
     unsigned nel = Bulk_mesh_pt->nregion_element(r);
     for (unsigned e = 0; e < nel; e++)
      {
       Bulk_mesh_pt->region_element_pt(r,e)->output(some_file,npts);
      }
     some_file.close();
    }
  }
 

 unsigned n_element=Bulk_mesh_pt->nelement();

 // Store a profile deflection slice through x=0
 //-----------------------------------------------------
 unsigned n_ppoint=101;
 Vector<double> ppoint(2,0.0), ppoint_s(2,0.0), z(8,0.0);
 GeomObject* ppoint_element_pt;
 double h = Params::L2/(n_ppoint-1);
  
 tempstring = Params::output_dir + "/profile" + std::to_string(Doc_info.number()) + ".dat";
 strcpy(filename, tempstring.c_str());
 some_file.open(filename);

 ppoint[1]=-0.5*Params::L2;
 for(unsigned j=0; j<n_ppoint; j++)
  {
   ppoint[1]+=h;
   for(unsigned i=0; i<n_element; i++)
    {
     dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i))
      ->locate_zeta(ppoint, ppoint_element_pt, ppoint_s);
     if(ppoint_element_pt!=NULL)
      {
       z = dynamic_cast<ELEMENT*>(ppoint_element_pt)
 	->interpolated_u_foeppl_von_karman(ppoint_s);
       some_file << ppoint[0]+z[6] << " " << ppoint[1]+z[7] << " " << z[0] << std::endl;
       break;
      }
    }  
  }
 some_file.close();
 

 // Write the pressure, degree of swelling and
 // deflection to the trace file
 //-----------------------------------------------------
 // Get the centre deflection first
 Vector<double> origin(2,0.0), s(2,0.0);
 GeomObject* centre_element_pt;
 Vector<double> w_centre(1,0.0);
  
 for(unsigned i=0; i<n_element; i++)
  {
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i))
    ->locate_zeta(origin, centre_element_pt, s);
   if(centre_element_pt!=NULL)
    {
     w_centre = dynamic_cast<ELEMENT*>(centre_element_pt)
      ->interpolated_u_foeppl_von_karman(s);
     break;
    }
  }
 
 Trace_file_nondim << Doc_info.number() << " "
		   << Params::p_mag   << " "
		   << Params::c_swell << " "
		   << w_centre[0] << endl;
 Trace_file_dim    << Doc_info.number()           << " "
		   << Params::p_mag*Params::P_dim << " "
		   << Params::c_swell             << " "
		   << w_centre[0]*Params::L_dim << endl;

 // Increment the doc_info number
 Doc_info.number()++;
 
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
 Params::output_dir="RESLT";
 CommandLineArgs::specify_command_line_flag("--dir", &Params::output_dir);

 // Poisson Ratio
 CommandLineArgs::specify_command_line_flag("--nu", &Params::nu);
 
 // Eta for the plate
 CommandLineArgs::specify_command_line_flag("--eta", &Params::eta);

 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--p", &Params::p_mag);

 // Swelling increase
 CommandLineArgs::specify_command_line_flag("--c", &Params::c_swell);

 // Element Area (no larger element than)
 double element_area=0.08;
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);

 // Parse command line
 CommandLineArgs::parse_and_assign();

 // How many nodes do we want along the long boundary (roughly length*width/thickness)
 Params::n_edge_nodes = ceil(0.2*Params::L2/Params::thickness);
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
 UnstructuredFvKProblem<ThermoFoepplVonKarmanC1CurvedBellElement<4>>
  problem(element_area);

 // Set up some problem paramters
 problem.newton_solver_tolerance()=1e-10;
 problem.max_residuals()=1e4;
 problem.max_newton_iterations()=30;
 //Params::output_dir+="_"+std::to_string(element_area);

 double p_inc = 10.0;
 // double c_inc = 0.001;

 // INITIAL SOLVE
 problem.newton_solve(); // SOLVE
 problem.doc_solution(); // AND DOCUMENT
 
 // INFLATION LOOP
 for(unsigned i=0; i<30; i++)
  {
   if(i%10==0)
    {p_inc*=10;}
   Params::p_mag += p_inc;
   problem.newton_solve(); // SOLVE
   problem.doc_solution(); // AND DOCUMENT
  }

 // // SWELLING LOOP
 // for(unsigned i=0; i<35; i++)
 //  {
 //   Params::c_swell += c_inc;
 //   problem.newton_solve(); // SOLVE
 //   problem.doc_solution(); // AND DOCUMENT
 //  }
 // c_inc/=5.0;
 // for(unsigned i=0; i<400; i++)
 //  {
 //   Params::c_swell += c_inc;
 //   problem.newton_solve(); // SOLVE
 //   problem.doc_solution(); // AND DOCUMENT
 //  }

 // Print success
 oomph_info<<"Exiting Normally\n";
} //End of main
