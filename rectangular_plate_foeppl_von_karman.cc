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
 double L1 = 5.0;
 double L2 = 1.0;
 double thickness = 0.0054;
 double nu = 0.495;
 double mu = 1000.0;
 double eta = 12*(1-nu*nu) / (thickness*thickness);
 double relevant_energy_factor = 20.0;
 // Control parameters
 double p_mag = 0.0;
 double c_swell = 0.0;
 Data *c_swell_data_pt;

 // Dimensions
 double L_dim = 5.0e-3;
 double E_dim = 1.2139e6;
 double P_dim = E_dim * thickness / eta;

 // Mesh parameters
 double element_area=0.1;
 unsigned n_long_edge_nodes=0;
 unsigned n_short_edge_nodes=0;

 // Output directory
 string output_dir="RESLT";

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
  swelling = c_swell_data_pt->value(0); // * (1 - pow(x[0]/L1,10) - pow(x[1]/L2,10)
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
 UnstructuredFvKProblem();

 /// Destructor
 ~UnstructuredFvKProblem()
 {
  Trace_file_nondim.close();
  Trace_file_dim.close();
  delete (Surface_mesh_pt);
  delete (Bulk_mesh_pt);
  // Clean up memory
  delete Boundary_pt;
  //  delete Inner_open_boundaries_pt[0];
  //  delete Inner_open_boundaries_pt[1];
  delete Boundary0_pt;
  delete Boundary1_pt;
  delete Boundary2_pt;
  delete Boundary3_pt;
  //  delete InnerBoudary0_pt;
  //  delete InnerBoudary1_pt;
 };

 /// Setup and build the mesh
 void build_mesh();

 /// Return the number of newton iterations in the last solve attempt
 unsigned nnewton_iter_taken()
  {
   return Nnewton_iter_taken;
  }
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
 
 /// Update after solve (empty)
 void actions_after_newton_solve()
 {
  /* No actions before newton solve */
 }

 /// Pin the in-plane displacements and set to zero at centre
 void pin_in_plane_displacements_at_centre_node();

 /// Print information about the parameters we are trying to solve for.
 void actions_before_newton_solve()
 {
  // Print some gubbins about parameters
  oomph_info << "-------------------------------------------------------" << std::endl;
  oomph_info << "Solving for p = " << Params::p_mag
	     << "  (" << Params::p_mag*Params::P_dim << "Pa)" << std::endl
	     << "      c_swell = " << Params::c_swell_data_pt->value(0) << std::endl
	     << "     Doc_info = " << Doc_steady_info.number()
	     << ", " << Doc_unsteady_info.number() << std::endl;
  oomph_info << "-------------------------------------------------------" << std::endl;
 }

 /// Adaptively try to swell the membrane until it becomes unfeasible.
 void adaptive_swell_solve(double c_inc);

 /// Use oomph-lib's pseudo arclength continuation to try swell the membrane.
 void arclength_swell_solve(double c_inc);

 /// Attempt an ordinary steady solve, but failing that solve an unsteady damped
 /// version of the equations with a run of adaptive unsteady solves.
 void damped_solve(double dt, double epsilon, bool doc_unsteady=false);

 /// Doc the solution
 void doc_solution(const bool steady, const std::string& comment="");

 /// \short Overloaded version of the problem's access function to
 /// the mesh. Recasts the pointer to the base Mesh object to
 /// the actual mesh type.
 TriangleMesh<ELEMENT>* mesh_pt()
 {
  return Bulk_mesh_pt;
 }

 /// Doc info object for labeling all output
 DocInfo Doc_steady_info;
 DocInfo Doc_unsteady_info;
 
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

 /// Helper function to apply boundary conditions
 void apply_boundary_conditions();

 /// \short Helper function to (re-)set boundary condition
 /// and complete the build of all elements
 void complete_problem_setup();

 /// Trace file to document norm of solution
 ofstream Trace_file_dim, Trace_file_nondim;

 // Keep track of boundary ids
 enum
  {
   Outer_boundary0 = 0,
   Outer_boundary1 = 1,
   Outer_boundary2 = 2,
   Outer_boundary3 = 3
  }; 

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
 // Build the mesh
 build_mesh();

 // Store number of bulk elements
 complete_problem_setup();

 char filename[100];
 ofstream Param_file;
 strcpy(filename, (Params::output_dir + "/parameters.dat").c_str());
 Param_file.open(filename);

 // Output plate parameters
 Param_file << "L1           " << Params::L1        << std::endl
	    << "L2           " << Params::L2        << std::endl
	    << "thickness    " << Params::thickness << std::endl
	    << "nu           " << Params::nu        << std::endl
	    << "eta          " << Params::eta       << std::endl
	    << std::endl
	    << "L_dim        " << Params::L_dim     << std::endl
	    << "E_dim        " << Params::E_dim     << std::endl
	    << "P_dim        " << Params::P_dim     << std::endl
	    << std::endl
	    << "Element area " << Params::element_area << std::endl
	    << "N_xedgenodes " << Params::n_long_edge_nodes << std::endl
	    << "N_yedgenodes " << Params::n_short_edge_nodes << std::endl;
 
 strcpy(filename, (Params::output_dir + "/trace_nondim.dat").c_str());
 Trace_file_nondim.open(filename);
 strcpy(filename, (Params::output_dir + "/trace_dim.dat").c_str());
 Trace_file_dim.open(filename);

 oomph_info << "Number of equations: "
	    << assign_eqn_numbers() << '\n';

 // Output the mesh
 strcpy(filename, (Params::output_dir + "/mesh.dat").c_str());
 ofstream meshfile(filename);
 Problem::mesh_pt()->output(meshfile);
 meshfile.close();
 
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

 
 unsigned Nl = Params::n_long_edge_nodes;
 unsigned Ns = Params::n_short_edge_nodes;
 double hl = L1 / (Nl-1);
 double hs = L2 / (Ns-1);
 Vector<Vector<double>> E0(Nl,Vector<double>(2,0.0)),
  E1(Ns,Vector<double>(2,0.0)),
  E2(Nl,Vector<double>(2,0.0)),
  E3(Ns,Vector<double>(2,0.0));
 for(unsigned i=0; i<Nl; i++)
  {
   E0[i][0] = -0.5*L1 + i*hl;
   E0[i][1] = -0.5*L2;
   E2[i][0] =  0.5*L1 - i*hl;
   E2[i][1] =  0.5*L2;
  }
 for(unsigned i=0; i<Ns; i++)
  {
   E1[i][0] =  0.5*L1;
   E1[i][1] = -0.5*L2 + i*hs;
   E3[i][0] = -0.5*L1;
   E3[i][1] =  0.5*L2 - i*hs;
  }
  
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

 TriangleMeshParameters Triangle_mesh_parameters(Boundary_pt);
 
 // Set the maximum element area
 Triangle_mesh_parameters.element_area()=Element_area;
 
 // Build an assign bulk mesh
 Bulk_mesh_pt=new TriangleMesh<ELEMENT>(Triangle_mesh_parameters, time_stepper_pt());
 
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
 // Create a new Data object whose one-and-only value contains the
 // (in principle) adjustible load
 Params::c_swell_data_pt = new Data(1);

 // And pin the pressure as it is our control.
 Params::c_swell_data_pt->pin(0);
 
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
     el_pt->fix_in_plane_displacement_dof(1, b, Params::get_null_fct);
    

     // Out-of-plane dofs:
     //   0 |  1  |  2  |  3  |  4  |  5
     //   w | w_x | w_y | w_xx| w_xy| w_yy
     for(unsigned idof=0; idof<6; ++idof)
      {
       // Clamp x-normal edges (left and right)
       if(b%2==1)
	{
	 // Pin all but the second x-derivative. (x derivativeS if pinned)
	 if(idof!=1&&idof!=3)
	  {
	   el_pt->fix_out_of_plane_displacement_dof(idof, b,
						    Params::get_null_fct);
	  }
	}
       // Clamp y-normal edges (top and bottom)
       if(b%2==0)
	{
	 // Pin all but the second y-derivative. (y derivativeS if pinned)
	 if(idof!=2&&idof!=5)
	  {
	   el_pt->fix_out_of_plane_displacement_dof(idof, b,
						    Params::get_null_fct);
	  }
	}
      }
    }
  }
} // end set bc




//==start_of_adaptive_swell_solve=========================================
/// Try to swell the membrane adaptively
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::adaptive_swell_solve(double c_inc)
{
    
 // Max number of newton iterations before shrinking step size
 unsigned max_desired_iterations = 20;
   
 // Min number of newton iterations before increasing step size
 unsigned min_desired_iterations = 5;

 // Step size relaxation multiplier
 double relaxation_const = 3.0/2.0;

 // Step size refinement multiplier
 double refinement_const = relaxation_const*relaxation_const;

 // Set absolute largest allowed value for dh
 double c_inc_max = 3.0*c_inc;
   
 // Set absolute smallest allowed value for dh
 double c_inc_min = c_inc / 1000.0;

 // Are things getting silly?
 bool getting_silly = false;

 while(!getting_silly)
  {
   // Loop for trying to solve a single displcement step
   
   // Make a copy of the current displacement incase the step fails
   double old_c = Params::c_swell_data_pt->value(0);
      
   // Make a copy of initial values incase the step fails
   DoubleVector old_dofs;
   get_dofs(old_dofs);
   
   // The number of newton iterations taken by the last solve
   unsigned n_newton_iterations = 0;
   
   // Has this step failed to converge
   bool step_rejected = false;
 
   do
    {
     cout << endl << endl
	  << "Trying to 'do' with c_inc = " << c_inc << endl;

       
     // Reset the initial values for this attempt at the step
     *(Params::c_swell_data_pt->value_pt(0)) = old_c;
     set_dofs(old_dofs);

     // Check we aren't swelling too slowly to be meaningful
     if (c_inc < c_inc_min)
      {
       oomph_info << std::endl << std::endl
	    << "This is just getting silly... c_inc is " << c_inc
	    << "! Which is smaller than " << c_inc_min
	    << std::endl << "It's time to give up..." << std::endl << std::endl;
       getting_silly=true;
       break;
      }
       
     // Assume for now this step wont be rejected
     step_rejected = false;

     // Increment the pinned value of c_swell
     *(Params::c_swell_data_pt->value_pt(0)) += c_inc;

     // Time to solve.
     // Actually do the newton solve stage
     try
      {
       newton_solve();
       n_newton_iterations = nnewton_iter_taken();
      }
   
     // Catch any exceptions thrown in the Newton solver
     catch (NewtonSolverError& error)
      {
       // Check whether it's the linear solver
       if (error.linear_solver_error)
	{
	 std::ostringstream error_stream;
	 error_stream << std::endl
		      << "USER-DEFINED ERROR IN NEWTON SOLVER " << std::endl;
	 oomph_info << "ERROR IN THE LINEAR SOLVER" << std::endl;
	 throw OomphLibError(error_stream.str(),
			     OOMPH_CURRENT_FUNCTION,
			     OOMPH_EXCEPTION_LOCATION);
	}
       // Otherwise mark the step as having failed
       else
	{
	 oomph_info << "STEP REJECTED --- TRYING AGAIN" << std::endl;
	 step_rejected = true;
	 // Let's take a smaller step
	 c_inc *= refinement_const;
	}
      } // End of error catch
    } while (step_rejected); // End of displacement increment

   if(getting_silly)
    {
     oomph_info << "Stopping the adaptive swell solve now..."
		<< std::endl << std::endl;
     break;
    }
   
   // Output solution
   doc_solution();
     
   // Time to reconsider our increment size dh...
   if (n_newton_iterations > max_desired_iterations)
    {
     c_inc *= refinement_const;
    }
   else if (n_newton_iterations < min_desired_iterations
	    && c_inc* relaxation_const < c_inc_max)
    {
     c_inc *= relaxation_const;
    }
  }
} // End of adaptive swell solve


//==start of arclength_swell_solve=========================================
/// Try to swell the membrane using pseudo-arclength continuation.
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::arclength_swell_solve(double c_inc)
{
 // Now we pin the pressure data
 Params::c_swell_data_pt->pin(0);
 
 //Re-assign the equation numbers
 cout << std::endl;
 cout << "-----------------ARC-LENGTH CONTINUATION --------------" 
      << std::endl;
 cout << "# of dofs " << assign_eqn_numbers() << std::endl;
 cout << "-------------------------------------------------------" 
      << std::endl;
 cout << std::endl;

 
 // Set the maximum number of Newton iterations
 Max_newton_iterations=7;

 // Set the minimum desired number of iterations before increasing arclength
 Desired_newton_iterations_ds=4;

 // Set the proportion of arclength??
 Desired_proportion_of_arc_length = 0.1;

 // Check for bifurcations using the sign of the determinant of the Jacobian
 Bifurcation_detection = true;
 Bisect_to_find_bifurcation =
  CommandLineArgs::command_line_flag_has_been_set("--Find_membrane_collapse");

 // Set an initial value for the step size
 double ds = 0.01;

 // Take ~so-many~ steps
 for(unsigned i=0;i<200;i++)
   {
    // Solve
    ds = arc_length_step_solve(Params::c_swell_data_pt, 0, ds);

    oomph_info << std::endl
	       << "An arclength step has been successfully completed :)"
	       << std::endl;

    // Now document the result!
    doc_solution();
    
    if(First_jacobian_sign_change)
     {
      oomph_info << std::endl << std::endl
		 << "We have passed a bifurcation!"
		 << std::endl << std::endl;
      First_jacobian_sign_change=false;
     }
   }
} // End of arclength_swell_solve()

template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::damped_solve(double dt,
						   double epsilon,
						   bool doc_unsteady)
{
 bool STEADY = false;
 while(!STEADY)
  {
   oomph_info << "NEW DAMPED PSEUDO-TIME STEP" << std::endl;
   // Try get us close to a steady solution
   double dt_next =
    adaptive_unsteady_newton_solve(dt,epsilon);
   dt = dt_next;
   if(doc_unsteady)
    {
     doc_solution(false);
    }
   
   // Now compare energies to check whether we are almost steady or not.
   double kinetic_energy = 0.0;
   double elastic_energy = 0.0;
   unsigned n_el = mesh_pt()->nelement();
   for(unsigned e=0; e<n_el; e++)
    {
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
     Vector<double> element_energies = el_pt->element_elastic_and_kinetic_energy();
     elastic_energy += element_energies[0];
     kinetic_energy += element_energies[1];
    }
   oomph_info << "Elastic energy: " << elastic_energy << std::endl
	      << "Kinetic energy: " << kinetic_energy << std::endl;
   double R = Params::relevant_energy_factor;
   if(R*abs(kinetic_energy) < abs(elastic_energy))
    {
     oomph_info << "ALMOST STEADY SO ATTEMPT A STEADY SOLVE" << std::endl;
     store_current_dof_values();
     try
      {
       steady_newton_solve();
       if(doc_unsteady)
	{
	 doc_solution(false);
	}
       time()=0.0;
       doc_solution(true);
       STEADY=true;
      }
     catch(OomphLibError& error)
      {
       oomph_info << "NOT STEADY ENOUGH. REDUCING THE KINETIC ENERGY PROPORTION FROM "
		  << R << " TO " << 10.0*R << " FROM NOW ON" << std::endl;
       R*=0.1;
       restore_dof_values();
       for (unsigned i=0; i<ntime_stepper(); i++)
	{
	 time_stepper_pt(i)->undo_make_steady();
	}
       assign_initial_values_impulsive();
       error.disable_error_message();
      }
    }
  } // End of while(UNSTEADY)
}



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(bool steady,
						   const std::string& comment)
{
 ofstream some_file;
 char filename[100];
 string tempstring;
  
 // Number of plot points
 unsigned npts = 2;

 if(steady)
  {
   tempstring = Params::output_dir + "/coarse_soln_"
    + std::to_string(Doc_steady_info.number()) + ".dat";
  }
 else
  {
   tempstring = Params::output_dir + "/coarse_soln_"
    + std::to_string(Doc_steady_info.number()) + "_"
    + std::to_string(Doc_unsteady_info.number()) + ".dat";
  }
 strcpy(filename, tempstring.c_str());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
	   << comment << "\"\n";
 some_file.close();

 // Number of plot points
 npts = 10;

 if(steady)
  {
   tempstring = Params::output_dir + "/soln_"
    + std::to_string(Doc_steady_info.number()) + ".dat";
  }
 else
  {
   tempstring = Params::output_dir + "/soln_"
    + std::to_string(Doc_steady_info.number()) + "_"
    + std::to_string(Doc_unsteady_info.number()) + ".dat";
  }
 strcpy(filename, tempstring.c_str());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
	   << comment << "\"\n";
 some_file.close();


 // Write the pressure, degree of swelling and
 // deflection to the trace file
 //-----------------------------------------------------
 // Get the centre deflection first
 Vector<double> origin(2,0.0), s(2,0.0);
 GeomObject* centre_element_pt;
 Vector<double> w_centre(1,0.0);
  
 unsigned n_element=Bulk_mesh_pt->nelement();
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
 
 Trace_file_nondim << Doc_steady_info.number()          << " "
		   << Doc_unsteady_info.number()        << " "
		   << time()                            << " "
		   << Params::p_mag                     << " "
		   << Params::c_swell_data_pt->value(0) << " "
		   << w_centre[0]                       << endl;
 Trace_file_dim    << Doc_steady_info.number()          << " "
		   << Doc_unsteady_info.number()        << " "
		   << time()                            << " "
		   << Params::p_mag*Params::P_dim       << " "
		   << Params::c_swell_data_pt->value(0) << " "
		   << w_centre[0]*Params::L_dim         << endl;

 // Increment the doc_info numbers
 if(steady)
  {
   Doc_steady_info.number()++;
   Doc_unsteady_info.number()=0;
  }
 else
  {
   Doc_unsteady_info.number()++;
  }
 
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
 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 // Directory for solution
 Params::output_dir="RESLT";
 CommandLineArgs::specify_command_line_flag("--dir",
					    &Params::output_dir);

 // Poisson Ratio
 CommandLineArgs::specify_command_line_flag("--nu",
					    &Params::nu);
 
 // Dampening coefficient
 CommandLineArgs::specify_command_line_flag("--mu",
					    &Params::mu);
 
 // Eta for the plate
 CommandLineArgs::specify_command_line_flag("--eta",
					    &Params::eta);

 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--p",
					    &Params::p_mag);

 // Temperature increase
 CommandLineArgs::specify_command_line_flag("--t",
					    &Params::c_swell);

 // Element Area (no larger element than)
 Params::element_area=0.009;
 CommandLineArgs::specify_command_line_flag("--element_area",
					    &Params::element_area);
 
 // Pin u_alpha everywhere
 CommandLineArgs::specify_command_line_flag("--pininplane");


 // Parse command line
 CommandLineArgs::parse_and_assign();

 // How many nodes do we want to manually place along the boundaries
 // (roughly length*width/thickness)
 Params::n_long_edge_nodes  = ceil(0.04*Params::L1/Params::thickness) + 2;
 Params::n_short_edge_nodes = ceil(0.04*Params::L2/Params::thickness) + 2;
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // CREATE THE PROBLEM
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

 // Set up some problem paramters
 problem.newton_solver_tolerance()=1e-10;
 problem.max_residuals()=1e4;
 problem.max_newton_iterations()=10;

 double dt = 1.0;
 double epsilon = 1.0e-3;

 problem.assign_initial_values_impulsive();
 problem.initialise_dt(dt);
 
 double p_inc = 24000.0;
 double c_inc = 0.002;

 oomph_info << "DO AN INITIAL STATE SOLVE" << std::endl;
 // INITIAL SOLVE
 problem.steady_newton_solve(); // SOLVE
 problem.doc_solution(true); // AND DOCUMENT

 // INFLATION
 Params::p_mag+=p_inc;   // 10
 
 // Pre-inflate the membrane
 oomph_info << "INFLATION STAGE" << std::endl;
 problem.damped_solve(dt, epsilon, true);

 // Swell the membrane
 oomph_info << "SWELLING STAGE" << std::endl;
 Params::mu*=1.0;
 for(unsigned i=0; i<200; i++)
  {
   Params::c_swell += c_inc;
   Params::c_swell_data_pt->set_value(0, Params::c_swell);
   oomph_info << "c_swell = " << Params::c_swell << std::endl;
   problem.damped_solve(dt, epsilon, true);
  } // End of swelling loop
 


// SWELLING LOOP
 //problem.arclength_swell_solve(c_inc);


 // Print success
 oomph_info<<"Exiting Normally\n";


} //End of main
