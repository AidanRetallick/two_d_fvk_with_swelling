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


template<class ELEMENT>
class MyMesh : public virtual TriangleMesh<ELEMENT>
{
public:
 MyMesh<ELEMENT>(TriangleMeshParameters& triangle_mesh_parameters,
		 TimeStepper *time_stepper_pt=&Mesh::Default_TimeStepper)
  : TriangleMesh<ELEMENT>(triangle_mesh_parameters, time_stepper_pt)
 {}

 void update_polyline_representation_from_restart()
 {
  oomph_info << "Called MyMesh::update_polyline_representation_from_restart()"
	     << std::endl
	     << "           NOT DOING ANYTHING IN THIS FUNCTION" << std::endl;
 }
};
 
namespace Params
{
 // Plate parameters
 double L1 = 5.0;
 double L2 = 1.0;
 double thickness = 0.0054;
 double nu = 0.495;
 double mu = 1000.0;
 double eta = 12*(1-nu*nu) / (thickness*thickness);
 double relevant_energy_factor = 100.0;

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

 // Swelling wrapper so we can output the degree of swelling function
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
  delete Boundary_pt;
  //  delete Inner_open_boundaries_pt[0];
  //  delete Inner_open_boundaries_pt[1];
  delete Boundary0_pt;
  delete Boundary1_pt;
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
	     << "  (" << Params::p_mag*Params::P_dim << "Pa)" << std::endl;
  oomph_info << "      c_swell = " << Params::c_swell_data_pt->value(0) << std::endl;
  oomph_info << "     Doc_info = " << Doc_steady_info.number()
	     << ", " << Doc_unsteady_info.number() << std::endl;
  oomph_info << "-------------------------------------------------------" << std::endl;
 }

 /// Remove surface mesh before reading
 void actions_before_read_unstructured_meshes()
 {
  // How many surface elements are in the surface mesh
  unsigned n_element = Surface_mesh_pt->nelement();

  // Loop over the surface elements
  for(unsigned e=0;e<n_element;e++)
   {
    // Kill surface element
    delete Surface_mesh_pt->element_pt(e);
   }

  Surface_mesh_pt->flush_element_and_node_storage();

  rebuild_global_mesh();
 }

 /// Add surface mesh back after reading
 void actions_after_read_unstructured_meshes()
 {
  rebuild_global_mesh();
  apply_boundary_conditions();
  complete_problem_setup();
 }

 /// Attempt an ordinary steady solve, but failing that solve an unsteady damped
 /// version of the equations with a run of adaptive unsteady solves.
 void damped_solve(double dt, double epsilon, bool doc_unsteady=false);

 /// Doc the solution
 void doc_solution(const bool steady, const std::string& comment="");

 /// Dump problem to disk to allow for restart.
 void dump_it(ofstream& dump_file)
 {
  dump_file << Doc_steady_info.number() << " # steady step number" << std::endl
	    << Doc_unsteady_info.number() << " # unsteady step number" << std::endl
	    << Params::p_mag << " # pressure" << std::endl
	    << Params::c_swell_data_pt->value(0) << " # swelling" << std::endl
	    << Next_dt << " # suggested next timestep" << std::endl;
  
  // Call generic dump()
  Problem::dump(dump_file); 
 }
 
 /// Read problem for restart from specified restart file.
 void restart(ifstream& restart_file)
 {
  string input_string;
  
  // Read line up to termination sign
  getline(restart_file,input_string,'#');
  // Ignore rest of line
  restart_file.ignore(80,'\n');
  // Read in steady step number
  Doc_steady_info.number()=unsigned(atof(input_string.c_str()));
  Doc_steady_info.number()++;
  
  getline(restart_file,input_string,'#');
  // Ignore rest of line
  restart_file.ignore(80,'\n');
  // Read in unsteady step number
  Doc_unsteady_info.number()=unsigned(atof(input_string.c_str()));
  Doc_unsteady_info.number()=0; //for now [witnessme]
  
  getline(restart_file,input_string,'#');
  // Ignore rest of line
  restart_file.ignore(80,'\n');
  // Read in pressure value
  Params::p_mag=double(atof(input_string.c_str()));

  getline(restart_file,input_string,'#');
  // Ignore rest of line
  restart_file.ignore(80,'\n');
  // Read in steady step number
  Params::c_swell=double(atof(input_string.c_str()));
  
  // Read line up to termination sign
  getline(restart_file,input_string,'#');
  // Ignore rest of line
  restart_file.ignore(80,'\n'); 
  // Read suggested next timestep
  Next_dt=double(atof(input_string.c_str()));
 
  // Refine the mesh and read in the generic problem data
  Problem::read(restart_file);
 
 } // end of restart
  
 /// \short Overloaded version of the problem's access function to
 /// the mesh. Recasts the pointer to the base Mesh object to
 /// the actual mesh type.
 MyMesh<ELEMENT>* mesh_pt()
 {
  return Bulk_mesh_pt;
 }

 /// Access function for Next_dt
 double next_dt()
 {
  return Next_dt;
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
 MyMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to "surface" mesh
 Mesh* Surface_mesh_pt;
 
 /// Suggestion for next timestep (stored to allow it to be written
 /// to (or read from) restart file)
 double Next_dt=0.1;

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
 
#ifdef OOMPH_HAS_MUMPS
 // Let everyone know we are going to use MUMPS
 oomph_info << std::endl << "Using MUMPS solver" << std::endl << std::endl;
 
 // Change solver
 linear_solver_pt()=new MumpsSolver;

 // Shut up
 dynamic_cast<MumpsSolver*>(linear_solver_pt())->
  enable_suppress_warning_about_MPI_COMM_WORLD();

#endif
 
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
 Bulk_mesh_pt=new MyMesh<ELEMENT>(Triangle_mesh_parameters, time_stepper_pt());
 
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
 DTSF_max_increase=2.0;
 DTSF_min_decrease=0.5;
 
 // Create a new Data object whose one-and-only value contains the
 // (in principle) adjustible load
 Params::c_swell_data_pt = new Data(1);

 // And pin the pressure as it is our control.
 Params::c_swell_data_pt->pin(0);
 Params::c_swell_data_pt->set_value(0,Params::c_swell);
 
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
   el_pt->nu_pt() = &Params::nu;
   el_pt->mu_pt() = &Params::mu;
   el_pt->eta_pt() = &Params::eta;
  }

 // Do we want to pin in-plane displacement?
 if(CommandLineArgs::command_line_flag_has_been_set("--pininplane"))
  {
   cout << "gonna pin em" << endl;
   // Pin the in-plane displacements
   unsigned nnode = mesh_pt()->nnode();
   for(unsigned inode=0; inode<nnode; inode++)
    {
     mesh_pt()->node_pt(inode)->pin(0);
     mesh_pt()->node_pt(inode)->pin(1);
    }
   cout << "successfully pinned" << endl;
   assign_eqn_numbers();
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

template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::damped_solve(double dt,
						   double epsilon,
						   bool doc_unsteady)
{
 bool STEADY = false;
 double R = Params::relevant_energy_factor;
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
   else
    {
     Doc_unsteady_info.number()++;
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
       if(R>1.0e20)
	{
	 oomph_info << "WARNING";
	 oomph_info << "R is now " << R << " which is smaller than it has any";
	 oomph_info << "right to be. Giving up on damped solves..." << std::endl;
	 throw error;
	}
       else
	{
	 oomph_info << "NOT STEADY ENOUGH. INCREASING THE KINETIC ENERGY PROPORTION FROM "
		    << R << " TO " << 100.0*R << " FROM NOW ON" << std::endl;
	 R*=100.0;
	 restore_dof_values();
	 for (unsigned i=0; i<ntime_stepper(); i++)
	 {
	  time_stepper_pt(i)->undo_make_steady();
	 }
	 assign_initial_values_impulsive(dt);
	 error.disable_error_message();
	}
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
 if(MPI_Helpers::communicator_pt()->my_rank()==0)
  {
   ofstream some_file;
   char filename[100];

   // Dump to a restart file if we are steady
   if(steady)
    {
     // Write restart file
     sprintf(filename,"%s/restart%i.dat",
   	     Params::output_dir.c_str(),
   	     Doc_steady_info.number());
     some_file.open(filename);
     dump_it(some_file);
     some_file.close();
    }
 
   unsigned npts;
   // // Number of plot points for coarse output
   // npts = 2;
   // if(steady)
   //  {
   //   sprintf(filename, "%s/coarse_soln_%i.dat",
   // 	     Params::output_dir.c_str(),
   // 	     Doc_steady_info.number());
   //  }
   // else
   //  {
   //   sprintf(filename, "%s/coarse_soln_%i_%i.dat",
   // 	     Params::output_dir.c_str(),
   // 	     Doc_steady_info.number(),
   // 	     Doc_unsteady_info.number());
   //  }
   // some_file.open(filename);
   // Bulk_mesh_pt->output(some_file,npts);
   // some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
   // 	     << comment << "\"\n";
   // some_file.close();

   // Number of plot points for fine outpout
   npts = 5;
   if(steady)
    {
     sprintf(filename, "%s/soln_%i.dat",
	     Params::output_dir.c_str(),
	     Doc_steady_info.number());
    }
   else
    {
     sprintf(filename, "%s/soln_%i_%i.dat",
	     Params::output_dir.c_str(),
	     Doc_steady_info.number(),
	     Doc_unsteady_info.number());
    }
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

  }
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
 // Initialise oomph-lib's MPI       
 MPI_Helpers::init(argc,argv);     
 
 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 // Directory for solution
 Params::output_dir="RESLT";
 CommandLineArgs::specify_command_line_flag("--dir",
					    &Params::output_dir);

 string restart_file_name;
 // Are we restarting from a dumped file?
 CommandLineArgs::specify_command_line_flag("--restart",
					    &restart_file_name);
 
 // Channel length Ratio
 CommandLineArgs::specify_command_line_flag("--L",
					    &Params::L1);
 
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
 double p_inc = 24000.0;
 CommandLineArgs::specify_command_line_flag("--p",
					    &p_inc);

 // Temperature increase
 CommandLineArgs::specify_command_line_flag("--t",
					    &Params::c_swell);

 // Element Area (no larger element than)
 Params::element_area=0.001;
 CommandLineArgs::specify_command_line_flag("--element_area",
					    &Params::element_area);
 
 // Pin u_alpha everywhere
 CommandLineArgs::specify_command_line_flag("--pininplane");


 // Parse command line
 CommandLineArgs::parse_and_assign();

 // How many nodes do we want to manually place along the boundaries
 // (roughly length*width/thickness)
 Params::n_long_edge_nodes  = ceil(0.00*Params::L1/Params::thickness) + 2;
 Params::n_short_edge_nodes = ceil(0.00*Params::L2/Params::thickness) + 2;
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // CREATE THE PROBLEM
 UnstructuredFvKProblem<ThermoFoepplVonKarmanC1CurvedBellElement<4>>
  problem;

 // Set up some problem paramters
 problem.newton_solver_tolerance()=1e-8;
 problem.max_residuals()=1e4;
 problem.max_newton_iterations()=10;

 if(CommandLineArgs::command_line_flag_has_been_set("--restart"))
  {
   oomph_info << "We are restarting from: " << restart_file_name << std::endl;
   ifstream restart_stream;
   restart_stream.open(restart_file_name);
   problem.restart(restart_stream);
   restart_stream.close();
  }

 double dt = problem.next_dt();
 double epsilon = 1.0e-3;
 problem.assign_initial_values_impulsive();
 problem.initialise_dt(dt);
 
 double c_inc = 0.002;

 if(CommandLineArgs::command_line_flag_has_been_set("--restart"))
  {
   problem.Doc_steady_info.number()-=1;
  }
 
 oomph_info << "DO AN INITIAL STATE SOLVE" << std::endl;
 // INITIAL SOLVE
 problem.steady_newton_solve(); // SOLVE
 problem.doc_solution(true); // AND DOCUMENT

 if(!CommandLineArgs::command_line_flag_has_been_set("--restart"))
  {
   //   while(Params::p_mag<p_max)
    {
     // INFLATION
     Params::p_mag+=p_inc;
     
     // Pre-inflate the membrane
     oomph_info << "INFLATION STAGE" << std::endl;
     problem.damped_solve(dt, epsilon, false);
    }  
  }
 
 // Swell the membrane
 oomph_info << "SWELLING STAGE" << std::endl;
 double c_swell_max=0.30;
 while(Params::c_swell<c_swell_max)
  {
   Params::c_swell += c_inc;
   Params::c_swell_data_pt->set_value(0, Params::c_swell);
   oomph_info << "c_swell = " << Params::c_swell << std::endl;
   // bool argument to specify whether we want to doc unsteady time steps.
   problem.damped_solve(dt, epsilon, false);
  } // End of swelling loop

 // Print success
 oomph_info<<"Exiting Normally\n";

 // Shut down oomph-lib's MPI
 MPI_Helpers::finalize();
} //End of main
