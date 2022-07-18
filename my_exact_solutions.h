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
namespace oomph {

// Class to contain exact solutions to Foppl von Karman systems
class FoepplVonKarmanExactSolution
{
 public:
 /// Constructor
 FoepplVonKarmanExactSolution() : Nu_pt(0), Pressure_mag_pt(0), Eta_pt(0) {};
 /// Destructor
 ~FoepplVonKarmanExactSolution(){};
 /// Exact solution
 virtual void exact_solution(const Vector<double>& X, Vector<double>& 
   displacements) const= 0;
 /// Get Nu
 double nu() const {return (Nu_pt!=0 ? *Nu_pt:/*Default*/0.5);};
 /// Get Nu_pt
 double*& nu_pt() {return Nu_pt;};
 /// Get Nu_pt (const version)
 const double* nu_pt() const {return Nu_pt;};
 /// Get Eta
 double eta() const {return (Eta_pt!=0 ? *Eta_pt:/*Default*/1.0);};
 /// Get Nu_pt
 double*& eta_pt() {return Eta_pt;};
 /// Get Nu_pt (const version)
 const double* eta_pt() const {return Eta_pt;};
 /// Get pressure_mag pt
 double*& pressure_mag_pt() {return Pressure_mag_pt;};
 /// Get pressure_mag pt (const version)
 const double* pressure_mag_pt() const {return Pressure_mag_pt;};
 // Out--of--plane forcing
 virtual void get_pressure(const Vector<double>& X, double& pressure) const
   {
   throw OomphLibError(
    "You must pressure forcing for your own object! \n",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
   };
 /// In plane forcing
 virtual void get_in_plane_force(const Vector<double>& X, Vector<double>& grad)const
   {
   throw OomphLibError(
    "You must in plane forcing for your own object! \n",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
   };
 /// Exact radial solution (broken)
 virtual void exact_radial_w_solution(const Vector<double>& X, Vector<double>& displacements)const
   {
   throw OomphLibError(
    "You must radial solution for your own object! \n",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
   };
 protected:
 /// Pointer to value of Poisson ratio
 double* Nu_pt;
 /// Pointer to pressure magnitude
 double* Pressure_mag_pt;
 /// Pointer to eta
 double* Eta_pt;
};

///\short Manufactured solution that has fourth order (axisymmetric) deflection
// and linear displacements 
class ManufacturedSolutionWithLinearDisplacements : protected FoepplVonKarmanExactSolution
{
 public:
 /// Constructor
 ManufacturedSolutionWithLinearDisplacements() : FoepplVonKarmanExactSolution(){};
 /// Constructor with nu_pt
 ManufacturedSolutionWithLinearDisplacements(double* poisson_pt, double* the_eta_pt) 
  : FoepplVonKarmanExactSolution()
  {
   nu_pt() = poisson_pt;
   eta_pt() = the_eta_pt;
  };

 // Assigns the value of pressure depending on the position (x,y)
 void get_pressure(const Vector<double>& X, double& pressure)const
 {
  // Convenient definitions
  double x=X[0], y=X[1];
  // Return the pressure
  pressure =-(-8*pow(-1 + pow(nu(), 2), -1)*(-2 + 3*eta()*(1 + nu())*(-1 + 2*pow(x, 2)
    + 2*pow(y, 2)) + 24*pow(eta(), 2)*(pow(x, 2) + pow(y, 2))*(-2 + 5*pow(x, 2) +
   5*pow(y, 2))*pow(-1 + pow(x, 2) + pow(y, 2), 2)))/3.;
  }

 // Assigns the value of in plane forcing depending on the position (x,y)
 void get_in_plane_force(const Vector<double>& X, Vector<double>& grad)const
  {
   grad.resize(2);
   // Convenient definitions
   double x=X[0], y=X[1];
   // Return the (constant) pressure
   grad[0] = 8*eta()*x*(-1 + pow(x, 2) + pow(y, 2))*(3 - nu() + (-7 + nu())*pow(x, 2) + 
   (-7 + nu())*pow(y, 2))*pow(-1 + pow(nu(), 2), -1); 
   grad[1] = 8*eta()*y*(-1 + pow(x, 2) + pow(y, 2))*(3 - nu() + (-7 + nu())*pow(x, 2) +
   (-7 + nu())*pow(y, 2))*pow(-1 + pow(nu(), 2), -1);
  }  
 
 void exact_radial_w_solution(const Vector<double>& X, 
  Vector<double>& exact_w)const
 {
  // Convenient definitions
  double x=X[0], y=X[1];
  exact_w.resize(6);
  exact_w[0] = pow(x*x+y*y-1,2);
  exact_w[1] = 4*sqrt(x*x+y*y)*(x*x+y*y-1);
  exact_w[2] = 0.0;
  exact_w[3] = 4*(3*(x*x+y*y)-1);
  exact_w[4] = 0.0;
  exact_w[5] = 0.0;
 }
 
 // Get the exact solution(s)
 void exact_solution(const Vector<double>& X, Vector<double>& exact_w)const
 {
  // Convenient definitions
  double x=X[0], y=X[1];
  // u_z and derivatives
  exact_w.resize(8);
  exact_w[0] = pow(x*x+y*y-1,2);
  exact_w[1] = 4*x*(x*x+y*y-1);
  exact_w[2] = 4*y*(x*x+y*y-1);
  exact_w[3] = 4*(3*x*x+y*y-1);
  exact_w[4] = 8*x*y;
  exact_w[5] = 4*(3*y*y+x*x-1);
  // u_x
  exact_w[6] = x;
  // u_y
  exact_w[7] = y;
 }
};

///\short Manufactured solution that extends the linearised limit (i.e. decoupled
// membrane forced by linear bending) to the whole range of p_mag
class ManufacturedSolutionDecoupledExtension : protected FoepplVonKarmanExactSolution
{
 public:
 /// Constructor
 ManufacturedSolutionDecoupledExtension() : FoepplVonKarmanExactSolution(){};
 /// Constructor with nu_pt
 ManufacturedSolutionDecoupledExtension(double* poisson_pt, 
  double* pressure_pt) 
  : FoepplVonKarmanExactSolution()
  {
   nu_pt() = poisson_pt;
   pressure_mag_pt() = pressure_pt;
  };

 // Assigns the value of pressure depending on the position (x,y)
 void get_pressure(const Vector<double>& X, double& pressure) const
 {
  // Convenient definitions
  double x=X[0], y=X[1], p = *pressure_mag_pt();
  // Return the pressure
  pressure = p + (-9*pow(-1 + pow(nu(),2),3)*pow(p,3)*(-1 + pow(x,2) + pow(y,2))*
   (-3 + 15*pow(y,2) + 5*(pow(x,6) + pow(y,4)*(-3 + pow(y,2)) + 3*pow(x,4)*(-1 + pow(y,2)) + 
        3*pow(x,2)*pow(-1 + pow(y,2),2))))/1024.;
  }

 // Assigns the value of in plane forcing depending on the position (x,y)
 void get_in_plane_force(const Vector<double>& X, Vector<double>& grad) const
  {
  grad.resize(2);
  grad[0]=0.0;
  grad[1]=0.0;
  }  
 
 void exact_radial_w_solution(const Vector<double>& X, 
  Vector<double>& exact_w)const
 {
  // Convenient definitions
  double x=X[0], y=X[1];
  // Exact Solution in limit p -> 0
  // We extend this by modifying loading
  double p = (12.*(1-nu()*nu()))*(*pressure_mag_pt());
  exact_w.resize(6);
  exact_w[0] = p*pow(x*x+y*y-1,2)/64.;
  exact_w[1] = p*4*sqrt(x*x+y*y)*(x*x+y*y-1)/64.;
  exact_w[2] = 0.0;
  exact_w[3] = p*4*(3*(x*x+y*y)-1)/64.;
  exact_w[4] = 0.0;
  exact_w[5] = 0.0;
 }
 
 // Get the exact solution(s)
 void exact_solution(const Vector<double>& X, Vector<double>& exact_w) const
 {
  // Convenient definitions
  double x=X[0], y=X[1];
  double p = *pressure_mag_pt(), r2 = x*x+y*y;
  exact_w.resize(8);
  // Out of plane disp.
  exact_w[0] = 12.*(1-nu()*nu())*p*pow(r2-1,2)/64.;
  exact_w[1] = 12.*(1-nu()*nu())*p*4*x*(r2-1)/64.;
  exact_w[2] = 12.*(1-nu()*nu())*p*4*y*(r2-1)/64.;
  exact_w[3] = 12.*(1-nu()*nu())*p*4*(3*x*x+y*y-1)/64.;
  exact_w[4] = 12.*(1-nu()*nu())*p*8*x*y/64.;
  exact_w[5] = 12.*(1-nu()*nu())*p*4*(3*y*y+x*x-1)/64.;
  // In plane disp.  
  exact_w[6] = (3*x*pow(-1 + pow(nu(),2),2)*pow(p,2)*(3 - 18*r2 + 20*pow(r2,2) - 
       7*pow(r2,3) + nu()*(-3 + 6*r2 - 4*pow(r2,2) + pow(r2,3))))/512.;
  
  exact_w[7] = (3*y*pow(-1 + pow(nu(),2),2)*pow(p,2)*(3 - 18*r2 + 20*pow(r2,2) - 
       7*pow(r2,3) + nu()*(-3 + 6*r2 - 4*pow(r2,2) + pow(r2,3))))/512.;
 }
};
}
