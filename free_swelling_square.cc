// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
// LIC//
// LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
#include <fenv.h>

// Generic routines
#include "generic.h"

// The equations
#include "c1_koiter_steigmann.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;


//===========================================================================
/// Namespace for problem parameters
//===========================================================================
namespace Parameters
{
  /// The plate thickness
  double Thickness = 0.01;

  /// Poisson ratio
  double Nu = 0.5;

  /// Coefficient of damping
  double Mu = 1.0e-1;

  /// What is this?
  double Eta_u = 1.0; // 12.0 * (1.0 - Nu * Nu) / (Thickness * Thickness);

  /// What is thissss?
  double Eta_sigma = 1.0;

  /// Pressure scale
  double P_scale = Thickness*Thickness / (12.0*(1.0-Nu*Nu));

  /// Magnitude of swelling
  double C_mag = 0.0;

  /// Magnitude of pressure
  double P_mag = 0.0;

  /// Element size
  double Element_area = 0.2;


  /// Pressure depending on the position (x,y) and deformation of the sheet
  void get_pressure(const Vector<double>& x,
		    const Vector<double>& u,
		    const DenseMatrix<double>& grad_u,
		    const Vector<double>& n,
		    Vector<double>& pressure)
  {
    // Metric tensor of deformed surface
    DenseMatrix<double> G(2,2,0.0);
    for(unsigned alpha = 0; alpha < 2; alpha++)
    {
      for (unsigned beta = 0; beta < 2; beta++)
      {
	for (unsigned i = 0; i < 3; i++)
	{
	  G(alpha,beta) += grad_u(i,0)*grad_u(i,1);
	}
      }
    }
    // Find the pressure per undeformed area in terms of the pressure per
    // deformed area
    double p = sqrt(G(0,0)*G(1,1) - G(1,0)*G(0,1)) * P_mag;
    // Assign pressure
    pressure.resize(3);
    pressure[0] = p * n[0];
    pressure[1] = p * n[1];
    pressure[2] = p * n[2];
  }


  /// Swelling induced prestrain
  void get_swelling_prestrain(const Vector<double>& x,
			      DenseMatrix<double>& prestrain)
  {
    double isostrain = -(C_mag + 0.5*C_mag*C_mag);
    prestrain(0,0) = isostrain;
    prestrain(0,1) = 0.0;
    prestrain(1,0) = 0.0;
    prestrain(1,1) = isostrain;
  }



  //-------- Boundary conditions -----------------------------------------------
  /// Function to specify boundary conditions (here all homogeneous) as a
  /// function of both coordinates (This is convenient for this problem; other
  /// interfaces that specify boundary conditions in terms of boundary
  /// coordinate exist).
  void get_null_fct(const Vector<double>& x, double& value)
  {
    value = 0.0;
  }

} // namespace Parameters


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredKSProblem : public virtual Problem
{
public:
  /// Constructor
  UnstructuredKSProblem();

  /// Destructor
  ~UnstructuredKSProblem()
  {
    // Close the trace file as we are done with the problem
    Trace_file.close();

    // Clean up memory
    delete Bulk_mesh_pt;
    delete Boundary_pt;
    delete Boundary0_pt;
    delete Boundary1_pt;
    delete Boundary2_pt;
    delete Boundary3_pt;
  };


  /// Actions to complete before each Newton solve
  void actions_before_newton_solve()
  {
    // Print a solve header to output
    oomph_info << "-------------------------------------------------------"
               << std::endl;
    oomph_info << "Solving for P = " << Parameters::P_mag << std::endl;
    oomph_info << "         step = " << Doc_info.number() << std::endl;
    oomph_info << "-------------------------------------------------------"
               << std::endl;
  }


  // /// Make the problem linear (biharmonic) by pinning all in-plane dofs and
  // /// setting eta=0
  // void make_linear()
  // {
  //   // Remove stretching coupling
  //     Parameters::Eta_u = 0.0;

  //   // Pin all in-plane displacements
  //   unsigned n_node = Bulk_mesh_pt->nnode();
  //   for (unsigned i_node = 0; i_node < n_node; i_node++)
  //   {
  //     Bulk_mesh_pt->node_pt(i_node)->pin(0);
  //     Bulk_mesh_pt->node_pt(i_node)->set_value(0, 0.0);
  //     Bulk_mesh_pt->node_pt(i_node)->pin(1);
  //     Bulk_mesh_pt->node_pt(i_node)->set_value(1, 0.0);
  //   }
  // } // End make_linear()

  /// Doc the solution
  void doc_solution(const std::string& comment = "");

  /// Overloaded version of the problem's access function to
  /// the mesh. Recasts the pointer to the base Mesh object to
  /// the actual mesh type.
  TriangleMesh<ELEMENT>* mesh_pt()
  {
    return dynamic_cast<TriangleMesh<ELEMENT>*>(Problem::mesh_pt());
  }

private:
  /// Setup and build the mesh
  void build_mesh();

  /// Helper function to (re-)set boundary condition
  /// and complete the build of all elements
  void complete_problem_setup();

  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();

  /// Triangle mesh parameters
  TriangleMeshParameters* Triangle_mesh_parameters_pt;

  /// Closed outer boundary
  TriangleMeshClosedCurve* Boundary_pt;

  /// Polyline defining boundary 0
  TriangleMeshPolyLine* Boundary0_pt;

  /// Polyline defining boundary 1
  TriangleMeshPolyLine* Boundary1_pt;

  /// Polyline defining boundary 2
  TriangleMeshPolyLine* Boundary2_pt;

  /// Polyline defining boundary 3
  TriangleMeshPolyLine* Boundary3_pt;

  /// Doc info object for labeling output
  DocInfo Doc_info;

  /// Trace file to document norm of solution
  ofstream Trace_file;

  /// Pointer to "bulk" mesh
  TriangleMesh<ELEMENT>* Bulk_mesh_pt;
}; // end_of_problem_class


//==start_of_problem_constructor===========================================
/// Constructor
//=========================================================================
template<class ELEMENT>
UnstructuredKSProblem<ELEMENT>::UnstructuredKSProblem()
{
  // Set output directory
  Doc_info.set_directory("RESLT");

  // Step number
  Doc_info.number() = 0;

  // Build the mesh
  build_mesh();

  // Complete problem setup
  complete_problem_setup();


  // Output parameters
  oomph_info << "Problem parameters:\n"
             << "thickness    " << Parameters::Thickness << std::endl
             << "nu           " << Parameters::Nu << std::endl
             << "eta_u        " << Parameters::Eta_u << std::endl
             << "eta_sigma    " << Parameters::Eta_sigma << std::endl
             << "Element area " << Parameters::Element_area << std::endl;


  // Open trace file
  char filename[100];
  strcpy(filename, (Doc_info.directory() + "/trace.dat").c_str());
  Trace_file.open(filename);


  // Assign equation numbers
  oomph_info << "Number of equations: " << assign_eqn_numbers() << '\n';

} // end Constructor


//==start_of_build_mesh====================================================
/// Build the rectangular mesh
//=========================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::build_mesh()
{
  /*================================
    Rectangular mesh boundary

    V3       E2       V2
    O-----------------O     ^
    |                 |     |
    |        (0,0)    |     |
 E3 |        x        | E1  1
    |                 |     |
    |                 |     |
    O-----------------O     v
    V0       E0       V1
    <--------L-------->
    ================================*/


  // Declare the vertices...
  Vector<double> vertex0(2, 0.0), vertex1(2, 0.0), vertex2(2, 0.0),
    vertex3(2, 0.0);

  // ...and place them at the corners of the rectangle
  vertex0[0] = 0.0;
  vertex0[1] = 0.0;
  vertex1[0] = 1.0;
  vertex1[1] = 0.0;
  vertex2[0] = 1.0;
  vertex2[1] = 1.0;
  vertex3[0] = 0.0;
  vertex3[1] = 1.0;

  // Declare the edges...
  Vector<Vector<double>> edge0(2, Vector<double>(2, 0.0)),
    edge1(2, Vector<double>(2, 0.0)), edge2(2, Vector<double>(2, 0.0)),
    edge3(2, Vector<double>(2, 0.0));

  // ...and assign their endpoints
  edge0[0] = vertex0;
  edge0[1] = vertex1;
  edge1[0] = vertex1;
  edge1[1] = vertex2;
  edge2[0] = vertex2;
  edge2[1] = vertex3;
  edge3[0] = vertex3;
  edge3[1] = vertex0;

  // Define boundaries from edges
  Boundary0_pt = new TriangleMeshPolyLine(edge0, 0);
  Boundary1_pt = new TriangleMeshPolyLine(edge1, 1);
  Boundary2_pt = new TriangleMeshPolyLine(edge2, 2);
  Boundary3_pt = new TriangleMeshPolyLine(edge3, 3);

  // Create closed outer boundary
  Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);
  boundary_polyline_pt[0] = Boundary0_pt;
  boundary_polyline_pt[1] = Boundary1_pt;
  boundary_polyline_pt[2] = Boundary2_pt;
  boundary_polyline_pt[3] = Boundary3_pt;
  Boundary_pt = new TriangleMeshClosedCurve(boundary_polyline_pt);


  // Define mesh parameters
  TriangleMeshParameters Triangle_mesh_parameters(Boundary_pt);

  // Set the maximum element area
  Triangle_mesh_parameters.element_area() = Parameters::Element_area;

  // Build  bulk mesh
  Bulk_mesh_pt = new TriangleMesh<ELEMENT>(Triangle_mesh_parameters);

  // Add submesh to problem
  add_sub_mesh(Bulk_mesh_pt);

  // Combine submeshes into a single Mesh (bit over the top here; could
  // have assigned bulk mesh to mesh_pt() directly).
  build_global_mesh();

} // end build_mesh



//==start_of_complete======================================================
/// Set boundary conditions and complete the build of
/// all elements
//=========================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::complete_problem_setup()
{
  // Complete the build of all elements so they are fully functional
  unsigned n_element = Bulk_mesh_pt->nelement();
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // Set the pressure & prestrain function pointers
    el_pt->pressure_fct_pt() = &Parameters::get_pressure;
    el_pt->prestrain_fct_pt() = &Parameters::get_swelling_prestrain;

    // Assign the parameter pointers for the element
    el_pt->thickness_pt() = &Parameters::Thickness;
    el_pt->nu_pt() = &Parameters::Nu;
    el_pt->mu_pt() = &Parameters::Mu;
    el_pt->eta_u_pt() = &Parameters::Eta_u;
    el_pt->eta_sigma_pt() = &Parameters::Eta_sigma;

    // Use the fd jacobian
    // el_pt->enable_finite_difference_jacobian();
  }

  // Set the boundary conditions
  apply_boundary_conditions();

  // Create a MeshAsGeomObject from the Mesh:
  MeshAsGeomObject mesh_as_geom_object(Bulk_mesh_pt);

} // end of complete


//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::apply_boundary_conditions()
{
  // [zdec] TODO REDO COMMENTS
  //------------------------------------------------------------------
  // Boundary conditions for KS elements are complicated and we provide an
  // illustration of how to apply all physically relevant boundary conditions
  // for problems with axis-aligned boundaries here. Other tutorials/driver
  // codes explain what to do about (i) straight boundaries that are not aligned
  // with the coordinate axes and (ii) curvilinear boundaries.
  //
  // FvK elements have two different types of degrees of freedom:
  // (1) The ones associated with the two in-plane displacements, u_x and u_y.
  //     These are interpolated with standard C0 continuous Lagrange
  //     interpolation between all the nodes in the underlying
  //     TElement<2,NNODE_1D>. Each node stores the values
  //     of the two in-plane displacements. We enumerate these dofs
  //     as 0 for the x displacement, and 1 for the y displacement.
  // (2) The dofs associated with the out-of-plane displacement, w.
  //     This is interpolated with Bell (Hermite) interpolants
  //     which involve six different types of degree of freedom,
  //     enumerated from 0 to 5 in the order w, w_x, w_y, w_xx, w_xy,
  //     w_yy. These values are only stored at three vertices of the element.
  //
  // Given that the book-keeping for which node stores which type of degree of
  // freedom is complicated, we let the element do the relevant assigments for
  // us. For this purpose the FvK elements provide two member functions:
  //
  //     fix_in_plane_displacement_dof(idof, b, fct_pt)
  //
  // assigns boundary conditions for the in-plane displacements via
  // the specification of
  //
  // idof   : the enumeration of the dof in the scheme listed above,
  //          so idof can take values 0 or 1.
  // b      : the mesh boundary along which the boundary condition is
  //          to be applied
  // fct_pt : a function pointer to a global function with arguments
  //          (const Vector<double> x, double& value) which computes
  //          the value for the relevant in-plane displacement as a
  //          function of the coordinate, x, a 2D vector. Note that,
  //          since the boundary is assumed to be aligned with the
  //          coordinate axes, one of the two coordinates will be
  //          irrelevant, but it is still passed to the function.
  //
  // So, if the function fix_in_plane_displacement_dof(idof, b, fct_pt) is
  // called with idof=1 and b=3, say, the y-in-plane displacement is pinned for
  // all the element's nodes (if any) that are located on mesh boundary b. The
  // value of the y-in-plane displacement is set to whatever the function
  // pointed to by fct_pt computes when evaluated at the nodal coordinate.
  //
  // Similarly,
  //
  //      fix_out_of_plane_displacement_dof(idof, b, fct_pt);
  //
  // hierher complete once Aidan has signed off the explanation above.
  // [zdec] "This is all good." -- Aidan
  //
  // Using the conventions introduced above, the following vectors identify the
  // in-plane and out-of-plane degrees of freedom to be pinned for various
  // physically meaningful boundary conditions:


  // Possible boundary conditions for out-of-plane displacements: Given that the
  // out-of-plane displacements feature in the fourth-order biharmonic operator,
  // we can apply boundary conditions on w and dw/dn, where n is the coordinate
  // direction normal to the (assumed to be axis aligned!) boundary. However if
  // w is given along the entire boundary (parametrised by the tangential
  // coordinate, t) we also know what dw/dt and d^2w/dt^2 are. Likewise if dw/dn
  // is known along the whole boundary we also know d^2w/dndt. In the various
  // cases below we identify physical scenarios of a pinned edge (w given, dw/dn
  // left free); a vertically sliding edge (w left free; dw/dn given) and fully
  // clamped (w and dw/dn given). Together with the two possible orientations of
  // the axis aligned boundaries (x aligned or y aligned) we get six different
  // cases:

  // Out-of-plane dofs:
  //-------------------
  // |   0   |   1   |   2   |   3   |   4   |   5   |
  // |  u_i  | u_i_x | u_i_y | u_i_xx| u_i_xy| u_i_yy|

  // Case: The plate is pinned (w given, dw/dn left free) along a boundary where
  // the outer unit normal points in the postive or negative x direction, so x
  // is constant and y varies along the boundary.  We therefore have to pin (and
  // assign values for) w, dw/dy and d^2w/dy^2
  const Vector<unsigned> pinned_edge_xn_dof{0, 2, 5};

  // Case: The plate is pinned (w given, dw/dn left free) along a boundary where
  // the outer unit normal points in the postive or negative y direction, so y
  // is constant and x varies along the boundary.  We therefore have to pin (and
  // assign values for) w, dw/dx and d^2w/dx^2
  const Vector<unsigned> pinned_edge_yn_dof{0, 1, 3};

  // Case: The plate is sliding (w left free, dw/dn given) along a boundary
  // where the outer unit normal points in the postive or negative x direction,
  // so x is constant and y varies along the boundary.  We therefore have to pin
  // (and assign values for) dw/dx and d^2w/dxdy
  const Vector<unsigned> sliding_clamp_xn_dof{1, 4};

  // Case: The plate is sliding (w left free, dw/dn given) along a boundary
  // where the outer unit normal points in the postive or negative y direction,
  // so y is constant and x varies along the boundary.  We therefore have to pin
  // (and assign values for) dw/dy and d^2w/dxdy
  const Vector<unsigned> sliding_clamp_yn_dof{2, 4};

  // Case: The plate is clamped (w given, dw/dn given) along a boundary where
  // the outer unit normal points in the postive or negative x direction, so x
  // is constant and y varies along the boundary.  We therefore have to pin (and
  // assign values for) w, dw/dx, dw/dy, d^2w/dxdy and d^2w/dy^2
  const Vector<unsigned> fully_clamped_xn_dof{0, 1, 2, 4, 5};

  // Case: The plate is clamped (w given, dw/dn given) along a boundary where
  // the outer unit normal points in the postive or negative y direction, so y
  // is constant and x varies along the boundary.  We therefore have to pin (and
  // assign values for) w, dw/dx, dw/dy, d^2w/dx^2 and d^2w/dxdy
  const Vector<unsigned> fully_clamped_yn_dof{0, 1, 2, 3, 4};

  // Free edge has no constraints
  const Vector<unsigned> free{};

  // [zdec] NONPHYSICAL - used for debugging
  const Vector<unsigned> uber_clamped{0, 1, 2, 3, 4, 5};

  //------------------------------------------------------------------
  //------------------------------------------------------------------

  unsigned n_field = 3;

  // Vector containers to store which boundary conditions we are applying to
  // each edge. (outlined above)
  Vector<Vector<Vector<unsigned>>> pinned_u_dofs(4, Vector<Vector<unsigned>>(3,Vector<unsigned>(free)));

  // Constrain the bottom against the x-axis
  pinned_u_dofs[0][1] = pinned_edge_yn_dof;

  // Constrain the left against the y-axis
  pinned_u_dofs[3][0] = pinned_edge_xn_dof;

  // Clamp all four edges flat
  pinned_u_dofs[0][2] = fully_clamped_yn_dof;
  pinned_u_dofs[1][2] = fully_clamped_xn_dof;
  pinned_u_dofs[2][2] = fully_clamped_yn_dof;
  pinned_u_dofs[3][2] = fully_clamped_xn_dof;


  // Loop over all the boundaries in our bulk mesh
  unsigned n_bound = Bulk_mesh_pt->nboundary();
  for (unsigned b = 0; b < n_bound; b++)
  {
    // Number of elements on b
    const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);

    // Loop over the elements on boundary b
    for (unsigned e = 0; e < nb_element; e++)
    {
      // Get pointer to bulk element adjacent to b
      ELEMENT* el_pt =
	dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b, e));

      // Loop over each of the three fields
      for(unsigned i_field = 0; i_field < n_field; i_field++)
      {
	// Number of dofs we are pinning on boundary b
	const unsigned n_pinned_u_dofs = pinned_u_dofs[b][i_field].size();

	// Pin in-plane dofs (enumerated as explained above) for all nodes on
	// boundary b. Here we're applying homogeneous BCs so all pinned values
	// are simply set to zero using Parameters::get_null_fct.
	for (unsigned k = 0; k < n_pinned_u_dofs; k++)
	{
	  unsigned k_type = pinned_u_dofs[b][i_field][k];
	  std::cout << "On boundary " << b
		    << " pinning deflection " << i_field
		    << " type " << k_type << std::endl;
	  el_pt->set_boundary_condition(
	    i_field, k_type, b, Parameters::get_null_fct);
	} // end for loop over types that need to be pinned [k]
      } // end for loop over displacements [i_field]
    } // end for loop over elements on b [e]
  } // end for loop over boundaries [b]

} // end set bc


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::doc_solution(const std::string& comment)
{
  ofstream some_file;
  char filename[100];


  // Number of plot points for coarse output (just showing the
  // element outline)
  unsigned npts = 2;
  sprintf(filename,
          "%s/coarse_soln_%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file, npts);
  some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" << comment << "\"\n";
  some_file.close();

  // Number of plot points for fine (full) output. Lots of plot points so
  // we see the goodness of the high-order Hermite/Bell
  // interpolation.
  npts = 10;
  sprintf(filename,
          "%s/soln_%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file, npts);
  some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" << comment << "\"\n";
  some_file.close();

  // Bump
  Doc_info.number()++;

} // end of doc


//=======start_of_main========================================
/// Driver code for demo of unstructured C1 Koiter-Steigmann
/// elements
//============================================================
int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

  // Create the problem, using KS elements derived from TElement<2,2>
  // elements.
  UnstructuredKSProblem<KoiterSteigmannC1CurvableBellElement> problem;

  problem.max_residuals() = 1.0e3;
  problem.max_newton_iterations() = 30;

  // Set pressure
  Parameters::C_mag = 1.0e-1;
  // Set the Poisson ratio
  Parameters::Nu = 0.5;

  // Document the initial state
  problem.doc_solution();
  // Solve the system
  problem.newton_solve();
  // Document the current solution
  problem.doc_solution();

} // End of main
