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
  /// Rectangle edge length
  double L = 1.0;

  /// The plate thickness
  double Thickness = 0.01;

  /// Poisson ratio
  double Nu = 0.5;

  /// Rubbish name for length scale
  double Eta_u = 1.0;

  /// Rubbish name for stress scale
  double Eta_sigma = 1.0;

  /// Magnitude of pressure
  double P_mag = 0.0;

  /// Magnitude of shear stress
  double T_mag = 0.0;

  /// Element size
  double Element_area = 0.01;

  /// Output directory
  std::string Output_dir = "RESLT";

  /// Pressure depending on the position (x,y)
  void
    get_pressure(const Vector<double>& x,
		 const Vector<double>& u,
		 const DenseMatrix<double>& grad_u,
		 const Vector<double>& n,
		 Vector<double>& pressure)
  {
    pressure.resize(3);
    pressure[0] = 0.0;
    pressure[1] = 0.0;
    pressure[2] = P_mag;
  }



  //----------------------------------------------------------------------
  // Mooney-Rivlin stress
  //----------------------------------------------------------------------
  // According to Li & Healey (2016)
  //   C1 + C2 = E / 6.0
  // Therefore, non-dimensionalising by E, we get that
  //   C2 = 1.0 / 6.0 - C1

  /// First dimensionless Mooney-Rivlin constant
  double C1 = 1.0 / 6.6;

  // Calculate in stress to ensure both don't fall out of sync
  // /// Second dimensionless Mooney-Rivlin constant
  // double C2 = 1.0 / 6.0 - C1;

  /// Mooney Rivlin stress function
  void mooney_rivlin_stress(const Vector<double>& x,
			    const Vector<double>& u,
			    const DenseMatrix<double>& e,
			    const DenseMatrix<double>& g,
			    DenseMatrix<double>& stress)
  {
    // Constants
    const double c1 = Parameters::C1;
    const double c2 = 1.0 / 6.0 - c1;

    // Matrix of cofactors of strain tensor
    DenseMatrix<double> cof_e(2, 2);
    cof_e(0, 0) = e(1, 1);
    cof_e(1, 1) = e(0, 0);
    cof_e(0, 1) = -e(0, 1);
    cof_e(1, 0) = -e(1, 0);
    // Matrix of cofactors of metric tensor
    DenseMatrix<double> cof_g(2, 2);
    cof_g(0, 0) = g(1, 1);
    cof_g(1, 1) = g(0, 0);
    cof_g(0, 1) = -g(0, 1);
    cof_g(1, 0) = -g(1, 0);

    // Determinants
    const double det_e = e(0, 0) * e(1, 1) - e(0, 1) * e(1, 0);
    const double det_g = g(0, 0) * g(1, 1) - g(0, 1) * g(1, 0);
    // Traces
    const double tr_e = e(0, 0) + e(1, 1);
    // const double tr_g = g(0,0)+g(1,1);
    // NB det(g) = 4 det(e) + 2 Tr(e) +1
    // Determinant of g squared minus one
    const double det_g_2_m1 =
      (4 * det_e + 2 * tr_e) * (4 * det_e + 2 * tr_e + 2);

    // Now fill in the stress
    // Loop over indices
    DenseMatrix<double> i2(2, 2, 0.0);
    i2(0, 0) = 1.0;
    i2(1, 1) = 1.0;

    // Now Fill in the Stress
    for (unsigned alpha = 0; alpha < 2; ++alpha)
    {
      for (unsigned beta = 0; beta < 2; ++beta)
      {
        // 2nd Piola Kirchhoff (Membrane) Stress for Mooney Rivlin model
        //  stress(alpha,beta)=2*(c1+ c2 / det_g) * kronecker(alpha,beta)
        //       + 2*((- c1 - tr_g *c2) / pow(det_g,2) + c2)*cof_g(alpha,beta);
        // For 2D:
        // Cof g = I + 2 Cof e
        // tr(g) = 2 + 2 tr(e)
        // Det(g) = 4 Det(e) + 2 Tr(e) + 1
        // 2nd Piola Kirchhoff (Membrane) Stress for Mooney Rivlin model
        // ( Modified so that all terms are order epsilon if possible *)
        stress(alpha, beta) =
	  2 * (c2 * (-4 * det_e - 2 * tr_e) / det_g) * i2(alpha, beta) -
	  4 * ((c1 + 2 * c2 + 2 * tr_e * c2) / pow(det_g, 2) - c2) *
	    cof_e(alpha, beta) -
          2 *
            ((c1 + 2 * c2 + 2 * tr_e * c2) * (-det_g_2_m1) / pow(det_g, 2) +
             2 * tr_e * c2) *
            i2(alpha, beta);
      }
    }
  }

  /// Mooney Rivlin stiffness tensor (only fills in dstrain, not du)
  void d_mooney_rivlin_stress_d_strain(const Vector<double>& x,
				       const Vector<double>& u,
				       const DenseMatrix<double>& strain,
				       const DenseMatrix<double>& g,
                                       RankThreeTensor<double>& d_stress_du,
                                       RankFourTensor<double>& d_stress_dstrain)
  {
    // Constants
    const double c1 = Parameters::C1;
    const double c2 = 1.0 / 6.0 - c1;

    // Matrix of cofactors of metric tensor
    DenseMatrix<double> cof_g(2, 2);
    cof_g(0, 0) = g(1, 1);
    cof_g(1, 1) = g(0, 0);
    cof_g(0, 1) = -g(0, 1);
    cof_g(1, 0) = -g(1, 0);

    // Fill in determinants
    const double det_g = g(0, 0) * g(1, 1) - g(0, 1) * g(1, 0);
    const double tr_g = g(0, 0) + g(1, 1);

    // Identity matrix
    DenseMatrix<double> i2(2, 2, 0.0);
    i2(0, 0) = 1.0;
    i2(1, 1) = 1.0;

    // Now Fill in the Stress
    for (unsigned alpha = 0; alpha < 2; ++alpha)
    {
      for (unsigned beta = 0; beta < 2; ++beta)
      {
        // 2nd Piola Kirchhoff (Membrane) Stress for Mooney Rivlin model
        // stress(alpha,beta)=2*(c1+ c2 / det_g) * i2(alpha,beta)
        //      + 2*((- c1 - tr_g *c2) / pow(det_g,2) + c2)*cof_g(alpha,beta);
        for (unsigned gamma = 0; gamma < 2; ++gamma)
        {
          for (unsigned delta = 0; delta < 2; ++delta)
          {
            // 2nd Piola Kirchhoff (Membrane) Stress for Mooney Rivlin model
            // d (cof_g) dg
            d_stress_dstrain(alpha, beta, gamma, delta) =
              2 * (-2 * (c2 / pow(det_g, 2)) * i2(alpha, beta) *
                     cof_g(gamma, delta) +
                   2 * (c2 - (c1 + tr_g * c2) / pow(det_g, 2)) *
                     // dcof(g)/dg = (cof(g) \otimes cof(g) - cof(g) . dg / dg .
                     // cof(g))/det(g)
                     (cof_g(alpha, beta) * cof_g(gamma, delta) -
                      cof_g(alpha, gamma) * cof_g(delta, beta)) /
                     det_g +
                   4 * ((c1 + tr_g * c2) / pow(det_g, 3)) * cof_g(alpha, beta) *
                     cof_g(gamma, delta) -
                   2 * (c2 / pow(det_g, 2)) * cof_g(alpha, beta) *
                     i2(gamma, delta));
          }
        }
      }
    }
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
    oomph_info << "Solving for T = " << Parameters::T_mag << std::endl;
    oomph_info << "         step = " << Doc_info.number() << std::endl;
    oomph_info << "-------------------------------------------------------"
               << std::endl;
  }


  /// Make the problem linear (biharmonic) by pinning all in-plane dofs and
  /// setting eta=0
  void make_linear()
  {
    // Remove stretching coupling
      Parameters::Eta_u = 0.0;

    // Pin all in-plane displacements
    unsigned n_node = Bulk_mesh_pt->nnode();
    for (unsigned i_node = 0; i_node < n_node; i_node++)
    {
      Bulk_mesh_pt->node_pt(i_node)->pin(0);
      Bulk_mesh_pt->node_pt(i_node)->set_value(0, 0.0);
      Bulk_mesh_pt->node_pt(i_node)->pin(1);
      Bulk_mesh_pt->node_pt(i_node)->set_value(1, 0.0);
    }
  } // End make_linear()


  /// Pin all out of plane dofs
  void pin_out_of_plane()
  {
    unsigned n_node = Bulk_mesh_pt->nnode();
    for(unsigned i_node = 0; i_node < n_node; i_node++)
    {
      for(unsigned j_type = 0; j_type < 6; j_type++)
      {
        Bulk_mesh_pt->node_pt(i_node)->pin(13+j_type);
      }
    }
  }


  /// Doc the solution
  void doc_solution(const std::string& comment = "");

  /// Overloaded version of the problem's access function to
  /// the mesh. Recasts the pointer to the base Mesh object to
  /// the actual mesh type.
  TriangleMesh<ELEMENT>* mesh_pt()
  {
    return dynamic_cast<TriangleMesh<ELEMENT>*>(Problem::mesh_pt());
  }

  /// Bulk mesh pointer
  TriangleMesh<ELEMENT>* bulk_mesh_pt()
  {
    return Bulk_mesh_pt;
  }

private:
  /// Setup and build the mesh
  void build_mesh();

  /// Helper function to (re-)set boundary condition
  /// and complete the build of all elements
  void complete_problem_setup();

  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();

  /// Helper function to output the stress
  void output_stress(std::ostream& outfile,
		     const unsigned& npts);

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

  /// Pointer to element that contains the central point
  GeomObject* Central_element_geom_obj_pt;

  /// Local coordinate in element pointed to by Central_element_pt
  /// that contains central point
  Vector<double> Central_point_local_coord;

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
             << "L            " << Parameters::L << std::endl
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


  double length = Parameters::L;


  // Declare the vertices...
  Vector<double> vertex0(2, 0.0), vertex1(2, 0.0), vertex2(2, 0.0),
    vertex3(2, 0.0);

  // ...and place them at the corners of the rectangle
  vertex0[0] = -0.5 * length;
  vertex0[1] = -0.5;
  vertex1[0] = 0.5 * length;
  vertex1[1] = -0.5;
  vertex2[0] = 0.5 * length;
  vertex2[1] = 0.5;
  vertex3[0] = -0.5 * length;
  vertex3[1] = 0.5;

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

    // Set the pressure & temperature function pointers and the physical
    // constants
    el_pt->pressure_fct_pt() = &Parameters::get_pressure;

    // Assign the parameter pointers for the element
    el_pt->thickness_pt() = &Parameters::Thickness;
    el_pt->nu_pt() = &Parameters::Nu;
    el_pt->eta_u_pt() = &Parameters::Eta_u;
    el_pt->eta_sigma_pt() = &Parameters::Eta_sigma;
    if(!CommandLineArgs::command_line_flag_has_been_set("--use_linear_stress"))
    {
      el_pt->stress_fct_pt() = &Parameters::mooney_rivlin_stress;
      el_pt->d_stress_fct_pt() = &Parameters::d_mooney_rivlin_stress_d_strain;
    }
    if(CommandLineArgs::command_line_flag_has_been_set("--use_fd_jacobian"))
    {
      el_pt->enable_finite_difference_jacobian();
    }
  }

  // Set the boundary conditions
  apply_boundary_conditions();

  // Find element (and local coordinate within it), that contains the
  // central point
  Vector<double> origin(2, 0.0);

  // Create a MeshAsGeomObject from the Mesh:
  MeshAsGeomObject mesh_as_geom_object(Bulk_mesh_pt);

  // Make space for local coordinate in element pointed to by Central_element_pt
  // that contains central point
  Central_point_local_coord.resize(2);

  // Find it!
  mesh_as_geom_object.locate_zeta(
    origin, Central_element_geom_obj_pt, Central_point_local_coord);

  // Did that work out?
  if (Central_element_geom_obj_pt == 0)
  {
    throw OomphLibError("couldn't find central element",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }

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

  // [zdec] NONPHYSICAL - used for debugging
  const Vector<unsigned> uber_clamped{0, 1, 2, 3, 4, 5};

  //------------------------------------------------------------------
  //------------------------------------------------------------------

  unsigned n_field = 3;

  // Vector containers to store which boundary conditions we are applying to
  // each edge. (outlined above)
  Vector<Vector<Vector<unsigned>>> pinned_u_dofs(4, Vector<Vector<unsigned>>(3));

  // We are doing uniaxial stress.
  // Totally left and right are free to slide in y, bottom can slide in x,
  // top can slide in x and y to allow for Poisson effect
  for(unsigned i_field = 0; i_field < n_field; i_field++)
  {
    // Clamp all but u_x
    if(i_field != 0)
    {
      pinned_u_dofs[0][i_field] = pinned_edge_yn_dof;
    }
    // Pin all but u_y
    if(i_field != 1)
    {
      pinned_u_dofs[1][i_field] = pinned_edge_xn_dof;
      pinned_u_dofs[3][i_field] = pinned_edge_xn_dof;
    }
    // Pin only u_z
    if(i_field == 2)
    {
      pinned_u_dofs[2][i_field] = pinned_edge_yn_dof;
    }
  }

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


//========================================================
/// Output helper function for stresses
///
/// Loop over all elements and plot (i.e. execute
/// the element's own output() function)
//========================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::output_stress(std::ostream& outfile,
						   const unsigned& npts)
{
  // Loop over the elements and call their output functions
  // Assign Element_pt_range
  unsigned long Element_pt_range = Bulk_mesh_pt->nelement();
  for (unsigned long e = 0; e < Element_pt_range; e++)
  {
    // Try to cast to ELEMENT
    ELEMENT* el_pt =
      dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
    if (el_pt == 0)
    {
      oomph_info << "Can't execute output_stress(...) for non"
		 << " KoiterSteigmanEquation based elements"
                 << std::endl;
    }
    else
    {
      el_pt->output_stress(outfile, npts);
    }
  }
}




//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::doc_solution(const std::string& comment)
{
  ofstream some_file;
  char filename[100];


  // Number of plot points for fine (full) output. Lots of plot points so
  // we see the goodness of the high-order Hermite/Bell
  // interpolation.
  unsigned npts = 10;
  sprintf(filename,
          "%s/soln_%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  output_stress(some_file, npts);
  some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" << comment << "\"\n";
  some_file.close();

  // Write the load magnitude and centrepoint displacements to trace file
  Vector<Vector<double>> u_centre(3, Vector<double>(6,0.0));

  dynamic_cast<ELEMENT*>(Central_element_geom_obj_pt)
    ->interpolated_koiter_steigmann_disp(Central_point_local_coord, u_centre);

  Trace_file << Parameters::P_mag << " "
	     << Parameters::T_mag << " "
             << u_centre[0][0] << " "
	     << u_centre[1][0] << " "
	     << u_centre[2][0] << " "
             << Doc_info.number() << endl;

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

  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  // Directory for solution
  CommandLineArgs::specify_command_line_flag("--dir", &Parameters::Output_dir);

  // Use linear stress rather than MR
  CommandLineArgs::specify_command_line_flag("--use_linear_stress");

  // Use the finite difference jacobian
  CommandLineArgs::specify_command_line_flag("--use_fd_jacobian");

  // Parse command line
  CommandLineArgs::parse_and_assign();


  // Create the problem, using KS elements derived from TElement<2,2>
  // elements.
  UnstructuredKSProblem<KoiterSteigmannC1CurvableBellElement> problem;

  problem.max_residuals() = 1.0e3;
  problem.max_newton_iterations() = 30;

  // // Document the initial state
  // problem.doc_solution();


  // Store a pointer to the mesh
  Mesh* mesh_pt = problem.bulk_mesh_pt();
  // Store the number of boundary nodes on the right boundary
  unsigned n_bnode = mesh_pt->nboundary_node(1);
  // Store a sample element pointer
  KoiterSteigmannC1CurvableBellElement* el_pt =
    dynamic_cast<KoiterSteigmannC1CurvableBellElement*>(mesh_pt->element_pt(0));

  // Displacement of the right boundary
  double u_imposed = 0.0;
  // Maximum displacement of the right boundary
  double u_max = 0.5;
  // Minimum displacement of the right boundary
  double u_min = -0.1;
  // Increment of the displacement each solve
  double u_inc = 0.05;

  // Pin all out of plane dofs
  problem.pin_out_of_plane();
  
  // Get to maximum compression
  for(unsigned i = 0; u_imposed > u_min; i++)
  {
    // Increment the stretch
    u_imposed -= u_inc/10.0;
    // Stretch the sheet by increasing the displacement of the right hand side
    for(unsigned j_node = 0; j_node < n_bnode; j_node++)
    {
      Node* node_pt = mesh_pt->boundary_node_pt(1, j_node);
      node_pt->set_value(0, u_imposed);
    }
    // Solve the system
    problem.newton_solve();
  }

  // Document the current solution
  problem.doc_solution();

  // Begin stretching loop
  for(unsigned i = 0; u_imposed < u_max; i++)
  {
    // Increment the stretch
    u_imposed += u_inc;
    // Stretch the sheet by increasing the displacement of the right hand side
    for(unsigned j_node = 0; j_node < n_bnode; j_node++)
    {
      Node* node_pt = mesh_pt->boundary_node_pt(1, j_node);
      node_pt->set_value(0, u_imposed);
    }
    // Solve the system
    problem.newton_solve();
    // Document the current solution
    problem.doc_solution();
  }

} // End of main
