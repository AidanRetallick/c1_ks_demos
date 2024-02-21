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
#include <fstream>
#include <ios>

//Generic routines
#include "generic.h"

// The equations
#include "c1_koiter_steigmann.h"

// The mesh
#include "meshes/triangle_mesh.h"


using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

//                      OUTLINE OF PROBLEM CONSTRUCTION
// The basic constuction is much the same as the usual order of things in a
// problem. Underneath is the order of actions (with stars next to non actions
// that are unique to these types of problems).
// 1.  Setup mesh parameters
// 2.  Build the mesh

// hierher: combine 3 and 4?

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

//                            REQUIRED DEFINITIONS
// Per Curve Section we will need:
// 1.  A parametric function defining the curve section.
// 2.  The tangential derivative of the parametric function defining
//     the curve section.


// hierher: we should probably insist on this being done in any case. Is it
//          obvious to the user that they have a fifth order boundary
//          representation?

// 3.* (For 5 order boundary representation) The second tangential derivative
//     of the parametric function defining the curve section.
// 4.  A unit normal and tangent to each curve section and corresponding
//     derivatives, to allow the rotation of boundary coordinates.

// hierher: doesn't really make sense and how is it convenient? It's
//          not used here. I think

// It also convenient to define:
// 1.  An inverse function (x,y) -> s (the arc coordinate) to help in setting
//     up the nodal positions in terms of this parametric coordinate.




//========================================================================
/// Namespace for problem parameters
//========================================================================
namespace Parameters
{

  // // Enumeration of cases
  // enum{
  //   Clamped_validation,
  //   Axisymmetric_shear_buckling,
  //   Nonaxisymmetric_shear_buckling
  // };

  /// Which case are we doing
  //  unsigned Problem_case = ; // Nonaxisymmetric_shear_buckling; // Axisymmetric_shear_buckling;

  /// Upper ellipse x span
  double A1 = 0.5;

  /// Upper ellipse y span
  double B1 = 1.0;

  /// Lower ellipse x span
  double A2 = 0.55;

  /// Lower ellipse y span
  double B2 = 0.5;

  /// Char array containing the condition on each boundary. Character array
  /// index corresponds to boundary enumeration and the entry to the contition
  /// type: // [zdec] use this at some point
  ///   - 'c' for clamped
  ///   - 'p' for pinned
  ///   - 's' for sliding
  ///   - 'f' for free
  char Bc_char[2][4] = {"ppp","ppp"};

  /// The plate thickness
  double Thickness = 0.01;

  /// Poisson ratio
  double Nu = 0.5;

  /// What is this?
  double Eta_u = 1.0; // 12.0 * (1.0 - Nu * Nu) / (Thickness * Thickness);

  /// What is thissss?
  double Eta_sigma = 1.0;

  /// Magnitude of pressure
  double P_mag = 0.0;

  /// Element size
  double Element_area = 0.2;

  /// Order of the polynomial interpolation of the boundary
  unsigned Boundary_order = 5;


  //----------------------------------------------------------------------------
  // [zdec] DEPENDENT VARIABLES -- need some sort of update hook

  /// x-component of intersection (positive value)
  double X_intersect = sqrt(A1*A1*A2*A2*(B1*B1-B2*B2)
                            / (A2*A2*B1*B1 - A1*A1*B2*B2));

  /// y-component of boundary intersection
  double Y_intersect = sqrt( B1*B1
                             * ( 1.0 - X_intersect*X_intersect/(A1*A1) ) );

  /// Theta of intersect 1
  // (shift by -pi/2 as the ellipse uses the angle about positive y)
  double Theta1 = atan2(Y_intersect/B1,X_intersect/A1) - Pi/2.0;

  /// Theta of intersect 2
  // (shift by -pi/2 as the ellipse uses the angle about positive y)
  double Theta2 = atan2(Y_intersect/B2,X_intersect/A2) - Pi/2.0;


  // Boundary info
  //                       __
  //                     -    -
  //                   -        -   *Upper ellipse arc*
  //                 /            \
  //               /                \
  //             /                    \
  //           /                        \
  //          /         ________         \
  //         /       --         --        \
  //       /       --              --       \
  //      /      -                    -      \
  //     /    /  *Lower ellipse arc*     \    \
  //    /  /                                \  \
  //   / /                                    \ \
  //   X(Theta2)                                X(Theta1)

  /// Parametric curve for the upper elliptical boundary arc
  // (true->anticlockwise parametrisation)
  CurvilineEllipseTop Upper_parametric_elliptical_curve(A1,B1,false);

  /// Parametric curve for the lower elliptical boundary arc
  // (true->clockwise parametrisation)
  CurvilineEllipseTop Lower_parametric_elliptical_curve(A2,B2,true);

  /// Vector of parametric boundaries
  Vector<CurvilineGeomObject*> Parametric_curve_pt =
  {
    &Upper_parametric_elliptical_curve,
    &Lower_parametric_elliptical_curve
  };

  /// Pressure depending on the position (x,y)
  void get_pressure(const Vector<double>& x, double& pressure)
  {
    pressure = P_mag;
  }

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



  // hierher: kill but check with Aidan first

  // // This metric will flag up any non--axisymmetric parts
  // void axiasymmetry_metric(const Vector<double>& x,
  //                      const Vector<double>& u,
  //                      const Vector<double>& u_exact,
  //                      double& error,
  //                      double& norm)
  // {
  //  // We use the theta derivative of the out of plane deflection
  //  error = pow((-x[1]*u[1] + x[0]*u[2])/sqrt(x[0]*x[0]+x[1]*x[1]),2);
  //  norm  = pow(( x[0]*u[1] + x[1]*u[2])/sqrt(x[0]*x[0]+x[1]*x[1]),2);
  // }

  // Get the null function for applying homogenous BCs
  void get_null_fct(const Vector<double>& X, double& exact_w)
  {
    exact_w = 0.0;
  }

  // Get the unit function for applying homogenous BCs
  void get_unit_fct(const Vector<double>& X, double& exact_w)
  {
    exact_w = 1.0;
  }

}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Problem definition
//====================================================================
using Parameters::Parametric_curve_pt;

template<class ELEMENT>
class UnstructuredKSProblem : public virtual Problem
{

public:

  /// Constructor
  UnstructuredKSProblem(double const& element_area = 0.09);

  /// Destructor
  ~UnstructuredKSProblem()
  {
    // Close trace file
    Trace_file.close();

    // Clean up memory
    delete Bulk_mesh_pt;
    delete Constraint_mesh_pt;
    delete Constraint_mesh_pt;
    delete Outer_boundary_pt;
    delete Outer_curvilinear_boundary_pt[0];
    delete Outer_curvilinear_boundary_pt[1];
  };

  /// Public function to reapply boundary conditions
  void reapply_boundary_conditions()
  {
    apply_boundary_conditions();
  }

  /// Update after solve (empty)
  void actions_after_newton_solve() {}

  /// Update the problem specs before solve: empty
  void actions_before_newton_solve(){}

  /// Doc the solution
  void doc_solution(const std::string& comment="");

  /// Pointer to the bulk mesh
  TriangleMesh<ELEMENT>* bulk_mesh_pt()
  {
    return Bulk_mesh_pt; //dynamic_cast<TriangleMesh<ELEMENT>*> (Problem::mesh_pt());
  }

  /// Pointer to the constraint mesh
  Mesh* constraint_mesh_pt()
  {
    return Constraint_mesh_pt;
  }

private:

  // [zdec] this is not necessarily inside the domain for this geometry, needs
  // fixing
  //
  // /// Pin all displacements and rotation (dofs 0-4) at the centre
  // void pin_all_displacements_and_rotation_at_centre_node();

  /// Setup and build the mesh
  void build_mesh();

  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();

  /// Helper function to (re-)set boundary condition
  /// and complete the build of  all elements
  void complete_problem_setup();

  /// Loop over all curved edges, then loop over elements and upgrade
  /// them to be curved elements
  void upgrade_edge_elements_to_curve(const unsigned &b);

  /// Duplicate nodes at corners in order to properly apply boundary
  /// conditions from each edge. Also adds (8) Lagrange multiplier dofs to the
  /// problem in order to constrain continuous interpolation here across its (8)
  /// vertex dofs. (Note "corner" here refers to the meeting point of any two
  /// sub-boundaries in the closed external boundary)
  void duplicate_corner_nodes();

  /// Loop over all edge elements and rotate the Hermite degrees of freedom
  /// to be in the directions of the two in-plane vectors specified in Parameters
  void rotate_edge_degrees_of_freedom();

  /// Delete traction elements and wipe the surface mesh
  void delete_traction_elements(Mesh* const &surface_mesh_pt);

  /// Trace file to document norm of solution
  ofstream Trace_file;

  /// Pointer to "bulk" mesh
  TriangleMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to mesh containing constraint elements
  Mesh* Constraint_mesh_pt;

  /// Enumeration to keep track of boundary ids
  enum
  {
    Outer_boundary0 = 0,
    Outer_boundary1 = 1
  };

  /// Target element area
  double Element_area;

  /// Doc info object for labeling output
  DocInfo Doc_info;

  // /// Outer boundary Geom Object
  // Ellipse* Outer_boundary_ellipse_pt;

  /// The outer curves
  Vector<TriangleMeshCurveSection*> Outer_curvilinear_boundary_pt;

  /// The close outer boundary
  TriangleMeshClosedCurve* Outer_boundary_pt;

}; // end_of_problem_class



//======================================================================
/// Constructor definition
//======================================================================
template<class ELEMENT>
UnstructuredKSProblem<ELEMENT>::UnstructuredKSProblem(const double& element_area)
  :
  Element_area(element_area)
{
  // Build the mesh
  build_mesh();

  // Curved Edge upgrade
  upgrade_edge_elements_to_curve(Outer_boundary0);
  upgrade_edge_elements_to_curve(Outer_boundary1);

  // Loop over boundary nodes and check they correspond to the nodes as they
  // understand
  unsigned n_bound = 2;
  for(unsigned i_bound = 0; i_bound < n_bound; i_bound++)
  {
    unsigned n_b_node = Bulk_mesh_pt->nboundary_node(i_bound);
    for(unsigned i_b_node = 0; i_b_node < n_b_node; i_b_node++)
    {
      oomph_info << "Node " << i_b_node << " is on boundary " << i_bound << std::endl;
    }
  }

  // Rotate degrees of freedom
  rotate_edge_degrees_of_freedom();

  // Store number of bulk elements
  complete_problem_setup();

  // Set directory
  Doc_info.set_directory("RESLT");

  // Open trace file
  char filename[100];
  sprintf(filename, "RESLT/trace.dat");
  Trace_file.open(filename);

  // Assign equation numbers
  oomph_info << "Number of equations: "
             << assign_eqn_numbers() << '\n';

  ofstream debug;
  debug.open("boundary_nodes0.txt");
  for(unsigned b = 0; b < 2; b++)
  {
    for(unsigned n = 0; n<Bulk_mesh_pt->nnode(); n++)
    {
      if(Bulk_mesh_pt->node_pt(n)->is_on_boundary(b))
      {
	debug << b << ", " << Bulk_mesh_pt->node_pt(n)->x(0)
	      << ", " << Bulk_mesh_pt->node_pt(n)->x(1) << std::endl;
      }
    }
  }
  debug.close();
  debug.open("boundary_nodes1.txt");
  for(unsigned b = 0; b < 2; b++)
  {
    for(unsigned n = 0; n<Bulk_mesh_pt->nboundary_node(b); n++)
    {
      debug << b << ", " << Bulk_mesh_pt->boundary_node_pt(b,n)->x(0)
	    << ", " << Bulk_mesh_pt->boundary_node_pt(b,n)->x(1) << std::endl;
    }
  }
  debug.close();

  // Doc the equations
  describe_dofs();

} // end Constructor




//======================================================================
/// Set up and build the mesh
//======================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::build_mesh()
{
  // Allow for a slightly larger mismatch between vertex positions
  ToleranceForVertexMismatchInPolygons::Tolerable_error = 1.0e-14;


  Vector<double> zeta(1);
  Vector<double> posn(2);

  //Outer boundary
  //--------------
  double theta1=Parameters::Theta1;
  double theta2=Parameters::Theta2;

  // [zdec] debug
  oomph_info << theta1 << " " << theta2 << std::endl;

  //First bit
  double zeta_start = theta1;
  double zeta_end = -theta1;
  unsigned nsegment = (int)(0.5*(theta2-theta1)/sqrt(Element_area))+2;

  Outer_curvilinear_boundary_pt.resize(2);
  Outer_curvilinear_boundary_pt[0] =
    new TriangleMeshCurviLine(Parameters::Parametric_curve_pt[0], zeta_start,
                              zeta_end, nsegment, Outer_boundary0);


  zeta_start = theta2;
  zeta_end = -theta2;
  Outer_curvilinear_boundary_pt[1] =
    new TriangleMeshCurviLine(Parameters::Parametric_curve_pt[1], zeta_start,
                              zeta_end, nsegment, Outer_boundary1);

  // Combine
  Outer_boundary_pt =
    new TriangleMeshClosedCurve(Outer_curvilinear_boundary_pt);

  // [zdec] debug: output mesh boundary
  ofstream mesh_debug;
  mesh_debug.open("boundary_file.dat");
  Outer_boundary_pt->output(mesh_debug, 200);
  mesh_debug.close();

  //Create mesh parameters object
  TriangleMeshParameters mesh_parameters(Outer_boundary_pt);

  // Element area
  mesh_parameters.element_area() = Element_area;

  // Build an assign bulk mesh
  Bulk_mesh_pt=new TriangleMesh<ELEMENT>(mesh_parameters);
  Bulk_mesh_pt->setup_boundary_element_info();

  // Build mesh to contain constraint elements
  Constraint_mesh_pt = new Mesh();

  // Reset the non-vertex node positions
  for(unsigned e=0; e<Bulk_mesh_pt->nelement(); e++)
  {
    // [zdec] debug
    oomph_info << "Repairing for element " << e << std::endl;
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
    el_pt->repair_lagrange_node_positions();
  }

  // Split elements that have two boundary edges
  TimeStepper* time_stepper_pt = Bulk_mesh_pt->Time_stepper_pt;
  Bulk_mesh_pt->
    template split_elements_with_multiple_boundary_edges<ELEMENT>(time_stepper_pt);

  // [debug] Print the number of boundaries each node is on (there should only
  // be 2 nodes on 2 boundaries)
  unsigned n_node = Bulk_mesh_pt->nnode();
  for(unsigned i_node = 0; i_node < n_node; i_node++)
  {
    unsigned tot = 0;
    Node* node_pt = Bulk_mesh_pt->node_pt(i_node);
    unsigned n_bound = Bulk_mesh_pt->nboundary();
    for(unsigned i_bound = 0; i_bound < n_bound; i_bound++)
    {
      tot += node_pt->is_on_boundary(i_bound);
    }
    cout << "Node " << i_node
      << " at (" << node_pt->x(0) << "," << node_pt->x(1)
      << ") is on " << tot << " boundaries: ";
    std::set<unsigned>* boundaries_pt;
    BoundaryNodeBase* bnode_pt = dynamic_cast<BoundaryNodeBase*>(node_pt);
    if(bnode_pt)
    {
      bnode_pt->get_boundaries_pt(boundaries_pt);
      for(unsigned b: *boundaries_pt)
      {
     oomph_info << b << ", ";
      }
    }
    oomph_info << std::endl;
  }

  // Add extra nodes at boundaries and constrain the dofs there.
  duplicate_corner_nodes();

  // [debug] Print the number of boundaries each node is on (there should only
  // be 2 nodes on 2 boundaries)
  n_node = Bulk_mesh_pt->nnode();
  for(unsigned i_node = 0; i_node < n_node; i_node++)
  {
    unsigned tot = 0;
    Node* node_pt = Bulk_mesh_pt->node_pt(i_node);
    unsigned n_bound = Bulk_mesh_pt->nboundary();
    for(unsigned i_bound = 0; i_bound < n_bound; i_bound++)
    {
      tot += node_pt->is_on_boundary(i_bound);
    }
    cout << "Node " << i_node
      << " at (" << node_pt->x(0) << "," << node_pt->x(1)
      << ") is on " << tot << " boundaries: ";
    std::set<unsigned>* boundaries_pt;
    BoundaryNodeBase* bnode_pt = dynamic_cast<BoundaryNodeBase*>(node_pt);
    if(bnode_pt)
    {
      bnode_pt->get_boundaries_pt(boundaries_pt);
      for(unsigned b: *boundaries_pt)
      {
     oomph_info << b << ", ";
      }
    }
    oomph_info << std::endl;
  }

  //Add submesh to problem
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Constraint_mesh_pt);

  // Combine submeshes into a single Mesh (over the top; could just have
  // assigned bulk mesh directly.
  build_global_mesh();

}// end build_mesh




// //==start_of_pin_all_displacements_and_rotation_at_centre_node======================
// /// pin all displacements and rotations in the centre
// //==============================================================================
// template<class ELEMENT>
// void UnstructuredKSProblem<ELEMENT>::pin_all_displacements_and_rotation_at_centre_node()
// {


//   // Choose non-centre node on which we'll supress
//   // the rigid body rotation around the z axis.
//   double max_x_potentially_pinned_node=-DBL_MAX;
//   Node* pinned_rotation_node_pt=0;

//   // Pin the node that is at the centre in the domain
//   unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(2);
//   for (unsigned inod=0;inod<num_int_nod;inod++)
//   {
//     // Get node point
//     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(2,inod);

//     // Check which coordinate increases along this boundary
//     // oomph_info << "node: "
//     //            << nod_pt->x(0) << " "
//     //            << nod_pt->x(1) << " "
//     //            << std::endl;

//     // Find the node with the largest x coordinate
//     if (fabs(nod_pt->x(0))>max_x_potentially_pinned_node)
//     {
//       max_x_potentially_pinned_node=fabs(nod_pt->x(0));
//       pinned_rotation_node_pt=nod_pt;
//     }


//     // If the node is on the other internal boundary too
//     if( nod_pt->is_on_boundary(3))
//     {
//       // Pin it! It's the centre of the domain!
//       // In-plane dofs are always 0 and 1
//       // Out of plane displacement is 2, x and y derivatives are 3 and 4.
//       nod_pt->pin(0);
//       nod_pt->set_value(0,0.0);
//       nod_pt->pin(1);
//       nod_pt->set_value(1,0.0);
//       nod_pt->pin(2);
//       nod_pt->set_value(2,0.0);
//       nod_pt->pin(3);
//       nod_pt->set_value(3,0.0);
//       nod_pt->pin(4);
//       nod_pt->set_value(4,0.0);


//     }
//   }


//   oomph_info << "rotation pinning node: "
//   << pinned_rotation_node_pt->x(0) << " "
//   << pinned_rotation_node_pt->x(1) << " "
//   << std::endl;
//   // Pin y displacement
//   pinned_rotation_node_pt->pin(1);

// }



//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of
/// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::complete_problem_setup()
{
  // Complete the build of all elements so they are fully functional
  unsigned n_element = Bulk_mesh_pt->nelement();
  for(unsigned e=0;e<n_element;e++)
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

    // Use the fd jacobian
    // el_pt->enable_finite_difference_jacobian();
  }

  // Set the boundary conditions
  apply_boundary_conditions();
}



//==Start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::apply_boundary_conditions()
{
  // Out-of-plane dofs:
  //-------------------
  // |  0  |  1  |  2  |  3  |  4  |  5  |
  // |  w  | w_n | w_t | w_nn| w_nt| w_tt|

  // Case: The plate edge is completely free
  const Vector<unsigned> free{};

  // Case: The plate is pinned (w given, dw/dn left free) along a boundary.
  // We therefore have to pin (and assign values for) w, dw/dt and d^2w/dt^2
  const Vector<unsigned> resting_pin_dofs{0, 2, 5};

  // Case: The plate is sliding (w left free, dw/dn given) along a boundary.
  // We therefore have to pin (and assign values for) dw/dn and d^2w/dndt
  const Vector<unsigned> sliding_clamp_dofs{1, 4};

  // Case: The plate is clamped (w given, dw/dn given) along a boundary.
  // We therefore have to pin (and assign values for) w, dw/dn,
  // dw/dt, d^2w/dndt and d^2w/dt^2
  const Vector<unsigned> true_clamp_dofs{0, 1, 2, 4, 5};

  //------------------------------------------------------------------
  //------------------------------------------------------------------

  // Number of fields we are interpolating
  unsigned n_field = 3;
  // Number of external boundaries to apply conditions to
  unsigned n_bound = 2;

  // Vector containers to store which boundary conditions we are applying to
  // each edge. (outlined above)
  Vector<Vector<Vector<unsigned>>>
    pinned_u_dofs(n_bound, Vector<Vector<unsigned>>(n_field));

  for(unsigned i_bound = 0; i_bound < n_bound; i_bound++)
  {
    for (unsigned i_field = 0; i_field < n_field; i_field++)
    {
      char cond = Parameters::Bc_char[i_bound][i_field];
      cout << "Edge " << i_bound
	   << " displacement " << i_field << " is " << cond << std::endl;
      switch (cond)
      {
        case 'f':
          pinned_u_dofs[i_bound][i_field] = free;
          break;
        case 'p':
          pinned_u_dofs[i_bound][i_field] = resting_pin_dofs;
          break;
        case 's':
          pinned_u_dofs[i_bound][i_field] = sliding_clamp_dofs;
          break;
        case 'c':
          pinned_u_dofs[i_bound][i_field] = true_clamp_dofs;
          break;
      }
    }
  }

  // Reset by unpinning all the nodal dofs in the mesh
  unsigned n_node = Bulk_mesh_pt->nnode();
  for (unsigned i_node = 0; i_node < n_node; i_node++)
  {
    Bulk_mesh_pt->node_pt(i_node)->unpin_all();
  }

  // Loop over all the boundaries in our bulk mesh
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


  // Update the corner constraintes based on boundary conditions
  unsigned n_el = Constraint_mesh_pt->nelement();
  for(unsigned i_el = 0; i_el < n_el; i_el++)
  {
    dynamic_cast<DuplicateNodeConstraintElement*>
      (Constraint_mesh_pt->element_pt(i_el))
      ->validate_and_pin_redundant_constraints();
  }

  // Assign the equation numbers
  assign_eqn_numbers();
} // end set bc



//==============================================================================
/// A function that upgrades straight sided elements to be curved. This involves
/// Setting up the parametric boundary, F(s) and the first derivative F'(s)
/// We also need to set the edge number of the upgraded element and the positions
/// of the nodes j and k (defined below) and set which edge (k) is to be exterior
/*            @ k                                                             */
/*           /(                                                               */
/*          /. \                                                              */
/*         /._._)                                                             */
/*      i @     @ j                                                           */
/// For RESTING or FREE boundaries we need to have a C2 CONTINUOUS boundary
/// representation. That is we need to have a continuous 2nd derivative defined
/// too. This is well discussed in by [Zenisek 1981] (Aplikace matematiky ,
/// Vol. 26 (1981), No. 2, 121--141). This results in the necessity for F''(s)
/// as well.
//==start_of_upgrade_edge_elements==============================================
template <class ELEMENT>
void UnstructuredKSProblem<ELEMENT >::
upgrade_edge_elements_to_curve(const unsigned &ibound)
{

  // [zdec] debug
  ofstream debug_stream;
  debug_stream.open("mesh_debug.dat");
  Bulk_mesh_pt->output(debug_stream);
  debug_stream.close();

  // Parametric curve describing the boundary
  // hierher code share with thing that defines the
  // boundaries of the mesh!
  CurvilineGeomObject* parametric_curve_pt=0;
  switch (ibound)
  {
  case 0:
    parametric_curve_pt = &Parameters::Upper_parametric_elliptical_curve;
    break;

  case 1:
    parametric_curve_pt = &Parameters::Lower_parametric_elliptical_curve;
    break;

  default:
    throw OomphLibError("Unexpected boundary number.",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
    break;
  } // end parametric curve switch


  // Loop over the bulk elements adjacent to boundary ibound
  const unsigned n_els=Bulk_mesh_pt->nboundary_element(ibound);
  for(unsigned e=0; e<n_els; e++)
  {
    // Get pointer to bulk element adjacent to b
    ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(ibound,e));

    // hierher what is that? why "My"?
    // Initialise enum for the curved edge
    MyC1CurvedElements::Edge edge(MyC1CurvedElements::none);

    // Loop over all (three) nodes of the element and record boundary nodes
    unsigned index_of_interior_node = 3;
    unsigned nnode_on_neither_boundary = 0;
    const unsigned nnode = 3;


    // hierher what does this comment mean?
    // Fill in vertices' positions (this step should be moved inside the curveable
    // Bell element)
    Vector<Vector<double> > xn(nnode,Vector<double>(2,0.0));
    for(unsigned n=0;n<nnode;++n)
    {
      Node* nod_pt = bulk_el_pt->node_pt(n);
      xn[n][0]=nod_pt->x(0);
      xn[n][1]=nod_pt->x(1);

      // Check if it is on the outer boundaries
      if(!(nod_pt->is_on_boundary(ibound)))
      {
        index_of_interior_node = n;
        ++nnode_on_neither_boundary;
      }
    }// end record boundary nodes


    // hierher: ouch! This seems to map (x,y) to zeta! This is at best possible to within
    // a tolerance. Needs a redesign!

    // s at the next (cyclic) node after interior
    const double s_ubar = parametric_curve_pt->get_zeta(xn[(index_of_interior_node+1) % 3]);

    // s at the previous (cyclic) node before interior
    const double s_obar = parametric_curve_pt->get_zeta(xn[(index_of_interior_node+2) % 3]);

    // Assign edge case
    edge = static_cast<MyC1CurvedElements::Edge>(index_of_interior_node);

    // Check nnode_on_neither_boundary
    if(nnode_on_neither_boundary == 0)
    {
      throw OomphLibError(
        "No interior nodes. One node per CurvedElement must be interior.",
        OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    else if (nnode_on_neither_boundary > 1)
    {
      throw OomphLibError(
        "Multiple interior nodes. Only one node per CurvedElement can be interior.",
        OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check for inverted elements
    if (s_ubar>s_obar)
    {
      throw OomphLibError(
        "Decreasing parametric coordinate. Parametric coordinate must increase as the edge is traversed anti-clockwise.",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    } // end checks

    // Upgrade it
    bulk_el_pt->upgrade_element_to_curved(edge,s_ubar,s_obar,parametric_curve_pt,5);
  }
}// end_upgrade_elements



//==============================================================================
/// Duplicate nodes at corners in order to properly apply boundary
/// conditions from each edge. Also adds (8) Lagrange multiplier dofs to the
/// problem in order to constrain continuous interpolation here across its (8)
/// vertex dofs. (Note "corner" here refers to the meeting point of any two
/// sub-boundaries in the closed external boundary)
//==============================================================================
template <class ELEMENT>
void UnstructuredKSProblem<ELEMENT >::duplicate_corner_nodes()
{
  // Loop over the sections of the external boundary
  unsigned n_bound = Bulk_mesh_pt->nboundary();
  for(unsigned i_bound = 0; i_bound < n_bound; i_bound++)
  {
    // Store the index of the next boundary
    unsigned ip1_bound = (i_bound+1)%n_bound;
    // Storage for node and el pts at the boundary vertex
    Node* old_node_pt = 0;
    Node* new_node_pt = 0;
    // FiniteElement* left_element_pt = 0;
    FiniteElement* right_element_pt = 0;

    // To find the node between boundaries i and i+1, we loop over all nodes on
    // boundary i until we find the one that is also on i+1, we then loop over
    // all boundary elements on i and i+1 until we find the elements that sit
    // either side of the corner. (there might be a smarter way but this was the
    // first idea I had -- Aidan)

    //----------------------------------------------------------------------
    // First, find corner the node
    unsigned n_b_node = Bulk_mesh_pt->nboundary_node(i_bound);
    for(unsigned i_b_node = 0; i_b_node < n_b_node; i_b_node++)
    {
      // Store the node we are checking
      Node* node_pt = Bulk_mesh_pt->boundary_node_pt(i_bound,i_b_node);

      // If it is on the next boundary we have found the corner node
      if(node_pt->is_on_boundary(ip1_bound))
      {
        // [zdec] debug
        oomph_info << "Found a corner node at " << std::endl << "  ("
                   << node_pt->position(0) << "," << node_pt->position(1) << ")"
                   << std::endl;
        old_node_pt = node_pt;
	break;
      }
    }

    //----------------------------------------------------------------------
    // Find the right (i+1th boundary) side element
    unsigned n_b_el = Bulk_mesh_pt->nboundary_element(ip1_bound);
    for (unsigned i_b_el = 0; i_b_el < n_b_el; i_b_el++)
    {
      // Get the element pointer
      FiniteElement* el_pt = Bulk_mesh_pt->boundary_element_pt(ip1_bound, i_b_el);
      // If the corner node pt is in the element we have found the right
      // element
      if (el_pt->get_node_number(old_node_pt) != -1)
      {
        right_element_pt = el_pt;
        break;
      }
    }

    //----------------------------------------------------------------------
    // Now we need to create a new node and substitute the right elements
    // old corner node for this new one
    new_node_pt = right_element_pt->construct_boundary_node(
      right_element_pt->get_node_number(old_node_pt));
    // Copy the position and other info from the old node into the new node
    // [debug]
    oomph_info << "About to copy node data" << std::endl;
    new_node_pt->x(0)=old_node_pt->x(0);
    new_node_pt->x(1)=old_node_pt->x(1);
    oomph_info << "Copied node data" << std::endl;
    // Then we add this node to the mesh
    Bulk_mesh_pt->add_node_pt(new_node_pt);
    // Then replace the old node for the new one on the right boundary
    Bulk_mesh_pt->remove_boundary_node(ip1_bound,old_node_pt);
    Bulk_mesh_pt->add_boundary_node(ip1_bound,new_node_pt);

    //----------------------------------------------------------------------
    // The final job is to constrain this duplication using the specialised
    // Lagrange multiplier elements which enforce equality of displacement and
    // its derivatives either side of this corner.
    CurvilineGeomObject* left_parametrisation_pt =
      Parameters::Parametric_curve_pt[i_bound];
    CurvilineGeomObject* right_parametrisation_pt =
      Parameters::Parametric_curve_pt[ip1_bound];

    // Get the coordinates on each node on their respective boundaries
    Vector<double> left_boundary_coordinate =
      {left_parametrisation_pt->get_zeta(old_node_pt->position())};
    Vector<double> right_boundary_coordinate =
      {right_parametrisation_pt->get_zeta(new_node_pt->position())};

    // Create the constraining element
    DuplicateNodeConstraintElement* constraint_element_pt =
      new DuplicateNodeConstraintElement(old_node_pt,
                                         new_node_pt,
                                         left_parametrisation_pt,
                                         right_parametrisation_pt,
                                         left_boundary_coordinate,
                                         right_boundary_coordinate);

    // Add the constraining element to the mesh
    Constraint_mesh_pt->add_element_pt(constraint_element_pt);
  }
}



//======================================================================
/// Function to set up rotated nodes on the boundary: necessary if we want to set
/// up physical boundary conditions on a curved boundary with Hermite type dofs.
/// For example if we know w(n,t) = f(t) (where n and t are the
/// normal and tangent to a boundary) we ALSO know dw/dt and d2w/dt2.
/// NB no rotation is needed if the edges are completely free!
//======================================================================
template <class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::rotate_edge_degrees_of_freedom()
{

  std::ofstream debugstream;
  debugstream.open("test_file");
  debugstream << "x y nx ny tx ty nx_x nx_y ny_x ny_y tx_x tx_y ty_x ty_y"
              << std::endl;
  debugstream.close();

  // Get the number of boundaries
  unsigned n_bound = 2;
  // Loop over the bulk elements
  unsigned n_element = Bulk_mesh_pt-> nelement();
  for(unsigned e=0; e<n_element; e++)
  {
    // Get pointer to bulk element adjacent
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // [zdec] debug
    oomph_info << "In element " << e << " (" << el_pt << ")" << std::endl;

    // Loop over each boundary and add the boundary parametrisation to the
    // relevant nodes' boundary data
    for(unsigned b=0; b<n_bound; b++)
    {

      // Calculate nodes on the relevant boundaries
      const unsigned nnode=3;
      // Count the number of boundary nodes on external boundaries
      Vector<unsigned> boundary_node;
      // Store the boundary coordinates of nodes on the boundaries
      Vector<double> boundary_coordinate_of_node;
      for (unsigned n=0; n<nnode;++n)
      {
        // If on external boundary b
        if (el_pt->node_pt(n)->is_on_boundary(b))
        {
          boundary_node.push_back(n);
          double coord = Parametric_curve_pt[b]
            ->get_zeta(el_pt->node_pt(n)->position());
          boundary_coordinate_of_node.push_back(coord);
        }
      }

      // If the element has nodes on the boundary, rotate the Hermite dofs
      if(!boundary_node.empty())
      {
        // [zdec] debug
        oomph_info << " Nodes ";
        for(unsigned n: boundary_node)
        {oomph_info << n << " "; }
        oomph_info << " are on boundary " << b << std::endl;
        // Rotate the nodes by passing the index of the nodes and the
        // normal / tangent vectors to the element
        el_pt->
          rotated_boundary_helper_pt()->
          set_nodal_boundary_parametrisation(boundary_node,
                                             boundary_coordinate_of_node,
                                             Parametric_curve_pt[b]);
      }
    }
  }
}// end rotate_edge_degrees_of_freedom



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::doc_solution(const
                                                   std::string& comment)
{
  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts = 30;

  sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file,npts);
  some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
            << comment << "\"\n";
  some_file.close();


  // // Find the solution at r=0
  // // ----------------------

  // // hierher precompute
  // MeshAsGeomObject Mesh_as_geom_obj(Bulk_mesh_pt);
  // Vector<double> s(2);
  // GeomObject* geom_obj_pt=0;
  // Vector<double> r(2,0.0);
  // Mesh_as_geom_obj.locate_zeta(r,geom_obj_pt,s);

  // // Compute the interpolated displacement vector
  // Vector<double> u_0(12,0.0);
  // u_0=dynamic_cast<ELEMENT*>(geom_obj_pt)->interpolated_u_foeppl_von_karman(s);

  // oomph_info << "w in the middle: " <<std::setprecision(15) << u_0[0] << std::endl;

  // Trace_file << Parameters::P_mag << " " << u_0[0] << '\n';

  // Increment the doc_info number
  Doc_info.number()++;

} // end of doc



//=======start_of_main========================================
///Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

  // Create the problem, using FvK elements derived from TElement<2,4>
  // elements (with 4 nodes per element edge and 10 nodes overall).
  UnstructuredKSProblem<KoiterSteigmannC1CurvableBellElement>
    problem(Parameters::Element_area);

  problem.max_residuals() = 1.0e3;
  problem.max_newton_iterations() = 30;
  problem.newton_solver_tolerance() = 1.0e-11;

  // Set pressure
  Parameters::P_mag = 1.0e-3;
  // Set the Poisson ratio
  Parameters::Nu = 0.5;

  // Document the initial state
  problem.doc_solution();
  // Solve the system
  problem.newton_solve();
  // Document the current solution
  problem.doc_solution();

  // Change the bcs to be clamped-free
  strcpy(Parameters::Bc_char[0], "ppc");
  strcpy(Parameters::Bc_char[1], "fff");
  problem.reapply_boundary_conditions();
  // Solve the system
  problem.newton_solve();
  // Document the current solution
  problem.doc_solution();

} //End of main
