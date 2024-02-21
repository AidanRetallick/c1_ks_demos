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
using MathematicalConstants::Pi;


//===========================================================================
/// Namespace for problem parameters
//===========================================================================
namespace Parameters
{
  /// Square edge length
  double L = 1.0;

  /// Angle of rotation (radians obviously)
  double Alpha = Pi / 5.0;

  /// Plate vertices. For the unrotated square, these coincide with:
  ///     L/2*{(-1,-1), (1,-1), (1,1), (-1,1)}
  Vector<Vector<double>> Vertices = {
    {-L / 2.0 * cos(Alpha) + L / 2.0 * sin(Alpha),
     -L / 2.0 * sin(Alpha) - L / 2.0 * cos(Alpha)},
    {L / 2.0 * cos(Alpha) + L / 2.0 * sin(Alpha),
     L / 2.0 * sin(Alpha) - L / 2.0 * cos(Alpha)},
    {L / 2.0 * cos(Alpha) - L / 2.0 * sin(Alpha),
     L / 2.0 * sin(Alpha) + L / 2.0 * cos(Alpha)},
    {-L / 2.0 * cos(Alpha) - L / 2.0 * sin(Alpha),
     -L / 2.0 * sin(Alpha) + L / 2.0 * cos(Alpha)}
  };

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

  /// Magnitude of shear stress
  double T_mag = 0.0;

  /// Element size
  double Element_area = 0.2;

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



  // ---- Parametric boundaries ------------------------------------------------

  /// Straight edge 0
  CurvilineLine Edge_0(Vertices[0], Vertices[1]);

  /// Straight edge 1
  CurvilineLine Edge_1(Vertices[1], Vertices[2]);

  /// Straight edge 2
  CurvilineLine Edge_2(Vertices[2], Vertices[3]);

  /// Straight edge 3
  CurvilineLine Edge_3(Vertices[3], Vertices[0]);

  /// Vector container of addresses for iterating over the edges
  Vector<CurvilineGeomObject*> Curviline_edge_pt =
  {
    &Edge_0,
    &Edge_1,
    &Edge_2,
    &Edge_3
  };


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

  /// Things to repeat after every newton iteration
  void actions_before_newton_step()
  {
    // File prefix and suffix strings
    std::string prefix = Doc_info.directory();
    std::string suffix = std::to_string(Doc_info.number()) + ".txt";

    // Get the residual and jacobian
    LinearAlgebraDistribution* dist = dof_distribution_pt();
    DoubleVector residual(dist);
    CRDoubleMatrix jacobian(dist);
    get_jacobian(residual, jacobian);
    residual.output("residual" + suffix);
    jacobian.sparse_indexed_output("jacobian" + suffix);

    // [zdec] debug
    // Output the solution after every step
    doc_solution();
  }

  /// Print information about the parameters we are trying to solve for.
  void actions_before_newton_solve()
  {
    oomph_info << "-------------------------------------------------------"
               << std::endl;
    oomph_info << "Solving for P = " << Parameters::P_mag << std::endl;
    oomph_info << "         step = " << Doc_info.number() << std::endl;
    oomph_info << "-------------------------------------------------------"
               << std::endl;
  }

  /// Doc the solution
  void doc_solution(const std::string& comment = "");

  // [zdec] This doesn't work, dynamic cast always fails -- returns 0
  // /// Overloaded version of the problem's access function to
  // /// the mesh. Recasts the pointer to the base Mesh object to
  // /// the actual mesh type.
  // TriangleMesh<ELEMENT>* mesh_pt()
  // {
  //   oomph_info << Problem::mesh_pt() << std::endl;
  //   oomph_info << dynamic_cast<TriangleMesh<ELEMENT>*> (Problem::mesh_pt())
  //   << std::endl; return dynamic_cast<TriangleMesh<ELEMENT>*>
  //   (Problem::mesh_pt());
  // }

  /// Overloaded version of the problem's access function to
  /// the mesh. Recasts the pointer to the base Mesh object to
  /// the actual mesh type.
  TriangleMesh<ELEMENT>* mesh_pt()
  {
    return Bulk_mesh_pt;
  }

private:
  /// Setup and build the mesh
  void build_mesh();

  /// Duplicate corner nodes and create constraint elements at those corners
  void duplicate_corner_nodes();

  /// Helper function to (re-)set boundary condition
  /// and complete the build of all elements
  void complete_problem_setup();

  /// Loop over all edge elements and rotate the Hermite degrees of freedom
  /// to be in the directions of the two in-plane vectors specified in
  /// Parameters
  void rotate_edge_degrees_of_freedom(Mesh* const& bulk_mesh_pt);

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

  /// Keep track of boundary ids, (b)ottom, (r)ight, (t)op, (l)eft
  // (slightly redundant in this example)
  // ((after rotation this naming convention is unhelpful))
  enum
  {
    Boundary_b_bnum = 0,
    Boundary_r_bnum = 1,
    Boundary_t_bnum = 2,
    Boundary_l_bnum = 3
  };

  /// Pointer to "bulk" mesh
  TriangleMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to constraint mesh
  Mesh* Constraint_mesh_pt;

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

  // Rotate degrees of freedom
  rotate_edge_degrees_of_freedom(Bulk_mesh_pt);

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
  // Get the vertices from the parameters
  Vector<double> vertex0 = Parameters::Vertices[0];
  Vector<double> vertex1 = Parameters::Vertices[1];
  Vector<double> vertex2 = Parameters::Vertices[2];
  Vector<double> vertex3 = Parameters::Vertices[3];

  // Declare the edges...
  Vector<Vector<double>> edge0(2, Vector<double>(2, 0.0));
  Vector<Vector<double>> edge1(2, Vector<double>(2, 0.0));
  Vector<Vector<double>> edge2(2, Vector<double>(2, 0.0));
  Vector<Vector<double>> edge3(2, Vector<double>(2, 0.0));

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

  // Split elements that have two boundary edges
  TimeStepper* time_stepper_pt = Bulk_mesh_pt->Time_stepper_pt;
  Bulk_mesh_pt->
    template split_elements_with_multiple_boundary_edges<ELEMENT>(time_stepper_pt);

  // Create the empty constraint element mesh
  Constraint_mesh_pt = new Mesh();

  // Add extra nodes at boundaries and constrain the dofs there.
  duplicate_corner_nodes();

  //Add submesh to problem
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Constraint_mesh_pt);

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

    // Use the fd jacobian
    // el_pt->enable_finite_difference_jacobian();
  }
  // Set the boundary conditions
  apply_boundary_conditions();

} // end of complete



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
  unsigned n_bound = 4;
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
      Parameters::Curviline_edge_pt[i_bound];
    CurvilineGeomObject* right_parametrisation_pt =
      Parameters::Curviline_edge_pt[ip1_bound];

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



//==============================================================================
/// Function to set up rotated nodes on the boundary: necessary if we want to
/// set up physical boundary conditions on a curved boundary with Hermite type
/// dofs. For example if we know w(n,t) = f(t) (where n and t are the normal and
/// tangent to a boundary) we ALSO know dw/dt and d2w/dt2. NB no rotation is
/// needed if the edges are completely free! begin
/// rotate_edge_degrees_of_freedom
//==============================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::rotate_edge_degrees_of_freedom(
  Mesh* const& bulk_mesh_pt)
{
  // Get the number of boundaries
  unsigned n_bound = 4;

  // Loop over the bulk elements
  unsigned n_element = Bulk_mesh_pt->nelement();
  for (unsigned e = 0; e < n_element; e++)
  {
    // Get pointer to bulk element adjacent
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // Loop over each boundary and add the boundary parametrisation to the
    // relevant nodes' boundary data
    for (unsigned b = 0; b < n_bound; b++)
    {
      // Calculate nodes on the relevant boundaries
      const unsigned nnode = 3;
      // Count the number of boundary nodes on external boundaries
      Vector<unsigned> boundary_node;
      // Store the boundary coordinates of nodes on the boundaries
      Vector<double> boundary_coordinate_of_node;
      for (unsigned n = 0; n < nnode; ++n)
      {
        // If on external boundary b
        if (el_pt->node_pt(n)->is_on_boundary(b))
        {
          boundary_node.push_back(n);
          double coord = Parameters::Curviline_edge_pt[b]->get_zeta(
            el_pt->node_pt(n)->position());
          boundary_coordinate_of_node.push_back(coord);
        }
      }

      // If the element has nodes on the boundary, rotate the Hermite dofs
      if (!boundary_node.empty())
      {
        // Rotate the nodes by passing the index of the nodes and the
        // normal / tangent vectors to the element
        el_pt->rotated_boundary_helper_pt()->set_nodal_boundary_parametrisation(
          boundary_node,
          boundary_coordinate_of_node,
          Parameters::Curviline_edge_pt[b]);
      }
    }
  }
} // end rotate_edge_degrees_of_freedom



//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredKSProblem<ELEMENT>::apply_boundary_conditions()
{
  //------------------------------------------------------------------
  // Boundary conditions for FvK elements are complicated and we provide an
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
  //     enumerated from 0 to 5 in the order w, w_n, w_t, w_nn, w_nt,
  //     w_tt. These values are only stored at three vertices of the element.
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
  //          function of the coordinate, x, a 2D vector.
  //
  // So, if the function fix_in_plane_displacement_dof(idof, b, fct_pt) is
  // called with idof=1 and b=3, say, the y-in-plane displacement is pinned for
  // all the element's nodes (if any) that are located on mesh boundary 3. The
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
  // Using the conventions introduced above, the following vectors identify
  // the in-plane and out-of-plane degrees of freedom to be pinned for
  // various physically meaningful boundary conditions:


  // Possible boundary conditions for out-of-plane displacements: Given that the
  // out-of-plane displacements feature in the fourth-order biharmonic operator,
  // we can apply boundary conditions on w and dw/dn, where n is the coordinate
  // direction normal to the (assumed to be axis aligned!) boundary. However if
  // w is given along the entire boundary (parametrised by the tangential
  // coordinate, t) we also know what dw/dt and d^2w/dt^2 are. Likewise if dw/dn
  // is known along the whole boundary we also know d^2w/dndt. In the various
  // cases below we identify physical scenarios of a pinned edge (w given, dw/dn
  // left free); a vertically sliding edge (w left free; dw/dn given) and fully
  // clamped (w and dw/dn given). The four possible combinations of boundary
  // condition are: fully free -- nothing given, pinned edge -- only w(t) given,
  // sliding clamp -- only dwdt(t) given, fully clamped -- both w(t) and dwdt(t)
  // given.

  // Out-of-plane dofs:
  //-------------------
  // |  0  |  1  |  2  |  3  |  4  |  5  |
  // |  w  | w_n | w_t | w_nn| w_nt| w_tt|

  // Case: The plate is pinned (w given, dw/dn left free) along a boundary.
  // We therefore have to pin (and assign values for) w, dw/dt and d^2w/dt^2
  const Vector<unsigned> pinned_edge_pinned_dof{0, 2, 5};

  // Case: The plate is sliding (w left free, dw/dn given) along a boundary.
  // We therefore have to pin (and assign values for) dw/dn and d^2w/dndt
  const Vector<unsigned> sliding_clamp_pinned_dof{1, 4};

  // Case: The plate is clamped (w given, dw/dn given) along a boundary.
  // We therefore have to pin (and assign values for) w, dw/dn,
  // dw/dt, d^2w/dndt and d^2w/dt^2
  const Vector<unsigned> fully_clamped_pinned_dof{0, 1, 2, 4, 5};


  //------------------------------------------------------------------
  //------------------------------------------------------------------

  unsigned n_field = 3;

  // Vector containers to store which boundary conditions we are applying to
  // each edge. (outlined above)
  Vector<Vector<Vector<unsigned>>> pinned_u_dofs(4, Vector<Vector<unsigned>>(3));

  // Pin all three displacements everywhere
  for(unsigned i_field = 0; i_field < n_field; i_field++)
  {
    pinned_u_dofs[0][i_field] = pinned_edge_pinned_dof;
    pinned_u_dofs[1][i_field] = pinned_edge_pinned_dof;
    pinned_u_dofs[2][i_field] = pinned_edge_pinned_dof;
    pinned_u_dofs[3][i_field] = pinned_edge_pinned_dof;
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


  // Update the corner constraintes based on boundary conditions
  unsigned n_el = Constraint_mesh_pt->nelement();
  for(unsigned i_el = 0; i_el < n_el; i_el++)
  {
    dynamic_cast<DuplicateNodeConstraintElement*>
      (Constraint_mesh_pt->element_pt(i_el))
      ->validate_and_pin_redundant_constraints();
  }

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

  Trace_file << Parameters::P_mag << " " << Parameters::T_mag << " "
             << Doc_info.number() << endl;

  // Bump
  Doc_info.number()++;

} // end of doc


//=======start_of_main========================================
/// Driver code for demo of unstructured C1 Foeppl-von Karman
/// elements
//============================================================
int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

  // Create the problem, using FvK elements derived from TElement<2,4>
  // elements (with 4 nodes per element edge and 10 nodes overall).
  UnstructuredKSProblem<KoiterSteigmannC1CurvableBellElement> problem;

  problem.max_residuals() = 1.0e3;
  problem.max_newton_iterations() = 30;
  problem.newton_solver_tolerance() = 1.0e-11;

  // Set pressure
  Parameters::P_mag = 1.0e-3;
  // Set the Poisson ratio
  Parameters::Nu = 0.5;

  // Document the initial state
  problem.doc_solution();

  // Describe the dofs
  ofstream dofstream("RESLT/dofs.txt");
  problem.describe_dofs(dofstream);
  dofstream.close();

  // Get the jacobian
  LinearAlgebraDistribution* dist = problem.dof_distribution_pt();
  DoubleVector residual(dist);
  CRDoubleMatrix jacobian(dist);
  problem.get_jacobian(residual, jacobian);
  jacobian.sparse_indexed_output("RESLT/jacobian.txt");

  // Solve the system
  problem.newton_solve();
  // Document the current solution
  problem.doc_solution();

} // End of main
