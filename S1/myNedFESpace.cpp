#include <comp.hpp>  
#include "myNedElement.hpp"
#include "myNedFESpace.hpp"
#include "myDiffOpHCurl.hpp"


namespace ngcomp
{

  MyNedFESpace :: MyNedFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    type = "mynedfespace";

    order = int(flags.GetNumFlag("order", 0));

    if (ma->GetDimension() == 2)
      {
	evaluator[VOL] = make_shared<T_DifferentialOperator<MyNedDiffOpId<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<MyNedDiffOpCurl<2>>>();
      }
    else
      {
	evaluator[VOL] = make_shared<T_DifferentialOperator<MyNedDiffOpId<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<MyNedDiffOpCurl<3>>>();
      }
  }
    
  void MyNedFESpace :: Update()
  {
    // some global update:

    int n_edge = ma->GetNEdges();

    first_edge_dof.SetSize(n_edge + 1);
    int ii = 0;
    for (int i = 0; i < n_edge; i++, ii += (order+1))
      first_edge_dof[i] = ii;
    first_edge_dof[n_edge] = ii;

//         first_cell_dof.SetSize(n_cell + 1);
//         for (int i = 0; i < n_cell; i++, ii += 6)
//             first_cell_dof[i] = ii;
//         first_cell_dof[n_cell] = ii;

    SetNDof (ii);
  }

  void MyNedFESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    // returns dofs of element number elnr
    dnums.SetSize(0);

    auto ngel = ma->GetElement(ei);

    for (auto e : ngel.Edges())
        for (auto j: Range(first_edge_dof[e], first_edge_dof[e + 1]))
    	  dnums.Append(j);
  }

  
  FiniteElement & MyNedFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    auto ngel = ma->GetElement (ei);
    switch (ngel.GetType())
      {
      case ET_TRIG:
        {
          auto trig = new (alloc) MyNedElementTrig(order);
          trig->SetVertexNumbers (ngel.vertices);
          return *trig;
        }
      case ET_TET:
        {
          auto tet = new (alloc) MyNedElementTet(order);
          tet->SetVertexNumbers (ngel.vertices);
          return *tet;
        }
      default:
        throw Exception (string("Element type ")+ToString(ngel.GetType())+" not supported");
      }
  }
  
}


void ExportNedFESpace(py::module m) {
    using namespace ngcomp;

    cout << "called ExportNedFESpace" << endl;

    ExportFESpace<MyNedFESpace>(m, "MyNedFESpace", true);
}