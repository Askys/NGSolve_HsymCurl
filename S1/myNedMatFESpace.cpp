#include <comp.hpp>
#include "myNedElement.hpp"
#include "myNedMatFESpace.hpp"
#include "myDiffOpHCurl.hpp"


namespace ngcomp {

    MyNedMatFESpace::MyNedMatFESpace(shared_ptr <MeshAccess> ama, const Flags &flags)
            : FESpace(ama, flags) {
        type = "mynedmatfespace";

        order = int(flags.GetNumFlag("order", 0));

        if (ma->GetDimension() == 2) {
            evaluator[VOL] = make_shared < T_DifferentialOperator < MyNedMatDiffOpId < 2>>>();
            flux_evaluator[VOL] = make_shared < T_DifferentialOperator < MyNedMatDiffOpCurl < 2>>>();
        } else {
            evaluator[VOL] = make_shared < T_DifferentialOperator < MyNedMatDiffOpId < 3>>>();
            flux_evaluator[VOL] = make_shared < T_DifferentialOperator < MyNedMatDiffOpCurl < 3>>>();
        }
    }

    void MyNedMatFESpace::Update() {
        // some global update:

        int n_edge = ma->GetNEdges();
        int n_cell = ma->GetNE();

        first_edge_dof.SetSize(3 * n_edge + 1);
        int ii = 0;
        for (int i = 0; i < n_edge; i++, ii += 3 * (order + 1))
            first_edge_dof[i] = ii;
        first_edge_dof[n_edge] = ii;

        first_cell_dof.SetSize(n_cell + 1);
        for (int i = 0; i < n_cell; i++, ii += 16)
            first_cell_dof[i] = ii;
        first_cell_dof[n_cell] = ii;

        SetNDof(ii);
    }

    void MyNedMatFESpace::GetDofNrs(ElementId ei, Array <DofId> &dnums) const {
        // returns dofs of element number elnr
        dnums.SetSize(0);

        auto ngel = ma->GetElement(ei);

        // edge dofs
        for (auto e: ngel.Edges())
            for (auto j: Range(first_edge_dof[e], first_edge_dof[e + 1]))
                dnums.Append(j);

        // inner dofs
        if (ei.IsVolume())
            for (auto j: Range(first_cell_dof[ei.Nr()], first_cell_dof[ei.Nr() + 1]))
                dnums.Append(j);
    }


    FiniteElement &MyNedMatFESpace::GetFE(ElementId ei, Allocator &alloc) const {
        auto ngel = ma->GetElement(ei);
        switch (ngel.GetType()) {
            case ET_TET: {
                auto tet = new(alloc) MyNedElementMatTet(order);
                tet->SetVertexNumbers(ngel.vertices);
                return *tet;
            }
            default:
                throw Exception(string("Element type ") + ToString(ngel.GetType()) + " not supported");
        }
    }

}


void ExportNedMatFESpace(py::module m) {
    using namespace ngcomp;

    cout << "called ExportNedMatFESpace" << endl;

    ExportFESpace<MyNedMatFESpace>(m, "MyNedMatFESpace", true);
}