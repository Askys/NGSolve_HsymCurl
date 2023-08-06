#include <comp.hpp>
#include "hsymcurl.hpp"
#include "hsymcurlfespace.hpp"
#include "myDiffOp.hpp"


namespace ngcomp {

    HsymCurlFESpace::HsymCurlFESpace(shared_ptr <MeshAccess> ama, const Flags &flags)
            : FESpace(ama, flags) {
        type = "HsymCurlFESpace";


        // default = 0
        order = int(flags.GetNumFlag("order", 0));
        evaluator[VOL] = make_shared < T_DifferentialOperator < DiffOpIdHsymCurl >> ();
        flux_evaluator[VOL] = make_shared < T_DifferentialOperator < DiffOpsymCurlHsymCurl >> ();
    }

    void HsymCurlFESpace::Update() {
        int n_edge = ma->GetNEdges();
        int n_cell = ma->GetNE();

        first_edge_dof.SetSize(n_edge + 1);
        int ii = 0;
        for (int i = 0; i < n_edge; i++, ii += 3)
            first_edge_dof[i] = ii;
        first_edge_dof[n_edge] = ii;

        first_cell_dof.SetSize(n_cell + 1);
        for (int i = 0; i < n_cell; i++, ii += 9)
            first_cell_dof[i] = ii;
        first_cell_dof[n_cell] = ii;

        SetNDof(ii);
    }

    void HsymCurlFESpace::GetDofNrs(ElementId ei, Array <DofId> &dnums) const {
        // returns dofs of element number elnr
        dnums.SetSize(0);
        auto ngel = ma->GetElement(ei);

        // edge dofs
        for (auto e: ngel.Edges())
            for (auto j: Range(first_edge_dof[e], first_edge_dof[e + 1]))
                dnums.Append(j);
        
         // face dofs
        //for (auto f: ngel.Faces())
        //    for (auto j: Range(first_face_dof[f], first_face_dof[f + 1]))
        //        dnums.Append(j);

        // inner dofs
        if (ei.IsVolume())
            for (auto j: Range(first_cell_dof[ei.Nr()], first_cell_dof[ei.Nr() + 1]))
                dnums.Append(j);
    }


    FiniteElement &HsymCurlFESpace::GetFE(ElementId ei, Allocator &alloc) const {
        auto ngel = ma->GetElement(ei);
        switch (ngel.GetType()) {
            case ET_TET: {
                auto trig = new(alloc) HsymCurlFE(order);
                trig->SetVertexNumbers(ngel.vertices);
                return *trig;
            }
            default:
                throw Exception(string("Element type ") + ToString(ngel.GetType()) + " not supported");
        }
    }


}

void ExportHsymCurlFESpace(py::module m) {
    using namespace ngcomp;

    cout << "called ExportHsymCurlFESpace" << endl;

    ExportFESpace<HsymCurlFESpace>(m, "HsymCurlFESpace", true);
}
