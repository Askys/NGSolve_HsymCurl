#ifndef FILE_HsymCurlFESPACE_HPP
#define FILE_HsymCurlFESPACE_HPP


namespace ngcomp
{
    class HsymCurlFESpace : public FESpace
    {
        int order;

        Array<int> first_edge_dof;
        Array<int> first_face_dof;
        Array<int> first_cell_dof;
    public:
        /*
          constructor. 
          Arguments are the access to the mesh data structure,
          and the flags from the define command in the pde-file
        */
        HsymCurlFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);

        // a name for our new fe-space
        string GetClassName () const override
        {
            return "HsymCurlFESpace";
        }

        void Update() override;

        void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
        FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;
    };
}

void ExportHsymCurlFESpace(py::module m);
#endif
