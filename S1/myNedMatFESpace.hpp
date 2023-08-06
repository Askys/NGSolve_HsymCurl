#ifndef FILE_MYNEDMATFESPACE_HPP
#define FILE_MYNEDMATFESPACE_HPP


namespace ngcomp
{
  class MyNedMatFESpace : public FESpace
  {
    int order;
    Array<int> first_edge_dof;
//    Array<int> first_face_dof;
   Array<int> first_cell_dof;
  public:
    /*
      constructor. 
      Arguments are the access to the mesh data structure,
      and the flags from the define command in the pde-file
    */
    MyNedMatFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);

    // a name for our new fe-space
    string GetClassName () const override
    {
      return "MyNedMatFESpace";
    }

    void Update() override;

    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;
  };
}  

void ExportNedMatFESpace(py::module m);

#endif
