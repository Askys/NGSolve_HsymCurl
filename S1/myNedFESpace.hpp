#ifndef FILE_MYNEDFESPACE_HPP
#define FILE_MYNEDFESPACE_HPP


namespace ngcomp
{
  class MyNedFESpace : public FESpace
  {
    int order;
    Array<int> first_edge_dof;
//    Array<int> first_face_dof;
//    Array<int> first_cell_dof;
  public:
    /*
      constructor. 
      Arguments are the access to the mesh data structure,
      and the flags from the define command in the pde-file
    */
    MyNedFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);

    // a name for our new fe-space
    string GetClassName () const override
    {
      return "MyNedFESpace";
    }

    void Update() override;

    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;
  };
}    

void ExportNedFESpace(py::module m);

#endif
