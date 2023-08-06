#ifndef FILE_MYDIFFOPHC_HPP
#define FILE_MYDIFFOPHC_HPP

#include <fem.hpp>

#include "myNedElement.hpp"

namespace ngfem
{
  // Our implementation of the identity operator, in two space dimensions.
  template <int D>
  class MyNedDiffOpId : public DiffOp<MyNedDiffOpId<D>>
  {
  public:
    // some constants for the diffop:
    
    static constexpr int DIM = 1;       // dimension of the input
    static constexpr int DIM_SPACE = D; // dimension of the space
    static constexpr int DIM_ELEMENT = D; // spatial dimension of the element
    static constexpr int DIM_DMAT = D;  // dimension of the output
    static constexpr int DIFFORDER = 0; // order of differentiation

    static string Name() { return "id"; }

    static const MyNedElement& Cast(const FiniteElement & fel)
    { return static_cast<const MyNedElement&> (fel); }
    

    // fill evaluation matrix of dimension 1 times fel.ndof 
    // the input mip is a mapped integration point, which has also
    // access to the integration point on the reference element
    template<typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement & fel, const MIP & mip,
                               MAT & mat, LocalHeap & lh)
    {
      // HeapReset hr(lh);
      
      // //shape \in R^(ndof x D)
      // FlatMatrix<double> shape(fel.GetNDof(), D, lh);
      
      // Cast(fel).CalcShape(mip.IP(), shape);
      // //shape^T \in R^(D x ndof)
      // // F^{-T} shape^T \in  R^(D x ndof)
      // mat = Trans(mip.GetJacobianInverse())*Trans(shape);

      Cast(fel).CalcMappedShape(mip, Trans(mat));
   
    }

    // can overload more functionality for performance optimization,
    // like evaluation in the whole integration rule
  };
    
  template <int D>
  class MyNedDiffOpCurl : public DiffOp<MyNedDiffOpCurl<D>>
  {
  public:
    // some constants for the diffop:
    
    static constexpr int DIM = 1;       // dimension of the input
    static constexpr int DIM_SPACE = D; // dimension of the space
    static constexpr int DIM_ELEMENT = D; // spatial dimension of the element
    static constexpr int DIM_DMAT = 2*D-3;  // dimension of the output
    static constexpr int DIFFORDER = 1; // order of differentiation


    static string Name() { return "curl"; }

    static const MyNedElement& Cast(const FiniteElement & fel)
    { return static_cast<const MyNedElement&> (fel); }
    

    // fill evaluation matrix of dimension 1 times fel.ndof 
    // the input mip is a mapped integration point, which has also
    // access to the integration point on the reference element
    template<typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement & fel, const MIP & mip,
                               MAT & mat, LocalHeap & lh)
    {
      // HeapReset hr(lh);
      // FlatVector<double> dshape(fel.GetNDof(), lh);
      
      // Cast(fel).CalcDShape(mip.IP(), dshape);
      // mat.Row(0) = 1/mip.GetJacobiDet()*dshape;
   
      Cast(fel).CalcMappedDShape(mip, Trans(mat));
    }
  };


  template <int D>
  class MyNedMatDiffOpId : public DiffOp<MyNedMatDiffOpId<D>>
  {
  public:
    // some constants for the diffop:
    
    static constexpr int DIM = 1;       // dimension of the input
    static constexpr int DIM_SPACE = D; // dimension of the space
    static constexpr int DIM_ELEMENT = D; // spatial dimension of the element
    static constexpr int DIM_DMAT = D*D;  // dimension of the output
    static constexpr int DIFFORDER = 0; // order of differentiation

    static string Name() { return "id"; }

    static Array<int> GetDimensions() { return Array<int> ({D,D}); }

    static const MyNedElement& Cast(const FiniteElement & fel)
    { return static_cast<const MyNedElement&> (fel); }
    

    // fill evaluation matrix of dimension 1 times fel.ndof 
    // the input mip is a mapped integration point, which has also
    // access to the integration point on the reference element
    template<typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement & fel, const MIP & mip,
                               MAT & mat, LocalHeap & lh)
    {

      Cast(fel).CalcMappedShape(mip, Trans(mat));
   
    }

    // can overload more functionality for performance optimization,
    // like evaluation in the whole integration rule
  };


  template <int D>
  class MyNedMatDiffOpCurl : public DiffOp<MyNedMatDiffOpCurl<D>>
  {
  public:
    // some constants for the diffop:
    
    static constexpr int DIM = 1;       // dimension of the input
    static constexpr int DIM_SPACE = D; // dimension of the space
    static constexpr int DIM_ELEMENT = D; // spatial dimension of the element
    static constexpr int DIM_DMAT = (2*D-3)*(2*D-3);  // dimension of the output
    static constexpr int DIFFORDER = 1; // order of differentiation


    static string Name() { return "curl"; }

    static Array<int> GetDimensions() { return Array<int> ({2*D-3,2*D-3}); }

    static const MyNedElement& Cast(const FiniteElement & fel)
    { return static_cast<const MyNedElement&> (fel); }
    

    // fill evaluation matrix of dimension 1 times fel.ndof 
    // the input mip is a mapped integration point, which has also
    // access to the integration point on the reference element
    template<typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement & fel, const MIP & mip,
                               MAT & mat, LocalHeap & lh)
    {
      Cast(fel).CalcMappedDShape(mip, Trans(mat));
    }
  };
    
}
#endif // FILE_MYDIFFOP_HPP

