#ifndef FILE_MYDIFFOP_HPP
#define FILE_MYDIFFOP_HPP

#include <fem.hpp>

#include "hsymcurl.hpp"

namespace ngfem
{
    /*
       DiffOps provide the link between function evaluation, and finite elements.
       Different DiffOps evaluate either shape functions, or derivatives.
       Typically, a FiniteElementSpace creates a DiffOp to for function evaluation,
       and for the canonical derivative, as well as DiffOps for evaluation at the
       boundary. These DiffOps are used when evaluating GridFunctions, and setting
       up element-matrices from trial- and test-functions.
       DiffOps use static polymorphism, aka Curiously Recurring Template Pattern (CRTP).
     */

    //

    class DiffOpIdHsymCurl : public DiffOp<DiffOpIdHsymCurl>
    {
    public:
        // some constants for the diffop:

        enum { DIM = 1 };
        enum { DIM_SPACE = 3 };
        enum { DIM_ELEMENT = 3 };
        enum { DIM_DMAT = 3*3 };
        enum { DIFFORDER = 0 };
        enum { DIM_STRESS = 3*3 };

//     static string Name() { return "id"; }

        static Array<int> GetDimensions() { return Array<int> ({3,3}); }
        
        //static const MyNedElement &Cast(const FiniteElement &fel) { return static_cast<const MyNedElement &> (fel); }

        // fill evaluation matrix of dimension 1 times fel.ndof
        // the input mip is a mapped integration point, which has also
        // access to the integration point on the reference element
         template<typename MIP, typename MAT>
         static void GenerateMatrix(const FiniteElement & fel, const MIP & mip,
                                    MAT & mat, LocalHeap & lh)
         {
             dynamic_cast<const HsymCurlFE&> (fel).CalcMappedShape(mip, Trans(mat));
         }

//        template<typename MIP, typename MAT>
//        static void GenerateMatrix(const FiniteElement &fel, const MIP &mip,
//                                   MAT &mat, LocalHeap &lh) {
//            HeapReset hr(lh);
//
//          //shape \in R^(ndof x D)
//            FlatMatrix<double> shape(fel.GetNDof(), 3 * 3, lh);
//
//            Cast(fel).CalcShape(mip, shape);
//            //shape^T \in R^(D x ndof)
//            // F^{-T} shape^T \in  R^(D x ndof)
//            mat = Trans(shape);
//
//        }
//
//        // can overload more functionality for performance optimization,
        // like evaluation in the whole integration rule
    };

    class DiffOpsymCurlHsymCurl : public DiffOp<DiffOpsymCurlHsymCurl>
    {
    public:
        // some constants for the diffop:
        enum { DIM = 1 };
        enum { DIM_SPACE = 3 };
        enum { DIM_ELEMENT = 3 };
        enum { DIM_DMAT = 3*3 };
        enum { DIFFORDER = 1 };
        enum { DIM_STRESS = 3*3 };

        static string Name() { return "curl"; }
        
        static Array<int> GetDimensions() { return Array<int> ({3,3}); }
        
        //static const MyNedElement &Cast(const FiniteElement &fel) { return static_cast<const MyNedElement &> (fel); }

        //static Array<int> GetDimensions() { return Array<int> ({2,1}); }

        // fill evaluation matrix of dimension 1 times fel.ndof
        // the input mip is a mapped integration point, which has also
        // access to the integration point on the reference element
         template<typename MIP, typename MAT>
         static void GenerateMatrix(const FiniteElement & fel, const MIP & mip,
                                    MAT & mat, LocalHeap & lh)
         {
             dynamic_cast<const HsymCurlFE&> (fel).CalcMappedDShape(mip, Trans(mat));
         }
        
//                template<typename MIP, typename MAT>
//        static void GenerateMatrix(const FiniteElement &fel, const MIP &mip,
//                                   MAT &mat, LocalHeap &lh) {
//            HeapReset hr(lh);
//            FlatMatrix<double> dshape(fel.GetNDof(), 3 * 3, lh);
//
//            Cast(fel).CalcDShape(mip, dshape);
//             mat.Row(0) = 1 / mip.GetJacobiDet() * dshape;
//            mat = Trans(dshape);
//
//        }

        // can overload more functionality for performance optimization,
        // like evaluation in the whole integration rule
    };
}
#endif // FILE_MYDIFFOP_HPP
