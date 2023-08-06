#ifndef FILE_HsymCurl_HPP
#define FILE_HsymCurl_HPP


namespace ngfem {

    class HsymCurlFE : public FiniteElement, public VertexOrientedFE<ET_TET> {
    public:
        enum {
            DIM = 3
        };
        enum {
            DIM_STRESS = (DIM * DIM)
        };

    public:
        using VertexOrientedFE<ET_TET>::SetVertexNumbers;

        HsymCurlFE(int order);

        virtual ELEMENT_TYPE ElementType() const { return ET_TET; }

        virtual void CalcShape(const IntegrationPoint &ip,
                               BareSliceMatrix<> shape) const {
            throw Exception("HsymCurlFE: can't evaluate from ip");
        }

        virtual void CalcDShape(const IntegrationPoint &ip,
                                BareSliceVector<> dshape) const {
            throw Exception("HsymCurlFE: can't evaluate from ip");
        }

        virtual void CalcMappedShape(const BaseMappedIntegrationPoint &bmip,
                                     BareSliceMatrix<> shape) const;

        virtual void CalcMappedDShape(const MappedIntegrationPoint<3, 3> &mip,
                                      BareSliceMatrix<> dshape) const;

        template<class T>
        void
        T_CalcShape(const T &x, const T &y, const T &z, const Mat<3, 3, T> &F_invt, BareSliceMatrix <T> shape) const;
    };


}

#endif

