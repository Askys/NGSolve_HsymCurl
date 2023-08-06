#ifndef FILE_MYNEDELEMENT_HPP
#define FILE_MYNEDELEMENT_HPP


namespace ngfem {

    class MyNedElement : public FiniteElement {
    public:
        int p;
    public:
        MyNedElement(int order);

        virtual void CalcShape(const IntegrationPoint &ip,
                               BareSliceMatrix<> shape) const = 0;

        virtual void CalcDShape(const IntegrationPoint &ip,
                                BareSliceVector<> dshape) const = 0;

        virtual void CalcMappedShape(const BaseMappedIntegrationPoint &bmip,
                                     BareSliceMatrix<> shape) const = 0;

        virtual void CalcMappedDShape(const BaseMappedIntegrationPoint &bmip,
                                      BareSliceMatrix<> dshape) const = 0;

    };


    class MyNedElementTrig : public MyNedElement, public VertexOrientedFE<ET_TRIG> {
    public:
        enum {
            DIM = 2
        };

    public:
        using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;

        MyNedElementTrig(int order);

        virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

        virtual void CalcShape(const IntegrationPoint &ip,
                               BareSliceMatrix<> shape) const {
            throw Exception("MyNedElementTrig: can't evaluate from ip");
        }

        virtual void CalcDShape(const IntegrationPoint &ip,
                                BareSliceVector<> dshape) const {
            throw Exception("MyNedElementTrig: can't evaluate from ip");
        }

        virtual void CalcMappedShape(const BaseMappedIntegrationPoint &bmip,
                                     BareSliceMatrix<> shape) const;

        virtual void CalcMappedDShape(const BaseMappedIntegrationPoint &bmip,
                                      BareSliceMatrix<> dshape) const;
    };


    class MyNedElementTet : public MyNedElement, public VertexOrientedFE<ET_TET> {
    public:
        enum {
            DIM = 3
        };
        enum {
            DIM_STRESS = DIM
        };

    public:
        using VertexOrientedFE<ET_TET>::SetVertexNumbers;

        MyNedElementTet(int order);

        virtual ELEMENT_TYPE ElementType() const { return ET_TET; }

        // virtual void CalcShape (const IntegrationPoint & ip, 
        //                         BareSliceMatrix<> shape) const;

        // virtual void CalcDShape (const IntegrationPoint & ip, 
        //                          BareSliceVector<> dshape) const;

        virtual void CalcShape(const IntegrationPoint &ip,
                               BareSliceMatrix<> shape) const {
            throw Exception("MyNedElementTet: can't evaluate from ip");
        }

        virtual void CalcDShape(const IntegrationPoint &ip,
                                BareSliceVector<> dshape) const {
            throw Exception("MyNedElementTet: can't evaluate from ip");
        }

        virtual void CalcMappedShape(const BaseMappedIntegrationPoint &bmip,
                                     BareSliceMatrix<> shape) const;

        virtual void CalcMappedDShape(const BaseMappedIntegrationPoint &bmip,
                                      BareSliceMatrix<> dshape) const;
    };


    class MyNedElementMatTet : public MyNedElement, public VertexOrientedFE<ET_TET> {
    public:
        enum {
            DIM = 3
        };
        enum {
            DIM_STRESS = DIM * DIM
        };

    public:
        using VertexOrientedFE<ET_TET>::SetVertexNumbers;

        MyNedElementMatTet(int order);

        virtual ELEMENT_TYPE ElementType() const { return ET_TET; }

        virtual void CalcShape(const IntegrationPoint &ip,
                               BareSliceMatrix<> shape) const {
            throw Exception("MyNedElementMatTet: can't evaluate from ip");
        }

        virtual void CalcDShape(const IntegrationPoint &ip,
                                BareSliceVector<> dshape) const {
            throw Exception("MyNedElementMatTet: can't evaluate from ip");
        }

        virtual void CalcMappedShape(const BaseMappedIntegrationPoint &bmip,
                                     BareSliceMatrix<> shape) const;

        virtual void CalcMappedDShape(const BaseMappedIntegrationPoint &bmip,
                                      BareSliceMatrix<> dshape) const;
    };

}

#endif

