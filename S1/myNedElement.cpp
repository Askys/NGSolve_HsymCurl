#include <fem.hpp>
#include "myNedElement.hpp"


namespace ngfem {

    // dyadic product of two vectors
    template<int H, int W, typename T>
    Mat <H, W, T> Dyad(Vec <H, T> a, Vec <W, T> b) {
        Mat <H, W, T> m;
        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++) {
                m(i, j) = a(i) * b(j);
            }
        }
        return m;
    }

    // symmetry operator
    template<int D, typename T>
    Mat <D, D, T> Sym(Mat <D, D, T> a) {
        Mat <D, D, T> m;
        for (int i = 0; i < D; i++) {
            for (int j = 0; j < D; j++) {
                m(i, j) = 0.5 * (a(i, j) + a(j, i));
            }
        }
        return m;
    }


    MyNedElement::MyNedElement(int order)
            : FiniteElement(3 * (order + 1), max(order, 1)), p(order) { ; }

    MyNedElementTrig::MyNedElementTrig(int order)
            : MyNedElement(order) {
        ndof = 3 * (order + 1);
    }


    void MyNedElementTrig::CalcMappedShape(const BaseMappedIntegrationPoint &bmip,
                                           BareSliceMatrix<> shape) const {
        auto mip = static_cast< const MappedIntegrationPoint<2, 2> &>(bmip);
        auto ip = mip.IP();
        AutoDiff<2> x(ip(0), 0);
        AutoDiff<2> y(ip(1), 1);

        Mat<2, 2> F_invt = Trans(mip.GetJacobianInverse());


        AutoDiff<2> lam[3] = {x, y, 1 - x - y};

        int ii = 0;

        for (int i = 0; i < 3; i++) {
            auto edge = GetVertexOrientedEdge(i);
            AutoDiff<2> ls = lam[edge[0]];
            AutoDiff<2> le = lam[edge[1]];
            Vec<2> grad_ls = GetGradient(ls);
            Vec<2> grad_le = GetGradient(le);

            Vec<2> tshape = ls.Value() * grad_le - le.Value() * grad_ls;
            shape.Row(ii++).Range(0, 2) = F_invt * tshape;
            if (p > 0) {
                shape.Row(ii++).Range(0, 2) = F_invt * GetGradient(ls * le);
            }
        }
    }

    void MyNedElementTrig::CalcMappedDShape(const BaseMappedIntegrationPoint &mip,
                                            BareSliceMatrix<> dshape) const {
        //TODO
    }


    MyNedElementTet::MyNedElementTet(int order)
            : MyNedElement(order) {
        ndof = 6 * (order + 1);
    }


    void MyNedElementTet::CalcMappedShape(const BaseMappedIntegrationPoint &bmip,
                                          BareSliceMatrix<> shape) const {
        auto mip = static_cast< const MappedIntegrationPoint<3, 3> &>(bmip);
        auto ip = mip.IP();
        AutoDiff<3> x(ip(0), 0);
        AutoDiff<3> y(ip(1), 1);
        AutoDiff<3> z(ip(2), 2);

        Mat<3, 3> F_invt = Trans(mip.GetJacobianInverse());


        AutoDiff<3> lam[4] = {x, y, z, 1 - x - y - z};

        int ii = 0;

        for (int i = 0; i < 6; i++) {
            auto edge = GetVertexOrientedEdge(i);
            AutoDiff<3> ls = lam[edge[0]];
            AutoDiff<3> le = lam[edge[1]];
            Vec<3> grad_ls = GetGradient(ls);
            Vec<3> grad_le = GetGradient(le);

            Vec<3> tshape = ls.Value() * grad_le - le.Value() * grad_ls;
            shape.Row(ii++).Range(0, 3) = F_invt * tshape;
            if (p > 0) {
                shape.Row(ii++).Range(0, 3) = F_invt * GetGradient(ls * le);
            }
        }
    }


    void MyNedElementTet::CalcMappedDShape(const BaseMappedIntegrationPoint &bmip,
                                           BareSliceMatrix<> dshape) const {
        auto mip = static_cast< const MappedIntegrationPoint<3, 3> &>(bmip);
        auto ip = mip.IP();
        AutoDiff<3> x(ip(0), 0);
        AutoDiff<3> y(ip(1), 1);
        AutoDiff<3> z(ip(2), 2);

        Mat<3, 3> F = mip.GetJacobian();
        auto det = 1.0 / mip.GetJacobiDet();

        AutoDiff<3> lam[4] = {x, y, z, 1 - x - y - z};

        int ii = 0;

        for (int i = 0; i < 6; i++) {
            auto edge = GetVertexOrientedEdge(i);
            AutoDiff<3> ls = lam[edge[0]];
            AutoDiff<3> le = lam[edge[1]];
            Vec<3> grad_ls = GetGradient(ls);
            Vec<3> grad_le = GetGradient(le);

            Vec<3> tshape = Cross(grad_ls, grad_le) - Cross(grad_le, grad_ls);
            dshape.Row(ii++).Range(0, 3) = det * F * tshape;
            if (p > 0) {
                dshape.Row(ii++).Range(0, 3) = 0;
            }
        }
    }


    MyNedElementMatTet::MyNedElementMatTet(int order)
            : MyNedElement(order + 2) {
        ndof = 36 + 6 * 2 + 4;
    }


    void MyNedElementMatTet::CalcMappedShape(const BaseMappedIntegrationPoint &bmip,
                                             BareSliceMatrix<> shape) const {
        auto mip = static_cast< const MappedIntegrationPoint<3, 3> &>(bmip);
        auto ip = mip.IP();
        AutoDiff<3> x(ip(0), 0);
        AutoDiff<3> y(ip(1), 1);
        AutoDiff<3> z(ip(2), 2);

        Mat<3, 3> F_invt = Trans(mip.GetJacobianInverse());


        Vec<3> e1 = {1, 0, 0};
        Vec<3> e2 = {0, 1, 0};
        Vec<3> e3 = {0, 0, 1};
        Vec<3> e[3] = {e1, e2, e3};

        Mat<3, 3> I = {{1, 0, 0},
                       {0, 1, 0},
                       {0, 0, 1}};

        AutoDiff<3> lam[4] = {x, y, z, 1 - x - y - z};

        int ii = 0;

        ArrayMem<AutoDiff<3>, 20> polx(4);

        // edge base functions
        for (int i = 0; i < 6; i++) {
            auto edge = GetVertexOrientedEdge(i);
            AutoDiff<3> ls = lam[edge[0]];
            AutoDiff<3> le = lam[edge[1]];
            Vec<3> grad_ls = GetGradient(ls);
            Vec<3> grad_le = GetGradient(le);
            ScaledIntegratedLegendrePolynomial(2, le - ls, le + ls, polx);

            Vec<3> tshape = ls.Value() * grad_le - le.Value() * grad_ls;
            for (int j = 0; j < 3; j++) {
                shape.Row(6 * i + j).Range(0, DIM_STRESS) = Dyad(e[j], F_invt * tshape).AsVector();
                shape.Row(6 * i + 3 + j).Range(0, DIM_STRESS) = Dyad(e[j], F_invt * GetGradient(polx[2])).AsVector();
                ii += 2;
            }
        }

        // cell edge base functions
        for (int i = 0; i < 6; i++, ii += 2) {
            auto edge = GetVertexOrientedEdge(i);
            AutoDiff<3> ls = lam[edge[0]];
            AutoDiff<3> le = lam[edge[1]];
            ScaledIntegratedLegendrePolynomial(3, le - ls, le + ls, polx);
            shape.Row(ii).Range(0, DIM_STRESS) = (polx[2].Value() * I).AsVector();
            shape.Row(ii + 1).Range(0, DIM_STRESS) = (polx[3].Value() * I).AsVector();
        }

        // cell face base functions
        for (int i = 0; i < 4; i++, ii++) {
            auto face = GetVertexOrientedFace(i);
            AutoDiff<3> l0 = lam[face[0]];
            AutoDiff<3> l1 = lam[face[1]];
            AutoDiff<3> l2 = lam[face[2]];
            shape.Row(ii).Range(0, DIM_STRESS) = (l0.Value() * l1.Value() * l2.Value() * I).AsVector();
        }
    }


    void MyNedElementMatTet::CalcMappedDShape(const BaseMappedIntegrationPoint &bmip,
                                              BareSliceMatrix<> dshape) const {
        auto mip = static_cast< const MappedIntegrationPoint<3, 3> &>(bmip);
        auto ip = mip.IP();
        AutoDiff<3> x(ip(0), 0);
        AutoDiff<3> y(ip(1), 1);
        AutoDiff<3> z(ip(2), 2);

        Mat<3, 3> F = mip.GetJacobian();
        auto det = 1.0 / mip.GetJacobiDet();

        AutoDiff<3> lam[4] = {x, y, z, 1 - x - y - z};
        Vec<3> e1 = {1, 0, 0};
        Vec<3> e2 = {0, 1, 0};
        Vec<3> e3 = {0, 0, 1};
        Vec<3> e[3] = {e1, e2, e3};

        int ii = 0;

        // edge base functions
        for (int i = 0; i < 6; i++) {
            auto edge = GetVertexOrientedEdge(i);
            AutoDiff<3> ls = lam[edge[0]];
            AutoDiff<3> le = lam[edge[1]];
            Vec<3> grad_ls = GetGradient(ls);
            Vec<3> grad_le = GetGradient(le);

            Vec<3> tshape = Cross(grad_ls, grad_le) - Cross(grad_le, grad_ls);
            for (int j = 0; j < 3; j++) {
                dshape.Row(6 * i + j).Range(0, DIM_STRESS) = Sym(Dyad(e[j], det * F * tshape)).AsVector();
                dshape.Row(6 * i + 3 + j).Range(0, DIM_STRESS) = 0;
                ii += 2;
            }
        }

        // identity base functions
        for (int i = 0; i < 16; i++, ii++) {
            dshape.Row(ii).Range(0, DIM_STRESS) = 0;
        }
    }
}
