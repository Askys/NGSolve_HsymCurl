#include <fem.hpp>
#include "hsymcurl.hpp"


namespace ngfem {

    // Cross product
    template<typename T>
    Vec<3, T> Cross(Vec<3, T> a, Vec<3, T> b) {
        Vec<3, T> m = {a(1) * b(2) - b(1) * a(2), a(2) * b(0) - b(2) * a(0), a(0) * b(1) - b(0) * a(1)};
        return m;
    }

    // dyadic product of two vectors
    template<int H, int W, typename T>
    Mat <H, W, T> Dyad(Vec <H, T> a, Vec <W, T> b) {
        Mat <H, W, T> m;
        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++)
                m(i, j) = a(i) * b(j);
        return m;
    }

    // Transposition operator
    template<int D, typename T>
    Mat <D, D, T> Trans(Mat <D, D, T> a) {
        Mat <D, D, T> m;
        for (int i = 0; i < D; i++)
            for (int j = 0; j < D; j++)
                m(i, j) = a(j, i);
        return m;
    }

    // Symmetry operator
    template<int D, typename T>
    Mat <D, D, T> Sym(Mat <D, D, T> a) {
        Mat <D, D, T> m = 0.5 * a + 0.5 * Trans(a);
        return m;
    }

    HsymCurlFE::HsymCurlFE(int order)
    // (order+1)*(order+2)/2-> #dof/trig
        //order+1 needed as full linear polynomials are used due to identity!!!!
            : FiniteElement(21, order+1) {
        if (order != 0)
            throw Exception("In HsymCurlFE: order must be 0");
    }

    template<class T>
    void HsymCurlFE::T_CalcShape(const T &x, const T &y, const T &z, const Mat<3, 3, T> &F_invt,
                                 BareSliceMatrix <T> shape) const {

        // Barycentric base functions
        AutoDiff<3> xx(x, 0);
        AutoDiff<3> yy(y, 1);
        AutoDiff<3> zz(z, 2);
        AutoDiff<3> lam[4] = {xx, yy, zz, 1 - xx - yy - zz};

        Vec<3> e1 = {1, 0, 0};
        Vec<3> e2 = {0, 1, 0};
        Vec<3> e3 = {0, 0, 1};
        Vec<3> e[3] = {e1, e2, e3};

        Mat<3, 3> I = {{1, 0, 0},
                       {0, 1, 0},
                       {0, 0, 1}};

        int c = 0;
        // Nedelec base functions
        for (int i = 0; i < 6; i++) {
            auto edge = GetVertexOrientedEdge(i);
            AutoDiff<3> ls = lam[edge[0]];
            AutoDiff<3> le = lam[edge[1]];
            Vec<3> grad_ls = F_invt * GetGradient(ls);
            Vec<3> grad_le = F_invt * GetGradient(le);
            Vec<3> ned = ls.Value() * grad_le - le.Value() * grad_ls;

            for (int j = 0; j < 3; j++) {
                shape.Row(3 * i + j).Range(0, DIM_STRESS) = Dyad(e[j], ned);
                c++;
            }
        }

        // Identity base functions
        shape.Row(c).Range(0, DIM_STRESS) = (2*x-1) * I;
        shape.Row(c + 1).Range(0, DIM_STRESS) = (2*y-1) * I;
        shape.Row(c + 2).Range(0, DIM_STRESS) = (2*z-1) * I;
    }

    void HsymCurlFE::CalcMappedShape(const BaseMappedIntegrationPoint &bmip,
                                     BareSliceMatrix<> shape) const {
        auto mip = static_cast< const MappedIntegrationPoint<3, 3> &>(bmip);
        IntegrationPoint ip = mip.IP();
        double x = mip(0);
        double y = mip(1);
        double z = mip(2);

        // Jacobian from mapping evaluated at integration point
        Mat<3, 3> F_invt = Trans(mip.GetJacobianInverse());

        T_CalcShape(x, y, z, F_invt, shape);
    }


    void HsymCurlFE::CalcMappedDShape(const MappedIntegrationPoint<3, 3> &mip,
                                      BareSliceMatrix<> dshape) const {

        Mat<3, 3> F_invt = Trans(mip.GetJacobianInverse());

        AutoDiff<3> x(mip(0), 0);
        AutoDiff<3> y(mip(1), 1);
        AutoDiff<3> z(mip(2), 2);

        Vec<3> e1 = {1, 0, 0};
        Vec<3> e2 = {0, 1, 0};
        Vec<3> e3 = {0, 0, 1};
        Vec<3> e[3] = {e1, e2, e3};

        Mat<3, 3> O = {{0, 0, 0},
                       {0, 0, 0},
                       {0, 0, 0}};

        AutoDiff<3> lam[4] = {x, y, z, 1 - x - y - z};

        int c = 0;
        // Nedelec base functions
        for (int i = 0; i < 6; i++) {
            auto edge = GetVertexOrientedEdge(i);
            AutoDiff<3> ls = lam[edge[0]];
            AutoDiff<3> le = lam[edge[1]];
            Vec<3> grad_ls = F_invt * GetGradient(ls);
            Vec<3> grad_le = F_invt * GetGradient(le);
            Vec<3> dned = 2 * Cross(grad_ls, grad_le);

            for (int j = 0; j < 3; j++) {
                dshape.Row(3 * i + j).Range(0, DIM_STRESS) = Sym(Dyad(e[j], dned));
                c++;
            }
        }

        // Identity base functions
        dshape.Row(c).Range(0, DIM_STRESS) = O;
        dshape.Row(c + 1).Range(0, DIM_STRESS) = O;
        dshape.Row(c + 2).Range(0, DIM_STRESS) = O;
    }
}
