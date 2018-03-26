#include <src/jaz/svd_helper.h>
#include <src/jaz/index_sort.h>

void SvdHelper::decompose(
        const Matrix2D<RFLOAT>& A,
        Matrix2D<RFLOAT>& U,
        Matrix1D<RFLOAT>& S,
        Matrix2D<RFLOAT>& Vt)
{
    Matrix2D<RFLOAT> U0, Vt0;
    Matrix1D<RFLOAT> S0;

    svdcmp(A, U0, S0, Vt0);

    const int rc = A.mdimy;
    const int cc = A.mdimx;

    std::vector<RFLOAT> Svec(cc);

    for (int i = 0; i < cc; i++)
    {
        Svec[i] = S0(i);
    }

    std::vector<int> order = IndexSort<RFLOAT>::sortIndices(Svec);

    U = Matrix2D<RFLOAT>(rc,cc);
    S = Matrix1D<RFLOAT>(cc);
    Vt = Matrix2D<RFLOAT>(cc,cc);

    for (int i = 0; i < cc; i++)
    {
        const int j = order[cc - i - 1];

        for (int c = 0; c < cc; c++)
        {
            Vt(c,i) = Vt0(c,j);
        }

        S(i) = S0(j);

        for (int r = 0; r < rc; r++)
        {
            U(r,i) = U0(r,j);
        }
    }
}
