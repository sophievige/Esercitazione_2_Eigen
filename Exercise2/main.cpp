#include <iostream>
#include "Eigen/Eigen"
using namespace std;
using namespace Eigen;

int main()
{
    Vector2d x = Vector2d::Constant(2,-1);
    Matrix2d A1;
    A1 << 5.547001962252291e-01, -3.770900990025203e-02,
        8.320502943378437e-01, -9.992887623566787e-01;
    Vector2d b1;
    b1 << -5.169911863249772e-01,
        1.672384680188350e-01;
    Vector2d x1_lu;
    x1_lu = A1.partialPivLu().solve(b1);
    Vector2d x1_qr;
    x1_qr = A1.householderQr().solve(b1);
    double rel_err1_lu, rel_err1_qr;
    rel_err1_lu =(x1_lu-x).norm()/x.norm();
    rel_err1_qr=(x1_qr - x).norm()/x.norm();
    cout<< "LU primo sistema:" << rel_err1_lu <<endl;
    cout << "QR primo sistema:" << rel_err1_qr << endl;

    Matrix2d A2;
    A2 << 5.547001962252291e-01, -5.540607316466765e-01,
        8.320502943378437e-01, -8.324762492991313e-01;
    Vector2d b2;
    b2 << -6.394645785530173e-04,
        4.259549612877223e-04;
    Vector2d x2_lu;
    x2_lu=A2.partialPivLu().solve(b2);
    Vector2d x2_qr;
    x2_qr = A2.householderQr().solve(b2);
    double rel_err2_lu, rel_err2_qr;
    rel_err2_lu = (x2_lu-x).norm()/x.norm();
    rel_err2_qr = (x2_qr - x).norm()/x.norm();
    cout<< "LU secondo sistema:" << rel_err2_lu <<endl;
    cout << "QR secondo sistema:" << rel_err2_qr << endl;

    Matrix2d A3;
    A3 << 5.547001962252291e-01, -5.547001955851905e-01,
        8.320502943378437e-01, -8.320502947645361e-01;
    Vector2d b3;
    b3 << -6.400391328043042e-10,
        4.266924591433963e-10;
    Vector2d x3_lu;
    x3_lu = A3.partialPivLu().solve(b3);
    Vector2d x3_qr;
    x3_qr = A3.householderQr().solve(b3);
    double rel_err3_lu, rel_err3_qr;
    rel_err3_lu = (x3_lu-x).norm()/x.norm();
    rel_err3_qr = (x3_qr - x).norm()/x.norm();
    cout<< "LU terzo sistema:" << rel_err3_lu <<endl;
    cout << "QR terzo sistema:" << rel_err3_qr << endl;

    return 0;
}

