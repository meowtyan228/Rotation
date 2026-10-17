#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

Matrix3d buildK(Vector3d& n) {
    Matrix3d K;
    K << 0, -n.z(),  n.y(),
         n.z(),  0, -n.x(),
        -n.y(), n.x(),  0;
    return K;
}

VectorXd residual(double teta, Vector3d& n, VectorXd& x, VectorXd& y, int N) {
    VectorXd r(3 * N);
    Matrix3d K = buildK(n);
    Matrix3d R = Matrix3d::Identity() + sin(teta) * K + (1 - cos(teta)) * K * K;
    for (int i = 0; i < N; ++i) {
        Vector3d xi = x.segment<3>(3 * i);
        Vector3d yi = y.segment<3>(3 * i);
        r.segment<3>(3 * i) = yi - R * xi;
    }
    return r;
}

VectorXd J(double teta, Vector3d& n, VectorXd& x, int N) {
    VectorXd Jk(3 * N);
    Matrix3d K = buildK(n);
    for (int i = 0; i < N; ++i) {
        Vector3d xi = x.segment<3>(3 * i);
        Vector3d xpl = n * n.dot(xi);
        Vector3d xp = xi - xpl;
        Jk.segment<3>(3 * i) = xp * sin(teta) - (K * xi) * cos(teta);
    }
    return Jk;
}

struct Result {
    double theta;
    int iters;
    double S;
    bool converged;
};

Result Newton(VectorXd& x, VectorXd& y, Vector3d& n,
              double theta0, int N,
              const string& fname,
              double eps = 1e-8, int max_iter = 100) {

    ofstream log(fname);

    double theta = theta0;
    double S_prev = numeric_limits<double>::max();
    double S = 0;
    bool converged = false;
    int iters = 0;

    for (int k = 0; k < max_iter; ++k) {
        VectorXd r = residual(theta, n, x, y, N);
        VectorXd Jk = J(theta, n, x, N);

        double JtJ = Jk.dot(Jk);
        double Jtr = Jk.dot(r);
        S = r.dot(r);

        iters = k + 1;
        log << k << " " << S << endl;

        if (k > 0 && abs(S - S_prev) / (abs(S) + 1e-12) < eps) {
            converged = true;
            break;
        }
        S_prev = S;

        if (abs(JtJ) < 1e-12) break;

        double delta = -Jtr / JtJ;
        theta += delta;
    }

    log.close();
    return {theta, iters, S, converged};
}

void gen(VectorXd& x, VectorXd& y, Vector3d& n, int N) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> nd(0.0, 1.0);

    Vector3d a(2.0, 1.0, 0.5);

    for (int i = 0; i < N; ++i) {
        Vector3d p;
        do {
            p = Vector3d(nd(gen), nd(gen), nd(gen));
        } while (p.norm() < 1e-6);

        p.normalize();
        p = p.cwiseProduct(a);

        x.segment<3>(3 * i) = p;
    }

    double th = 0.7;
    double sig = 0.1;

    Matrix3d K = buildK(n);
    Matrix3d R = Matrix3d::Identity()
               + sin(th) * K
               + (1 - cos(th)) * K * K;

    for (int i = 0; i < N; ++i) {
        Vector3d xi = x.segment<3>(3 * i);
        y.segment<3>(3 * i) = R * xi + sig * Vector3d::Random();
    }
}

void save(const string& fname, const VectorXd& v, int N) {
    ofstream f(fname);
    for (int i = 0; i < N; ++i) {
        f << v(3*i) << " " << v(3*i+1) << " " << v(3*i+2) << "\n";
    }
}
int main() {
    int N = 2000;

    Vector3d n(1, 1, 1);
    n.normalize();

    VectorXd x(3 * N);
    VectorXd y(3 * N);
    gen(x, y, n, N);

    vector<double> theta0 = {0.3, 0.7, 1.5, 3.0, 5.0};

    cout << "\nConvergence comparison" << endl;
    cout << "theta0\t\ttheta*\t\tS*\t\titers\tconv" << endl;

    for (double theta0 : theta0) {
        string fname = "S_theta0_" + to_string(theta0) + ".txt";
        Result res = Newton(x, y, n, theta0, N, fname);
        cout << theta0 << "\t\t"
             << res.theta << "\t\t"
             << res.S << "\t\t"
             << res.iters << "\t"
             << (res.converged ? "yes" : "no") << endl;
    }

    double theta_est = theta0[1];

    Matrix3d K = buildK(n);
    Matrix3d R_est = Matrix3d::Identity()
                   + sin(theta_est) * K
                   + (1 - cos(theta_est)) * K * K;

    VectorXd y_fit(3 * N);
    for (int i = 0; i < N; ++i) {
        Vector3d xi = x.segment<3>(3 * i);
        y_fit.segment<3>(3 * i) = R_est * xi;
    }

    save("x.txt", x, N);
    save("y.txt", y, N);
    save("y_fit.txt", y_fit, N);

    return 0;
}
