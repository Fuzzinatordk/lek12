#include <iostream>
#include "nr3.h"
#include "tridag.h"
string vec2string(VecDoub a)
{
    string s = "[";
    for (int i = 0; i < a.size(); i++)
    {
        s += to_string(a[i]);
        if (i != a.size() - 1)
            s += ", ";
    }
    s += "]";
    return s;
}
VecDoub operator*(double scale, VecDoub &b)
{
    VecDoub result(b.size());
    for (int i = 0; i < b.size(); i++)
    {
        result[i] = scale * b[i];
    }
    return result;
}
VecDoub operator+(VecDoub &a, VecDoub &b)
{
    VecDoub result(a.size());
    for (int i = 0; i < a.size(); i++)
    {
        result[i] = a[i] + b[i];
    }
    return result;
}
double F(double dy, double x, double y)
{
    return (2 * x + sin(dy) - std::cos(y));
}
double Fy(double dy, double x, double y)
{
    return std::sin(y);
}
double Fdy(double dy, double x, double y)
{
    return std::cos(dy);
}
VecDoub JCalcA(double alpha, double beta, VecDoub x, VecDoub y, int n, double h, bool firstRun)
{
    VecDoub triDagA(n);
    if (firstRun)
        triDagA[0] = 0;
    else
    {
        triDagA[0] = 0;
        for (int j = 1; j < n - 1; j++)
        {
            double dy = (y[j + 1] - y[j - 1]) / (2 * h);
            triDagA[j] = -1 - (h / 2) * Fdy(dy, x[j], y[j]);
        }
    }
    triDagA[n - 1] = -1 - (h / 2) * Fdy(x[n - 1], y[n - 1], (beta - y[n - 2]) / (2 * h));

    return triDagA;
}
VecDoub JCalcB(double alpha, double beta, VecDoub x, VecDoub y, int n, double h, bool firstRun)
{
    VecDoub triDagB(n);
    if (firstRun)
        triDagB[0] = 2 + std::pow(h, 2) * Fy((beta - alpha) / (2 * h), x[0], y[0]);
    else
    {
        triDagB[0] = 2 + std::pow(h, 2) * Fy((y[1] - alpha) / (2 * h), x[0], y[0]);
        for (int j = 1; j < n - 1; j++)
        {
            double dy = (y[j + 1] - y[j - 1]) / (2 * h);
            triDagB[j] = 2 + std::pow(h, 2) * Fy(dy, x[j], y[j]);
        }
    }
    triDagB[n - 1] = 2 + std::pow(h, 2) * Fy((beta - y[n - 2]) / 2 * h, x[n - 1], y[n - 1]);
    return triDagB;
}
VecDoub JCalcC(double alpha, double beta, VecDoub x, VecDoub y, int n, double h, bool firstRun)
{
    VecDoub triDagC(n);
    if (firstRun)
        triDagC[0] = 0;
    else
    {
        triDagC[0] = -1 + (h / 2) * Fdy((y[1] - alpha) / (2 * h), x[0], y[0]);
        for (int j = 1; j < n - 1; j++)
        {
            double dy = (y[j + 1] - y[j - 1]) / (2 * h);
            triDagC[j] = -1.0 + (h / 2.0) * Fdy(dy, x[j], y[j]);
        }
    }
    triDagC[n - 1] = 0;
    return triDagC;
}
VecDoub JCalcR(double alpha, double beta, VecDoub x, VecDoub y, int n, double h, bool firstRun)
{
    VecDoub triDagR(n);
    if (firstRun)
    {
        triDagR[0] = alpha - 2 * y[0] + beta - std::pow(h, 2) * F((beta - alpha) / (2 * h), x[0], y[0]);
        return triDagR;
    }
    triDagR[0] = alpha - 2 * y[0] + y[1] - std::pow(h, 2) * F((y[1] - alpha) / (2 * h), x[0], y[0]);
    for (int j = 1; j < n - 1; j++)
    {
        double dy = (y[j + 1] - y[j - 1]) / (2 * h);
        triDagR[j] = y[j - 1] - 2 * y[j] + y[j + 1] - std::pow(h, 2) * F(dy, x[j], y[j]);
    }
    triDagR[n - 1] = y[n - 2] - 2 * y[n - 1] + beta - std::pow(h, 2) * F((beta - y[n - 2]) / (2 * h), x[n - 1], y[n - 1]);
    return triDagR;
}
VecDoub yUpdater(double a, double b, double alpha, double beta, VecDoub y, VecDoub x, int n, double h)
{
    VecDoub yPrev = y;
    y.resize(n);
    if (y.size() == 1)
    {
        return yPrev;
    }
    for (Int i = 0; i < n; i++)
    {
        if (i == 0)
        {
            y[i] = (yPrev[0] + alpha) / 2;
        }
        else if (i == n - 1)
        {
            y[i] = (beta + yPrev[i / 2 - 1]) / 2;
        }
        else if (i % 2 == 0)
        {
            y[i] = (yPrev[i / 2] + yPrev[i / 2]) / 2;
        }
        else
        {
            y[i] = yPrev[i / 2];
        }
    }
    return y;
}

VecDoub xUpdater(double a, double b, double alpha, double beta, VecDoub y, VecDoub x, int n, double h)
{
    x.resize(n);
    for (Int i = 0; i < n; i++)
    {
        x[i] = a + (i + 1) * h;
    }
    return x;
}
void finiteDiff(double h, double alpha, double beta, double x0, double xEnd, VecDoub y, VecDoub x, double eps)
{
    VecDoub yOld(1), yOldOld(1);
    double xWanted = 0, xOldWanted = 0, xOldOldWanted = 0;
    yOld[0] = alpha, yOld[1] = beta;
    // Checking if its the first iteration that is being run, which is a condition for the calculation of yn and matrix elements
    bool firstRun = true;
    for (int i = 0; i <= 20; i++)
    {
        int N = (xEnd - x0) / h + 1;
        int n = N - 2;
        y = yUpdater(x0, xEnd, alpha, beta, y, x, n, h);
        x = xUpdater(x0, xEnd, alpha, beta, y, x, n, h);
        double error = 1;
        while (error > eps)
        {
            //  Calculating the tridiagonal matrix elements a from a2 to a_n
            VecDoub a = JCalcA(alpha, beta, x, y, n, h, firstRun);
            //  Calculating the tridiagonal matrix elements b from b1 to b_n-1
            VecDoub b = JCalcB(alpha, beta, x, y, n, h, firstRun);
            //  Calculating the tridiagonal matrix elements c from c1 to c_n-1
            VecDoub c = JCalcC(alpha, beta, x, y, n, h, firstRun);
            //  Calculating the right hand side of the equation
            VecDoub r = JCalcR(alpha, beta, x, y, n, h, firstRun);
            VecDoub u(n);
            // Solving the tridiagonal matrix
            tridag(a, b, c, r, u);
            // Adding the solution to the previous solution
            y = y + u;
            error = 0;
            for (int i = 0; i < n; i++)
            {
                error += std::pow(u[i], 2);
            }
            error = std::sqrt(error);
        }
        Doub alpha_k = 0, richardson = 0;
        if (i > 2)
        {
            alpha_k = (yOldOld[xOldOldWanted] - yOld[xOldWanted]) / (yOld[xOldWanted] - y[xWanted]);
            richardson = (yOld[xOldWanted] - y[xWanted]) / (alpha_k - 1);
        }
        

        std::cout << std::left << std::setw(15) << i;
        std::cout << std::left << std::setw(15) << y[xWanted];
        std::cout << std::left << std::setw(25) << yOld[xOldWanted] - y[xWanted];
        std::cout << std::left << std::setw(15) << alpha_k;
        std::cout << std::left << std::setw(15) << richardson << std::endl;
        yOldOld = yOld;
        yOld = y;
        xOldOldWanted = xOldWanted;
        xOldWanted = xWanted;
        xWanted = 2 * xWanted + 1;
        firstRun = false;
        h = h / 2;
    }
}
int main(int, char **)
{
    std::string spacer = "-";
    for (int i = 0; i < 95; i++)
        spacer += "-";
    std::cout << spacer << std::endl;
    double lastA = 0.0;
    std::cout << "solving y(1)" << ":" << std::endl;
    std::cout << std::left << std::setw(15) << "Iteration";
    std::cout << std::left << std::setw(15) << "A(hi)";
    std::cout << std::left << std::setw(25) << "A(hi-1) - A(hi)";
    std::cout << std::left << std::setw(15) << "alpha_k";
    std::cout << std::left << std::setw(15) << "richardson" << std::endl;
    std::cout << spacer << std::endl;
    double a = 0, b = 2;
    double alpha = 0, beta = 1;
    double h = (b - a) / 2;
    VecDoub y(1);
    VecDoub x(1);
    double eps = 1e-6;
    y[0] = (1.0 / 2.0) * (alpha + beta);
    x[0] = a + h;
    finiteDiff(h, alpha, beta, a, b, y, x, eps);
}
