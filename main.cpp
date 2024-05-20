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
VecDoub operator*(double scale, VecDoub b)
{
    VecDoub result(b.size());
    for (int i = 0; i < b.size(); i++)
    {
        result[i] = scale * b[i];
    }
    return result;
}
VecDoub operator+(VecDoub a, VecDoub b)
{
    VecDoub result(a.size());
    for (int i = 0; i < a.size(); i++)
    {
        result[i] = a[i] + b[i];
    }
    return result;
}
double F(double dy, double x, Doub y)
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
    {
        triDagA[0] = 0;
    }
    else
    {
        triDagA[0] = -1 - h / 2 * Fdy((y[1] - alpha) / (2 * h), x[0], y[0]);
        for (int j = 1; j < n - 1; j++)
        {
            double dy = (y[j + 1] - y[j - 1]) / (2 * h);
            triDagA[j] = -1 - h / 2 * Fdy(dy, x[j], y[j]);
        }
        triDagA[n - 1] = 2 + std::pow(h, 2) * Fy((beta - y[n - 2]) / (2 * h), x[n - 1], y[n - 1]);
    }
    return triDagA;
}
VecDoub JCalcB(double alpha, double beta, VecDoub x, VecDoub y, int n, double h, bool firstRun)
{
    VecDoub triDagB(n);
    if (firstRun)
    {
        triDagB[0] = 2 + std::pow(h, 2) * Fy((beta - alpha) / (2 * h), x[0], y[0]);
    }
    else
    {
        triDagB[0] = 2 + std::pow(h, 2) * Fy((y[1] - alpha) / (2 * h), x[0], y[0]);
        for (int j = 1; j < n - 1; j++)
        {
            double dy = (y[j + 1] - y[j - 1]) / (2 * h);
            triDagB[j] = 2 + std::pow(h, 2) * Fy(dy, x[j], y[j]);
        }
    }
    triDagB[n - 1] = 2 + std::pow(h, 2) * Fy((beta - y[n - 2]) / (2 * h), x[n - 1], y[n - 1]);
    return triDagB;
}
VecDoub JCalcC(double alpha, double beta, VecDoub x, VecDoub y, int n, double h, bool firstRun)
{
    VecDoub triDagC(n);
    if (firstRun)
    {
        triDagC[0] = -1 + h / 2 * Fdy((beta - alpha) / (2 * h), x[0], y[0]);
    }
    else
    {
        triDagC[0] = -1 + h / 2 * Fdy((y[1] - alpha) / (2 * h), x[0], y[0]);
        for (int j = 1; j < n - 1; j++)
        {
            double dy = (y[j + 1] - y[j - 1]) / (2 * h);
            triDagC[j] = -1 + h / 2 * Fdy(dy, x[j], y[j]);
        }
    }
    triDagC[n - 1] = 0;
    return triDagC;
}
VecDoub JCalcR(double alpha, double beta, VecDoub x, VecDoub y, int n, double h, bool firstRun)
{
    VecDoub triDagR(n);
    triDagR[0] = y[0];
    triDagR[1] = -alpha + 2 * y[1] - y[2] + h * h * F((y[2] - alpha) / (2 * h), x[1], y[1]);
    for (int j = 2; j < n - 2; j++)
    {
        double dy = (y[j + 1] - y[j - 1]) / (2 * h);
        triDagR[j] = -y[j - 1] + 2 * y[j] - y[j + 1] + h * h * F(dy, x[j], y[j]);
    }
    triDagR[n - 2] = -y[n - 2] + 2 * y[n - 1] - beta + h * h * F((beta - y[n - 2]) / (2 * h), x[n - 1], y[n - 1]);
    triDagR[n - 1] = y[n];
    return triDagR;
}
void finiteDiff(double alpha, double beta, double x0, double xEnd)
{
    VecDoub yOld(2), yOldOld(2);
    yOld[0] = alpha, yOld[1] = beta;
    // Checking if its the first iteration that is being run, which is a condition for the calculation of yn and matrix elements
    bool firstRun = true;
    for (int i = 1; i <= 20; i++)
    {
        int n = pow(2, i);
        VecDoub y(n);
        y[0] = alpha, y[n - 1] = beta;
        Doub h = (xEnd - x0) / n;
        if (firstRun)
        {
            for (int j = 1; j < n - 1; j++)
            {
                y[j] = alpha + j/n * (beta - alpha);
            }
        }
        else
        {
            for (int j = 0; j < n + 1; j++)
            {
                if (j % 2 == 0)
                {
                    y[j] = yOld[j / 2];
                }
                else
                    y[j] = (yOld[(j - 1) / 2] + yOld[(j + 1) / 2]) / 2;
            }
        }
        VecDoub x(n);
        for (int j = 0; j < n; j++)
        {
            x[j] = x0 + j * h;
        }
        // Calculating the tridiagonal matrix elements a from a1 to a_n-1
        VecDoub a = JCalcA(alpha, beta, x, y, n, h, firstRun);
        // Calculating the tridiagonal matrix elements b from b1 to b_n-1
        VecDoub b = JCalcB(alpha, beta, x, y, n, h, firstRun);
        // Calculating the tridiagonal matrix elements c from c1 to c_n-1
        VecDoub c = JCalcC(alpha, beta, x, y, n, h, firstRun);
        // Calculating the right hand side of the equation
        VecDoub r = JCalcR(alpha, beta, x, y, n, h, firstRun);
        VecDoub u(n);
        r = -1 * r;
        // Solving the tridiagonal matrix
        tridag(a, b, c, r, u);
        // Adding the solution to the previous solution
        y = y + u;
        Doub alpha_k = 0, richardson = 0;
        if (i >= 3)
        {
            alpha_k = (yOldOld[1] - yOld[1]) / (yOld[1] - y[1]);
            richardson = (yOld[1] - y[1]) / (alpha_k - 1);
        }

        yOldOld = yOld;
        yOld = y;
        std::cout << i << " " << y[n/2] << " " << alpha_k << " " << richardson << std::endl;
        firstRun = false;
    }
}
int main(int, char **)
{
    finiteDiff(0, 1, 0, 2);


}
