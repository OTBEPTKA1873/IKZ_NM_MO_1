#include <iostream>
#include <cmath>

using namespace std;

double J(double x, int choice) // Функция
{
    switch (choice)
    {
    case 1:
        return pow(x, 2) - sin(x);
    case 2:
        return x * cos(2 * x) + 1;
    }
}

double HalfDivision(int Jchoice, int text, double A, double B, double EPS)
{
    double x1 = A, x2 = B, delta = EPS / 2; // Динамическая граница метода деления пополам
    cout << "You using Half Division method!\n\n";
    int iter = 0;
    if (text == 1)
    {
        cout << "Iter=0:  x1=" << x1 << "  x2=" << x2 << endl;
    }
    do
    {
        iter++;
        if (J((x1 + x2 - delta) / 2, Jchoice) <= J((x1 + x2 + delta) / 2, Jchoice))
        {
            x2 = (x1 + x2 + delta) / 2;
        }
        else
        {
            x1 = (x1 + x2 - delta) / 2;
        }
        if (text == 1)
        {
            cout << "Iter=" << iter << "  A=" << x1 << "  B=" << x2 << endl;
        }
    } while (abs(x1 - x2) > EPS);
    cout << "\n---------------------------\n\nMinimum=" << (x1 + x2) / 2 << "  Jmin=" << J((x1 + x2) / 2, Jchoice) << "  iter=" << iter << endl;
    return 1;
}

double GoldenRatio(int Jchoice, int text, double A, double B, double EPS)
{
    double x1 = A, x2 = B; // Динамическая граница метода деления пополам
    cout << "You using Golden ratio method!\n\n";
    int iter = 0;
    if (text == 1)
    {
        cout << "Iter=0:  x1=" << x1 << "  x2=" << x2 << endl;
    }
    do
    {
        iter++;
        double A = x1 + (3 - sqrt(5)) * (x2 - x1) / 2;
        double B = x1 + (sqrt(5) - 1) * (x2 - x1) / 2;
        if (J(A, Jchoice) <= J(B, Jchoice))
        {
            x2 = B;
        }
        else
        {
            x1 = A;
        }
        if (text == 1)
        {
            cout << "Iter=" << iter << "  A=" << x1 << "  B=" << x2 << endl;
        }
    } while (pow((sqrt(5) - 1) / 2, iter) * (x2 - x1) > EPS);
    cout << "\n---------------------------\n\nMinimum=" << (x1 + x2) / 2 << "  Jmin=" << J((x1 + x2) / 2, Jchoice) << "  iter=" << iter << endl;
    return 1;
}

double Fn(int n) // Нахождение числа Фибоначи через индукцию
{
    return (pow((1 + sqrt(5)) / 2, n) - pow((1 - sqrt(5)) / 2, n)) / sqrt(5);
}

double Fibonachi(int Jchoice, int text, double A, double B, double EPS)
{
    double an = A, bn = B, x1, x2, xn;
    int n = (int)(log(EPS / ((bn - an) * sqrt(5))) / log(2 / (sqrt(5) + 1)));
    cout << "You using Fibonachi method!\n\n";
    int iter = 0;
    if (text == 1)
    {
        cout << "Iter=0:  x1=" << an << "  x2=" << bn << endl;
    }
    for (int i = 1; i < n; i++)
    {
        iter++;
        double x1 = an + (bn - an) * Fn(n - i) / Fn(n - i + 2);
        double x2 = an + (bn - an) * Fn(n - i + 1) / Fn(n - i + 2);
        if (J(x1, Jchoice) <= J(x2, Jchoice))
        {
            bn = x2;
            xn = x1;
        }
        else
        {
            an = x1;
            xn = x2;
        }
        if (text == 1)
        {
            cout << "Iter=" << iter << "  A=" << an << "  B=" << bn << endl;
        }
    }
    double eps_fin = (B - A) / (10 * Fn(n + 1));
    if (J(xn - eps_fin, Jchoice) <= J(xn + eps_fin, Jchoice))
    {
        cout << "\n---------------------------\n\nMinimum=" << xn - eps_fin << "  Jmin=" << J(xn - eps_fin, Jchoice) << "  iter=" << iter << endl;
    }
    else
    {
        cout << "\n---------------------------\n\nMinimum=" << xn + eps_fin << "  Jmin=" << J(xn + eps_fin, Jchoice) << "  iter=" << iter << endl;
    }
    return 1;
}

double dJ(double x, int choice)
{
    switch (choice)
    {
    case 1:
        return 2 * x - cos(x);
    case 2:
        return cos(2 * x) - 2 * x * sin(2 * x);
    }
}

double L(int choice, double a, double b)
{
    double max = -pow(10, 10);
    double N = abs((b - a) / 1000);
    for (double i = a; i <= b; i+=N)
    {
        if (abs(dJ(i , choice)) > max)
        {
            max = abs(dJ(i, choice));
        }
    }
    return max;
}

void BrokenLine(int Jchoice, double x1, double p1, double& Lip, double EPS, double& Xmin, double& Jmin, int& Branch)
{
    double delta_n = 1. / (2 * Lip) * (J(x1, Jchoice) - p1);
    if (2 * Lip * delta_n > EPS)
    {
        p1 = (1. / 2) * (J(x1, Jchoice) + p1);
        BrokenLine(Jchoice, x1 - delta_n, p1, Lip, EPS, Xmin, Jmin, Branch);
        BrokenLine(Jchoice, x1 + delta_n, p1, Lip, EPS, Xmin, Jmin, Branch);
        Branch += 2;
    }
    else
    {
        if (Jmin > J(x1, Jchoice))
        {
            Jmin = J(x1, Jchoice);
            Xmin = x1;
        }
    }
}

int main()
{
    int Jchoice, Mchoice, text;
    cout << "You using programm for find minimum J1 or J2 different various\n\nPlease choose function J\n\n";
    cout << "1) x^2 + sin(x)\n2) x * cos(2x) + 1\n\nYour choice: ";
    cin >> Jchoice;
    cout << "\n---------------------------\n\n";
    cout << "Choose method:\n1) Half division\n2) Golden Ratio\n3) Fibonachi\n4) Broken line\n\nYour choice: ";
    cin >> Mchoice;
    cout << "\n---------------------------\n\n";
    if (Mchoice != 4)
    {
        cout << "Do you want see all calculations? 1-yes  2-no\n\nYour choice ";
        cin >> text;
        cout << "\n---------------------------\n\n";
    }
    double a, b, eps;
    cout << "Please enter left edge: ";
    cin >> a;
    cout << "\n---------------------------\n\n";
    cout << "Please enter right edge: ";
    cin >> b;
    cout << "\n---------------------------\n\n";
    cout << "Please enter eps: ";
    cin >> eps;
    cout << "\n---------------------------\n\n";
    switch (Mchoice)
    {
    case 1:
        HalfDivision(Jchoice, text, a, b, eps);
        break;
    case 2:
        GoldenRatio(Jchoice, text, a, b, eps);
        break;
    case 3:
        Fibonachi(Jchoice, text, a, b, eps);
        break;
    case 4:
        double Xmin = a, Jmin = J(a, Jchoice);
        double Lip = L(Jchoice, a, b);
        double x1 = 1. / (2 * Lip) * (J(a, Jchoice) - J(b, Jchoice) + Lip * (a + b));
        double p1 = 1. / 2 * (J(a, Jchoice) + J(b, Jchoice) + Lip * (a - b));
        int branch = 1;
        BrokenLine(Jchoice, x1, p1, Lip, eps, Xmin, Jmin, branch);
        cout << "Minimum=" << Xmin << "  Jmin=" << Jmin << "  Total branch=" << branch << endl;
    }
    return 100;
}