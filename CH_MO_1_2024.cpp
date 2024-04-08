#include <iostream>
#include <cmath>

using namespace std;

// Функции
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

// Первая производная функции
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

// Вторая производная функции
double ddJ(double x, int choice)
{
    switch (choice)
    {
    case 1:
        return 2 + sin(x);
    case 2:
        return -4 * sin(2 * x) - 4 * x * cos(2 * x);
    }
}

// Нахождение n-го исла Фибоначи
double Fn(int n)
{
    return (pow((1 + sqrt(5)) / 2, n) - pow((1 - sqrt(5)) / 2, n)) / sqrt(5);
}

// Метод Половинного деления
double HalfDivision(int Jchoice, int text, double A, double B, double EPS, int& iter)
{
    double x1 = A, x2 = B, delta = EPS / 2; // Динамическая граница метода деления пополам
    cout << "You using Half Division method!\n\n";
    if (text == 1)
    {
        cout << "Iter=0  x1=" << x1 << "  x2=" << x2 << endl;
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
    } while ((x2 - x1 - delta) / pow(2, iter + 1) > EPS);
    return (x1 + x2) / 2;
}

// Метод Золотого сечения
double GoldenRatio(int Jchoice, int text, double A, double B, double EPS, int& iter)
{
    double x1 = A, x2 = B;
    cout << "You using Golden ratio method!\n\n";
    if (text == 1)
    {
        cout << "Iter=0  x1=" << x1 << "  x2=" << x2 << endl;
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
    return (x1 + x2) / 2;
}

// Метод Фибоначи
double Fibonachi(int Jchoice, int text, double A, double B, double EPS, int& iter)
{
    double an = A, bn = B, x1, x2, xn;
    int n = (int)(log(EPS / ((bn - an) * sqrt(5))) / log(2 / (sqrt(5) + 1)));
    cout << "You using Fibonachi method!\n\n";
    if (text == 1)
    {
        cout << "Iter=0  x1=" << an << "  x2=" << bn << endl;
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
        return xn - eps_fin;
    }
    return xn + eps_fin;
}

// Нахождение константы Липищица
double L(int choice, double a, double b)
{
    double max = -1;
    double N = abs((b - a) / 1000);
    for (double i = a; i <= b; i += N)
    {
        if (abs(dJ(i, choice)) > max)
        {
            max = abs(dJ(i, choice));
        }
    }
    return max;
}

// Метод Ломанных
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
            cout << "Xmin=" << Xmin << "   J=" << J(x1, Jchoice) << endl;
            Jmin = J(x1, Jchoice);
            Xmin = x1;
        }
    }
}

// Метод касательных
double Tangents(int Jchoice, int text, double A, double B, double EPS, int& iter)
{
    cout << "You using Tangents method!\n";
    double an = A, bn = B, xn = A;
    if (text == 1)
    {
        cout << "\nIter=0:  x1=" << an << "  x2=" << bn;
    }
    while (abs(dJ(xn, Jchoice)) > EPS && dJ(an, Jchoice) * dJ(bn, Jchoice) < 0)
    {
        iter++;
        xn = (bn * dJ(bn, Jchoice) - an * dJ(an, Jchoice) + J(an, Jchoice) - J(bn, Jchoice)) / (dJ(bn, Jchoice) - dJ(an, Jchoice));
        if (dJ(xn, Jchoice) >= 0)
        {
            bn = xn;
        }
        else
        {
            an = xn;
        }
        if (text == 1)
        {
            cout << "\nIter=" << iter << "  A=" << an << "  B=" << bn;
        }
    }
    if (dJ(an, Jchoice) == 0 || (dJ(an, Jchoice) > 0 && dJ(bn, Jchoice) > 0))
    {
        return an;
    }
    else
    {
        return bn;
    }
}

// Метод Ньютона
double Newton(int Jchoice, int text, double A, double B, double EPS, int& iter)
{
    cout << "You using Newton method!\n";
    double xn = (A + B) / 2;
    if (text == 1)
    {
        cout << "\nIter=0  xk=" << xn;
    }
    while (abs(dJ(xn, Jchoice)) > EPS)
    {
        iter++;
        xn -= dJ(xn, Jchoice) / ddJ(xn, Jchoice);
        if (text == 1)
        {
            cout << "\nIter=" << iter << "   xk=" << xn;
        }
    }
    return xn;
}

int main()
{
    int Jchoice, Mchoice, text, branch = 1, Iter = 0;
    double a, b, eps, Xmin, Jmin, Lip, x1, p1;
    cout << "You using programm for find minimum J1 or J2 different various\n\nPlease choose function J\n\n";
    cout << "1) x^2 + sin(x)\n2) x * cos(2x) + 1\n\nYour choice: ";
    cin >> Jchoice;
    cout << "\n---------------------------\n\n";
    cout << "Choose method:\n1) Half division\n2) Golden Ratio\n3) Fibonachi\n4) Broken line\n5) Tangents\n6) Newton\n\nYour choice: ";
    cin >> Mchoice;
    cout << "\n---------------------------\n\n";
    if (Mchoice != 4)
    {
        cout << "Do you want see all calculations? 1-yes  2-no\n\nYour choice ";
        cin >> text;
        cout << "\n---------------------------\n\n";
    }
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
        Xmin = HalfDivision(Jchoice, text, a, b, eps, Iter);
        break;
    case 2:
        Xmin = GoldenRatio(Jchoice, text, a, b, eps, Iter);
        break;
    case 3:
        Xmin = Fibonachi(Jchoice, text, a, b, eps, Iter);
        break;
    case 4:
        cout << "You using Broken line method!\n\n";
        Xmin = a;
        Jmin = J(a, Jchoice);
        Lip = L(Jchoice, a, b);
        x1 = 1. / (2 * Lip) * (J(a, Jchoice) - J(b, Jchoice) + Lip * (a + b));
        p1 = 1. / 2 * (J(a, Jchoice) + J(b, Jchoice) + Lip * (a - b));
        BrokenLine(Jchoice, x1, p1, Lip, eps, Xmin, Jmin, branch);
        cout << "Minimum=" << Xmin << "  Jmin=" << Jmin << "  Total branch=" << branch << endl;
        break;
    case 5:
        Xmin = Tangents(Jchoice, text, a, b, eps, Iter);
        break;
    case 6:
        Xmin = Newton(Jchoice, text, a, b, eps, Iter);
        break;
    }
    if (Mchoice != 4)
    {
        cout << "\n---------------------------\n\nIter=" << Iter << "   Minimum = " << Xmin << "  Jmin = " << J(Xmin, Jchoice) << endl;
    }
    return 100;
}