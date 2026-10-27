#include "ode_solvers.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <sstream>

namespace Task41
{
    double ExactY(double x)
    {
        return (1.0 + x) * std::exp(x * x);
    }

    double ExactZ(double x)
    {
        return (1.0 + 2.0 * x + 2.0 * x * x) * std::exp(x * x);
    }

    double f1(double, double, double z)
    {
        return z;
    }

    double f2(double x, double y, double z)
    {
        return 4.0 * x * z - (4.0 * x * x - 2.0) * y;
    }

    struct SolutionGrid
    {
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
    };

    // Метод Эйлера
    SolutionGrid SolveEuler(double a, double b, double y0, double z0, double h)
    {
        int N = static_cast<int>(std::round((b - a) / h));
        SolutionGrid grid;
        grid.x.resize(N + 1);
        grid.y.resize(N + 1);
        grid.z.resize(N + 1);

        grid.x[0] = a;
        grid.y[0] = y0;
        grid.z[0] = z0;

        for (int i = 0; i < N; ++i)
        {
            double xi = grid.x[i];
            double yi = grid.y[i];
            double zi = grid.z[i];

            grid.x[i + 1] = a + (i + 1) * h;
            grid.y[i + 1] = yi + h * f1(xi, yi, zi);
            grid.z[i + 1] = zi + h * f2(xi, yi, zi);
        }
        return grid;
    }

    // Рунге-Кутта 4 порядка
    SolutionGrid SolveRK4(double a, double b, double y0, double z0, double h)
    {
        int N = static_cast<int>(std::round((b - a) / h));
        SolutionGrid grid;
        grid.x.resize(N + 1);
        grid.y.resize(N + 1);
        grid.z.resize(N + 1);

        grid.x[0] = a;
        grid.y[0] = y0;
        grid.z[0] = z0;

        for (int i = 0; i < N; ++i)
        {
            double xi = grid.x[i];
            double yi = grid.y[i];
            double zi = grid.z[i];

            double k1_y = f1(xi, yi, zi);
            double k1_z = f2(xi, yi, zi);

            double k2_y = f1(xi + 0.5 * h, yi + 0.5 * h * k1_y, zi + 0.5 * h * k1_z);
            double k2_z = f2(xi + 0.5 * h, yi + 0.5 * h * k1_y, zi + 0.5 * h * k1_z);

            double k3_y = f1(xi + 0.5 * h, yi + 0.5 * h * k2_y, zi + 0.5 * h * k2_z);
            double k3_z = f2(xi + 0.5 * h, yi + 0.5 * h * k2_y, zi + 0.5 * h * k2_z);

            double k4_y = f1(xi + h, yi + h * k3_y, zi + h * k3_z);
            double k4_z = f2(xi + h, yi + h * k3_y, zi + h * k3_z);

            grid.x[i + 1] = a + (i + 1) * h;
            grid.y[i + 1] = yi + (h / 6.0) * (k1_y + 2.0 * k2_y + 2.0 * k3_y + k4_y);
            grid.z[i + 1] = zi + (h / 6.0) * (k1_z + 2.0 * k2_z + 2.0 * k3_z + k4_z);
        }
        return grid;
    }

    // Метод Адамса 4 порядка
    SolutionGrid SolveAdams4(double a, double b, double y0, double z0, double h)
    {
        int N = static_cast<int>(std::round((b - a) / h));
        SolutionGrid rk = SolveRK4(a, b, y0, z0, h);

        SolutionGrid grid;
        grid.x = rk.x;
        grid.y.resize(N + 1);
        grid.z.resize(N + 1);

        // первые 4 точки из РК4
        for (int i = 0; i < 4 && i <= N; ++i)
        {
            grid.y[i] = rk.y[i];
            grid.z[i] = rk.z[i];
        }

        for (int i = 3; i < N; ++i)
        {
            double f1_i = f1(grid.x[i], grid.y[i], grid.z[i]);
            double f1_im1 = f1(grid.x[i - 1], grid.y[i - 1], grid.z[i - 1]);
            double f1_im2 = f1(grid.x[i - 2], grid.y[i - 2], grid.z[i - 2]);
            double f1_im3 = f1(grid.x[i - 3], grid.y[i - 3], grid.z[i - 3]);

            double f2_i = f2(grid.x[i], grid.y[i], grid.z[i]);
            double f2_im1 = f2(grid.x[i - 1], grid.y[i - 1], grid.z[i - 1]);
            double f2_im2 = f2(grid.x[i - 2], grid.y[i - 2], grid.z[i - 2]);
            double f2_im3 = f2(grid.x[i - 3], grid.y[i - 3], grid.z[i - 3]);

            // предиктор
            double y_pred = grid.y[i] + (h / 24.0) * (55.0 * f1_i - 59.0 * f1_im1 + 37.0 * f1_im2 - 9.0 * f1_im3);
            double z_pred = grid.z[i] + (h / 24.0) * (55.0 * f2_i - 59.0 * f2_im1 + 37.0 * f2_im2 - 9.0 * f2_im3);

            // корректор
            double f1_next = f1(grid.x[i + 1], y_pred, z_pred);
            double f2_next = f2(grid.x[i + 1], y_pred, z_pred);

            grid.y[i + 1] = grid.y[i] + (h / 24.0) * (9.0 * f1_next + 19.0 * f1_i - 5.0 * f1_im1 + f1_im2);
            grid.z[i + 1] = grid.z[i] + (h / 24.0) * (9.0 * f2_next + 19.0 * f2_i - 5.0 * f2_im1 + f2_im2);
        }
        return grid;
    }
}

namespace Task42
{
    double ExactY(double x)
    {
        return x + std::exp(-2.0 * x);
    }

    double p_func(double x)
    {
        return 4.0 * x / (2.0 * x + 1.0);
    }

    double q_func(double x)
    {
        return -4.0 / (2.0 * x + 1.0);
    }

    double f_func(double)
    {
        return 0.0;
    }

    struct SolutionGrid
    {
        std::vector<double> x;
        std::vector<double> y;
    };

    void RK4_CauchyStep(double x, double y, double z, double h, double &ynext, double &znext)
    {
        auto f1 = [](double, double, double z_val) { return z_val; };
        auto f2 = [](double x_val, double y_val, double z_val) {
            return (4.0 * y_val - 4.0 * x_val * z_val) / (2.0 * x_val + 1.0);
        };

        double k1_y = f1(x, y, z);
        double k1_z = f2(x, y, z);

        double k2_y = f1(x + 0.5 * h, y + 0.5 * h * k1_y, z + 0.5 * h * k1_z);
        double k2_z = f2(x + 0.5 * h, y + 0.5 * h * k1_y, z + 0.5 * h * k1_z);

        double k3_y = f1(x + 0.5 * h, y + 0.5 * h * k2_y, z + 0.5 * h * k2_z);
        double k3_z = f2(x + 0.5 * h, y + 0.5 * h * k2_y, z + 0.5 * h * k2_z);

        double k4_y = f1(x + h, y + h * k3_y, z + h * k3_z);
        double k4_z = f2(x + h, y + h * k3_y, z + h * k3_z);

        ynext = y + (h / 6.0) * (k1_y + 2.0 * k2_y + 2.0 * k3_y + k4_y);
        znext = z + (h / 6.0) * (k1_z + 2.0 * k2_z + 2.0 * k3_z + k4_z);
    }

    void SolveCauchyForShooting(double a, double b, double eta, double h, std::vector<double> &xs, std::vector<double> &ys, std::vector<double> &zs)
    {
        int N = static_cast<int>(std::round((b - a) / h));
        xs.resize(N + 1);
        ys.resize(N + 1);
        zs.resize(N + 1);

        xs[0] = a;
        ys[0] = eta;
        zs[0] = -1.0;

        for (int i = 0; i < N; ++i)
        {
            xs[i + 1] = a + (i + 1) * h;
            RK4_CauchyStep(xs[i], ys[i], zs[i], h, ys[i + 1], zs[i + 1]);
        }
    }

    SolutionGrid SolveShooting(double a, double b, double h, double &found_eta)
    {
        auto Discrepancy = [a, b, h](double eta, std::vector<double> &xs, std::vector<double> &ys, std::vector<double> &zs) {
            SolveCauchyForShooting(a, b, eta, h, xs, ys, zs);
            return zs.back() + 2.0 * ys.back() - 3.0;
        };

        std::vector<double> xs, ys, zs;
        double eta0 = 0.0;
        double phi0 = Discrepancy(eta0, xs, ys, zs);

        double eta1 = 1.0;
        double phi1 = Discrepancy(eta1, xs, ys, zs);

        // метод секущих для пристрелки
        double eta_curr = eta1;
        while (std::abs(phi1) > 1e-8)
        {
            eta_curr = eta1 - phi1 * (eta1 - eta0) / (phi1 - phi0);
            eta0 = eta1;
            phi0 = phi1;
            eta1 = eta_curr;
            phi1 = Discrepancy(eta1, xs, ys, zs);
        }

        found_eta = eta_curr;
        SolveCauchyForShooting(a, b, found_eta, h, xs, ys, zs);

        SolutionGrid result;
        result.x = xs;
        result.y = ys;
        return result;
    }

    // Конечно-разностный метод
    SolutionGrid SolveFiniteDifference(double a, double b, double h)
    {
        int N = static_cast<int>(std::round((b - a) / h));
        SolutionGrid result;
        result.x.resize(N + 1);
        result.y.resize(N + 1);

        for (int i = 0; i <= N; ++i)
        {
            result.x[i] = a + i * h;
        }

        std::vector<double> A(N + 1, 0.0);
        std::vector<double> B(N + 1, 0.0);
        std::vector<double> C(N + 1, 0.0);
        std::vector<double> D(N + 1, 0.0);

        // левое краевое условие
        double p0 = p_func(a);
        double q0 = q_func(a);
        B[0] = -2.0 + h * h * q0;
        C[0] = 2.0;
        D[0] = -2.0 * h * (1.0 - 0.5 * h * p0);

        // внутренние узлы
        for (int i = 1; i < N; ++i)
        {
            double xi = result.x[i];
            double pi = p_func(xi);
            double qi = q_func(xi);

            A[i] = 1.0 - 0.5 * h * pi;
            B[i] = -2.0 + h * h * qi;
            C[i] = 1.0 + 0.5 * h * pi;
            D[i] = h * h * f_func(xi);
        }

        // правое краевое условие
        double pN = p_func(b);
        double qN = q_func(b);
        A[N] = 2.0;
        B[N] = -2.0 + h * h * qN - 4.0 * h * (1.0 + 0.5 * h * pN);
        D[N] = -6.0 * h * (1.0 + 0.5 * h * pN);

        // прогонка
        std::vector<double> P(N + 1, 0.0);
        std::vector<double> Q(N + 1, 0.0);

        P[0] = -C[0] / B[0];
        Q[0] = D[0] / B[0];

        for (int i = 1; i <= N; ++i)
        {
            double denom = B[i] + A[i] * P[i - 1];
            if (i < N)
            {
                P[i] = -C[i] / denom;
            }
            Q[i] = (D[i] - A[i] * Q[i - 1]) / denom;
        }

        result.y[N] = Q[N];
        for (int i = N - 1; i >= 0; --i)
        {
            result.y[i] = P[i] * result.y[i + 1] + Q[i];
        }

        return result;
    }
}

void Solve_Task_4_1(const std::string &filename)
{
    std::cout << "\n--- Задача 4.1: Задача Коши ---\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    std::string line;
    double a = 0.0, b = 1.0, y0 = 1.0, z0 = 1.0, h = 0.1;
    std::vector<double> values;

    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::stringstream ss(line);
        double val;
        while (ss >> val)
            values.push_back(val);
    }
    file.close();

    if (values.size() >= 5)
    {
        a = values[0];
        b = values[1];
        y0 = values[2];
        z0 = values[3];
        h = values[4];
    }

    std::cout << "Отрезок: [" << a << ", " << b << "], y(" << a << ") = " << y0 << ", y'(" << a << ") = " << z0 << ", h = " << h << "\n\n";

    auto euler_h = Task41::SolveEuler(a, b, y0, z0, h);
    auto rk4_h = Task41::SolveRK4(a, b, y0, z0, h);
    auto adams_h = Task41::SolveAdams4(a, b, y0, z0, h);

    double h2 = h / 2.0;
    auto euler_h2 = Task41::SolveEuler(a, b, y0, z0, h2);
    auto rk4_h2 = Task41::SolveRK4(a, b, y0, z0, h2);
    auto adams_h2 = Task41::SolveAdams4(a, b, y0, z0, h2);

    int N = static_cast<int>(euler_h.x.size());

    std::cout << ">>> МЕТОД ЭЙЛЕРА <<<\n";
    std::cout << std::string(75, '-') << "\n";
    std::cout << std::setw(4) << "k" << " | "
              << std::setw(6) << "x_k" << " | "
              << std::setw(12) << "y_Euler" << " | "
              << std::setw(12) << "y_exact" << " | "
              << std::setw(14) << "|y - y_ex|" << " | "
              << std::setw(14) << "Runge-Romberg\n";
    std::cout << std::string(75, '-') << "\n";

    double max_err_euler = 0.0;
    for (int i = 0; i < N; ++i)
    {
        double xk = euler_h.x[i];
        double yk = euler_h.y[i];
        double y_ex = Task41::ExactY(xk);
        double err_ex = std::abs(yk - y_ex);
        double yk_h2 = euler_h2.y[2 * i];
        double err_rr = std::abs(yk - yk_h2) / (std::pow(2.0, 1.0) - 1.0);
        max_err_euler = std::max(max_err_euler, err_ex);

        std::cout << std::setw(4) << i << " | "
                  << std::setw(6) << std::fixed << std::setprecision(2) << xk << " | "
                  << std::setw(12) << std::fixed << std::setprecision(6) << yk << " | "
                  << std::setw(12) << std::fixed << std::setprecision(6) << y_ex << " | "
                  << std::setw(14) << std::scientific << std::setprecision(4) << err_ex << " | "
                  << std::setw(14) << std::scientific << std::setprecision(4) << err_rr << "\n";
    }

    std::cout << "\n>>> МЕТОД РУНГЕ-КУТТЫ 4 <<<\n";
    std::cout << std::string(75, '-') << "\n";
    std::cout << std::setw(4) << "k" << " | "
              << std::setw(6) << "x_k" << " | "
              << std::setw(12) << "y_RK4" << " | "
              << std::setw(12) << "y_exact" << " | "
              << std::setw(14) << "|y - y_ex|" << " | "
              << std::setw(14) << "Runge-Romberg\n";
    std::cout << std::string(75, '-') << "\n";

    double max_err_rk4 = 0.0;
    for (int i = 0; i < N; ++i)
    {
        double xk = rk4_h.x[i];
        double yk = rk4_h.y[i];
        double y_ex = Task41::ExactY(xk);
        double err_ex = std::abs(yk - y_ex);
        double yk_h2 = rk4_h2.y[2 * i];
        double err_rr = std::abs(yk - yk_h2) / (std::pow(2.0, 4.0) - 1.0);
        max_err_rk4 = std::max(max_err_rk4, err_ex);

        std::cout << std::setw(4) << i << " | "
                  << std::setw(6) << std::fixed << std::setprecision(2) << xk << " | "
                  << std::setw(12) << std::fixed << std::setprecision(6) << yk << " | "
                  << std::setw(12) << std::fixed << std::setprecision(6) << y_ex << " | "
                  << std::setw(14) << std::scientific << std::setprecision(4) << err_ex << " | "
                  << std::setw(14) << std::scientific << std::setprecision(4) << err_rr << "\n";
    }

    std::cout << "\n>>> МЕТОД АДАМСА 4 <<<\n";
    std::cout << std::string(75, '-') << "\n";
    std::cout << std::setw(4) << "k" << " | "
              << std::setw(6) << "x_k" << " | "
              << std::setw(12) << "y_Adams" << " | "
              << std::setw(12) << "y_exact" << " | "
              << std::setw(14) << "|y - y_ex|" << " | "
              << std::setw(14) << "Runge-Romberg\n";
    std::cout << std::string(75, '-') << "\n";

    double max_err_adams = 0.0;
    for (int i = 0; i < N; ++i)
    {
        double xk = adams_h.x[i];
        double yk = adams_h.y[i];
        double y_ex = Task41::ExactY(xk);
        double err_ex = std::abs(yk - y_ex);
        double yk_h2 = adams_h2.y[2 * i];
        double err_rr = std::abs(yk - yk_h2) / (std::pow(2.0, 4.0) - 1.0);
        max_err_adams = std::max(max_err_adams, err_ex);

        std::cout << std::setw(4) << i << " | "
                  << std::setw(6) << std::fixed << std::setprecision(2) << xk << " | "
                  << std::setw(12) << std::fixed << std::setprecision(6) << yk << " | "
                  << std::setw(12) << std::fixed << std::setprecision(6) << y_ex << " | "
                  << std::setw(14) << std::scientific << std::setprecision(4) << err_ex << " | "
                  << std::setw(14) << std::scientific << std::setprecision(4) << err_rr << "\n";
    }

    std::cout << "\nМаксимальные погрешности:\n";
    std::cout << "Эйлер:   " << std::scientific << std::setprecision(4) << max_err_euler << "\n";
    std::cout << "РК4:     " << std::scientific << std::setprecision(4) << max_err_rk4 << "\n";
    std::cout << "Адамс 4: " << std::scientific << std::setprecision(4) << max_err_adams << "\n";

    std::ofstream csv("data/lab4/plot_4_1.csv");
    if (csv.is_open())
    {
        csv << "x,exact,euler,rk4,adams,euler_err,rk4_err,adams_err\n";
        for (int i = 0; i < N; ++i)
        {
            double xk = euler_h.x[i];
            double ex = Task41::ExactY(xk);
            double ye = euler_h.y[i];
            double yr = rk4_h.y[i];
            double ya = adams_h.y[i];
            csv << std::fixed << std::setprecision(6) << xk << ","
                << ex << "," << ye << "," << yr << "," << ya << ","
                << std::scientific << std::setprecision(6)
                << std::abs(ye - ex) << ","
                << std::abs(yr - ex) << ","
                << std::abs(ya - ex) << "\n";
        }
        csv.close();
    }
}

void Solve_Task_4_2(const std::string &filename)
{
    std::cout << "\n--- Задача 4.2: Краевая задача ---\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    std::string line;
    double a = 0.0, b = 1.0, h = 0.1;
    std::vector<double> values;

    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::stringstream ss(line);
        double val;
        while (ss >> val)
            values.push_back(val);
    }
    file.close();

    if (values.size() >= 3)
    {
        a = values[0];
        b = values[1];
        h = values.back();
    }

    std::cout << "Отрезок: [" << a << ", " << b << "], h = " << h << "\n\n";

    double found_eta = 0.0;
    auto shoot_h = Task42::SolveShooting(a, b, h, found_eta);
    double dummy_eta = 0.0;
    double h2 = h / 2.0;
    auto shoot_h2 = Task42::SolveShooting(a, b, h2, dummy_eta);

    std::cout << ">>> МЕТОД СТРЕЛЬБЫ <<<\n";
    std::cout << "Найденный y(0) = " << std::fixed << std::setprecision(6) << found_eta << "\n";
    std::cout << std::string(75, '-') << "\n";
    std::cout << std::setw(4) << "i" << " | "
              << std::setw(6) << "x_i" << " | "
              << std::setw(12) << "y_Shoot" << " | "
              << std::setw(12) << "y_exact" << " | "
              << std::setw(14) << "|y - y_ex|" << " | "
              << std::setw(14) << "Runge-Romberg\n";
    std::cout << std::string(75, '-') << "\n";

    int N = static_cast<int>(shoot_h.x.size());
    double max_err_shoot = 0.0;

    for (int i = 0; i < N; ++i)
    {
        double xi = shoot_h.x[i];
        double yi = shoot_h.y[i];
        double y_ex = Task42::ExactY(xi);
        double err_ex = std::abs(yi - y_ex);
        double yi_h2 = shoot_h2.y[2 * i];
        double err_rr = std::abs(yi - yi_h2) / (std::pow(2.0, 4.0) - 1.0);
        max_err_shoot = std::max(max_err_shoot, err_ex);

        std::cout << std::setw(4) << i << " | "
                  << std::setw(6) << std::fixed << std::setprecision(2) << xi << " | "
                  << std::setw(12) << std::fixed << std::setprecision(6) << yi << " | "
                  << std::setw(12) << std::fixed << std::setprecision(6) << y_ex << " | "
                  << std::setw(14) << std::scientific << std::setprecision(4) << err_ex << " | "
                  << std::setw(14) << std::scientific << std::setprecision(4) << err_rr << "\n";
    }

    auto fd_h = Task42::SolveFiniteDifference(a, b, h);
    auto fd_h2 = Task42::SolveFiniteDifference(a, b, h2);

    std::cout << "\n>>> КОНЕЧНО-РАЗНОСТНЫЙ МЕТОД <<<\n";
    std::cout << std::string(75, '-') << "\n";
    std::cout << std::setw(4) << "i" << " | "
              << std::setw(6) << "x_i" << " | "
              << std::setw(12) << "y_FinDiff" << " | "
              << std::setw(12) << "y_exact" << " | "
              << std::setw(14) << "|y - y_ex|" << " | "
              << std::setw(14) << "Runge-Romberg\n";
    std::cout << std::string(75, '-') << "\n";

    double max_err_fd = 0.0;
    for (int i = 0; i < N; ++i)
    {
        double xi = fd_h.x[i];
        double yi = fd_h.y[i];
        double y_ex = Task42::ExactY(xi);
        double err_ex = std::abs(yi - y_ex);
        double yi_h2 = fd_h2.y[2 * i];
        double err_rr = std::abs(yi - yi_h2) / (std::pow(2.0, 2.0) - 1.0);
        max_err_fd = std::max(max_err_fd, err_ex);

        std::cout << std::setw(4) << i << " | "
                  << std::setw(6) << std::fixed << std::setprecision(2) << xi << " | "
                  << std::setw(12) << std::fixed << std::setprecision(6) << yi << " | "
                  << std::setw(12) << std::fixed << std::setprecision(6) << y_ex << " | "
                  << std::setw(14) << std::scientific << std::setprecision(4) << err_ex << " | "
                  << std::setw(14) << std::scientific << std::setprecision(4) << err_rr << "\n";
    }

    std::cout << "\nМаксимальные погрешности:\n";
    std::cout << "Стрельба:  " << std::scientific << std::setprecision(4) << max_err_shoot << "\n";
    std::cout << "Разности:  " << std::scientific << std::setprecision(4) << max_err_fd << "\n";

    std::ofstream csv("data/lab4/plot_4_2.csv");
    if (csv.is_open())
    {
        csv << "x,exact,shoot,fd,shoot_err,fd_err\n";
        for (int i = 0; i < N; ++i)
        {
            double xi = shoot_h.x[i];
            double ex = Task42::ExactY(xi);
            double ys = shoot_h.y[i];
            double yf = fd_h.y[i];
            csv << std::fixed << std::setprecision(6) << xi << ","
                << ex << "," << ys << "," << yf << ","
                << std::scientific << std::setprecision(6)
                << std::abs(ys - ex) << ","
                << std::abs(yf - ex) << "\n";
        }
        csv.close();
    }
}
