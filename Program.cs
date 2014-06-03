using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Finite_element_method
{
    class TDM
    {
        // p(x)result_system'' + q(x)result_system' + r(x)result_system = f(x)

        // c11*result_system(a) + c12*result_system'(a) = c13

        // c21*result_system(b) + c22*result_system'(b) = c23

        double c11 = 0;
        double c12 = 1;
        double c13 = 0;

        double c21 = 1;
        double c22 = 0;
        double c23 = 0;

        double left = 1;
        double right = Math.Pow(Math.E, Math.PI);

        double h;

        double find_p(double x)
        {
            return 1;
        }

        double find_q(double x)
        {
            return 1 / x;
        }

        double find_r(double x)
        {
            return 4 / (x * x);
        }

        double find_f(double x)
        {
            return -10.0 / x + (8 * Math.Pow(Math.E, Math.PI)) / (x * x);
        }

        public double integral(int from_i_to_i1, int type)
        {
            double result;
            double x, xn, xk, dx;
            int i;
            int n1 = 1000;
            xn = left + from_i_to_i1 * h;               //предел интегрирования x(i-1)
            xk = left + (from_i_to_i1 + 1) * h;         //предел интегрирования x(i)
            dx = h / n1;
            x = xn + dx / 2;
            if (from_i_to_i1 <= 0)
                return 0;
            result = 0;
            if (type == 1)
            {
                for (i = 0; i < n1; i++)
                {
                    result += find_p(x);
                    x = x + dx;
                }
                result = result * dx;
                return result;
            }
            if (type == 4)
            {
                for (i = 0; i < n1; i++)
                {
                    result += find_q(x) * (x - xn) * (xk - x);
                    x = x + dx;
                }
                result = result * dx;
                return result;
            }
            if (type == 2)
            {
                for (i = 0; i < n1; i++)
                {
                    result += find_q(x) * Math.Pow((x - xn), 2);
                    x = x + dx;
                }
                result = result * dx;
                return result;
            }
            if (type == 3)
            {
                for (i = 0; i < n1; i++)
                {
                    result += find_q(x) * Math.Pow((xk - x), 2);
                    x = x + dx;
                }
                result = result * dx;
                return result;
            }
            if (type == 5)
            {
                for (i = 0; i < n1; i++)
                {
                    result += (x - xn) * find_f(x);
                    x = x + dx;
                }
                result = result * dx;
                return result;
            }
            if (type == 6)
            {
                for (i = 0; i < n1; i++)
                {
                    result += (xk - x) * find_f(x);
                    x = x + dx;
                }
                result = result * dx;
                return result;
            }
            return result;
        }

        private double fi(int i, int x, double h)
        {
            double x_i = left + (i-1) * h;               
            double xi = left + i * h;        
            double xi1 = left + (i+1) * h;       
            if (x_i<x||xi>x)
                return 0;
            else
            {
                if (x > x_i && x < xi)
                    return (x - x_i) / h;
                else
                    return -(x - xi1) / h;
            }
        }


        public void algoritm(int n, StreamWriter file_x, StreamWriter file_y) // tridiagonal matrix algorithm
        {
            double[] x, result_system, f, p, q, r;
            double[] a, b, c, fmk;
            double[] alpha, beta;
            double kapa1, kapa2, mu1, mu2;

            int i;

            x = new double[n + 1];
            result_system = new double[n + 1];  //решение прогонки
            f = new double[n + 1];
            p = new double[n + 1];
            q = new double[n + 1];
            r = new double[n + 1];

            a = new double[n + 1];
            b = new double[n + 1];
            c = new double[n + 1];
            fmk = new double[n + 1];  //коэфициент для МКЭ

            alpha = new double[n + 1];
            beta = new double[n + 1];

            h = (right - left) / n;

            double k = 1 / Math.Pow(h, 2);   // =1/h^2 замена для упрощения

            for (i = 0; i < n + 1; i++)
            {
                x[i] = left + i * h;
                f[i] = find_f(x[i]);
                p[i] = find_p(x[i]);
                q[i] = find_q(x[i]);
                r[i] = find_r(x[i]);

                a[i] = (-1) * k * integral(i - 1, 1) + k * integral(i - 1, 4);
                b[i] = k * integral(i, 1) + k * integral(i - 1, 1) + k * integral(i, 2) + k * integral(i, 3);
                c[i] = (-1) * k * integral(i, 1) + k * integral(i, 4);
                fmk[i] = (1 / h) * integral(i - 1, 5) + (1 / h) * integral(i, 6);
            }

            kapa1 = (-c12 / h) / (c11 - c12 / h);
            kapa2 = (-c22 * h) / (c21 - c22 / h);
            mu1 = c13 / (c11 - c12 / h);
            mu2 = c23 / (c21 - c22 / h);

            alpha[1] = kapa1;
            beta[1] = mu1;

            for (i = 1; i < n; i++)
            {
                alpha[i + 1] = a[i] / (b[i] - c[i] * alpha[i]);
                beta[i + 1] = (fmk[i] - c[i] * beta[i]) / (c[i] * alpha[i] - b[i]);
            }

            result_system[n] = (kapa2 * beta[n] + mu2) / (1 - kapa2 * alpha[n]);
            for (i = n - 1; i > 0; i--)
            {
                result_system[i] = (alpha[i + 1] * result_system[i + 1] + beta[i + 1]);
            }
            result_system[0] = kapa1 * result_system[1] + mu1;

            if (n == 80)
                for (i = 0; i < n + 1; i = i + n / 20)
                {
                    file_x.WriteLine(x[i]);
                    file_y.WriteLine(result_system[i]);
                }
            else
                for (i = 0; i < n + 1; i++)
                {
                    file_x.WriteLine(x[i]);
                    file_y.WriteLine(result_system[i]);
                }
            printMass(x, result_system);
        }

        void printMass(double[] m, double[] k)
        {
            int n = m.Length;
            Console.WriteLine("x          result_system ");
            for (int i = 0; i < n; i++)
                Console.WriteLine("{0,10:F4}{1,10:F4}", m[i], k[i]);
            Console.WriteLine();
        }


        //    void writeToFile(double* x, double n, FILE* f)
        //    {
        //        for (int i = 0; i < n; i++) fprintf(f, "%lf\n", x[i]);
        //        fclose(f);
        //    }
    }

    class Program
    {

        static void Main(string[] args)
        {
            TDM tdm = new TDM();
            string path = "F:\\Dropbox\\Visual Studio\\Projects\\Finite_element_method\\";
            StreamWriter fx = new System.IO.StreamWriter(@"" + path + "fx1.txt");
            StreamWriter fy = new System.IO.StreamWriter(@"" + path + "fy1.txt");
            StreamWriter fx1 = new System.IO.StreamWriter(@"" + path + "fx2.txt");
            StreamWriter fy1 = new System.IO.StreamWriter(@"" + path + "fy2.txt");

            tdm.algoritm(400, fx, fy);
            fy.Close();
            fx.Close();
            tdm.algoritm(200, fx1, fy1);
            fy1.Close();
            fx1.Close();

            Graphics e = new Graphics();
            e.ShowDialog();
            e.Select();
            Console.ReadKey();

        }
    }
}
