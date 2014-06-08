using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Finite_element_method
{
    class FEM
    {
        // p(x)result_system'' + q(x)result_system' + r(x)result_system = f(x)

        // c11*result_system(a) + c12*result_system'(a) = c13

        // c21*result_system(b) + c22*result_system'(b) = c23

        double c11 = 0;
        double c12 = -1;
        double c13 = -2;

        double c21 = 1;
        double c22 = 0;
        double c23 = 0;

        double left = 1;
        double right = Math.Pow(Math.E,Math.PI);

        double h;

        double find_p(double x)
        {
            return x;
        }

        double find_q(double x)
        {
            return 0;
        }

        double find_r(double x)
        {
            return -4 /  x;
        }

        double find_f(double x)
        {
            return 0;
        }

        public double integral(int low_index_integrate, int type)
        {
            double result;
            double x, xn, xk, dx;
            int i;
            int n1 = 1000;
            xn = left + low_index_integrate * h;               //предел интегрирования x(i-1)  x начало
            xk = left + (low_index_integrate + 1) * h;         //предел интегрирования x(i)   x конец
            dx = h / n1;
            x = xn + dx / 2;
            if (low_index_integrate < 0)
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
                    result += find_r(x) * (x - xn) * (xk - x);
                    x = x + dx;
                }
                result = result * dx;
                return result;
            }
            if (type == 2)
            {
                for (i = 0; i < n1; i++)
                {
                    result += find_r(x) * Math.Pow((x - xn), 2);
                    x = x + dx;
                }
                result = result * dx;
                return result;
            }
            if (type == 3)
            {
                for (i = 0; i < n1; i++)
                {
                    result += find_r(x) * Math.Pow((xk - x), 2);
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

        private double fi(int i, double x, double h)
        {
            double x_i = left + (i - 1) * h;
            double xi = left + i * h;
            double xi1 = left + (i + 1) * h;

            if (Math.Abs(x - xi) > h)
                return 0.0;
            else
            {
                if (x >= x_i && x <= xi)
                    return 1 + (x - x_i) / h;

                if ((x <= xi1 && x > xi))
                
                    return 1 - (x - xi1) / h;
                else return 0.0;
            } 
        }

        public void algoritm(int n_system, StreamWriter file_x, StreamWriter file_y) // tridiagonal matrix algorithm
        {
            double[] x, result_system, f, p, q, r, y;
            double[] a, b, c, fmk;
            double[] alpha, beta;

            int n_points = 1000;

            int i;

            x = new double[n_system + 1];
            result_system = new double[n_system + 1];  //решение прогонки
            f = new double[n_system + 1];
            p = new double[n_system + 1];
            q = new double[n_system + 1];
            r = new double[n_system + 1];
            y = new double[n_points + 1];

            a = new double[n_system + 1];
            b = new double[n_system + 1];
            c = new double[n_system + 1];
            fmk = new double[n_system + 1];  //коэфициент для МКЭ

            alpha = new double[n_system + 1];
            beta = new double[n_system + 1];

            h = (right - left) / n_system;

            double k = 1 / Math.Pow(h, 2);   // =1/h^2 замена для упрощения

            for (i = 0; i < n_system + 1; i++)
            {
                x[i] = left + i * h;
                f[i] = find_f(x[i]);
                p[i] = find_p(x[i]);
                q[i] = find_q(x[i]);
                r[i] = find_r(x[i]);

                a[i] = (-1) * k * integral(i - 1, 1) + k * integral(i - 1, 4);
                b[i] = k * integral(i, 1) + k * integral(i - 1, 1) + k * integral(i - 1, 2) + k * integral(i, 3);
                c[i] = (-1) * k * integral(i, 1) + k * integral(i, 4);
                fmk[i] = (1 / h) * integral(i - 1, 5) + (1 / h) * integral(i, 6);
            }

            double k0 = -c[0] / (b[0] + c11 * p[0]);
            //double k1 = -a[n_system] / (b[n_system-2] - c21 * p[n_system]);
            double k1 = 0;

            double n0 = (fmk[0] +p[0]*c13)/ (b[0] + c11 * p[0]);
            ///double n1 = (fmk[n_system] +p[n_system]*c23) / (b[n_system-2] + c21 * p[n_system]);
            double n1 = c23;
            alpha[1] = k0;
            beta[1] = n0;

            for (i = 1; i < n_system - 1; i++)
            {
                alpha[i + 1] = -c[i] / (b[i] + a[i] * alpha[i]);
                beta[i + 1] = (fmk[i] - a[i] * beta[i]) / (b[i] + alpha[i] * a[i]);
            }

            result_system[n_system] = (n1 + k1 * beta[n_system - 1]) / (1 - k1 * alpha[n_system - 1]);
            for (i = n_system - 1; i > 0; i--)
            {
                result_system[i] = (alpha[i] * result_system[i + 1] + beta[i]);
            }
            result_system[0] = k0 * result_system[1]+n0;

            y[0] = 0;
            for ( i = 0; i < n_system; i++)
            {
                for (int j = 0; j < n_system; j++)
                {
                    y[j] += result_system[i] * fi(i, left + j * h, h);
                }
            }


            if (n_system == 80)
                for (i = 0; i < n_system; i = i + n_system / 20)
                {
                    file_x.WriteLine(x[i]);
                    file_y.WriteLine(y[i]);
                }
            else
                for (i = 0; i < n_system; i++)
                {
                    file_x.WriteLine(left + i * h);
                    file_y.WriteLine(result_system[i]);
                }
           
        }

        void printMass(double[] m, double[] k)
        {
            int n = m.Length;
            Console.WriteLine("x          result_system ");
            for (int i = 0; i < n; i+=2)
                Console.WriteLine("{0,10:F4}{1,10:F4}", m[i], k[i]);
            Console.WriteLine();
        }

    }

    class Program
    {

        static void Main(string[] args)
        {
            FEM tdm = new FEM();
            string path = "E:\\Dropbox\\Visual Studio\\Projects\\Finite_element_method\\";
            StreamWriter fx = new System.IO.StreamWriter(@"" + path + "fx1.txt");
            StreamWriter fy = new System.IO.StreamWriter(@"" + path + "fy1.txt");
            StreamWriter fx1 = new System.IO.StreamWriter(@"" + path + "fx2.txt");
            StreamWriter fy1 = new System.IO.StreamWriter(@"" + path + "fy2.txt");

            tdm.algoritm(200, fx, fy);
            fy.Close();
            fx.Close();
             tdm.algoritm(100, fx1, fy1);
             fy1.Close();
             fx1.Close();

            Graphics e = new Graphics();
            e.ShowDialog();
            e.Select();
            Console.ReadKey();

        }
    }
}
