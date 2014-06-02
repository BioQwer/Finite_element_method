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
        // p(x)y'' + q(x)y' + r(x)y = f(x)

        // c11*y(a) + c12*y'(a) = c13

        // c21*y(b) + c22*y'(b) = c23

        double c11 = 0;
        double c12 = 1;
        double c13 = 2;

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
            return 0;
        }


        public double sigma(int t, int l)
        {
            double sigm;
            double x, xn, xk, h1;
            int i;
            int n1 = 1000;
            xn = left + t * h;
            xk = left + (t + 1) * h;
            h1 = h / n1;
            x = xn + h1 / 2;
            if (t <= 0)
                return 0;
            sigm = 0;
            if (l == 1)
            {
                for (i = 0; i < n1; i++)
                {
                    sigm = sigm + left/ x;
                    x = x + h1;
                }
                sigm = sigm * h1;
                return sigm;
            }
            sigm = 0;
            if (l == 4)
            {
                for (i = 0; i < n1; i++)
                {
                    sigm = sigm - 4 * left* x * (x - xn) * (xk - x);
                    x = x + h1;
                }
                sigm = sigm * h1;
                return sigm;
            }
            sigm = 0;
            if (l == 2)
            {
                for (i = 0; i < n1; i++)
                {
                    sigm = sigm - 4 * left* x * Math.Pow((x - xn), 2);
                    x = x + h1;
                }
                sigm = sigm * h1;
                return sigm;
            }
            sigm = 0;
            if (l == 3)
            {
                for (i = 0; i < n1; i++)
                {
                    sigm = sigm - 4 * left* x * Math.Pow((xk - x), 2);
                    x = x + h1;
                }
                sigm = sigm * h1;
                return sigm;
            }
            return 0.0;
        }

        public void algoritm(int n, StreamWriter file_x, StreamWriter file_y) // tridiagonal matrix algorithm
        {
            double[] x, y, f, p, q, r;
            double[] a, b, c;
            double[] alpha, beta;
            double kapa1, kapa2, mu1, mu2;
           
            int i;

            x = new double[n + 1];
            y = new double[n + 1];
            f = new double[n + 1];
            p = new double[n + 1];
            q = new double[n + 1];
            r = new double[n + 1];

            a = new double[n + 1];
            b = new double[n + 1];
            c = new double[n + 1];

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

                a[i] = (-1) * k * sigma(i - 1, 1) + k * sigma(i - 1, 4);
                b[i] = k * sigma(i, 1) + k * sigma(i - 1, 1) + k * sigma(i, 2) + k * sigma(i - 1, 3);
                c[i] = (-1) * k * sigma(i, 1) + k * sigma(i, 4);
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
                beta[i + 1] = (f[i] - c[i] * beta[i]) / (c[i] * alpha[i] - b[i]);
            }

            y[n] = (kapa2 * beta[n] + mu2) / (1 - kapa2 * alpha[n]);
            for (i = n - 1; i > 0; i--)
            {
                y[i] = (alpha[i + 1] * y[i + 1] + beta[i + 1]);
            }
            y[0] = kapa1 * y[1] + mu1;

            if (n == 80)
                for (i = 0; i < n + 1; i = i + n / 20)
                {
                    file_x.WriteLine(x[i]);
                    file_y.WriteLine(y[i]);
                }
            else
                for (i = 0; i < n + 1; i++)
                {
                    file_x.WriteLine(x[i]);
                    file_y.WriteLine(y[i]);
                }
            printMass(x, y);
        }

        void printMass(double[] m, double[] k)
        {
            int n = m.Length;
            Console.WriteLine("x          y ");
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

            tdm.algoritm(40, fx, fy);
            fy.Close();
            fx.Close();
            tdm.algoritm(80, fx1, fy1);
            fy1.Close();
            fx1.Close();

            Graphics e = new Graphics();
            e.ShowDialog();
            e.Select();
            Console.ReadKey();

        }
    }
}
