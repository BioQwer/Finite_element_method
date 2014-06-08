using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Finite_element_method
{
    public partial class Graphics : Form
    {
        public Graphics()
        {
            InitializeComponent();
            points();

        }

        private void points()
        {
            string fxline,fyline;
            int counter1 = -1, counter2 = -1;
            chart1.Visible = true;
            chart1.Series["Series1"].Points.Clear();
            chart1.Series["Series2"].Points.Clear();
            string path = "E:\\Dropbox\\Visual Studio\\Projects\\Finite_element_method\\";
            StreamReader fx = new System.IO.StreamReader(@"" + path + "fx1.txt");
            StreamReader fy = new System.IO.StreamReader(@"" + path + "fy1.txt");
            StreamReader fx1 = new System.IO.StreamReader(@"" + path + "fx2.txt");
            StreamReader fy1 = new System.IO.StreamReader(@"" + path + "fy2.txt");
            while (((fxline = fx.ReadLine()) != null) && ((fyline = fy.ReadLine()) != null))
            {
                chart1.Series["Series1"].Points.AddXY(Convert.ToDouble(fxline), Convert.ToDouble(fyline));
                counter1++;
            }
            while ((fxline = fx1.ReadLine()) != null && (fyline = fy1.ReadLine()) != null)
            {
                chart1.Series["Series2"].Points.AddXY(Convert.ToDouble(fxline), Convert.ToDouble(fyline));
                counter2++;
            }
            fy.Close();
            fx.Close();
            fy1.Close();
            fx1.Close();
            chart1.Series["Series1"].LegendText = "Приближение "+counter1.ToString();
            chart1.Series["Series2"].LegendText = "Приближение "+counter2.ToString();
        }

    }
}
