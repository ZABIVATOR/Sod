using CsvHelper;
using CsvHelper.Configuration.Attributes;
using System.Globalization;
using System.IO;
using System;
using Sod.CFL;
using Sod.Utility;

namespace Sod
{
    internal static class Program
    {
        static void Main()
        {

            int N = 600;
            string result_name_file = $"";
            double gamma = 1.4;
            double time = 0.1;

            // Sod, Modified Sod, 2 vacuum waves, Big gradient, 3 breaks
            double[,] left = { { 1, 0, 1 }, { 1.0, 0.75, 1.0 }, { 1.0, -2.0, 0.4 }, { 1.0, 0.0, 1000.0 }, { 5.99924, 19.5975, 460.894 } };
            double[,] right = { { 0.125, 0.0, 0.1 }, { 0.125, 0.0, 0.1 }, { 1.0, 2.0, 0.4 }, { 1.0, 0.0, 0.01 }, { 5.99242, -6.19633, 46.0950 } };

            for (int i = 0; i <= 2; i++)
            {
                double ro_left = left[i,0];
                double u_left = left[i,1];
                double p_left = left[i,2];
                
                double ro_right = right[i,0];
                double u_right = right[i, 1];
                double p_right = right[i, 2];

                var param = new Parameters(N, time, [ro_left, u_left, p_left], [ro_right, u_right, p_right], gamma);

                var CIR = new CIR(param);
                var k = CIR.CalculateWrite(result_name_file+"Test "+i.ToString()+ " ",false,3,0.2);

                var gg = new Godunov(param);
                var g = gg.CalculateWrite(result_name_file + "Test " + i.ToString() + " ",false, 3, 0.2);

                Plots.PlotBoth(param, k, g, "Test " + i.ToString() + " ");
            }
            return;

        }
    }
}