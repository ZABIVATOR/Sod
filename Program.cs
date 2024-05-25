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

            int N = 800;
            string result_name_file = $"Sod ";
            double gamma = 1.4;
            double time = 0.4;

            double ro_left = 1;
            double u_left = 0;
            double p_left = 1;

            double ro_right = 0.125;
            double u_right = 0;
            double p_right = 0.1;

            var param = new Parameters(N, time, [ro_left, u_left, p_left], [ro_right, u_right, p_right], gamma);

            var CIR = new CIR(param);
            var k = CIR.CalculateWrite(result_name_file,false);

            var gg = new Godunov(param);
            var g = gg.CalculateWrite(result_name_file,false);

            Plots.PlotBoth(param, k, g, result_name_file); 

            return;

        }
    }
}