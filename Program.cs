using CsvHelper;
using CsvHelper.Configuration.Attributes;
using System.Globalization;
using System.IO;
using System;

namespace Sod
{
    internal static class Program
    {
        /// <summary>
        ///  The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            // To customize application configuration such as set high DPI settings or default font,
            // see https://aka.ms/applicationconfiguration.

            string result_name_file = "sod_fix";
            double GAMMA = 1.4;
            double T_END = 1;

            double LEFT_R = 1;
            double LEFT_U = 0;
            double LEFT_P = 1;

            double RIGHT_R = 0.125;
            double RIGHT_U = 0;
            double RIGHT_P = 0.1;
            int N = 5000;
            var param = new Parameters(N, T_END, [LEFT_R, LEFT_U, LEFT_P], [RIGHT_R, RIGHT_U, RIGHT_P], GAMMA);

            var kir = new KIR(param);
            kir.Calculate_Write(result_name_file);


            var gg = new Godunov(param);
            //Console.WriteLine(gg.Calculate_Write(result_name_file)); return base from KIR


            return;

        }
    }
}