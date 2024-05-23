using CsvHelper;
using CsvHelper.Configuration.Attributes;
using System.Globalization;
using System.IO;
using System;
using static System.Windows.Forms.VisualStyles.VisualStyleElement.TaskbarClock;

namespace Graph
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

            string result_name_file = "скорость_";
            double GAMMA = 1.4;
            double T_END = 0.6;

            double LEFT_R = 1;
            double LEFT_U = 1;
            double LEFT_P = 1;

            double RIGHT_R = 1;
            double RIGHT_U = 0;
            double RIGHT_P = 1;
            int N = 1000;




            var param = new Parameters(N, T_END, [LEFT_R, LEFT_U, LEFT_P], [RIGHT_R, RIGHT_U, RIGHT_P], GAMMA);
            Sod.KIR_Calculate_Write(param, result_name_file);


            Plots.PLotAll();


            return;

        }
    }
}