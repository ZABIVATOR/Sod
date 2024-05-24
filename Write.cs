using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sod
{
    public abstract class Write : KernelOperations
    {
        public static void writeCSV_solution(double[] xc, double[,] v_cons, double time, int N, double GAMMA, string name)
        {
            StreamWriter f = new StreamWriter(name + Math.Round(time * 1000) + "ms" + ".csv");
            var csv = new CsvHelper.CsvWriter(f, CultureInfo.InvariantCulture);

            csv.WriteField("X");
            csv.WriteField("density");
            csv.WriteField("speed");
            csv.WriteField("pressure");
            csv.WriteField("energy");
            csv.NextRecord();
            double[] v_ncons = new double[3];
            int i, j;
            for (i = 0; i < N; i++)
            {
                csv.WriteField(xc[i]);
                double[] t = new double[3];
                for (j = 0; j < 3; j++)
                    t[j] = v_cons[i, j];
                v_ncons = ConvertConsToPrimitive(t, GAMMA);
                for (j = 0; j < 3; j++)
                {
                    csv.WriteField(v_ncons[j]);
                    f.Write(" ");
                }
                var te = v_ncons[2] / v_ncons[0] / (GAMMA - 1.0);
                csv.WriteField(te);

                csv.Configuration.Delimiter.Replace(',', ';');
                csv.NextRecord();
            }
            f.Close();
        }

        public static void write_solution(double[] xc, double[,] v_cons, double time, int N, double GAMMA)
        {
            StreamWriter f = new StreamWriter("sod" + Math.Round(time, 2) + ".txt");
            double[] v_ncons = new double[3];
            int i, j;
            for (i = 0; i < N; i++)
            {
                f.Write(xc[i]);
                f.Write(" ");
                double[] t = new double[3];
                for (j = 0; j < 3; j++)
                    t[j] = v_cons[i, j];
                v_ncons = ConvertConsToPrimitive(t, GAMMA);
                for (j = 0; j < 3; j++)
                {
                    f.Write(v_ncons[j]);
                    f.Write(" ");
                }
                var te = v_ncons[2] / v_ncons[0] / (GAMMA - 1.0);
                f.Write(te);
                f.Write(" ");
                f.Write("\n");
            }
            f.Close();
        }

    }
}
