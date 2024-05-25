using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CsvHelper;
using CsvHelper.Configuration.Attributes;
using System.IO;
using System;
using OxyPlot;
using OxyPlot.Series;
using OxyPlot.Axes;
using Sod.Utility;

namespace Sod.CFL
{
    public class CIR : Plots
    {
        Parameters param;
        public CIR(Parameters par) {
            param = par;
        }

        public double[,] CalculateWrite(string result_name_file, bool write = true, int presision = 3, double CFL = 0.2)
        {
            int boundary = 1;
            result_name_file += "CIR_";

            double GAMMA = param.g;
            double T_END = param.stop_time;

            double LEFT_R = param.left_params[(int)Vector_index.R];
            double LEFT_U = param.left_params[(int)Vector_index.V];
            double LEFT_P = param.left_params[(int)Vector_index.P];

            double RIGHT_R = param.right_params[(int)Vector_index.R];
            double RIGHT_U = param.right_params[(int)Vector_index.V];
            double RIGHT_P = param.right_params[(int)Vector_index.P];
            int N = param.cells_number;
            int M = (int)Vector_index.M;


            double[] x = new double[N + 1];
            double[] xc = new double[N];
            double[,] u_prev = new double[N, M];
            double[,] u_next = new double[N, M];
            double[] boun_v = new double[M];
            double[] flux_left = new double[M];
            double[] flux_right = new double[M];
            double dt = 1e-6;
            int steps_num = 0;
            double curr_t = 0.0;
            double h = 1.0 / N;

            for (int i = 0; i < N + 1; i++)
            {
                x[i] = i * h;
            }

            for (int i = 0; i < N; i++)
            {
                xc[i] = 0.5 * (x[i] + x[i + 1]);
            }


            double[] v = new double[M];

            for (int i = 0; i < N; i++)
            {
                if (i < Math.Round(N * 0.5))
                {
                    v[0] = LEFT_R;
                    v[1] = LEFT_U;
                    v[2] = LEFT_P;
                }
                else
                {
                    v[0] = RIGHT_R;
                    v[1] = RIGHT_U;
                    v[2] = RIGHT_P;
                }
                var temp = ConvertPrimToConservative(v, GAMMA);
                for (int j = 0; j < M; j++)
                {
                    u_prev[i, j] = temp[j];
                }
            }
            if (write)
                PlotSod(param, xc, u_prev, curr_t, N, GAMMA, result_name_file);
            while (T_END - curr_t > 0)
            {

                dt = CFL*CalculateTimeStep(x, u_prev, N, GAMMA);
                for (int i = 0; i < N; i++)
                {
                    if (i != 0)
                    {
                        double[] temp1 = new double[M];
                        double[] temp2 = new double[M];
                        for (int j = 0; j < M; j++)
                        {
                            temp1[j] = u_prev[i - 1, j];
                            temp2[j] = u_prev[i, j];
                        }
                        flux_left = CalculateFlux(temp1, temp2, GAMMA);
                    }
                    else
                    {
                        double[] temp = new double[M];
                        for (int j = 0; j < M; j++)
                        {
                            temp[j] = u_prev[0, j];
                        }
                        boun_v = KernelOperations.BoundaryCondition(temp, boundary);
                        flux_left = CalculateFlux(boun_v, temp, GAMMA);
                    }

                    if (i != N - 1)
                    {
                        double[] temp1 = new double[M];
                        double[] temp2 = new double[M];
                        for (int j = 0; j < M; j++)
                        {
                            temp1[j] = u_prev[i, j];
                            temp2[j] = u_prev[i + 1, j];
                        }

                        flux_right = CalculateFlux(temp1, temp2, GAMMA);
                    }
                    else
                    {
                        double[] temp = new double[M];
                        for (int j = 0; j < M; j++)
                        {
                            temp[j] = u_prev[N - 1, j];
                        }
                        boun_v = KernelOperations.BoundaryCondition(temp, boundary);
                        flux_right = CalculateFlux(temp, boun_v, GAMMA);
                    }

                    for (int j = 0; j < M; j++)
                        u_next[i, j] = u_prev[i, j] - dt * (flux_right[j] - flux_left[j]) / (x[i + 1] - x[i]);
                }

                for (int i = 0; i < N; i++)
                    for (int j = 0; j < M; j++)
                        u_prev[i, j] = u_next[i, j];

                curr_t += dt;
                steps_num += 1;
                if (write)
                    if (0.0001 < curr_t & curr_t < 0.1 & Math.Round(curr_t % 0.005, presision) == 0)
                    {
                        PlotSod(param, xc, u_prev, curr_t, N, GAMMA, result_name_file);
                        
                    }
                    else
                    {
                        if (curr_t >= 0.1 & Math.Round(curr_t % 0.1, presision) == 0)
                        {
                            PlotSod(param, xc, u_prev, curr_t, N, GAMMA, result_name_file);
                        }
                    }
            }
            return u_prev;
        }

        static double[] CalculateFlux(double[] left_params, double[] right_params, double GAMMA)
        {
            ConvertConsToPrimitive(left_params, GAMMA);
            ConvertConsToPrimitive(right_params, GAMMA);
            var leftspeed = CalculateSoundVelocity(left_params, GAMMA)+left_params[1];
            var rightspeed = CalculateSoundVelocity(right_params, GAMMA)+right_params[1];
            var maxspeed = Math.Max(leftspeed, rightspeed);
            ConvertPrimToConservative(left_params, GAMMA);
            ConvertPrimToConservative(right_params, GAMMA);

            int M = (int)Vector_index.M;
            int i, j;
            double[] flux = new double[M];

            var left_diff_flux = ConservativeDifference(left_params, GAMMA);
            var right_diff_flux = ConservativeDifference(right_params, GAMMA);

            double[,] omega_r_L_omega_l = new double[3,3];
            //rusanov
            for (i = 0; i < M; i++)
                for (j = 0; j < M; j++)
                {
                    omega_r_L_omega_l[i, j] = 0;
                    if (i==j)
                        omega_r_L_omega_l[i, j] = maxspeed;
                }

            for (i = 0; i < M; i++)
            {
                flux[i] = 0.5 * (left_diff_flux[i] + right_diff_flux[i]);
                for (j = 0; j < M; j++)
                    flux[i] += 0.5 * (omega_r_L_omega_l[i, j]) * (left_params[j] - right_params[j]);
            }

            return flux;
        }

    }
}
