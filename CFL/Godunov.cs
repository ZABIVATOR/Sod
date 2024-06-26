﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Sod.Utility;
using static System.Windows.Forms.VisualStyles.VisualStyleElement.TaskbarClock;


namespace Sod.CFL
{
    public class Godunov : Plots
    {

        Parameters param;
        public Godunov(Parameters par)
        {
            param = par;
        }

        public double[,] CalculateWrite(string result_name_file, bool write = true, int presision = 3, double CFL = 0.2, string path = ""  )
        {
            int boundary = 1;
            result_name_file += "Godunov_";

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
                PlotSod(param, xc, u_prev, curr_t, N, GAMMA, result_name_file, path);
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
                        boun_v = BoundaryCondition(temp, boundary);
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
                        boun_v = BoundaryCondition(temp, boundary);
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
                        PlotSod(param, xc, u_prev, curr_t, N, GAMMA, result_name_file, path);

                    }
                    else
                    {
                        if (curr_t >= 0.1 & Math.Round(curr_t % 0.1, presision) == 0)
                        {
                            PlotSod(param, xc, u_prev, curr_t, N, GAMMA, result_name_file, path);
                        }
                    }
            }
            return u_prev;
        }

        static double[] CalculateFlux(double[] left_params, double[] right_params, double GAMMA)
        {
            //int i;
            //double h;   // шаг сетки  
            //h = (right_boundary_x - left_boundary_x) / cells_num;
            //координаты узлов
            //for (i = 0; i < cells_num + 1; i++)
            //{
            //    x[i] = left_boundary_x + i * h;
            //}
            //координаты центров ячеек
            //for (i = 0; i < cells_num; i++)
            //{
            //    xc[i] = 0.5 * (x[i] + x[i + 1]);
            //}

            var parameters = new Parameters((int)Vector_index.M,1e-6,left_params,right_params,GAMMA);
            var exact = ExactSolution.Calculate(parameters, GAMMA);
            exact = ConvertConsToPrimitive(exact,GAMMA);

            //double r, v, p;                 // примитивные переменные  
            //double g;                       // показатель адиабаты  
            //double p_ratio, fg, q;          // вспомогательные переменные  

            //r = v_ncons[(int)Vector_index.R];
            //v = v_ncons[(int)Vector_index.V];
            //p = v_ncons[(int)Vector_index.P];
            //g = param.g;

            //p_ratio = curr_press / p;
            //if (curr_press <= p)
            //{
            //     волна разрежения  
            //    fg = 2.0 / (g - 1.0);
            //    F = fg * c * (Math.Pow(p_ratio, 1.0 / fg / g) - 1.0);
            //    DF = 1.0 / r / c * Math.Pow(p_ratio, -0.5 * (g + 1.0) / g);
            //}
            //else
            //{
            //     ударная волна  
            //    q = Math.Sqrt(0.5 * (g + 1.0) / g * p_ratio + 0.5 * (g - 1.0) / g);
            //    F = (curr_press - p) / c / r / q;
            //    DF = 0.25 * ((g + 1.0) * p_ratio + 3 * g - 1.0) / g / r / c / Math.Pow(q, 3.0);
            //}

            return res(left_params, right_params, GAMMA);
        }

    }
}
