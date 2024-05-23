using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Windows.Forms.VisualStyles.VisualStyleElement.TaskbarClock;

namespace Sod
{
    public class Godunov : Plots
    {
        
        Parameters param;
        public Godunov(Parameters par)
        {
            param = par;
        }

        
        public double[,] Calculate_Write(string result_name_file, bool write = true, int presision = 3)
        {
            int boundary = 0;
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
                var temp = convert_noncons_to_cons(v, GAMMA);
                for (int j = 0; j < M; j++)
                {
                    u_prev[i, j] = temp[j];
                }
            }
            if (write)
                PlotSodOxyPlot(param, xc, u_prev, curr_t, N, GAMMA, result_name_file);
            while (T_END - curr_t > 0)
            {

                dt = calc_time_step(x, u_prev, steps_num, N, GAMMA);
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
                        flux_left = calc_flux(temp1, temp2, GAMMA, dt);
                    }
                    else
                    {
                        double[] temp = new double[M];
                        for (int j = 0; j < M; j++)
                        {
                            temp[j] = u_prev[0, j];
                        }
                        boun_v = KernelOperations.boundary(temp, boundary);
                        flux_left = calc_flux(boun_v, temp, GAMMA, dt);
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

                        flux_right = calc_flux(temp1, temp2, GAMMA, dt);
                    }
                    else
                    {
                        double[] temp = new double[M];
                        for (int j = 0; j < M; j++)
                        {
                            temp[j] = u_prev[N - 1, j];
                        }
                        boun_v = KernelOperations.boundary(temp, boundary);
                        flux_right = calc_flux(temp, boun_v, GAMMA, dt);
                    }

                    for (int j = 0; j < M; j++)
                        u_next[i, j] = u_prev[i, j] - dt * (flux_right[j] - flux_left[j]) / (x[i + 1] - x[i]);
                }

                for (int i = 0; i < N; i++)
                    for (int j = 0; j < M; j++)
                        u_prev[i, j] = u_next[i, j];

                curr_t += dt;
                steps_num += 1;
                PlotSodOxyPlot(param, xc, u_prev, curr_t, N, GAMMA, result_name_file);
                if (write)
                    if (0.01 < curr_t & curr_t < 0.1 & Math.Round(curr_t % 0.02, presision) == 0)
                    {
                        PlotSodOxyPlot(param, xc, u_prev, curr_t, N, GAMMA, result_name_file);

                    }
                    else
                    {
                        if (curr_t >= 0.1 & Math.Round(curr_t % 0.1, presision) == 0)
                        {
                            PlotSodOxyPlot(param, xc, u_prev, curr_t, N, GAMMA, result_name_file);
                        }
                    }
            }
            return u_prev;
        }
        
        static double[] calc_flux(double[] left_params, double[] right_params, double GAMMA, double dt)
        {
            throw new NotImplementedException();
            var param = new Parameters(101, dt, left_params, right_params, GAMMA);
            var res = ExactSolution.Calculate(param);

            return [res[50,0], res[50, 1], res[50, 2]];
        }

    }

}
