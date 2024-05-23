using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sod
{
    public class KernelOperations
    {
        public static double[,] mult_matrixes(double[,] A, double[,] B)
        {
            double[,] C = new double[3, 3];
            int i, j, k;

            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    C[i, j] = 0.0;
                    for (k = 0; k < 3; k++)
                        C[i, j] += A[i, k] * B[k, j];
                }
            }
            return C;
        }

        public static double[] convert_noncons_to_cons(double[] v_ncons, double GAMMA)
        {
            double[] res = new double[(int)Vector_index.M];
            res[0] = v_ncons[0];                                                                     /* масса */
            res[1] = v_ncons[0] * v_ncons[1];                                                        /* импульс */
            res[2] = 0.5 * v_ncons[0] * Math.Pow(v_ncons[1], 2.0) + v_ncons[2] / (GAMMA - 1.0);        /* полная энергия */
            return res;
        }

        public static double[] convert_cons_to_noncons(double[] v_cons, double GAMMA)
        {
            double[] v_ncons = new double[(int)Vector_index.M];
            v_ncons[0] = v_cons[0];                                                                     /* плотность */
            v_ncons[1] = v_cons[1] / v_cons[0];                                                         /* скорость */
            v_ncons[2] = (GAMMA - 1) * (v_cons[2] - 0.5 * v_cons[0] * Math.Pow(v_ncons[1], 2.0));      /* давление */
            return v_ncons;
        }

        public static double[] boundary(double[] v_cons, int wall)
        {
            double[] boun_v = new double[(int)Vector_index.M];
            if (wall == 0)
            {
                boun_v[0] = v_cons[0];
                boun_v[1] = -v_cons[1];
                boun_v[2] = v_cons[2];
            }
            if (wall == 1)
            {
                boun_v[0] = v_cons[0];
                boun_v[1] = v_cons[1];
                boun_v[2] = v_cons[2];
            }
            return boun_v;
        }

        public static double calc_time_step(double[] x, double[,] v_cons, int time_step_number, int N, double GAMMA)
        {
            int i;
            double new_step = 1e10;
            double[] v_ncons = new double[(int)Vector_index.M];
            double c;
            double curr_step;

            for (i = 0; i < N; i++)
            {
                double[] t = new double[3];
                for (int j = 0; j < 3; j++)
                    t[j] = v_cons[i, j];
                v_ncons = convert_cons_to_noncons(t, GAMMA);
                c = calc_sound_velocity(v_ncons, GAMMA);
                curr_step = (x[i + 1] - x[i]) / (Math.Abs(v_ncons[1]) + c);
                if (time_step_number < 5)
                {
                    curr_step *= 0.2;
                }
                else
                {
                    curr_step *= 0.2;
                }
                if (curr_step < new_step)
                    new_step = curr_step;
            }

            return new_step;

        }

        public static double calc_sound_velocity(double[] v_ncons, double GAMMA)
        {
            return (Math.Sqrt(GAMMA * v_ncons[2] / v_ncons[0]));
        }


    }
}
