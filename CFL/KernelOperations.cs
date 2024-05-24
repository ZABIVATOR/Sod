using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sod.CFL
{
    public class KernelOperations
    {
        protected static double[] ConservativeDifference(double[] v_cons, double GAMMA)
        {
            int M = (int)Vector_index.M;
            double[] flux = new double[M];
            double[] v_ncons = new double[M];

            v_ncons = ConvertConsToPrimitive(v_cons, GAMMA);

            flux[0] = v_cons[1];                                    /* масса */
            flux[1] = v_cons[1] * v_ncons[1] + v_ncons[2];          /* импульс */
            flux[2] = (v_cons[2] + v_ncons[2]) * v_ncons[1];      /* полная энергия */
            return flux;
        }

        protected static double[,] MultuplyMatrixes(double[,] A, double[,] B)
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

        protected static double[] ConvertPrimToConservative(double[] v_ncons, double GAMMA)
        {
            double[] res = new double[(int)Vector_index.M];
            res[0] = v_ncons[0];                                                                     /* масса */
            res[1] = v_ncons[0] * v_ncons[1];                                                        /* импульс */
            res[2] = 0.5 * v_ncons[0] * Math.Pow(v_ncons[1], 2.0) + v_ncons[2] / (GAMMA - 1.0);        /* полная энергия */
            return res;
        }

        protected static double[,] ConvertPrimToConservative(double[,] v_ncons, double GAMMA)
        {
            double[,] res = new double[v_ncons.GetLength(0), v_ncons.GetLength(1)];
            for (int i = 0; i < v_ncons.GetLength(0); i++)
            {
                double[] t1 = new double[(int)Vector_index.M];
                for (int j = 0; j < (int)Vector_index.M; j++)
                {
                    t1[j] = v_ncons[i, j];
                }
                var temp = ConvertPrimToConservative(t1, GAMMA);
                for (int j = 0; j < (int)Vector_index.M; j++)
                {
                    res[i, j] = temp[j];
                }
            }
            return res;
        }

        protected static double[] ConvertConsToPrimitive(double[] v_cons, double GAMMA)
        {
            double[] v_ncons = new double[(int)Vector_index.M];
            v_ncons[0] = v_cons[0];                                                                     /* плотность */
            v_ncons[1] = v_cons[1] / v_cons[0];                                                         /* скорость */
            v_ncons[2] = (GAMMA - 1) * (v_cons[2] - 0.5 * v_cons[0] * Math.Pow(v_ncons[1], 2.0));      /* давление */
            return v_ncons;
        }

        protected static double[,] ConvertConsToPrimitive(double[,] v_cons, double GAMMA)
        {
            double[,] res = new double[v_cons.GetLength(0), v_cons.GetLength(1)];
            for (int i = 0; i < v_cons.GetLength(0); i++)
            {
                double[] t1 = new double[(int)Vector_index.M];
                for (int j = 0; j < (int)Vector_index.M; j++)
                {
                    t1[j] = v_cons[i, j];
                }
                var temp = ConvertConsToPrimitive(t1, GAMMA);
                for (int j = 0; j < (int)Vector_index.M; j++)
                {
                    res[i, j] = temp[j];
                }
            }
            return res;
        }

        protected static double[] BoundaryCondition(double[] v_cons, int wall)
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

        protected static double CalculateTimeStep(double[] x, double[,] v_cons, int time_step_number, int N, double GAMMA)
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
                v_ncons = ConvertConsToPrimitive(t, GAMMA);
                c = CalculateSoundVelocity(v_ncons, GAMMA);
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

        protected static double CalculateSoundVelocity(double[] v_ncons, double GAMMA)
        {
            return Math.Sqrt(GAMMA * v_ncons[(int)Vector_index.P] / v_ncons[(int)Vector_index.R]);
        }

        protected static double[,] calc_omega(double[] v_cons, double GAMMA)
        {
            int M = (int)Vector_index.M;
            double[,] omega = new double[M, M];
            double[] v_ncons = new double[M];
            double c;               /* скорость звука */
            double teta, b, h;

            v_ncons = ConvertConsToPrimitive(v_cons, GAMMA);

            c = CalculateSoundVelocity(v_ncons, GAMMA);

            teta = 0.5 * v_ncons[1] * v_ncons[1];
            b = GAMMA - 1.0;
            h = 0.5 * v_ncons[1] * v_ncons[1] + GAMMA * v_ncons[2] / v_ncons[0] / (GAMMA - 1.0);

            omega[0, 0] = 1.0;
            omega[0, 1] = 1.0;
            omega[0, 2] = 1.0;

            omega[1, 0] = v_ncons[1] - c;
            omega[1, 1] = v_ncons[1];
            omega[1, 2] = v_ncons[1] + c;

            omega[2, 0] = h - v_ncons[1] * c;
            omega[2, 1] = h - c * c / b;
            omega[2, 2] = h + v_ncons[1] * c;

            return omega;
        }

        protected static double[,] calc_omega_inverse(double[] v_cons, double GAMMA)
        {
            int M = (int)Vector_index.M;
            double[,] omega_inverse = new double[M, M];
            double[] v_ncons = new double[M];
            double c;               /* скорость звука */
            double teta, b, h;

            v_ncons = ConvertConsToPrimitive(v_cons, GAMMA);

            c = CalculateSoundVelocity(v_ncons, GAMMA);

            teta = 0.5 * v_ncons[1] * v_ncons[1];
            b = GAMMA - 1.0;
            h = 0.5 * v_ncons[1] * v_ncons[1] + GAMMA * v_ncons[2] / v_ncons[0] / (GAMMA - 1.0);

            omega_inverse[0, 0] = 0.5 * b * (teta + v_ncons[1] * c / b) / Math.Pow(c, 2.0);
            omega_inverse[0, 1] = 0.5 * b * (-v_ncons[1] - c / b) / Math.Pow(c, 2.0);
            omega_inverse[0, 2] = 0.5 * b / Math.Pow(c, 2.0);

            omega_inverse[1, 0] = 0.5 * b * (2.0 * h - 2.0 * Math.Pow(v_ncons[1], 2.0)) / Math.Pow(c, 2.0);
            omega_inverse[1, 1] = 0.5 * b * (2.0 * v_ncons[1]) / Math.Pow(c, 2.0);
            omega_inverse[1, 2] = 0.5 * b * -2.0 / Math.Pow(c, 2.0);

            omega_inverse[2, 0] = 0.5 * b * (teta - v_ncons[1] * c / b) / Math.Pow(c, 2.0);
            omega_inverse[2, 1] = 0.5 * b * (-v_ncons[1] + c / b) / Math.Pow(c, 2.0);
            omega_inverse[2, 2] = 0.5 * b / Math.Pow(c, 2.0);
            return omega_inverse;
        }

        protected static double[,] calc_lambda(double[] v_cons, double GAMMA)
        {
            int M = (int)Vector_index.M;
            double[,] lambda = new double[M, M];
            double[] v_ncons = new double[M];
            double c;
            v_ncons = ConvertConsToPrimitive(v_cons, GAMMA);
            c = CalculateSoundVelocity(v_ncons, GAMMA);

            lambda[0, 0] = Math.Abs(v_ncons[1] - c);
            lambda[0, 1] = 0.0;
            lambda[0, 2] = 0.0;

            lambda[1, 0] = 0.0;
            lambda[1, 1] = Math.Abs(v_ncons[1]);
            lambda[1, 2] = 0.0;

            lambda[2, 0] = 0.0;
            lambda[2, 1] = 0.0;
            lambda[2, 2] = Math.Abs(v_ncons[1] + c);
            return lambda;
        }
    }
}
