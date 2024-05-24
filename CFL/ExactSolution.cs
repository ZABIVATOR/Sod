using System;
using System.Collections.Generic;
using System.Linq;
using System.Linq.Expressions;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Documents;

namespace Sod.CFL
{
    public class ExactSolution
    {
        static int MaxIter = 20;
        static double EPS = 1e-6;
        static double P_MAX_RATIO = 2.0;           /* перепад давлений слева и справа, при котором начинает использоваться начальное приближение
                                                       по двум ударным волнам или двум волнам разрежения */
        public ExactSolution()
        {

        }

        public static void BuildGrid(ref double[] xc, ref double[] x, double left_boundary_x, double right_boundary_x, int cells_num)
        {
            int i;
            double h;   // шаг сетки  
            // пока реализовано только равномерное распределение узлов сетки  
            h = (right_boundary_x - left_boundary_x) / cells_num;
            // координаты узлов  
            for (i = 0; i < cells_num + 1; i++)
            {
                x[i] = left_boundary_x + i * h;
            }
            // координаты центров ячеек  
            for (i = 0; i < cells_num; i++)
            {
                xc[i] = 0.5 * (x[i] + x[i + 1]);
            }
        }

        static double CalculateSoundVelocity(Parameters param, double[] v_ncons)
        {
            return Math.Sqrt(param.g * v_ncons[(int)Vector_index.P] / v_ncons[(int)Vector_index.R]);
        }

        static void CalculateFdF(ref double F, ref double DF, Parameters param, double curr_press, double[] v_ncons, double c)
        {

            double r, v, p;                 // примитивные переменные  
            double g;                       // показатель адиабаты  
            double p_ratio, fg, q;          // вспомогательные переменные  

            r = v_ncons[(int)Vector_index.R];
            v = v_ncons[(int)Vector_index.V];
            p = v_ncons[(int)Vector_index.P];
            g = param.g;

            p_ratio = curr_press / p;
            if (curr_press <= p)
            {
                // волна разрежения  
                fg = 2.0 / (g - 1.0);
                F = fg * c * (Math.Pow(p_ratio, 1.0 / fg / g) - 1.0);
                DF = 1.0 / r / c * Math.Pow(p_ratio, -0.5 * (g + 1.0) / g);
            }
            else
            {
                // ударная волна  
                q = Math.Sqrt(0.5 * (g + 1.0) / g * p_ratio + 0.5 * (g - 1.0) / g);
                F = (curr_press - p) / c / r / q;
                DF = 0.25 * ((g + 1.0) * p_ratio + 3 * g - 1.0) / g / r / c / Math.Pow(q, 3.0);
            }
        }

        static double InitialPressure(Parameters param, double[] v_ncons_l, double[] v_ncons_r, double cl, double cr)
        {

            double rl, vl, pl;                  // примитивные переменные слева от разрыва  
            double rr, vr, pr;                  // примитивные переменные справа от разрыва  
            double g;                           // показатель адиабаты  
            // начальное приближение, рассчитанное на освановании рассмотрения линеаризованной системы в примитивных переменных
            double p_lin;
            double p_min, p_max;                // минимальное и максимальное давления слева и справа от разрыва  
            double p_ratio;                     // перепад по давлению слева и справа от разрыва  
            double p1, p2, g1, g2;              // вспомогательные переменные для промежуточных расчетов  

            rl = v_ncons_l[(int)Vector_index.R];
            vl = v_ncons_l[(int)Vector_index.V];
            pl = v_ncons_l[(int)Vector_index.P];
            rr = v_ncons_r[(int)Vector_index.R];
            vr = v_ncons_r[(int)Vector_index.V];
            pr = v_ncons_r[(int)Vector_index.P];
            g = param.g;

            //Начальное приближение из линейной задачи
            p_lin = Math.Max(0.0, 0.5 * (pl + pr) - 0.125 * (vr - vl) * (rl + rr) * (cl + cr));
            p_min = Math.Min(pl, pr);
            p_max = Math.Max(pl, pr);
            p_ratio = p_max / p_min;

            if (p_ratio <= P_MAX_RATIO &&
               (p_min < p_lin && p_lin < p_max || Math.Abs(p_min - p_lin) < EPS || Math.Abs(p_max - p_lin) < EPS))
            {
                // Начальное приближение из линеаризованной задачи  
                return p_lin;
            }
            else
            {
                if (p_lin < p_min)
                {
                    // Начальное приближение по двум волнам разрежения
                    g1 = 0.5 * (g - 1.0) / g;
                    return Math.Pow((cl + cr - 0.5 * (g - 1.0) * (vr - vl)) / (cl / Math.Pow(pl, g1) + cr / Math.Pow(pr, g1)), 1.0 / g1);
                }
                else
                {
                    // Начальное приближение по двум ударным волнам
                    g1 = 2.0 / (g + 1.0);
                    g2 = (g - 1.0) / (g + 1.0);
                    p1 = Math.Sqrt(g1 / rl / (g2 * pl + p_lin));
                    p2 = Math.Sqrt(g1 / rr / (g2 * pr + p_lin));
                    return (p1 * pl + p2 * pr - (vr - vl)) / (p1 + p2);
                }
            }
        }

        static void CalculateContactParameters(ref double p_cont, ref double v_cont, Parameters param,
        double[] v_ncons_l, double[] v_ncons_r, double cl, double cr)
        {

            double vl, vr;      // скорости слева и справа от разрыва  
            double p_old;       // значение давления на предыдущей итерации  
            double fl = 0, fr = 0;      // значения функций  
            double fld = 0, frd = 0;    // значения производных  
            int iter_num = 0;   // количество проведенных итераций  
            double criteria;    // переменная для определения сходимости  
            double g;           // показатель адиабаты  

            // введение обозначений для удобства записи формул  
            vl = v_ncons_l[(int)Vector_index.V];
            vr = v_ncons_r[(int)Vector_index.V];
            g = param.g;

            if (2.0 * (cl + cr) / (g - 1.0) <= vr - vl)
            {
                // случай возникновения вакуума  
                throw new Exception("vacuum is generated");
            }

            // расчет начального приближения для давления  
            p_old = InitialPressure(param, v_ncons_l, v_ncons_r, cl, cr);
            if (p_old < 0.0)
            {
                throw new Exception("initial pressure guess is negative");
            }

            // решение нелинейного уравнения для нахождения давления на контактном разрыве методом Ньютона-Рафсона  
            do
            {
                CalculateFdF(ref fl, ref fld, param, p_old, v_ncons_l, cl);
                CalculateFdF(ref fr, ref frd, param, p_old, v_ncons_r, cr);
                p_cont = p_old - (fl + fr + vr - vl) / (fld + frd);
                criteria = 2.0 * Math.Abs((p_cont - p_old) / (p_cont + p_old));
                iter_num++;
                if (iter_num > MaxIter)
                {
                    throw new Exception("number of iterations exceeds the maximum value.");
                }
                if (p_cont < 0.0)
                {
                    throw new Exception("pressure is negative.");
                }
                p_old = p_cont;
            } while (criteria > EPS);

            // скорость контактного разрыва  
            v_cont = 0.5 * (vl + vr + fr - fl);

        }

        static double[] SampleSolution(Parameters param, double[] v_ncons_l, double[] v_ncons_r, double cl, double cr,
            double p_cont, double v_cont, double s)
        {
            double[] v_ncons_res = new double[(int)Vector_index.M];

            double rl, vl, pl;                      // примитивные переменные слева от разрыва  
            double rr, vr, pr;                      // примитивные переменные справа от разрыва  
            double g1, g2, g3, g4, g5, g6, g7;      // вспомогательные переменные, производные от показателя адиабаты

            // скорости левых волн  
            double shl, stl;        // скорости "головы" и "хвоста" левой волны разрежения  
            double sl;              // скорость левой ударной волны  

            // скорости правых волн  
            double shr, str;        // скорости "головы" и "хвоста" правой волны разрежения  
            double sr;              // скорость правой ударной волны  

            double cml, cmr;        // скорости звука слева и справа от контактного разрыва  
            double c;               // локальная скорость звука внутри волны разрежения  
            double p_ratio;
            double r, v, p;         // отобранные значения объемной доли, плотности, скорости и давления  

            // вспомогательные переменные  
            // параметры слева от разрыва  
            rl = v_ncons_l[(int)Vector_index.R];
            vl = v_ncons_l[(int)Vector_index.V];
            pl = v_ncons_l[(int)Vector_index.P];
            // параметры справа от разрыва  
            rr = v_ncons_r[(int)Vector_index.R];
            vr = v_ncons_r[(int)Vector_index.V];
            pr = v_ncons_r[(int)Vector_index.P];
            // производные от показателя адиабаты  
            g1 = 0.5 * (param.g - 1.0) / param.g;
            g2 = 0.5 * (param.g + 1.0) / param.g;
            g3 = 2.0 * param.g / (param.g - 1.0);
            g4 = 2.0 / (param.g - 1.0);
            g5 = 2.0 / (param.g + 1.0);
            g6 = (param.g - 1.0) / (param.g + 1.0);
            g7 = 0.5 * (param.g - 1.0);

            if (s <= v_cont)
            {
                // рассматриваемая точка - слева от контактного разрыва  
                if (p_cont <= pl)
                {
                    // левая волна разрежения  
                    shl = vl - cl;
                    if (s <= shl)
                    {
                        // параметры слева от разрыва  
                        r = rl;
                        v = vl;
                        p = pl;
                    }
                    else
                    {
                        cml = cl * Math.Pow(p_cont / pl, g1);
                        stl = v_cont - cml;
                        if (s > stl)
                        {
                            // параметры слева от контактного разрыва  
                            r = rl * Math.Pow(p_cont / pl, 1.0 / param.g);
                            v = v_cont;
                            p = p_cont;
                        }
                        else
                        {
                            // параметры внутри левой волны разрежения  
                            v = g5 * (cl + g7 * vl + s);
                            c = g5 * (cl + g7 * (vl - s));
                            r = rl * Math.Pow(c / cl, g4);
                            p = pl * Math.Pow(c / cl, g3);
                        }
                    }
                }
                else
                {
                    // левая ударная волна  
                    p_ratio = p_cont / pl;
                    sl = vl - cl * Math.Sqrt(g2 * p_ratio + g1);
                    if (s <= sl)
                    {
                        // параметры слева от разрыва  
                        r = rl;
                        v = vl;
                        p = pl;
                    }
                    else
                    {
                        // параметры за левой ударной волной  
                        r = rl * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                        v = v_cont;
                        p = p_cont;
                    }
                }
            }
            else
            {
                // рассматриваемая точка - справа от контактного разрыва  
                if (p_cont > pr)
                {
                    // правая ударная волна  
                    p_ratio = p_cont / pr;
                    sr = vr + cr * Math.Sqrt(g2 * p_ratio + g1);
                    if (s >= sr)
                    {
                        // параметры справа от разрыва  
                        r = rr;
                        v = vr;
                        p = pr;
                    }
                    else
                    {
                        // параметры за правой ударной волной  
                        r = rr * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                        v = v_cont;
                        p = p_cont;
                    }
                }
                else
                {
                    // правая волна разрежения  
                    shr = vr + cr;
                    if (s >= shr)
                    {
                        // параметры справа от разрыва  
                        r = rr;
                        v = vr;
                        p = pr;
                    }
                    else
                    {
                        cmr = cr * Math.Pow(p_cont / pr, g1);
                        str = v_cont + cmr;
                        if (s <= str)
                        {
                            // параметры справа от контактного разрыва  
                            r = rr * Math.Pow(p_cont / pr, 1.0 / param.g);
                            v = v_cont;
                            p = p_cont;
                        }
                        else
                        {
                            // параметры внутри правой волны разрежения  
                            v = g5 * (-cr + g7 * vr + s);
                            c = g5 * (cr - g7 * (vr - s));
                            r = rr * Math.Pow(c / cr, g4);
                            p = pr * Math.Pow(c / cr, g3);
                        }
                    }
                }
            }

            // формирование выходного вектора с результатом  
            v_ncons_res[(int)Vector_index.R] = r;
            v_ncons_res[(int)Vector_index.V] = v;
            v_ncons_res[(int)Vector_index.P] = p;
            return v_ncons_res;
        }

        public static double[,] Calculate(Parameters SODparam, double curr_time = -1)//отрицательное время для расчета по моменту из структуры параметров
        {
            var param = new Parameters();
            if (curr_time < 0)
            {
                param = SODparam;
            }
            else
                param = new Parameters(SODparam.cells_number, curr_time, SODparam.left_params, SODparam.right_params, SODparam.g);

            double[,] v_ncons_res = new double[param.cells_number, (int)Vector_index.M];      // решение  

            double[] xc = new double[param.cells_number];                 // массив координат центров ячеек сетки  
            double[] x = new double[param.cells_number + 1];                 // массив координат узлов сетки  
            double cl, cr;              // скорости звука слева и справа от разрыва  
            double p_cont = 0, v_cont = 0;      // давление и скорость на контактном разрыве  
            double s;                   // текущее значение автомодельной переменной  

            // определение координат центров ячеек сетки  
            BuildGrid(ref xc, ref x, -0.5, 0.5, param.cells_number);

            // расчет скоростей звука слева и справа от разрыва  
            cl = CalculateSoundVelocity(param, param.left_params);
            cr = CalculateSoundVelocity(param, param.right_params);

            // расчет давления и скорости на контактном разрыве  
            CalculateContactParameters(ref p_cont, ref v_cont, param, param.left_params, param.right_params, cl, cr);

            // цикл по ячейкам  
            for (int i_cell = 0; i_cell < param.cells_number; i_cell++)
            {
                // изначально решение строится на отрезке [-0.5;0.5]  
                s = xc[i_cell] / param.stop_time;
                // отбор решения для заданного значения s  
                var temp = SampleSolution(param, param.left_params, param.right_params, cl, cr, p_cont, v_cont, s);
                for (int j = 0; j < (int)Vector_index.M; j++)
                {
                    v_ncons_res[i_cell, j] = temp[j];
                }
            }
            return v_ncons_res;
        }
    }
}
