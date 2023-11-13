using CsvHelper;
using CsvHelper.Configuration.Attributes;
using System.Globalization;

class Sod
{
    static void writeCSV_solution(double[] xc, double[,] v_cons, double time, int N, double GAMMA,string name)
    {
        StreamWriter f = new StreamWriter(name + Math.Round(time * 1000) +"ms"+ ".csv");
        CsvWriter csv = new CsvWriter(f,CultureInfo.InvariantCulture);

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
            v_ncons = convert_cons_to_noncons(t, GAMMA);
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


    static void write_solution(double[] xc, double[,] v_cons, double time, int N, double GAMMA)
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
            v_ncons = convert_cons_to_noncons(t, GAMMA);
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

    static double calc_time_step(double[] x, double[,] v_cons, int time_step_number, int N, double GAMMA)
    {
        int i;
        double new_step = 1e10;   
        double[] v_ncons = new double[3];               
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

    static double[] convert_noncons_to_cons(double[] v_ncons, double GAMMA)
    {
        double[] res = new double[3];
        res[0] = v_ncons[0];                                                                     /* масса */
        res[1] = v_ncons[0] * v_ncons[1];                                                        /* импульс */
        res[2] = 0.5 * v_ncons[0] * Math.Pow(v_ncons[1], 2.0) + v_ncons[2] / (GAMMA - 1.0);        /* полная энергия */
        return res;
    }

    static double[] convert_cons_to_noncons(double[] v_cons, double GAMMA)
    {
        double[] v_ncons = new double[3];
        v_ncons[0] = v_cons[0];                                                                     /* плотность */
        v_ncons[1] = v_cons[1] / v_cons[0];                                                         /* скорость */
        v_ncons[2] = (GAMMA - 1) * (v_cons[2] - 0.5 * v_cons[0] * Math.Pow(v_ncons[1], 2.0));      /* давление */
        return v_ncons;
    }

    static double calc_sound_velocity(double[] v_ncons, double GAMMA)
    {
        return (Math.Sqrt(GAMMA * v_ncons[2] / v_ncons[0]));
    }

    static double[] diff_flux_cons(double[] v_cons, double GAMMA)
    {
        int M = 3;
        double[] flux = new double[M];
        double[] v_ncons = new double[M];                                      

        v_ncons = convert_cons_to_noncons(v_cons, GAMMA);

        flux[0] = v_cons[1];                                    /* масса */
        flux[1] = v_cons[1] * v_ncons[1] + v_ncons[2];          /* импульс */
        flux[2] = (v_cons[2] + v_ncons[2]) * v_ncons[1];      /* полная энергия */
        return flux;
    }
     
    static double[] boundary(double[] v_cons, int wall)
    {
        double[] boun_v = new double[3];
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

    static double[,] calc_omega(double[] v_cons, double GAMMA)
    {
        int M = 3;
        double[,] omega = new double[M, M];
        double[] v_ncons = new double[M];
        double c;               /* скорость звука */
        double teta, b, h;

        v_ncons = convert_cons_to_noncons(v_cons, GAMMA);

        c = calc_sound_velocity(v_ncons, GAMMA);

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
     
    static double[,] calc_omega_inverse(double[] v_cons, double GAMMA)
    {
        int M = 3;
        double[,] omega_inverse = new double[M, M];
        double[] v_ncons = new double[M];
        double c;               /* скорость звука */
        double teta, b, h;

        v_ncons = convert_cons_to_noncons(v_cons, GAMMA);

        c = calc_sound_velocity(v_ncons, GAMMA);

        teta = 0.5 * v_ncons[1] * v_ncons[1];
        b = GAMMA - 1.0;
        h = 0.5 * v_ncons[1] * v_ncons[1] + GAMMA * v_ncons[2] / v_ncons[0] / (GAMMA - 1.0);

        omega_inverse[0, 0] = 0.5 * b * (teta + v_ncons[1] * c / b) / Math.Pow(c, 2.0);
        omega_inverse[0, 1] = 0.5 * b * (-v_ncons[1] - c / b) / Math.Pow(c, 2.0);
        omega_inverse[0, 2] = 0.5 * b / Math.Pow(c, 2.0);

        omega_inverse[1, 0] = 0.5 * b * (2.0 * h - 2.0 * Math.Pow(v_ncons[1], 2.0)) / Math.Pow(c, 2.0);
        omega_inverse[1, 1] = 0.5 * b * (2.0 * v_ncons[1]) / Math.Pow(c, 2.0);
        omega_inverse[1, 2] = 0.5 * b * (-2.0) / Math.Pow(c, 2.0);

        omega_inverse[2, 0] = 0.5 * b * (teta - v_ncons[1] * c / b) / Math.Pow(c, 2.0);
        omega_inverse[2, 1] = 0.5 * b * (-v_ncons[1] + c / b) / Math.Pow(c, 2.0);
        omega_inverse[2, 2] = 0.5 * b / Math.Pow(c, 2.0);
        return omega_inverse;
    }

    static double[,] calc_lambda(double[] v_cons, double GAMMA)
    {
        int M = 3;
        double[,] lambda = new double[M, M];
        double[] v_ncons = new double[M];
        double c;              
        v_ncons = convert_cons_to_noncons(v_cons, GAMMA);
        c = calc_sound_velocity(v_ncons, GAMMA);

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

    static double[] calc_flux(double[] left_params, double[] right_params, double GAMMA)
    {
        int M = 3;
        int i, j;
        double[] flux = new double[M];
        double[,] m_left = new double[M, M];
        double[,] m_right = new double[M, M];

        var left_diff_flux = diff_flux_cons(left_params, GAMMA);
        var right_diff_flux = diff_flux_cons(right_params, GAMMA);

        m_left = cir_util(left_params, GAMMA);
        m_right = cir_util(right_params, GAMMA);

        for (i = 0; i < M; i++)
        {
            flux[i] = 0.5 * (left_diff_flux[i] + right_diff_flux[i]);
            for (j = 0; j < M; j++)
                flux[i] += 0.5 * (0.5 * (m_left[i, j] + m_right[i, j])) * (left_params[j] - right_params[j]);
        }
        return flux;
    }

    static double[,] cir_util(double[] cons_params, double GAMMA)
    {
        var omega = calc_omega(cons_params, GAMMA);
        var omega_inverse = calc_omega_inverse(cons_params, GAMMA);
        var lambda = calc_lambda(cons_params, GAMMA);
        
        var m_tmp = mult_matrixes(omega, lambda);
        return mult_matrixes(m_tmp, omega_inverse);
    }

    static double[,] mult_matrixes(double[,] A, double[,] B)
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



    static void Main()
    {

        string result_name_file = "sod_t1000_";
        double GAMMA = 1.4;       
        double T_END = 0.01 ;        

        double LEFT_R = 1;
        double LEFT_U = 0;
        double LEFT_P = 1;

        double RIGHT_R = 0.1;
        double RIGHT_U = 0;
        double RIGHT_P = 0.1;
        int N = 1000;           
        int M = 3;

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
            if (i < Math.Round(N*0.5))
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
        writeCSV_solution(xc, u_prev, curr_t, N, GAMMA,result_name_file);
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
                    flux_left = calc_flux(temp1, temp2, GAMMA);
                }
                else
                {
                    double[] temp = new double[M];
                    for (int j = 0; j < M; j++)
                    {
                        temp[j] = u_prev[0, j];
                    }
                    boun_v = boundary(temp,1);
                    flux_left = calc_flux(boun_v, temp, GAMMA);
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

                    flux_right = calc_flux(temp1, temp2, GAMMA);
                }
                else
                {
                    double[] temp = new double[M];
                    for (int j = 0; j < M; j++)
                    {
                        temp[j] = u_prev[N - 1, j];
                    }
                    boun_v = boundary(temp,1);
                    flux_right = calc_flux(temp, boun_v, GAMMA);
                }

                for (int j = 0; j < M; j++)
                    u_next[i, j] = u_prev[i, j] - dt * (flux_right[j] - flux_left[j]) / (x[i + 1] - x[i]);
            }

            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
                    u_prev[i, j] = u_next[i, j];

            curr_t += dt;
            steps_num += 1;

            if (Math.Round(curr_t % 0.2, 3) == 0)
            {
                
                writeCSV_solution(xc, u_prev, curr_t, N, GAMMA, result_name_file);
                Console.WriteLine("stpep num: "+steps_num+" time: "+curr_t);
            }

        }

        Console.WriteLine("end: "+steps_num);

        write_solution(xc, u_prev, curr_t, N, GAMMA);

        return;

    }
}
