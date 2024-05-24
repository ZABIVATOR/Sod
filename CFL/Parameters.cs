using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sod.CFL
{
    public enum Vector_index
    {
        R = 0, // плотность газовой фазы  
        V = 1, // скорость газовой фазы  
        P = 2,  // давление газовой фазы  
        M = 3   //Размерность вектора решений
    };

    public struct Parameters
    {
        public int cells_number { get; set; }        // число ячеек  
        public double stop_time { get; set; }            // момент времени, для которого строится точное решение  
        public double[] left_params { get; set; }     // вектор примитивных переменных слева от разрыва  
        public double[] right_params { get; set; }    // вектор примитивных переменных справа от разрыва  
        public double g { get; set; }                     // показатель адибаты  

        public Parameters(int cells_number, double stop_time, double[] left_params, double[] right_params, double g)
        {
            this.cells_number = cells_number;
            this.stop_time = stop_time;
            this.left_params = left_params;
            this.right_params = right_params;
            this.g = g;
        }

    }

    public interface CFLmethod
    {
        public double[,] CalculateWrite(string result_name_file, bool write = true, int presision = 3);
    }
}
