using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sod
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
        // СЕТКА ДЛЯ ПОСТРОЕНИЯ ТОЧНОГО РЕШЕНИЯ  

        public int cells_number { get; set; }        // число ячеек  

        // ВРЕМЯ ВЫВОДА  

        public double stop_time { get; set; }            // момент времени, для которого строится точное решение  

        // НАЧАЛЬНЫЕ УСЛОВИЯ  

        public double[] left_params { get; set; }     // вектор примитивных переменных слева от разрыва  
        public double[] right_params { get; set; }    // вектор примитивных переменных справа от разрыва  

        // УРАВНЕНИе СОСТОЯНИЯ  

        public double g { get; set; }                     // показатель адибаты  

        public Parameters(int cells_number, double stop_time, double[] left_params, double[] right_params, double g)
        {
            this.cells_number = cells_number;
            this.stop_time = stop_time;
            this.left_params = left_params;
            this.right_params = right_params;
            this.g = g;
        }

    };
}
