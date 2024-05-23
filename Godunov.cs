using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Sod
{
    public class Godunov : KIR
    {
        public Godunov(Parameters par):
            base(par)
        {
        }

        public new double[,] Calculate_Write(string result_name_file, bool write = true)
        {
            int boundary = 1;
            result_name_file += "KIR_";

            return null;
        }
    }
}
