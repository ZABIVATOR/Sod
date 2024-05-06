using OxyPlot.Axes;
using OxyPlot.Series;
using OxyPlot;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Sod.CFL;
using OxyPlot.Legends;
using System.Windows.Media.Animation;

namespace Sod.Utility
{
    public abstract class Plots : Write
    {
        public static void PlotBoth(Parameters param, double[,] v_cons1, double[,] v_cons2, string method,string method1, string method2, string path)
        {
            method = method + "both_";
            var time = param.stop_time;
            var N = param.cells_number;
            var GAMMA = param.g;

            double[] x = new double[N + 1];
            double[] xc = new double[N];
            for (int i = 0; i < N + 1; i++)
            {
                x[i] = i * 1.0 / N;
            }

            for (int i = 0; i < N; i++)
            {
                xc[i] = 0.5 * (x[i] + x[i + 1]);
            }


            double[,] exact = ExactSolution.Calculate(param, time);


            var name = path+ method + Math.Round(time * 1000) + "ms" + ".png";

            var pm = new PlotModel()
            {
                Title = $" time = {Math.Round(time * 1000)} ms"
            };
            var ls1 = new LineSeries()
            {
                Title = method1+" density",
                MarkerStrokeThickness = 0,
                MarkerStroke = OxyColors.Transparent,
                LineStyle = LineStyle.None,
                MarkerSize = 3,
                MarkerType = MarkerType.Circle,
                MarkerFill = OxyColors.Green
            };
            var ls1_ = new LineSeries()
            {
                Title = method2+" density",
                MarkerStrokeThickness = 0,
                MarkerStroke = OxyColors.Transparent,
                LineStyle = LineStyle.None,
                MarkerSize = 3,
                MarkerType = MarkerType.Triangle,
                MarkerFill = OxyColors.Red,
            };
                var ls1_exact = new LineSeries() { Title = "Exact density", MarkerStrokeThickness = 0.5, Color = OxyColors.Black };


            //поменял местами для красоты
            var ls2 = new LineSeries()
            {
                Title = method1 +" speed",
                MarkerStrokeThickness = 0,
                MarkerStroke = OxyColors.Transparent,
                LineStyle = LineStyle.None,
                MarkerSize = 3,
                MarkerType = MarkerType.Circle,
                MarkerFill = OxyColors.GreenYellow
            };
            var ls2_ = new LineSeries()
            {
                Title = method2+ " speed",
                MarkerStrokeThickness = 0,
                MarkerStroke = OxyColors.Transparent,
                LineStyle = LineStyle.None,
                MarkerSize = 3,
                MarkerType = MarkerType.Triangle,
                MarkerFill = OxyColors.OrangeRed,
            };
            var ls2_exact = new LineSeries() { Title = "Exact speed", MarkerStrokeThickness = 0.5, Color = OxyColors.Black,  };


            var ls3 = new LineSeries()
            {
                Title = method1 + " pressure",
                MarkerStrokeThickness = 0,
                MarkerStroke = OxyColors.Transparent,
                LineStyle = LineStyle.None,
                MarkerSize = 3,
                MarkerType = MarkerType.Circle,
                MarkerFill = OxyColors.LawnGreen
            };
            var ls3_ = new LineSeries()
            {
                Title = method2 +" speed",
                MarkerStrokeThickness = 0,
                MarkerStroke = OxyColors.Transparent,
                LineStyle = LineStyle.None,
                MarkerSize = 3,
                MarkerType = MarkerType.Triangle,
                MarkerFill = OxyColors.Crimson
            };
            var ls3_exact = new LineSeries() { Title = "Exact pressure", MarkerStrokeThickness = 0.5, Color = OxyColors.Black,  };

            /*
            var ls4 = new LineSeries()
            {
                Title = "energy",
                MarkerStrokeThickness = 0,
                MarkerStroke = OxyColors.Transparent,
                LineStyle = LineStyle.None,
                MarkerSize = 3,
                MarkerType = MarkerType.Circle,
                MarkerFill = OxyColors.Red
            };
            var ls4_ = new LineSeries()
            {
                Title = "energy",
                MarkerStrokeThickness = 0,
                MarkerStroke = OxyColors.Transparent,
                LineStyle = LineStyle.None,
                MarkerSize = 3,
                MarkerType = MarkerType.Square,
                MarkerFill = OxyColors.Crimson
            };
            var ls4_exact = new LineSeries() { Title = "energy_exact", MarkerStroke = OxyColors.Black };
            */

            double[] density = new double[N];
            double[] speed = new double[N];
            double[] pressure = new double[N];
            double[] energy = new double[N];

            for (int i = 0; i < N; i++)
            {
                double[] t1 = new double[(int)Vector_index.M];
                double[] temp2 = new double[(int)Vector_index.M];
                for (int j = 0; j < (int)Vector_index.M; j++)
                {
                    t1[j] = v_cons1[i, j];
                    temp2[j] = v_cons2[i, j];
                }
                var v_ncons1 = ConvertConsToPrimitive(t1, GAMMA);
                var v_ncons2 = ConvertConsToPrimitive(temp2, GAMMA);


                var te1 = v_ncons1[2] / v_ncons1[0] / (GAMMA - 1.0);
                var te2 = v_ncons2[2] / v_ncons2[0] / (GAMMA - 1.0);

                ls1.Points.Add(new DataPoint(xc[i], v_ncons1[0]));
                ls1_exact.Points.Add(new DataPoint(xc[i], exact[i, 0]));

                ls2.Points.Add(new DataPoint(xc[i], v_ncons1[2]));
                ls2_exact.Points.Add(new DataPoint(xc[i], exact[i, 2]));

                ls3.Points.Add(new DataPoint(xc[i], v_ncons1[1]));
                ls3_exact.Points.Add(new DataPoint(xc[i], exact[i, 1]));

                //ls4.Points.Add(new OxyPlot.DataPoint(xc[i], te1));
                //ls4_exact.Points.Add(new OxyPlot.DataPoint(xc[i], exact[i, 2] / exact[i, 0] / (GAMMA - 1.0)));



                ls1_.Points.Add(new DataPoint(xc[i], v_ncons2[0]));

                ls2_.Points.Add(new DataPoint(xc[i], v_ncons2[2]));

                ls3_.Points.Add(new DataPoint(xc[i], v_ncons2[1]));

                //ls4_.Points.Add(new OxyPlot.DataPoint(xc[i], te2));


            }

            var linearAxis1 = new LinearAxis();
            linearAxis1.EndPosition = 0.99;
            linearAxis1.StartPosition = 0.67;
            linearAxis1.MajorGridlineStyle = LineStyle.Solid;
            linearAxis1.MinorGridlineStyle = LineStyle.Dot;
            linearAxis1.Title = "density";
            linearAxis1.Key = "Series 1";
            ls1.YAxisKey = "Series 1";
            ls1_.YAxisKey = "Series 1";
            ls1_exact.YAxisKey = "Series 1";
            

            pm.Axes.Add(linearAxis1);
            pm.Series.Add(ls1_exact);
            pm.Series.Add(ls1);
            pm.Series.Add(ls1_);


            var linearAxis2 = new LinearAxis();
            linearAxis2.EndPosition = 0.65;
            linearAxis2.StartPosition = 0.34;
            linearAxis2.MajorGridlineStyle = LineStyle.Solid;
            linearAxis2.MinorGridlineStyle = LineStyle.Dot;
            linearAxis2.Title = "pressure";
            linearAxis2.Key = "Series 2";
            ls2.YAxisKey = "Series 2";
            ls2_.YAxisKey = "Series 2";
            ls2_exact.YAxisKey = "Series 2";

            pm.Axes.Add(linearAxis2);
            pm.Series.Add(ls2_exact);
            pm.Series.Add(ls2);
            pm.Series.Add(ls2_);

            var linearAxis3 = new LinearAxis();
            linearAxis3.EndPosition = 0.32;
            linearAxis3.StartPosition = 0.01;
            linearAxis3.MajorGridlineStyle = LineStyle.Solid;
            linearAxis3.MinorGridlineStyle = LineStyle.Dot;
            linearAxis3.Title = "speed";
            linearAxis3.Key = "Series 3";
            ls3.YAxisKey = "Series 3";
            ls3_.YAxisKey = "Series 3";
            ls3_exact.YAxisKey = "Series 3";

            pm.Axes.Add(linearAxis3);
            pm.Series.Add(ls3_exact);
            pm.Series.Add(ls3);
            pm.Series.Add(ls3_);

            /*
            var linearAxis4 = new LinearAxis();
            linearAxis4.EndPosition = 0.23;
            linearAxis4.StartPosition = 0.01;
            linearAxis4.MajorGridlineStyle = LineStyle.Solid;
            linearAxis4.MinorGridlineStyle = LineStyle.Dot;
            linearAxis4.Title = "energy";
            linearAxis4.Key = "Series 4";
            ls4.YAxisKey = "Series 4";
            ls4_.YAxisKey = "Series 4";
            ls4_exact.YAxisKey = "Series 4";

            pm.Axes.Add(linearAxis4);
            pm.Series.Add(ls4);
            pm.Series.Add(ls4_);
            pm.Series.Add(ls4_exact);
            */

            pm.Background = OxyColors.White;

            pm.Legends.Add(new Legend
            {
                LegendPlacement = LegendPlacement.Outside,
                LegendPosition = LegendPosition.RightMiddle,
                LegendOrientation = LegendOrientation.Vertical,
            });

            var pngExporter = new OxyPlot.WindowsForms.PngExporter { Width = 1000, Height = 1000 };
            var stream = File.Create(name);

            pngExporter.Export(pm, stream);
            stream.Close();
        }

        protected static void PlotSod(Parameters param, double[] xc, double[,] v_cons, double time, int N, double GAMMA, string method, string path)
        {
            double[,] exact = ExactSolution.Calculate(param, time);



            var name = path + method + Math.Round(time * 1000) + "ms" + ".png";

            var pm = new PlotModel()
            {
                Title = $" time = {Math.Round(time * 1000)} ms",
            };
            var ls1 = new LineSeries()
            {
                Title = "Density",
                MarkerStrokeThickness = 0,
                MarkerStroke = OxyColors.Transparent,
                LineStyle = LineStyle.None,
                MarkerSize = 2,
                MarkerType = MarkerType.Circle,
                MarkerFill = OxyColors.Green
            };
            var ls1_exact = new LineSeries() { Title = "EXact density", Color = OxyColors.Black, };


            //поменял местами для красоты
            var ls2 = new LineSeries()
            {
                Title = "Speed",
                MarkerStrokeThickness = 0,
                MarkerStroke = OxyColors.Transparent,
                LineStyle = LineStyle.None,
                MarkerSize = 2,
                MarkerType = MarkerType.Circle,
                MarkerFill = OxyColors.Blue
            };
            var ls2_exact = new LineSeries() { Title = "Exact speed", Color = OxyColors.Black, };


            var ls3 = new LineSeries()
            {
                Title = "Pressure",
                MarkerStrokeThickness = 0,
                MarkerStroke = OxyColors.Transparent,
                LineStyle = LineStyle.None,
                MarkerSize = 2,
                MarkerType = MarkerType.Circle,
                MarkerFill = OxyColors.Red
            };
            var ls3_exact = new LineSeries() { Title = "Exact pressure", Color = OxyColors.Black, };

            /*
            var ls4 = new LineSeries() { Title = "energy",
                MarkerStrokeThickness = 0,
                MarkerStroke = OxyColors.Transparent,
                LineStyle = LineStyle.None,
                MarkerSize = 3,
                MarkerType = MarkerType.Circle,
                MarkerFill = OxyColors.Red
            };
            var ls4_exact = new LineSeries() { Title = "energy_exact", MarkerStroke = OxyColors.Black };
            */

            double[] density = new double[N];
            double[] speed = new double[N];
            double[] pressure = new double[N];
            double[] energy = new double[N];

            double[] v_ncons = new double[3];
            for (int i = 0; i < N; i++)
            {
                double[] t = new double[3];
                for (int j = 0; j < 3; j++)
                    t[j] = v_cons[i, j];
                v_ncons = ConvertConsToPrimitive(t, GAMMA);

                density[i] = v_ncons[0];
                speed[i] = v_ncons[1];
                pressure[i] = v_ncons[2];

                var te = v_ncons[2] / v_ncons[0] / (GAMMA - 1.0);
                energy[i] = te;

                ls1.Points.Add(new DataPoint(xc[i], v_ncons[0]));
                ls1_exact.Points.Add(new DataPoint(xc[i], exact[i, 0]));

                ls2.Points.Add(new DataPoint(xc[i], v_ncons[2]));
                ls2_exact.Points.Add(new DataPoint(xc[i], exact[i, 2]));

                ls3.Points.Add(new DataPoint(xc[i], v_ncons[1]));
                ls3_exact.Points.Add(new DataPoint(xc[i], exact[i, 1]));

                //ls4.Points.Add(new OxyPlot.DataPoint(xc[i], te));
                //ls4_exact.Points.Add(new OxyPlot.DataPoint(xc[i], exact[i, 2]/ exact[i, 0] / (GAMMA - 1.0)));


            }

            var linearAxis1 = new LinearAxis();
            linearAxis1.EndPosition = 0.99;
            linearAxis1.StartPosition = 0.67;
            linearAxis1.MajorGridlineStyle = LineStyle.Solid;
            linearAxis1.MinorGridlineStyle = LineStyle.Dot;
            
            linearAxis1.Title = "density";
            linearAxis1.Key = "Series 1";
            ls1.YAxisKey = "Series 1";
            ls1_exact.YAxisKey = "Series 1";

            pm.Axes.Add(linearAxis1);
            pm.Series.Add(ls1_exact);
            pm.Series.Add(ls1);

            var linearAxis2 = new LinearAxis();
            linearAxis2.EndPosition = 0.65;
            linearAxis2.StartPosition = 0.34;
            linearAxis2.MajorGridlineStyle = LineStyle.Solid;
            linearAxis2.MinorGridlineStyle = LineStyle.Dot;
            linearAxis2.Title = "pressure";
            linearAxis2.Key = "Series 2";
            ls2.YAxisKey = "Series 2";
            ls2_exact.YAxisKey = "Series 2";

            pm.Axes.Add(linearAxis2);
            pm.Series.Add(ls2_exact);
            pm.Series.Add(ls2);

            var linearAxis3 = new LinearAxis();
            linearAxis3.EndPosition = 0.32;
            linearAxis3.StartPosition = 0.01;
            linearAxis3.MajorGridlineStyle = LineStyle.Solid;
            linearAxis3.MinorGridlineStyle = LineStyle.Dot;
            linearAxis3.Title = "speed";
            linearAxis3.Key = "Series 3";
            ls3.YAxisKey = "Series 3";
            ls3_exact.YAxisKey = "Series 3";

            pm.Axes.Add(linearAxis3);
            pm.Series.Add(ls3_exact);
            pm.Series.Add(ls3);


            pm.Legends.Add(new Legend
             {
                 LegendPlacement = LegendPlacement.Outside,
                 LegendPosition = LegendPosition.BottomCenter,
                 LegendOrientation = LegendOrientation.Horizontal,
             });


            /*
            var linearAxis4 = new LinearAxis();
            linearAxis4.EndPosition = 0.23;
            linearAxis4.StartPosition = 0.01;
            linearAxis4.MajorGridlineStyle = LineStyle.Solid;
            linearAxis4.MinorGridlineStyle = LineStyle.Dot;
            linearAxis4.Title = "energy";
            linearAxis4.Key = "Series 4";
            //ls4.YAxisKey = "Series 4";
            //ls4_exact.YAxisKey = "Series 4";

            pm.Axes.Add(linearAxis4);
            //pm.Series.Add(ls4);
            //pm.Series.Add(ls4_exact);
            */

            pm.Background = OxyColors.White;

            var pngExporter = new OxyPlot.WindowsForms.PngExporter { Width = 1000, Height = 1000 };
            var stream = File.Create(name);

            pngExporter.Export(pm, stream);
            stream.Close();
        }

    }
}
