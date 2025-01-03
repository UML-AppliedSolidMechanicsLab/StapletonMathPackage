using System;
using System.Collections.Generic;
using System.Text;

namespace RandomMath.CubicSplineFit
{
    /// <summary>
	/// Creates and Solves for all of the variables for a cubic spline for a set of input data, X and Y.  Then, the Interpolation 
	/// method allows one to get the y, y', and y'' for a certain x as long as it's within the range of X.  Based on the routine
	/// in Numerical Methods for Engineers by Chapra and Canale, pg 504
	/// </summary>
	public class CubicSplineFit
    {
        #region Private members
        private double[] x;
        private double[] y;
        private int n;
        private double[] e;
        private double[] f;
        private double[] g;
        private double[] r;
        private double[] d2x;
        #endregion

        #region Constructor
        #region
        /// <summary>
        /// Fits a Cubic spline to a data set X and Y.  Interpolation method available to get intermediate y, y', and y'' for a given x
        /// </summary>
        #endregion
        public CubicSplineFit(double[] X, double[] Y)
        {
            x = X;
            y = Y;

            if (X.Length != Y.Length)
            {
                throw new System.ArgumentException("X and Y must be the same length");
            }

            n = x.Length;

            FormTriDiagonal();

            double[] d2xShort = TriDiagonalSolver.Solve(e, f, g, r);

            //This section is because the tridiagonal solver only gets the 2nd der at the interior nodes:
            //We have to copy that to a longer array so that we include the 0 2nd der at the boundary nodes.
            d2x = new double[n];

            Array.Copy(d2xShort, 0, d2x, 1, n - 2);

        }
        #endregion

        #region Private Methods
        private void FormTriDiagonal()
        {
            e = new double[n - 2];
            f = new double[n - 2];
            g = new double[n - 2];
            r = new double[n - 2];

            f[0] = 2 * (x[2] - x[0]);
            g[0] = x[2] - x[1];
            r[0] = 6 / (x[2] - x[1]) * (y[2] - y[1]);
            r[0] += 6 / (x[1] - x[0]) * (y[0] - y[1]);

            for (int i = 1; i < n - 3; i++)
            {

                e[i] = x[i + 1] - x[i];

                f[i] = 2 * (x[i + 2] - x[i]);

                g[i] = (x[i + 2] - x[i + 1]);

                r[i] = 6 / (x[i + 2] - x[i + 1]) * (y[i + 2] - y[i + 1]);
                r[i] += 6 / (x[i + 1] - x[i]) * (y[i] - y[i + 1]);

            }

            e[n - 3] = (x[n - 2] - x[n - 3]);
            f[n - 3] = 2 * (x[n - 1] - x[n - 3]);
            r[n - 3] = 6 / (x[n - 1] - x[n - 2]) * (y[n - 1] - y[n - 2]);
            r[n - 3] += 6 / (x[n - 2] - x[n - 3]) * (y[n - 3] - y[n - 2]);
        }
        #endregion

        #region Public Methods
        #region
        /// <summary>
        /// Given a value of x, returns the cubic spline interpolated value of y
        /// </summary>
        /// <param name="x">some value x (within sets X and Y)</param>
        /// <returns>the interpolated y</returns>
        #endregion
        public double Interpolate(double x)
        {
            double dy = 0.0;
            double d2y = 0.0;

            return Interpolate(x, ref dy, ref d2y);
        }
        #region
        /// <summary>
        /// Interpolates y based on some x
        /// </summary>
        /// <param name="x">value within X set</param>
        /// <param name="dy">returns: 1st derivative of y wrt x at given x</param>
        /// <param name="d2y">returns: 2st derivative of y wrt x at given x</param>
        /// <returns>interpolated value of y at given x</returns>
        #endregion
        public double Interpolate(double inX, ref double dy, ref double d2y)
        {
            double outY = 0.0;

            int flag = 0;
            int i = 1;

            double[] c = new double[4];
            double[] t = new double[4];



            while (flag == 0)
            {

                if ((inX >= x[i - 1] && inX <= x[i]) || (Math.Abs((inX - x[i]) / inX) < 0.0000001))
                {

                    c[0] = d2x[i - 1] / 6 / (x[i] - x[i - 1]);
                    c[1] = d2x[i] / 6 / (x[i] - x[i - 1]);
                    c[2] = y[i - 1] / (x[i] - x[i - 1]) - d2x[i - 1] * (x[i] - x[i - 1]) / 6;
                    c[3] = y[i] / (x[i] - x[i - 1]) - d2x[i] * (x[i] - x[i - 1]) / 6;

                    t[0] = c[0] * System.Math.Pow((x[i] - inX), 3);
                    t[1] = c[1] * System.Math.Pow((inX - x[i - 1]), 3);
                    t[2] = c[2] * (x[i] - inX);
                    t[3] = c[3] * (inX - x[i - 1]);

                    outY = t[0] + t[1] + t[2] + t[3];

                    t[0] = -3 * c[0] * System.Math.Pow((x[i] - inX), 2);
                    t[1] = 3 * c[1] * System.Math.Pow((inX - x[i - 1]), 2);
                    t[2] = -c[2];
                    t[3] = c[3];

                    dy = t[0] + t[1] + t[2] + t[3];

                    t[0] = 6 * c[0] * (x[i] - inX);
                    t[1] = 6 * c[1] * (inX - x[i - 1]);

                    d2y = t[0] + t[2];

                    flag = 1;
                }
                else
                {
                    i++;
                }

                if (i == n)
                {

                    throw new System.ArgumentException("x outside of the data range");
                }
            }

            return outY;
        }
        #endregion
    }
}
