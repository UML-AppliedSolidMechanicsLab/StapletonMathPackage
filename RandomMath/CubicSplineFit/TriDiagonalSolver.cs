using System;
using System.Collections.Generic;
using System.Text;

namespace RandomMath.CubicSplineFit
{
    public class TriDiagonalSolver
    {
		public TriDiagonalSolver()
        {
        }

        #region
        /// <summary>
        /// Solve the tridiagonal system Ax=r for x.  Based on Numberical Methods for Engineers page 289
        /// </summary>
        /// <param name="e">vector of values in A below the diagonal: from index 1 to n-1</param>
        /// <param name="f">vector of valuess in A at the diagonal: from index 0 to n-1</param>
        /// <param name="g">vector of valuess in A above the diagonal: from index 0 to n-2</param>
        /// <param name="r">vector of values on the other side of the equation (r): from index 0 to n-1</param>
        /// <returns>the answer (x): from index 0 to n-1</returns>
        #endregion
        public static double[] Solve(double[] e, double[] f, double[] g, double[] r)
        {
            int n = f.Length;

            Decompostion(ref e, ref f, g);

            double[] x = Substitution(e, f, g, r);

            return x;
        }

        #region
        /// <summary>
        /// Decomposition decomposes the tridiagonal system
        /// </summary>
        /// <param name="e">vector of values below the diagonal: from index 1 to n-1</param>
        /// <param name="f">vector of values at the diagonal: from index 0 to n-1</param>
        /// <param name="g">vector of values above the diagonal: from index 1 to n-2</param>
        #endregion
        private static void Decompostion(ref double[] e, ref double[] f, double[] g)
        {
            int n = f.Length;

            for (int k = 1; k < n; k++)
            {

                e[k] = e[k] / f[k - 1];

                f[k] += -e[k] * g[k - 1];
            }

        }

        private static double[] Substitution(double[] e, double[] f, double[] g, double[] r)
        {
            int n = f.Length;

            double[] x = new double[n];

            for (int k = 1; k < n; k++)
            {

                r[k] += -e[k] * r[k - 1];

            }

            x[n - 1] = r[n - 1] / f[n - 1];

            for (int k = n - 2; k > -1; k--)
            {

                x[k] = (r[k] - g[k] * x[k + 1]) / f[k];

            }

            return x;
        }
    }
}