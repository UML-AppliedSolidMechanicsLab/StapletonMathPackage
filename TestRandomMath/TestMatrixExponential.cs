using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;


namespace TestRandomMath
{
    [TestFixture]
    public class TestMatrixExponential
    {
        /// <summary>
        /// Purpose: Test matrix exponential
        /// Created By: Scott_Stapleton
        /// Created On: 12/17/2024 8:27:27 PM
        /// </summary>
        #region Protected Members

        [Test]
        public void TestMatrixExponentialStability_SLJ()
        {
            double[,] A = new double[12, 12]
        {
            {  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
            {  0.0255,  0.0,  0.0,  0.01275,  0.0,  0.0, -0.0255,  0.0,  0.0,  0.01275,  0.0,  0.0 },
            {  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
            {  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
            {  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
            {  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
            {  0.153, -0.8571428571428571,  0.0,  0.0765,  0.0,  0.0, -0.153,  0.8571428571428571,  0.0,  0.0765,  0.0,  0.0 },
            {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
            { -0.0255,  0.0,  0.0, -0.01275,  0.0,  0.0,  0.0255,  0.0,  0.0, -0.01275,  0.0,  0.0 },
            {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0 },
            {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0 },
            {  0.153,  0.8571428571428571,  0.0,  0.0765,  0.0,  0.0, -0.153, -0.8571428571428571,  0.0,  0.0765,  0.0,  0.0 },
        };

            RandomMath.ODESystemSolver.MatrixExponential me = new RandomMath.ODESystemSolver.MatrixExponential(A, 80.0, 100, 0.000000001);
            double[,] result = me.Solve(1.0);
            double[,] result2 = me.Solve(80.0);


        }   
        #endregion

    }
}