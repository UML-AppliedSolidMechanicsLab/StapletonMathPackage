using NUnit.Framework;
using RandomMath;

namespace TestRandomMath
{
    [TestFixture]

    public class testMath
    {

        double[,] coeffs = new double[3, 3];
        double[] multMatrix = new double[3];

        double[,] myA = new double[3, 3];
        double[,] myB = new double[3, 1];
        int myn;
        int myC;

        [Test]
        public void testGauss()
        {
            myn = 3;
            myC = 0;
            myA[0, 0] = 1.0;
            myA[0, 1] = 5.0;
            myA[0, 2] = 2.0;
            myA[1, 0] = 5.0;
            myA[1, 1] = 2.0;
            myA[1, 2] = 4.0;
            myA[2, 0] = 3.0;
            myA[2, 1] = 2.0;
            myA[2, 2] = 3.0;

            myB[0, 0] = 5.0;
            myB[1, 0] = 4.0;
            myB[2, 0] = 7.0;
            MatrixMath.Gauss(myA, myB, myn, myC);
            double acceptablePrecision = 0.000000000001;
            Assert.AreEqual(-6.44444444444444, myB[0, 0], acceptablePrecision);
            Assert.AreEqual(-1.66666666666667, myB[1, 0], acceptablePrecision);
            Assert.AreEqual(9.88888888888889, myB[2, 0], acceptablePrecision);
        }
        [Test]
        public void testMatrixMultiplication()
        {
            double[,] myA = new double[2, 2];
            double[,] myB = new double[2, 2];
            double[,] myC;

            double acceptablePrecision = 1.0E-9;

            myA[0, 0] = 0.5;
            myA[0, 1] = 0.25;
            myA[1, 0] = 0.75;
            myA[1, 1] = 0.5;

            myB[0, 0] = 0.25;
            myB[0, 1] = 0.75;
            myB[1, 0] = 0.5;
            myB[1, 1] = 0.5;

            myC = MatrixMath.Multiply(myA, myB);

            Assert.AreEqual(0.25, myC[0, 0], acceptablePrecision);
            Assert.AreEqual(0.5, myC[0, 1], acceptablePrecision);
            Assert.AreEqual(0.4375, myC[1, 0], acceptablePrecision);
            Assert.AreEqual(0.8125, myC[1, 1], acceptablePrecision);
        }
        [Test]
        public void testInvertMatrix()
        {
            double[,] myA = new double[3, 3];
            double[,] myC = new double[3, 3];
            double acceptablePrecision = 1.0E-9;

            myA[0, 0] = 0.5;
            myA[0, 1] = 0.25;
            myA[0, 2] = 2.0;
            myA[1, 0] = 5.0;
            myA[1, 1] = 2.0;
            myA[1, 2] = 4.0;
            myA[2, 0] = 3.0;
            myA[2, 1] = 2.0;
            myA[2, 2] = 3.0;

            myC = MatrixMath.InvertMatrix(myA);

            Assert.AreEqual(-0.32, myC[0, 0], acceptablePrecision);
            Assert.AreEqual(0.52, myC[0, 1], acceptablePrecision);
            Assert.AreEqual(-0.48, myC[0, 2], acceptablePrecision);
            Assert.AreEqual(-0.48, myC[1, 0], acceptablePrecision);
            Assert.AreEqual(-0.72, myC[1, 1], acceptablePrecision);
            Assert.AreEqual(1.28, myC[1, 2], acceptablePrecision);
            Assert.AreEqual(0.64, myC[2, 0], acceptablePrecision);
            Assert.AreEqual(-0.04, myC[2, 1], acceptablePrecision);
            Assert.AreEqual(-0.04, myC[2, 2], acceptablePrecision);
        }
        [Test]
        public void testSymmetricBandedSolver()
        {
            double[,] myA = new double[5, 2] { { 1, 2 }, { 3, 4 }, { 12, 1 }, { 2, 2 }, { 5, 0 } };
            double[] myB = new double[5] { 1, 0, 14, 76, -2 };
            double[] myC = MatrixMath.SymmetricBandedSolver(myA, myB);
            double acceptablePrecision = 1.0E-9;

            Assert.AreEqual(14.079754601226993865, myC[0], acceptablePrecision);
            Assert.AreEqual(-6.5398773006134969325, myC[1], acceptablePrecision);
            Assert.AreEqual(-2.1349693251533742331, myC[2], acceptablePrecision);
            Assert.AreEqual(65.779141104294478528, myC[3], acceptablePrecision);
            Assert.AreEqual(-26.711656441717791411, myC[4], acceptablePrecision);

            //Try another test, a little bit hairier!!

            double[,] myA1 = new double[8, 4]{{1,5,9,2},{5,1,9,-7},{4,7,-9,11},
                {1,12,1,0},{1,8,9,89},{4,1,6,0},{5,2,0,0},{7,0,0,0}};
            double[] myB1 = new double[8] { 6, 7, 8, 9, 0, 2, 3, 4 };
            double[] myC1 = MatrixMath.SymmetricBandedSolver(myA1, myB1);

            Assert.AreEqual(3.6792959471903019873, myC1[0], acceptablePrecision);
            Assert.AreEqual(-.11799136471442218435, myC1[1], acceptablePrecision);
            Assert.AreEqual(.58407892062480899170, myC1[2], acceptablePrecision);
            Assert.AreEqual(-1.1730247046207359954, myC1[3], acceptablePrecision);
            Assert.AreEqual(.11905421305965486397, myC1[4], acceptablePrecision);
            Assert.AreEqual(-1.6408479084195896295, myC1[5], acceptablePrecision);
            Assert.AreEqual(.59635467612136649060, myC1[6], acceptablePrecision);
            Assert.AreEqual(.29379330513793170041, myC1[7], acceptablePrecision);

        }
        [Test]
        public void testSwapRows()
        {
            int n = 4;
            double[,] myA = new double[n, n];
            double acceptablePrecision = 1.0E-9;

            double k = 1;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    myA[i, j] = k;
                    k++;
                }
            }

            double[,] myC = MatrixMath.SwapRows(myA, 1, 2);

            Assert.AreEqual(1, myC[0, 0], acceptablePrecision);
            Assert.AreEqual(2, myC[0, 1], acceptablePrecision);
            Assert.AreEqual(9, myC[1, 0], acceptablePrecision);
            Assert.AreEqual(10, myC[1, 1], acceptablePrecision);
            Assert.AreEqual(5, myC[2, 0], acceptablePrecision);
            Assert.AreEqual(6, myC[2, 1], acceptablePrecision);
            Assert.AreEqual(13, myC[3, 0], acceptablePrecision);
            Assert.AreEqual(14, myC[3, 1], acceptablePrecision);

            myC = MatrixMath.SwapCols(myC, 1, 2);

            Assert.AreEqual(1, myC[0, 0], acceptablePrecision);
            Assert.AreEqual(2, myC[0, 2], acceptablePrecision);
            Assert.AreEqual(9, myC[1, 0], acceptablePrecision);
            Assert.AreEqual(10, myC[1, 2], acceptablePrecision);
            Assert.AreEqual(5, myC[2, 0], acceptablePrecision);
            Assert.AreEqual(6, myC[2, 2], acceptablePrecision);
            Assert.AreEqual(13, myC[3, 0], acceptablePrecision);
            Assert.AreEqual(14, myC[3, 2], acceptablePrecision);

            myC = MatrixMath.SwapRows(myA, 0, 0, 2, 3);

            Assert.AreEqual(9, myC[0, 0], acceptablePrecision);
            Assert.AreEqual(10, myC[0, 1], acceptablePrecision);
            Assert.AreEqual(13, myC[1, 0], acceptablePrecision);
            Assert.AreEqual(14, myC[1, 1], acceptablePrecision);
            Assert.AreEqual(5, myC[2, 0], acceptablePrecision);
            Assert.AreEqual(6, myC[2, 1], acceptablePrecision);
            Assert.AreEqual(1, myC[3, 0], acceptablePrecision);
            Assert.AreEqual(2, myC[3, 1], acceptablePrecision);

            myC = MatrixMath.SwapCols(myA, 1, 1, 2, 3);

            Assert.AreEqual(1, myC[0, 0], acceptablePrecision);
            Assert.AreEqual(3, myC[0, 1], acceptablePrecision);
            Assert.AreEqual(4, myC[0, 2], acceptablePrecision);
            Assert.AreEqual(2, myC[0, 3], acceptablePrecision);
            Assert.AreEqual(5, myC[1, 0], acceptablePrecision);
            Assert.AreEqual(7, myC[1, 1], acceptablePrecision);
            Assert.AreEqual(8, myC[1, 2], acceptablePrecision);
            Assert.AreEqual(6, myC[1, 3], acceptablePrecision);
        }
    }
}
  