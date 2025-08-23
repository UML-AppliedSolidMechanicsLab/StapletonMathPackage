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
            Assert.That(myB[0, 0], Is.EqualTo(-6.44444444444444).Within(acceptablePrecision));
            Assert.That(myB[1, 0], Is.EqualTo(-1.66666666666667).Within(acceptablePrecision));
            Assert.That(myB[2, 0], Is.EqualTo(9.88888888888889).Within(acceptablePrecision));
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

            Assert.That(myC[0, 0], Is.EqualTo(0.25).Within(acceptablePrecision));
            Assert.That(myC[0, 1], Is.EqualTo(0.5).Within(acceptablePrecision));
            Assert.That(myC[1, 0], Is.EqualTo(0.4375).Within(acceptablePrecision));
            Assert.That(myC[1, 1], Is.EqualTo(0.8125).Within(acceptablePrecision));
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

            Assert.That(myC[0, 0], Is.EqualTo(-0.32).Within(acceptablePrecision));
            Assert.That(myC[0, 1], Is.EqualTo(0.52).Within(acceptablePrecision));
            Assert.That(myC[0, 2], Is.EqualTo(-0.48).Within(acceptablePrecision));
            Assert.That(myC[1, 0], Is.EqualTo(-0.48).Within(acceptablePrecision));
            Assert.That(myC[1, 1], Is.EqualTo(-0.72).Within(acceptablePrecision));
            Assert.That(myC[1, 2], Is.EqualTo(1.28).Within(acceptablePrecision));
            Assert.That(myC[2, 0], Is.EqualTo(0.64).Within(acceptablePrecision));
            Assert.That(myC[2, 1], Is.EqualTo(-0.04).Within(acceptablePrecision));
            Assert.That(myC[2, 2], Is.EqualTo(-0.04).Within(acceptablePrecision));
        }
        [Test]
        public void testSymmetricBandedSolver()
        {
            double[,] myA = new double[5, 2] { { 1, 2 }, { 3, 4 }, { 12, 1 }, { 2, 2 }, { 5, 0 } };
            double[] myB = new double[5] { 1, 0, 14, 76, -2 };
            double[] myC = MatrixMath.SymmetricBandedSolver(myA, myB);
            double acceptablePrecision = 1.0E-9;

            Assert.That(myC[0], Is.EqualTo(14.079754601226993865).Within(acceptablePrecision));
            Assert.That(myC[1], Is.EqualTo(-6.5398773006134969325).Within(acceptablePrecision));
            Assert.That(myC[2], Is.EqualTo(-2.1349693251533742331).Within(acceptablePrecision));
            Assert.That(myC[3], Is.EqualTo(65.779141104294478528).Within(acceptablePrecision));
            Assert.That(myC[4], Is.EqualTo(-26.711656441717791411).Within(acceptablePrecision));

            //Try another test, a little bit hairier!!

            double[,] myA1 = new double[8, 4]{{1,5,9,2},{5,1,9,-7},{4,7,-9,11},
                {1,12,1,0},{1,8,9,89},{4,1,6,0},{5,2,0,0},{7,0,0,0}};
            double[] myB1 = new double[8] { 6, 7, 8, 9, 0, 2, 3, 4 };
            double[] myC1 = MatrixMath.SymmetricBandedSolver(myA1, myB1);

            Assert.That(myC1[0], Is.EqualTo(3.6792959471903019873).Within(acceptablePrecision));
            Assert.That(myC1[1], Is.EqualTo(-.11799136471442218435).Within(acceptablePrecision));
            Assert.That(myC1[2], Is.EqualTo(.58407892062480899170).Within(acceptablePrecision));
            Assert.That(myC1[3], Is.EqualTo(-1.1730247046207359954).Within(acceptablePrecision));
            Assert.That(myC1[4], Is.EqualTo(.11905421305965486397).Within(acceptablePrecision));
            Assert.That(myC1[5], Is.EqualTo(-1.6408479084195896295).Within(acceptablePrecision));
            Assert.That(myC1[6], Is.EqualTo(.59635467612136649060).Within(acceptablePrecision));
            Assert.That(myC1[7], Is.EqualTo(.29379330513793170041).Within(acceptablePrecision));

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

            Assert.That(myC[0, 0], Is.EqualTo(1).Within(acceptablePrecision));
            Assert.That(myC[0, 1], Is.EqualTo(2).Within(acceptablePrecision));
            Assert.That(myC[1, 0], Is.EqualTo(9).Within(acceptablePrecision));
            Assert.That(myC[1, 1], Is.EqualTo(10).Within(acceptablePrecision));
            Assert.That(myC[2, 0], Is.EqualTo(5).Within(acceptablePrecision));
            Assert.That(myC[2, 1], Is.EqualTo(6).Within(acceptablePrecision));
            Assert.That(myC[3, 0], Is.EqualTo(13).Within(acceptablePrecision));
            Assert.That(myC[3, 1], Is.EqualTo(14).Within(acceptablePrecision));

            myC = MatrixMath.SwapCols(myC, 1, 2);

            Assert.That(myC[0, 0], Is.EqualTo(1).Within(acceptablePrecision));
            Assert.That(myC[0, 2], Is.EqualTo(2).Within(acceptablePrecision));
            Assert.That(myC[1, 0], Is.EqualTo(9).Within(acceptablePrecision));
            Assert.That(myC[1, 2], Is.EqualTo(10).Within(acceptablePrecision));
            Assert.That(myC[2, 0], Is.EqualTo(5).Within(acceptablePrecision));
            Assert.That(myC[2, 2], Is.EqualTo(6).Within(acceptablePrecision));
            Assert.That(myC[3, 0], Is.EqualTo(13).Within(acceptablePrecision));
            Assert.That(myC[3, 2], Is.EqualTo(14).Within(acceptablePrecision));

            myC = MatrixMath.SwapRows(myA, 0, 0, 2, 3);

            Assert.That(myC[0, 0], Is.EqualTo(9).Within(acceptablePrecision));
            Assert.That(myC[0, 1], Is.EqualTo(10).Within(acceptablePrecision));
            Assert.That(myC[1, 0], Is.EqualTo(13).Within(acceptablePrecision));
            Assert.That(myC[1, 1], Is.EqualTo(14).Within(acceptablePrecision));
            Assert.That(myC[2, 0], Is.EqualTo(5).Within(acceptablePrecision));
            Assert.That(myC[2, 1], Is.EqualTo(6).Within(acceptablePrecision));
            Assert.That(myC[3, 0], Is.EqualTo(1).Within(acceptablePrecision));
            Assert.That(myC[3, 1], Is.EqualTo(2).Within(acceptablePrecision));

            myC = MatrixMath.SwapCols(myA, 1, 1, 2, 3);

            Assert.That(myC[0, 0], Is.EqualTo(1).Within(acceptablePrecision));
            Assert.That(myC[0, 1], Is.EqualTo(3).Within(acceptablePrecision));
            Assert.That(myC[0, 2], Is.EqualTo(4).Within(acceptablePrecision));
            Assert.That(myC[0, 3], Is.EqualTo(2).Within(acceptablePrecision));
            Assert.That(myC[1, 0], Is.EqualTo(5).Within(acceptablePrecision));
            Assert.That(myC[1, 1], Is.EqualTo(7).Within(acceptablePrecision));
            Assert.That(myC[1, 2], Is.EqualTo(8).Within(acceptablePrecision));
            Assert.That(myC[1, 3], Is.EqualTo(6).Within(acceptablePrecision));
        }
    }
}
  