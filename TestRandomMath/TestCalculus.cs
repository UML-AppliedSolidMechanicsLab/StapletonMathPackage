using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using RandomMath;


namespace TestRandomMath
{
    [TestFixture]
    public class TestCalculus
    {

        [Test]
        public void TestTrapezoidRule()
        {

            TestFunction f = new TestFunction();
            double[,] integral = Calculus.TrapezoidMatrix(0, 1, 1000, f);
            double acceptablePrecision = 0.000001;

            Assert.That(integral[0, 0], Is.EqualTo(0.4596976941).Within(acceptablePrecision));
            Assert.That(integral[1, 0], Is.EqualTo(0.8414709848).Within(acceptablePrecision));
            Assert.That(integral[2, 0], Is.EqualTo(3.951615162).Within(acceptablePrecision));
        }

        [Test]
        public void TestRhomberg()
        {
            TestFunction f = new TestFunction();
            double acceptablePrecision = 0.00000001;
            double[,] integral = Calculus.RhombergMatrix(0, 1, 1000, acceptablePrecision, f);


            Assert.That(integral[0, 0], Is.EqualTo(0.4596976941).Within(acceptablePrecision));
            Assert.That(integral[1, 0], Is.EqualTo(0.8414709848).Within(acceptablePrecision));
            Assert.That(integral[2, 0], Is.EqualTo(3.951615162).Within(acceptablePrecision));

        }

        [Test]
        public void TestDoubleIntegration()
        {
            TestFunction f = new TestFunction();
            double acceptablePrecision = 0.00000001;
            double[,] integral = Calculus.DoubleIntegration(0.0, 1.0, 0.0, 1.0, 1000, acceptablePrecision, f);


            Assert.That(integral[0, 0], Is.EqualTo(1.142639664).Within(acceptablePrecision));
            Assert.That(integral[1, 0], Is.EqualTo(1.445884302).Within(acceptablePrecision));
            Assert.That(integral[2, 0], Is.EqualTo(5.019159109).Within(acceptablePrecision));

        }
    }

    public class TestFunction : IMatrixFunctionOneVariable, IMatrixFunctionTwoVariables
    {
        public TestFunction() { }
        public double[,] Eval(double x)
        {
            double[,] b = new double[3, 1];
            b[0, 0] = Math.Sin(x);
            b[1, 0] = Math.Cos(x);
            b[2, 0] = Math.Pow(x, 2.0) + 3.8 * x + Math.Exp(x);
            return b;
        }
        public double[,] Eval(double x, double z)
        {
            double[,] b = new double[3, 1];
            b[0, 0] = x * Math.Sin(x) + Math.Cos(z);
            b[1, 0] = Math.Exp(x) * Math.Cos(z);
            b[2, 0] = z * Math.Pow(x, 2.0) + 3.8 * x + Math.Exp(z) * Math.Exp(x);
            return b;
        }
    }
}
